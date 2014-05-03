data.sim<-function(m,distr,link.inv,beta=NULL,...){
  if(m<=0) 
    stop("\'m\' must be >= 1")
  
  if(!is.null(beta) & !all(dim(beta)==c(4,m))) 
    stop("dim(beta) must be 4 by m")
  
  if(!distr %in% c('poisson','binomial','gaussian','CP'))
    stop("Please specify either 'poisson', 'binomial', 'gaussian' or 'CP'")
  
  # Fixed and random effects factors
  period<-factor(1:4)
  ID<-factor(1:30)
  
  dat<-expand.grid(period=period,ID=ID)
  
  if(is.null(beta))
    beta<-matrix(runif(4*m,-2,2),4,m)
  X<-model.matrix(~-1+period,dat)
  
  # Function to generate a PD matrix for the RE's
  f.S<-function(n){
    S<-matrix(0,n,n)
    diag(S)<-rep(1,n)
    
    corr<-runif(n*(n-1)/2,-.9,.9)
    S[lower.tri(S)]<-S[upper.tri(S)]<-corr
    
    nearPD(S,corr=T)$mat
  }
  
  S<-f.S(m)
  u<-mvrnorm(length(levels(dat$ID)),rep(0,m),S)
  Z<-model.matrix(~-1+ID,dat)
  
  
  Eta<-X%*%beta+Z%*%u
  Mu<-link.inv(Eta)
  
  Mu<-data.frame(Mu)
  names(Mu)<-paste("C",1:m,sep="")
  
  dat<-cbind(dat,Mu)
  dat<-melt(dat,id=c('period','ID'))
  dat$value<-switch(distr,
                    poisson=rpois(nrow(dat),dat$value),
                    binomial=rbinom(nrow(dat),prob=dat$value,...),
                    gaussian=rnorm(nrow(dat),mean=dat$value,...),
                    CP=tweedie::rtweedie(nrow(dat),mu=dat$value,...)
  )
  
  list(Data=dat,beta=beta,S=S)
}


glmmMultiCP <-
function(formula, id, data, cl = NULL) {
  lev <- paste(unique(id))
  tes <- combn(lev, m = 2)

  glmm.fit <- function(x, data, formula) {
    ind <- id %in% x
    dat0 <- subset(data, ind)
    cpglmm(formula, data = dat0)
  }

  res <- foreach(i = 1:ncol(tes)) %dopar% {
           glmm.fit(tes[,i], data, formula)
         }

  return(res)
}


format0CP <-
function(mod) {
  betas <- fixef(mod)
  df.b <- melt(betas)
  df.b$Parametro <- rownames(df.b)
  rownames(df.b) <- 1:nrow(df.b)

  df.b <- df.b[,c(2,1)]

  S <- VarCorr(mod)
  uS <- unlist(S)
  triS <- unique(uS)

  nameS <- unlist(attr(S[[1]], "dimnames"))
  nameS <- combn(sort(nameS), 2)
  nameS <- unique(t(nameS))
  nameS <- apply(nameS, 1, paste, collapse = ":")

  df.s <- data.frame(Parametro = nameS, value = triS)

  phi <- as.numeric(attr(S, 'sc'))^2
  phi <- data.frame(Parametro = 'phi', 
                    value = ifelse(is.na(phi), 1, phi))

  p <- data.frame(Parametro = 'p', value = mod@p)

  rbind(df.b, df.s, phi, p)
}


format1CP <-
function(df.m, formula, data) {
  m1 <- lmer(formula, data = data, doFit = F)

  matn <- names(m1$fr$fixef)
  ind <- match(matn, df.m$Parameter)
  est <- df.m[ind,]

  ind.phi <- match('phi', df.m$Parameter)
  phi <- df.m[ind.phi,-1]

  ind.p <- match('p', df.m$Parameter)
  p <- df.m[ind.p,-1]

  SigHat <- m1$FL$trms[[1]]$ST
  indS <- !seq(1, nrow(df.m)) %in% c(ind, ind.phi, ind.p)
  SigHat[lower.tri(SigHat, diag = T)] <- df.m[indS,2]
  SigHat <- SigHat + t(lower.tri(SigHat)*SigHat)
  SigHat <- nearPD(SigHat)$mat

  list(m1 = m1, est = est, phi = phi,
       p = p, SigHat = SigHat)
}


llik.fim <-
function (mod, formula, beta, S, phi, p, B = 10000, cl = NULL) {
    mf.ind <- split(mod$fr$mf, mod$FL$fl[[1]])
    mf.r <- lapply(mf.ind, f0.ep, formula = formula)

    q <- ncol(S)
    U <- MASS::mvrnorm(B, rep(0, q), S)

    ll <- foreach(i = 1:length(mf.r), .combine = c) %dopar% {
        r <- mf.r[[i]]
        .Call("mllCPMC", r$y, beta, r$x, r$z, U, phi, p, PACKAGE = "pair.mglmm")
    }

    sum(log(ll))
}


f0.ep <-
function (data, formula) {
    formula0 <- eval(formula)
    fix.form <- lme4:::nobars(formula)
    fr <- model.frame(fix.form, data = data)

    Fr <- list()
    Fr$Y <- model.response(fr)
    Fr$X <- model.matrix(attr(fr, "terms"), data = fr)
    Fr$wts <- rep(1, nrow(Fr$X))
    Fr$off <- rep(0, nrow(Fr$X))
    Fr$mf <- data
    Fr$fixef <- rep(0, ncol(Fr$X))

    fl <- lme4:::lmerFactorList(formula0, Fr, 0L, 0L)

    x <- Fr$X; y <- Fr$Y
    z <- t(as.matrix(fl$trms[[1]]$Zt))

    list(y = y, x = x, z = z)
}


mglmmCP <-
function (formula, id, data, cl = NULL) {
  m0 <- pair.mglmm:::glmmMultiCP(formula, id, data)

  prm0 <- lapply(m0, pair.mglmm:::format0CP)
  prm0 <- do.call(rbind, prm0)
  prm0 <- aggregate(value~Parametro, prm0, mean)
  names(prm0) <- c('Parameter', 'Estimate')

  prm1 <- pair.mglmm:::format1CP(prm0, formula, data)

  beta <- as.vector(prm1$est[,2])
  S <- as.matrix(prm1$SigHat)
  phi <- as.vector(prm1$phi)
  p <- as.vector(prm1$p)

  re <- re.mglmm(prm1$m1, formula, beta, S, phi, p)
  fit <- fit.mglmm(prm1$m1, beta, as.vector(re))
  resid <- resid.mglmm(prm1$m1, fit, p)
  
  VC <- rcov(prm1$m1, formula, S, phi, p, fit)
  prm1$est$stdErr <- sqrt(diag(VC))

  LL <- pair.mglmm:::llik.fim(prm1$m1, formula, beta, S, phi, p)
 
  q0 <- ncol(S)
  df <- ncol(prm1$m1$fr$X) + (q0*(q0+1)/2) + 1

  list(fixef = prm1$est, VarCov = S, phi = phi, p = p, logLik = LL, df = df, 
       ranef = re, fitted = fit, residuals = resid, VC.fixef = VC, frames = prm1$m1,
       formula = formula)
}


cll<-function(u, r, beta, S, phi, p){
  lli <- .Call("cll_call", r, beta, u, S, phi, p, PACKAGE = "pair.mglmm") 
  -lli
}


uhat <- function(r, beta, S, phi, p){
  q <- ncol(S)
  ini <- rep(0, q)

  u.hat <- powell(par = ini, fn = cll, r = r, beta = beta, 
                  S = S, phi = phi, p = p)
  u.hat$par
}


re.mglmm <- function(mod, formula, beta, S, phi, p){
  mf.ind <- split(mod$fr$mf, mod$FL$fl[[1]])
  mf.r <- lapply(mf.ind, pair.mglmm:::f0.ep, formula = formula)
 
  re <- lapply(mf.r, function(r) uhat(r, beta, S, phi, p))
  do.call(rbind, re)
}


fit.mglmm <- function(mod, fixef, ranef){
  x <- mod$fr$X
  z <- t(mod$FL$trms[[1]]$Zt)

  y.hat <- exp(x %*% fixef + z %*% ranef)
  as.vector(y.hat)
}


resid.mglmm <- function(mod, fit, p){
  y <- mod$fr$Y
  (y-fit)/sqrt(fit^p)
}


rcov <- function(mod, formula, S, phi, p, fit) {
  q <- length(unique(mod$FL$fl[[1]]))

  Lambdat <- chol(S/phi)
  Lambdat <- lapply(1:q, function(i) Lambdat)
  Lambdat <- bdiag(Lambdat)


  fzt <- lme4:::findbars(formula[[3]])[[1]]
  formZt <- paste("~", deparse(fzt[[2]]), ":", fzt[[3]], collapse = "")
  Zt <- t(model.matrix(as.formula(formZt), mod$fr$mf))

  W <- diag(fit ^ (2 - p))

  LtZt <- Lambdat %*% Zt

  LtZtWZL <- LtZt %*% W %*% t(LtZt)
  q <- ncol(LtZtWZL)
 
  L_z <- t(chol(LtZtWZL + diag(1, q)))

  X <- mod$fr$X
  L_xz <- t(solve(L_z) %*% LtZt %*% W %*% X)

  L_x <- chol(t(X) %*% W %*% X - L_xz %*% t(L_xz))

  phi * chol2inv(L_x)
}


summary.cp <- function(mod) mod[1:6]
fitted.cp <- function(mod) mod[["fitted"]]
ranef.cp <- function(mod) mod[["ranef"]]
resid.cp <- function(mod) mod[["residuals"]]
logLik.cp <- function(mod) mod[5:6]