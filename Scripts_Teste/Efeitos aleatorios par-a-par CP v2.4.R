library(pair.mglmm)
library(powell)


# Simulação de dados
phi <- 1; p <- 1.6 

mydat <- data.sim(3, 'CP', exp, xi=p, phi=phi)
dat <- mydat$Data

mainForm0 <- value~-1+variable:period+(-1+variable|ID)
mainForm1 <- value~-1+variable+(-1+variable|ID)


# Modelos par-a-par
cl<-makeCluster(4)
registerDoParallel(cl)
clusterEvalQ(cl, library(pair.mglmm))

system.time(m0.1 <- mglmmCP2(mainForm0, dat$variable, dat))
system.time(m1.1 <- mglmmCP2(mainForm1, dat$variable, dat))

stopCluster(cl)


# Modelos multivariados
system.time(m0 <- cpglmm(mainForm0, data=dat))
system.time(m1 <- cpglmm(mainForm1, data=dat))



# Comparação das log-verossimilhanças
logLik(m0)-logLik(m1)
m0.1$logLik-m1.1$logLik


# Efeitos aleatórios - modas condicionais
m0.1$ranef
ranef(m0)

m1.1$ranef
ranef(m1)


# Valores ajustados
fit0.1 <- m0.1$fitted
fit0 <- fitted(m0)
#cbind(fit0, fit0.1)

fit1.1 <- m1.1$fitted
fit1 <- fitted(m1)
#cbind(fit1, fit1.1)


# Modelos com erro padrão
m0
summary.cp(m0.1)

m1
summary.cp(m1.1)


## Ajuste par-a-par 
mglmmCP2 <-
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

  prm1$est$stdErr <- rcov(prm1$m1, formula, S, phi, p, fit)

  LL <- pair.mglmm:::llik.fim(prm1$m1, formula, beta, S, phi, p)
 
  q0 <- ncol(S)
  df <- ncol(prm1$m1$fr$X) + (q0*(q0+1)/2) + 1

  list(fixef = prm1$est, VarCov = S, phi = phi, p = p, logLik = LL, df = df, 
       ranef = re, fitted = fit, residuals = resid)
}


## Efeitos aleatórios - substituir cll() em R pela versão em Rcpp
cll <- function(u, r, beta, S, phi, p){
  eta <- r$x %*% beta + r$z %*% u
  mu <- exp(eta)
  
  li.u <- mvtnorm::dmvnorm(u, sigma = S, log = T)
  lli <- sum(log(tweedie::dtweedie.series(r$y, p, mu, phi))) + li.u
  
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


## Valores ajustados às observações
fit.mglmm <- function(mod, fixef, ranef){
  x <- mod$fr$X
  z <- t(mod$FL$trms[[1]]$Zt)

  y.hat <- exp(x %*% fixef + z %*% ranef)
  as.vector(y.hat)
}


## Resíduos
resid.mglmm <- function(mod, fit, p){
  y <- mod$fr$Y
  (y-fit)/sqrt(fit^p)
}


## Erro padrão das estimativas de efeitos fixos
rcov <- function(mod, formula, S, phi, p, fit) {
  q <- length(unique(mod$FL$fl[[1]]))

  Lambdat <- chol(S)
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

  L_x <- t(chol(t(X) %*% W %*% X - L_xz %*% t(L_xz)))

  se <- phi * chol2inv(L_x)
  sqrt(diag(se))
}


## Funções convenientes
summary.cp <- function(mod) mod[1:6]
fitted.cp <- function(mod) mod[["fitted"]]
ranef.cp <- function(mod) mod[["ranef"]]
resid.cp <- function(mod) mod[["residuals"]]
logLik.cp <- function(mod) mod[5:6]