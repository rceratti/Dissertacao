library(lme4)

pirls <- function(y, Z, X, beta, Lambda) {
  ZL <- Z %*% Lambda
  new <- rep(0, ncol(ZL))
  err <- 100
  k <- 1

  while(err > 1e-8) {
    old <- new

    eta <- as.vector(X %*% beta + ZL %*% old)
    mu <- exp(eta)
    W <- diag(mu); I <- diag(1, ncol(ZL))
    H <- t(ZL) %*% W %*% ZL + I
    L <- chol(H)

    q <- (y-mu)/mu
    g <- (t(ZL) %*% W %*% q) - old

    new <- old + as.vector(chol2inv(L) %*% g)

    err <- max(abs(new - old))
    k <- k + 1
  }

  list(par = new, L = L)
}


cll <- function(u, y, x, z, beta, Lambda) {
  eta <- as.vector(x %*% beta + z %*% Lambda %*% u)
  mu <- exp(eta)
  cll.i <- dpois(y, mu, log = T)
  - sum(cll.i)
}


mll <- function(th, y, x, z, L0) {
  nx <- ncol(x); nz <- ncol(L0)
  nu <- sum(lower.tri(L0))
  q <- ncol(Z)/nz

  beta <- th[1:nx]
  Lu.d <- th[(nx + 1):(nx + nz)]
  if(nu > 0)
    Lu.o <- th[(nx + nz + 1):(nx + nz + nu)]

  diag(L0) <- Lu.d
  if(nu > 0) 
    L0[lower.tri(L0)] <- Lu.o
  Lambda <- lapply(1:q, function(i) L0)
  Lambda <- as.matrix(bdiag(Lambda))

  uhat <- pirls(y, z, x, beta, Lambda)

  L <- uhat$L
  uh <- uhat$par; utu <- as.vector(crossprod(uh))

  cll.i <- cll(uh, y, x, z, beta, Lambda)
  
  2 * cll.i + utu + log(det(uhat$L)^2)
}


glmmPo <- function(formula, data, start = NULL) {
  mc <- match.call()
  fr <- lme4:::lmerFrames(mc, formula, NULL)
  fl <- lme4:::lmerFactorList(formula, fr, 0L, 0L)

  y <- fr$Y; x <- fr$X; z <- t(as.matrix(fl$trms[[1]]$Zt))
  L0 <- fl$trms[[1]]$ST

  nx <- ncol(x); nl <- ncol(L0); nu <- sum(lower.tri(L0))

  if(is.null(start))
    ini <- c(rep(0, nx), rep(1, nl + nu))
  else
    ini <- c(start[1:nx], rep(1, nl + nu))

  inf <- c(rep(-Inf, nx), rep(1e-2, nl), rep(-Inf, nu))
  sup <- c(rep(Inf, nx + nl + nu))

  th.op <- nlminb(ini, mll, lower = inf, upper = sup, y = y, x = x, z = z, L0 = L0)

  return(th.op)
}



# Simulação de dados
period <- factor(1:4)
herd <- factor(1:10)

dat <- expand.grid(period = period, herd = herd)

beta<-c(-.3,1.7,2.5,3.4)                  # Vetor de efeitos fixos
X<-model.matrix(~period,dat)              # Matriz de delineamento de efeitos fixos

u<-rnorm(length(levels(dat$herd)),sd=.5)  # Vetor de efeitos aleatórios
Z<-model.matrix(~-1+herd,dat)             # Matriz de delineamento de ef. aleatórios

mu<-as.vector(exp(X%*%beta+Z%*%u))        # Vetor de efeitos médios

dat$resp<-rpois(nrow(dat),mu)  

# Ajuste do modelo
system.time(gm0 <- glmer(resp~period+(1|herd), data = dat, poisson))
system.time(gm1 <- glmmPo(resp~period+(1|herd), dat, fixef(gm0)))