library(cplm)
library(tweedie)


pirls <- function(y, Z, X, beta, Lambda, p) {
  ZL <- Z %*% Lambda
  new <- rep(0, ncol(ZL))
  err <- 100
  k <- 1

  while(err > 1e-8) {
    old <- new

    eta <- as.vector(X %*% beta + ZL %*% old)
    mu <- exp(eta)
    W <- diag(mu^(2-p)); I <- diag(1, ncol(ZL))
    H <- t(ZL) %*% W %*% ZL + I
    L <- chol(H)

    q <- (y-mu)/mu
    g <- (t(ZL) %*% W %*% q) - old

    new <- old + as.vector(chol2inv(L) %*% g)

    err <- max(abs(new - old))
    if(is.nan(err))
      stop('Problema na iteração ', k)
    k <- k + 1
  }

  list(par = new, L = L)
}


cll <- function(u, y, x, z, beta, Lambda, phi, p) {
  eta <- as.vector(x %*% beta + z %*% Lambda %*% u)
  mu <- exp(eta)
  cll.i <- log(dtweedie.series(y, p, mu, phi))
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
  phi <- exp(th[nx + nz + nu + 1])
  p <- th[nx + nz + nu + 2]

  diag(L0) <- Lu.d
  if(nu > 0) 
    L0[lower.tri(L0)] <- Lu.o
  Lambda <- lapply(1:q, function(i) L0)
  Lambda <- as.matrix(bdiag(Lambda))

  uhat <- pirls(y, z, x, beta, Lambda, p)

  L <- uhat$L
  uh <- uhat$par; utu <- as.vector(crossprod(uh))

  cll.i <- cll(uh, y, x, z, beta, Lambda, phi, p)
  
  2 * cll.i + utu/phi + log(det(uhat$L)^2)
}


glmmCP <- function(formula, data, start = NULL) {
  mc <- match.call()
  fr <- lme4:::lmerFrames(mc, formula, NULL)
  fl <- lme4:::lmerFactorList(formula, fr, 0L, 0L)

  y <- fr$Y; x <- fr$X; z <- t(as.matrix(fl$trms[[1]]$Zt))
  L0 <- fl$trms[[1]]$ST

  nx <- ncol(x); nl <- ncol(L0); nu <- sum(lower.tri(L0))

  if(is.null(start))
    ini <- c(rep(0, nx), rep(1, nl + nu), 1, 1.5)
  else
    ini <- c(start[1:nx], rep(1, nl + nu), start[nx + 1], start[nx + 2])

  inf <- c(rep(-Inf, nx), rep(1e-2, nl), rep(-Inf, nu + 1), 1.001)
  sup <- c(rep(Inf, nx + nl + nu + 1), 1.999)

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
phi<-1; p<-1.6                            # Parâmetros da dist. Poisson Composta

dat$resp<-rtweedie(nrow(dat),p,mu,phi)  

# Ajuste do modelo
system.time(gm0 <- cpglmm(resp~period+(1|herd), data = dat))
system.time(gm1 <- glmmCP(resp~period+(1|herd), dat, 
                          c(fixef(gm0), log(gm0@phi), gm0@p)))


gm0
gm1$par
(phi <- exp(gm1$par[6]))
sqrt(phi) * gm1$par[5]
gm1$obj


##
L1 <- diag(c(gm0@ST[[1]]), 10)
tes <- pirls(dat$resp, Z, X, fixef(gm0), L1, gm0@p)
tes$par
gm0@u

prm <- c(as.vector(fixef(gm0)), c(gm0@ST[[1]]), log(gm0@phi), gm0@p)
mll(prm, dat$resp, X, Z, gm0@ST[[1]])
deviance(gm0)

fitted(gm0)
(fit <- exp(as.vector(X %*% fixef(gm0) + Z %*% L1 %*% gm0@u)))


as.vector(gm0@sqrtXWt)
sqrt(fit^(2-gm0@p))


log(det(expand(gm0@L)$L)^2)
log(det(tes$L)^2)


#####
pirls <- function(y, Z, X, beta, Lambda, phi, p) {
  ZL <- Z %*% Lambda
  new <- rep(0, ncol(ZL))
  err <- 100
  k <- 1
  step <- 1

  while(err > 1e-6) {
    old <- new

    eta <- as.vector(X %*% beta + ZL %*% old)
    mu <- exp(eta)
    W <- diag(mu^(2-p)); I <- rep(1, ncol(ZL))
    H <- t(ZL) %*% W %*% ZL + I
    L <- chol(H)

    q <- (y-mu)/mu
    g <- (t(ZL) %*% W %*% q) - old

    step <- ot.uni(old, as.vector(g), y, X, Z, beta, Lambda, phi, p)
    new <- old + step * as.vector(chol2inv(L) %*% g)
    cat(step, '\n')
    
    err <- max(abs(new - old))
    if(is.nan(err))
      stop('Problema na iteração ', k)
    k <- k + 1
  }

  list(par = new, L = L)
}


ot.uni <- function(u, g, y, x, z, beta, Lambda, phi, p) {
  a <- 1e-3; b <- 1
  
  while(b - a > 1e-3){
    z.a <- b - .618 * (b - a)
    z.b <- a + .618 * (b - a)
    f.a <- cll(u + z.a * g, y, x, z, beta, Lambda, phi, p)
    f.b <- cll(u + z.b * g, y, x, z, beta, Lambda, phi, p)
    if(f.a < f.b) b <- z.b else a <- z.a
  }

  (a+b)/2
}