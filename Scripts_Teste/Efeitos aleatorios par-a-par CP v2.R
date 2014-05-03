library(pair.mglmm)

# Simulação de dados
phi <- 1; p <- 1.6 

mydat <- data.sim(3, 'CP', exp, xi=p, phi=phi)
dat <- mydat$Data


cl<-makeCluster(2)
registerDoParallel(cl)
clusterEvalQ(cl, library(pair.mglmm))

mainForm0 <- value~-1+variable:period+(-1+variable|ID)
mainForm1 <- value~-1+variable+(-1+variable|ID)

system.time(m0.1 <- mglmmCP(mainForm0, dat$variable, dat))
system.time(m1.1 <- mglmmCP(mainForm1, dat$variable, dat))

stopCluster(cl)


system.time(m0 <- cpglmm(mainForm0, data=dat))
system.time(m1 <- cpglmm(mainForm1, data=dat))


# Modas condicionais dos efeitos aleatórios
# Versao 1
cll<-function(r,beta,u,L,phi,p){
  eta<-r$x%*%beta+r$z%*%L%*%u
  mu<-exp(eta)
  
  li.u<-dnorm(u,sd=sqrt(phi),log=T)
  lli<-sum(log(tweedie::dtweedie.series(r$y,p,mu,phi)))+li.u
  
  -lli
}


uhat<-function(r,beta,L,phi,p,U){
  te <- apply(U, 1, cll, r=r, beta=beta, L=L, phi=phi, p=p)
  ini <- as.vector(U[which.min(te),])
  
  u.hat<-nlminb(ini,cll,r=r,beta=beta,L=L,phi=phi,p=p)
  u.hat$par
}


re.mglmm <- function(mod, formula, beta, S, phi, p){
  L <- t(chol(S))/sqrt(phi)
  
  mf.ind <- split(mod$fr$mf, mod$FL$fl[[1]])
  mf.r <- lapply(mf.ind, pair.mglmm:::f0.ep, formula = formula)
  
  q <- ncol(S)
  U <- MASS::mvrnorm(1e3, rep(0, q), diag(1, q))
  
  re <- lapply(mf.r, function(r) uhat(r, beta, L, phi, p, U))
  re <- do.call(rbind, re)
  
  
}


# Versao 2
cll<-function(r,beta,u,S,phi,p){
  eta<-r$x%*%beta+r$z%*%u
  mu<-exp(eta)
  
  li.u<-mvtnorm::dmvnorm(u,sigma=S,log=T)
  lli<-sum(log(tweedie::dtweedie.series(r$y,p,mu,phi)))+li.u
  
  -lli
}


uhat<-function(r,beta,S,phi,p){
  q <- ncol(S)
  ini <- rep(0, q)

  # Modas condicionais via Nelder-Mead
  u.hat <- optim(par = ini, fn = cll, r=r, beta=beta, S=S, phi=phi, p=p)
  u.hat$par
}


re.mglmm <- function(mod, formula, beta, S, phi, p){
  mf.ind <- split(mod$fr$mf, mod$FL$fl[[1]])
  mf.r <- lapply(mf.ind, pair.mglmm:::f0.ep, formula = formula)
  
  re <- lapply(mf.r, function(r) uhat(r, beta, S, phi, p))
  do.call(rbind, re)
}


# Versao 3
cll<-function(r,beta,u,S,phi,p){
  eta<-r$x%*%beta+r$z%*%u
  mu<-exp(eta)
  
  li.u<-mvtnorm::dmvnorm(u,sigma=S,log=T)
  lli<-sum(log(tweedie::dtweedie.series(r$y,p,mu,phi)))+li.u
  
  -lli
}


uhat<-function(r,beta,S,phi,p,U){
  te <- apply(U, 1, cll, r=r, beta=beta, S=S, phi=phi, p=p)
  ini <- as.vector(U[which.min(te),])
  
  u.hat<-nlminb(ini,cll,r=r,beta=beta,S=S,phi=phi,p=p)
  u.hat$par
}


re.mglmm <- function(mod, formula, beta, S, phi, p){
  mf.ind <- split(mod$fr$mf, mod$FL$fl[[1]])
  mf.r <- lapply(mf.ind, pair.mglmm:::f0.ep, formula = formula)
  
  q <- ncol(S)
  U <- MASS::mvrnorm(1e3, rep(0, q), S)
  
  re <- lapply(mf.r, function(r) uhat(r, beta, S, phi, p, U))
  do.call(rbind, re)
}


# Resultados e comparaçoes
m0.2 <- lmer(mainForm0, dat, doFit = F)
re0 <- re.mglmm(m0.2, mainForm0, as.vector(m0.1$fixef[,2]), 
                as.matrix(m0.1$VarCov), c(m0.1$phi[1])[[1]], 
                c(m0.1$p[1])[[1]])

re0
ranef(m0)



m1.2 <- lmer(mainForm1, dat, doFit = F)
re1 <- re.mglmm(m1.2, mainForm1, as.vector(m1.1$fixef[,2]), 
                as.matrix(m1.1$VarCov), c(m1.1$phi[1])[[1]], 
                c(m1.1$p[1])[[1]])

re1
ranef(m1)



#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

mll <- function(r, beta, u, S, phi, p){
  lli <- (-1)*cll(r, beta, u, S, phi, p)
  H.hat <- hessian(cll, u, r = r, beta = beta, S = S, phi = phi, p = p)

  lli - .5*log(abs(det(H.hat)))
}


llik.fim <- function(mod, formula, beta, u, S, phi, p){
  mf.ind <- split(mod$fr$mf, mod$FL$fl[[1]])
  mf.r <- lapply(mf.ind, pair.mglmm:::f0.ep, formula = formula)
  
  ll <- sapply(1:length(mf.r), function(i, r) {
          mll(r[[i]], beta, u[i,], S, phi, p)
        }, r = mf.r)
  
  sum(ll)
}

EP.hess <- hessian(llik.fim, as.vector(m0.1$fixef[,2]), mod = m0.2,
                   formula = mainForm0, u = re0, S = as.matrix(m0.1$VarCov),
                   phi = c(m0.1$phi[1])[[1]], p = c(m0.1$p[1])[[1]])



ll0 <- llik.fim(m0.2, mainForm0, as.vector(m0.1$fixef[,2]), re0,
         as.matrix(m0.1$VarCov), c(m0.1$phi[1])[[1]], c(m0.1$p[1])[[1]])

ll1 <- llik.fim(m1.2, mainForm1, as.vector(m1.1$fixef[,2]), re1,
         as.matrix(m1.1$VarCov), c(m1.1$phi[1])[[1]], c(m1.1$p[1])[[1]])


logLik(m0)-logLik(m1)
m0.1$logLik-m1.1$logLik
ll0-ll1
