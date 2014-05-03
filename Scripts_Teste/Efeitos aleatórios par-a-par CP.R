library(pair.mglmm)
library(numDeriv)


# Simulação de dados
beta.c1<-c(0.70,1.45,1.65,1.90)
beta.c2<-c(0.96,1.39,0.40,1.19)
beta.c3<-c(1.25,1.86,0.19,-0.39)
beta<-matrix(c(beta.c1,beta.c2,beta.c3),4,3)

phi<-1; p<-1.6 

mydat<-data.sim(6,'CP',exp,xi=p,phi=phi)
dat<-mydat$Data

detach('package:pair.mglmm')


m0<-cpglmm(value~-1+variable+variable:period+(-1+variable|ID),data=dat)
m0.1<-cpglmm(value~-1+variable+variable:period+(-1+variable|ID),data=dat,doFit=F)


ranef(m0)
re.mglmm(m0.1,fixef(m0),VarCorr(m0)[[1]],m0@phi,m0@p)


##
f0<-function(mod,data){
  formula<-eval(mod@call$formula)

  fix.form<-lme4:::nobars(mod@formula)
  fr<-model.frame(fix.form,data=data)
  Fr<-list()
  Fr$Y<-model.response(fr); Fr$X<-model.matrix(attr(fr,"terms"),data=fr)
  Fr$wts<-mod@pWt; Fr$off<-mod@offset
  Fr$mf<-data; Fr$fixef<-mod@fixef

  fl<-lme4:::lmerFactorList(formula,Fr,0L,0L)

  x<-Fr$X; y<-Fr$Y; z<-t(as.matrix(fl$trms[[1]]$Zt))

  list(y=y,x=x,z=z)
}


## Versão 1
cll<-function(r,beta,u,S,phi,p){
  eta<-r$x%*%beta+r$z%*%u
  mu<-exp(eta)

  li.u<-mvtnorm::dmvnorm(u,sigma=S,log=T)
  lli<-sum(log(tweedie::dtweedie.series(r$y,p,mu,phi)))+li.u

  -lli
}


uhat<-function(r,beta,S,phi,p){
  q<-ncol(S)
  u.hat<-nlminb(rep(0,q),cll,r=r,beta=beta,S=S,phi=phi,p=p)
  u.hat$par
}


re.mglmm<-function(mod,beta,S,phi,p){
  mf.ind<-split(mod@frame,mod@flist[[1]])
  mf.r<-lapply(mf.ind,f0,mod=mod)
  
  re<-lapply(mf.r,function(r) uhat(r,beta,S,phi,p))
  do.call(rbind,re)
}


## Versão 2
cll<-function(r,beta,u,L,phi,p){
  eta<-r$x%*%beta+r$z%*%L%*%u
  mu<-exp(eta)

  li.u<-dnorm(u,sd=sqrt(phi),log=T)
  lli<-sum(log(tweedie::dtweedie.series(r$y,p,mu,phi)))+li.u

  -lli
}


uhat<-function(r,beta,L,phi,p){
  q<-ncol(L)
  u.hat<-nlminb(rep(0,q),cll,r=r,beta=beta,L=L,phi=phi,p=p)
  u.hat$par
}


re.mglmm<-function(mod,beta,S,phi,p){
  eS<-eigen(S/phi,symmetric=TRUE,EISPACK=TRUE)
  ev<-eS$values
  L<-eS$vectors%*%diag(sqrt(pmax(ev,0)))

  mf.ind<-split(mod@frame,mod@flist[[1]])
  mf.r<-lapply(mf.ind,f0,mod=mod)
  
  re<-lapply(mf.r,function(r) uhat(r,beta,L,phi,p))
  do.call(rbind,re)
}