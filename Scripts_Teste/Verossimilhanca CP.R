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
m1<-cpglmm(value~-1+variable+(-1+variable|ID),data=dat)

m0.1<-cpglmm(value~-1+variable+variable:period+(-1+variable|ID),data=dat,doFit=F)
m1.1<-cpglmm(value~-1+variable+(-1+variable|ID),data=dat,doFit=F)

llm0<-llik.fim(m0.1,fixef(m0),VarCorr(m0)[[1]],m0@phi,m0@p)
llm1<-llik.fim(m1.1,fixef(m1),VarCorr(m1)[[1]],m1@phi,m1@p)

llm0-llm1
logLik(m0)-logLik(m1)


jkl0<-llik.ep(m0)
jkl1<-llik.ep(m1)


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


cll<-function(r,beta,u,S,phi,p){
  eta<-r$x%*%beta+r$z%*%u
  mu<-exp(eta)

  li.u<-mvtnorm::dmvnorm(u,sigma=S,log=T)
  lli<-sum(log(tweedie::dtweedie.series(r$y,p,mu,phi)))+li.u

  -lli
}


# P/ log-lik final
mll<-function(r,beta,S,phi,p){
  q<-ncol(S)
  u.hat<-nlminb(rep(0,q),cll,r=r,beta=beta,S=S,phi=phi,p=p)
  H.hat<-hessian(cll,u.hat$par,r=r,beta=beta,S=S,phi=phi,p=p)

  -u.hat$objective-.5*log(abs(det(H.hat)))
}


llik.fim<-function(mod,beta,S,phi,p){
  mf.ind<-split(mod@frame,mod@flist[[1]])
  mf.r<-lapply(mf.ind,f0,mod=mod)
  
  ll<-sapply(mf.r,function(r) mll(r,beta,S,phi,p))
  
  sum(ll)
}


# P/ EP:
mll2<-function(th,r,u){
  nx<-ncol(r$x); nz<-ncol(r$z)
  
  beta<-th[1:nx]
  
  S<-matrix(0,nz,nz); nSlt<-length(S[lower.tri(S,diag=T)])
  S[lower.tri(S,diag=T)]<-th[(nx+1):(nx+nSlt)]
  S<-S+t(lower.tri(S)*S)
  
  phi<-th[nx+nSlt+1]; p<-th[nx+nSlt+2]
  
  q<-ncol(S)
  obj<-cll(r,beta,u,S,phi,p)
  H.hat<-hessian(cll,u,r=r,beta=beta,S=S,phi=phi,p=p)
  
  -obj-.5*log(abs(det(H.hat)))
}


llik.ep<-function(mod){
  beta<-fixef(mod); phi<-mod@phi; p<-mod@p
  L<-t(.Call('mer_ST_chol',mod)[[1]])
  S<-phi*(L%*%t(L)); Slt<-S[lower.tri(S,diag=T)]
  th<-c(beta,Slt,phi,p)
  
  u<-ranef(mod)[[1]]
  
  mf.ind<-split(mod@frame,mod@flist[[1]])
  mf.r<-lapply(mf.ind,f0,mod=mod)

  res<-lapply(1:length(mf.r),function(i,r){
    g<-grad(mll2,th,r=r[[i]],u=unlist(u[i,]))
    h<-hessian(mll2,th,r=r[[i]],u=unlist(u[i,]))
    list(g=g,h=h)
  },r=mf.r)

  return(res)
}





rcov<-function(jkl){
  l.jkl<-length(jkl)
  
  # Matriz K
  gs0<-lapply(jkl,function(x){
    gi<-x[['g']]
    outer(gi,gi)
  })
  K<-(1/length(gs0))*Reduce('+',gs0)
  
  # Matriz J
  hs0<-lapply(jkl,'[[','h')
  J<-(1/length(hs0))*Reduce('+',hs0)
  Jinv<-solve(J)
  
  JKJ<-(Jinv%*%K%*%Jinv)/l.jkl
  sqrt(diag(JKJ))
}

rcov(jkl0)
rcov(jkl1)




tes<-hessian(llik.fim,fixef(m1.1),mod=m1,S=VarCorr(m1)[[1]],phi=m1@phi,p=m1@p)

