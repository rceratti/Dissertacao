library(lme4)
library(numDeriv)

period<-factor(1:4)
herd<-factor(1:10)

dat<-expand.grid(period=period,herd=herd)

beta<-c(-.3,1.7,2.5,3.4)
x<-model.matrix(~period,dat)

b<-rnorm(length(levels(dat$herd)),sd=.3)  
z<-model.matrix(~-1+herd,dat)             

mu<-as.vector(exp(x%*%beta+z%*%b))        

dat$resp<-rpois(nrow(dat),mu)  



glmm.pois(resp~period+(1|herd),dat)
glmer(resp~period+(1|herd),dat,poisson)



glmm.pois<-function(formula,data){
  contrasts<-NULL
  mc<-match.call()
  fr<-lme4:::lmerFrames(mc,formula,contrasts)
  fl<-lme4:::lmerFactorList(formula,fr,0L,0L)

  y<-fr$Y; x<-fr$X; z<-t(as.matrix(fl$trms[[1]]$Zt))


  ll<-function(th,y,X,z){
    beta<-th[1:ncol(x)]
    sig<-th[(ncol(x)+1)]

    l.u<-function(u,y,z,X,beta,sig){
      eta.u<-exp(X%*%beta+z%*%u)
      f.y<-dpois(y,as.vector(eta.u),log=T)
      l.yu<-sum(f.y)+sum(dnorm(u,sd=sig,log=T))
      return(-l.yu)
    }

    uhat<-nlminb(rep(0,ncol(z)),l.u,y=y,z=z,X=x,beta=beta,sig=sig)
    hess.u<-hessian(l.u,uhat$par,y=y,z=z,X=x,beta=beta,sig=sig)

    l.y<--uhat$objective-.5*log(abs(det(hess.u)))

    return(-l.y)
  }

  ini<-c(rep(0,ncol(x)),1)
  inf<-c(rep(-Inf,ncol(x)),1e-2)
  sup<-c(rep(Inf,ncol(x)+1))

  b<-nlminb(ini,ll,lower=inf,upper=sup,y=y,X=x,z=z)
  hess.b<-hessian(ll,b$par,y=y,X=x,z=z)
  
  coefficients<-b$par[1:ncol(x)]; vcov<-solve(hess.b)
  colnames(vcov)<-rownames(vcov)<-c(attr(fr$fixef,'names'),'sigma.u')

  se<-sqrt(diag(vcov))[1:ncol(x)]
  tval<-coefficients/se
  pval<-2*pnorm(-abs(tval))
  pval<-ifelse(pval<1e-4,paste("<1e-4"),round(pval,4))

  coef<-data.frame(Estimate=coefficients,StdErr=se,t.value=tval)
  coef<-round(coef,4)
  coef<-data.frame(coef,p.value=pval)

  list(par=coef,sigma.u=round(b$par[ncol(x)+1],4))
}