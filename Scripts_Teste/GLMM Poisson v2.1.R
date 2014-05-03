library(lme4)
library(MASS)
library(reshape)
library(mvtnorm)

f.S<-function(n){
  S<-matrix(0,n,n)
  diag(S)<-rep(1,n)

  corr<-runif(n*(n-1)/2,-.9,.9)
  S[lower.tri(S)]<-S[upper.tri(S)]<-corr

  nearPD(S,corr=T)$mat
}


period<-factor(1:2)
herd<-factor(1:10)

dat<-expand.grid(period=period,herd=herd)

beta<-c(-.3,1.7)
X<-model.matrix(~period,dat)

S<-f.S(2)
b<-mvrnorm(length(levels(dat$herd)),c(0,0),S)
b<-as.vector(b)
Z<-model.matrix(~-1+herd:period,dat)           


eta<-X%*%beta+Z%*%b
mu<-exp(eta)

dat$resp<-rpois(nrow(dat),mu)



m0<-glmm.pois(resp~period+(0+period|herd),dat)
m1<-glmer(resp~period+(0+period|herd),dat,poisson)



glmm.pois<-function(formula,data){
  contrasts<-NULL
  mc<-match.call()
  fr<-lme4:::lmerFrames(mc,formula,contrasts)
  fl<-lme4:::lmerFactorList(formula,fr,0L,0L)

  y<-fr$Y; x<-fr$X; z<-t(as.matrix(fl$trms[[1]]$Zt))
  S<-fl$trms[[1]]$ST

  ll<-function(th,y,x,z,S){
    nx<-ncol(x); nz<-ncol(S)
    nu<-sum(lower.tri(S))

    beta<-th[1:nx]
    Su.d<-th[(nx+1):(nx+nz)]
    Su.o<-th[(nx+nz+1):(nx+nz+nu)]

    diag(S)<-Su.d
    S[lower.tri(S)]<-S[upper.tri(S)]<-Su.o

    l.u<-function(u,y,z,x,beta,S){
      eta.u<-exp(x%*%beta+z%*%u)
      u<-matrix(u,nc=ncol(S))
      f.y<-dpois(y,as.vector(eta.u),log=T)
      l.yu<-sum(f.y)+sum(dmvnorm(u,sigma=S,log=T))
      return(l.yu)
    }

    uhat<-optim(rep(0,ncol(z)),l.u,method="BFGS",hessian=T,
                y=y,z=z,x=x,beta=beta,S=S,control=list(fnscale=-1))

    hess.u<-uhat$hessian

    l.y<-uhat$value-.5*log(abs(det(hess.u)))

    return(-l.y)
  }

  ini.S<-c(rep(1,ncol(S)),c(S[lower.tri(S)]))+1e-2
  nx<-ncol(x); ns<-ncol(S); nu<-sum(lower.tri(S))

  ini<-c(rep(0,nx),ini.S)
  inf<-c(rep(-Inf,nx),rep(1e-3,ns),rep(-Inf,nu))
  sup<-c(rep(Inf,nx+ns+nu))

  b<-optim(ini,ll,method="L-BFGS-B",lower=inf,upper=sup,hessian=T,
           y=y,x=x,z=z,S=S,control=list(trace=T))
  
  coefficients<-b$par[1:nx]; vcov<-solve(b$hessian[1:nx,1:nx])
  colnames(vcov)<-rownames(vcov)<-c(attr(fr$fixef,'names'))

  se<-sqrt(diag(vcov))[1:nx]
  tval<-coefficients/se

  coef<-cbind(Estimate=coefficients,StdErr=se,t.value=tval,
  	        p.value=2*pnorm(-abs(tval)))

  list(par=round(coef,4),sigma.u=round(b$par[(nx+1):(nx+ns+nu)],4))
}
