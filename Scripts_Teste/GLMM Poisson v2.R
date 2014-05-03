library(lme4)

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


  ll<-function(th,y,x,z){
    beta<-th[1:ncol(x)]
    sig<-th[(ncol(x)+1)]


    l.u<-function(u,y,z,x,beta,sig){
      eta.u<-exp(x%*%beta+z%*%u)
      f.y<-dpois(y,as.vector(eta.u),log=T)
      l.yu<-sum(f.y)+sum(dnorm(u,sd=sig,log=T))
      return(l.yu)
    }

    uhat<-optim(rep(0,ncol(z)),l.u,method="BFGS",hessian=T,
                y=y,z=z,x=x,beta=beta,sig=sig,control=list(fnscale=-1))

    hess.u<-uhat$hessian

    l.y<-uhat$value-.5*log(abs(det(hess.u)))

    return(-l.y)
  }

  m0<-glm.fit(x,y,family=poisson())
    
  ini<-c(coef(m0),1)  # rep(0,ncol(x))
  inf<-c(rep(-Inf,ncol(x)),1e-2)
  sup<-c(rep(Inf,ncol(x)+1))

  b<-optim(ini,ll,method="L-BFGS-B",lower=inf,upper=sup,hessian=T,
           y=y,x=x,z=z,control=list(trace=T))
  
  coefficients<-b$par[1:ncol(x)]; vcov<-solve(b$hessian)
  colnames(vcov)<-rownames(vcov)<-c(attr(fr$fixef,'names'),'sigma.u')

  se<-sqrt(diag(vcov))[1:ncol(x)]
  tval<-coefficients/se

  coef<-cbind(Estimate=coefficients,StdErr=se,t.value=tval,
  	  p.value=2*pnorm(-abs(tval)))

  list(par=round(coef,4),sigma.u=round(b$par[ncol(x)+1],4))
}

