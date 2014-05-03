period<-factor(1:4)
herd<-factor(1:10)

dat<-expand.grid(period=period,herd=herd)

beta<-c(-.3,1.7,2.5,3.4)
x<-model.matrix(~period,dat)

b<-rnorm(length(levels(dat$herd)),sd=.3)  
z<-model.matrix(~-1+herd,dat)             

mu<-as.vector(exp(x%*%beta+z%*%b))        

dat$resp<-rpois(nrow(dat),mu)  



gm1<-glmm.pois(resp~period,~-1+herd,dat)
lme4::glmer(resp~period+(1|herd),dat,poisson)



glmm.pois<-function(form.fix,rand,data){
  mfx<-model.frame(form.fix,data=data)
  x<-model.matrix(attr(mfx,"terms"),data=mfx)

  mfz<-model.frame(rand,data=data)
  z<-model.matrix(attr(mfz,"terms"),data=mfz)

  mr<-model.response(mfx)
  y<-as.vector(mr)


  ll<-function(th,y,x,z){
    beta<-th[1:ncol(x)]
    sig<-th[(ncol(x)+1)]

    z.ast<-sig*z

    l.u<-function(u,y,z,x,beta){
      eta.u<-exp(x%*%beta+z%*%u)
      f.y<-log(dpois(y,as.vector(eta.u)))
      l.yu<-sum(f.y)+sum(dnorm(u,log=T))
      return(-l.yu)
    }

    uhat<-optim(rep(0,ncol(z)),l.u,method="BFGS",hessian=T,
                y=y,z=z.ast,x=x,beta=beta)

    hess.u<-uhat$hessian
    uh<-uhat$par; utu<-as.vector(crossprod(uh))

    eta<-as.vector(x%*%beta+z%*%uh)
    mu<-exp(eta)
    f.y<-sum(dpois(y,mu,log=T))
    l.y<-f.y-.5*utu-.5*log(abs(det(hess.u)))

    return(-l.y)
  }

  ini<-c(rep(0,ncol(x)),1)
  inf<-c(rep(-Inf,ncol(x)),1e-2)
  sup<-c(rep(Inf,ncol(x)+1))

  th.op<-optim(ini,ll,method="L-BFGS-B",lower=inf,upper=sup,hessian=T,
               y=y,x=x,z=z,control=list(trace=T))

  return(th.op)
}


