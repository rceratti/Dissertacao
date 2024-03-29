library(tweedie)
library(cplm)


glmm.cp<-function(form.fix,rand,data,eps=1/6){
  mfx<-model.frame(form.fix,data=data)
  x<-model.matrix(attr(mfx,"terms"),data=mfx)

  mfz<-model.frame(rand,data=data)
  z<-model.matrix(attr(mfz,"terms"),data=mfz)

  mr<-model.response(mfx)
  y<-as.vector(mr)


  ll<-function(th,y,x,z){
    n<-ncol(x)

    beta<-th[1:n]
    sig<-th[n+1]
    phi<-th[n+2]
    p<-th[n+3]

    l.u<-function(u,y,x,z,beta,phi,p){
      eta.u<-as.vector(x%*%beta+z%*%u)
      f.y<-log(dtweedie.saddle(y,p,exp(eta.u),phi,eps=eps))
      l.yu<-sum(f.y)+sum(dnorm(u,sd=sig,log=T))
      return(l.yu)
    }

    uhat<-optim(par=rep(0,ncol(z)),fn=l.u,method="BFGS",hessian=T,
                y=y,x=x,z=z,beta=beta,phi=phi,p=p,control=list(fnscale=-1))

    hess.u<-uhat$hessian
    uh<-uhat$par; utu<-as.vector(crossprod(uh))

    f.y<-l.u(uh,y,x,z,beta,phi,p)
    l.y<-f.y-.5*log(abs(det(hess.u)))

    return(-l.y)
  }

  m0<-cpglm(form.fix,data=data)

  ini<-c(coef(m0),1,m0$phi,m0$p)
  inf<-c(rep(-Inf,ncol(x)),1e-2,1e-2,1.0001)
  sup<-c(rep(Inf,ncol(x)+2),1.9999)

  b<-optim(ini,ll,method="L-BFGS-B",lower=inf,upper=sup,hessian=T,
           y=y,x=x,z=z,control=list(trace=T))

  coefficients<-b$par[1:ncol(x)]; vcov<-solve(b$hessian)

  se<-sqrt(diag(vcov))[1:ncol(x)]
  tval<-coefficients/se

  coef<-cbind(Estimate=coefficients,StdErr=se,t.value=tval,
      p.value=2*pnorm(-abs(tval)))

  list(par=round(coef,4),sigma.u=round(b$par[ncol(x)+1],4),
       phi=round(b$par[ncol(x)+2],4),p=round(b$par[ncol(x)+3],4))
}




# Simula��o de dados
period<-factor(1:4)
herd<-factor(1:10)

dat<-expand.grid(period=period,herd=herd)

beta<-c(-.3,1.7,2.5,3.4)                  # Vetor de efeitos fixos
X<-model.matrix(~period,dat)              # Matriz de delineamento de efeitos fixos

u<-rnorm(length(levels(dat$herd)),sd=.3)  # Vetor de efeitos aleat�rios
Z<-model.matrix(~-1+herd,dat)             # Matriz de delineamento de ef. aleat�rios

mu<-as.vector(exp(X%*%beta+Z%*%u))        # Vetor de efeitos m�dios
phi<-1; p<-1.6                            # Par�metros da dist. Poisson Composta

dat$resp<-rtweedie(nrow(dat),p,mu,phi)  

# Ajuste do modelo
system.time(gm0<-cpglmm(resp~period+(1|herd),data=dat))
system.time(gm1<-glmm.cp(form.fix=resp~period,~-1+herd,dat))





##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
l.yu<-function(y,beta,u,sig) dpois(y,exp(beta+sig*u),log=T)-.5*u^2 #dnorm(u,log=T)
f.yu<-function(...) exp(l.yu(...))

I<-integrate(f.yu,-Inf,Inf,y=1,beta=-1,sig=.3)
log(I$value)

uh<-optim(0,l.yu,method="BFGS",hessian=T,y=1,beta=-1,sig=.3,control=list(fnscale=-1))
.5*log(2*pi)+l.yu(1,-1,uh$par,.3)-.5*log(-(uh$hessian))
.5*log(2*pi)+l.yu(1,-1,uh$par,.3)+.5*log(-solve(uh$hessian))