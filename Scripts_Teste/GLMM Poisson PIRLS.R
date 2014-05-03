pirls<-function(y,Z,X,beta){
  new<-rep(0,ncol(Z))
  k<-0
  e<-100

  while(e>1e-6){
    k<-k+1
    old<-new

    eta<-as.vector(X%*%beta+Z%*%old)
    mu<-exp(eta)
    q<-eta+diag(1/mu)%*%(y-mu)
    W<-solve(diag(1/mu))
    I<-diag(rep(1,length(new)))
    L<-chol((t(Z)%*%W%*%Z)+I)

    new<-chol2inv(L)%*%t(Z)%*%W%*%q
    e<-max(abs(new-old))
  }

  list(par=new,L=L)
}




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

    uhat<-pirls(y=y,Z=z.ast,X=x,beta=beta)

    L<-uhat$L
    uh<-uhat$par; utu<-as.vector(crossprod(uh))

    eta<-as.vector(x%*%beta+z%*%uh)
    mu<-exp(eta)
    f.y<-sum(dpois(y,mu,log=T))
    l.y<-2*f.y-utu-log(det(L)^2)

    return(-l.y)
  }

  ini<-c(rep(0,ncol(x)),1)
  inf<-c(rep(-Inf,ncol(x)),1e-2)
  sup<-c(rep(Inf,ncol(x)+1))

  th.op<-optim(ini,ll,method="L-BFGS-B",lower=inf,upper=sup,hessian=T,
               y=y,x=x,z=z,control=list(trace=T))

  return(th.op)
}

