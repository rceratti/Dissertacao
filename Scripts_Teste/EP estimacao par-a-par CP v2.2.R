library(pair.mglmm)
library(numDeriv)
library(SparseGrid)

# Simulação de dados
beta.c1<-c(0.70,1.45,1.65,1.90)
beta.c2<-c(0.96,1.39,0.40,1.19)
beta.c3<-c(1.25,1.86,0.19,-0.39)
beta<-matrix(c(beta.c1,beta.c2,beta.c3),4,3)

phi<-1; p<-1.6 

mydat<-data.sim(3,'CP',exp,beta,xi=p,phi=phi)
dat<-mydat$Data


##
f0<-function(mod,data,formula){
  formula0<-eval(formula)
  
  fix.form<-lme4:::nobars(formula)
  fr<-model.frame(fix.form,data=data)
  Fr<-list()
  Fr$Y<-model.response(fr); Fr$X<-model.matrix(attr(fr,"terms"),data=fr)
  Fr$wts<-mod@pWt; Fr$off<-mod@offset
  Fr$mf<-data; Fr$fixef<-mod@fixef
  
  fl<-lme4:::lmerFactorList(formula0,Fr,0L,0L)
  
  x<-Fr$X; y<-Fr$Y; z<-t(as.matrix(fl$trms[[1]]$Zt))
  
  list(y=y,x=x,z=z)
}


# P/ log-lik final
cll<-function(u,r,beta,L,phi,p){
  eta<-r$x%*%beta+r$z%*%L%*%u
  mu<-exp(eta)
  
  p1<-tweedie::dtweedie.series(r$y,p,mu,phi)
  p2<-dnorm(u,0,sqrt(phi))
  prod(p1)*prod(p2)*prod(exp(u^2/2))*sqrt(2*pi)^length(u)
}


mll2<-function(th,r){
  nx<-ncol(r$x); nz<-ncol(r$z)
  
  beta<-th[1:nx]
  
  L<-matrix(0,nz,nz); nLlt<-length(L[lower.tri(L,diag=T)])
  L[lower.tri(L,diag=T)]<-th[(nx+1):(nx+nLlt)]
  
  phi<-th[nx+nLlt+1]; p<-th[nx+nLlt+2]
  
  q<-ncol(L)
  pr.grid<-createProductRuleGrid('GQN',2,11,sym=FALSE)
  
  cll.pr<-apply(pr.grid$nodes,1,cll,r=r,beta=beta,L=L,phi=phi,p=p)
  res<-log(cll.pr%*%pr.grid$weights)
  
  as.vector(res)
}


llik.ep<-function(mod,formula){
  beta<-fixef(mod); phi<-mod@phi; p<-mod@p
  L<-t(.Call('mer_ST_chol',mod)[[1]])
  Llt<-L[lower.tri(L,diag=T)]
  th<-c(beta,Llt,phi,p)
  
  mf.ind<-split(mod@frame,mod@flist[[1]])
  mf.r<-lapply(mf.ind,f0,mod=mod,formula=formula)
  
  res<-foreach(i=1:length(mf.r)) %dopar% {
    g<-grad(mll2,th,r=mf.r[[i]])
    h<-hessian(mll2,th,r=mf.r[[i]])
    list(g=g,h=h)
  }
  
  return(res)
}


# Para lista de modelos par-a-par
rcov2<-function(lmods,formula){
  jkl<-lapply(lmods,llik.ep,formula=formula)
  l.jkl<-length(jkl[[1]])
  
  # Matriz K
  gs0<-lapply(1:l.jkl,function(i,x){
    p0<-lapply(x,'[[',i)
    p1<-lapply(p0,'[[','g')
    p2<-lapply(1:length(p1),function(i) lapply(p1,outer,p1[[i]]))
    p3<-lapply(p2,do.call,what=cbind)
    do.call(rbind,p3)
  },x=jkl)
  
  K<-(1/length(gs0))*Reduce('+',gs0)
  
  # Matriz J
  hs0<-lapply(jkl,function(x){
    p0<-lapply(x,'[[','h')
    (-1/length(x))*Reduce('+',p0)
  })
  
  J<-bdiag(hs0)
  Jinv<-solve(J)
  
  (Jinv%*%K%*%Jinv)/l.jkl
}


cl<-makeCluster(4)
registerDoParallel(cl)
clusterEvalQ(cl,{library(testpack);
                 library(SparseGrid);
                 library(numDeriv);
                 library(cplm)})
clusterExport(cl,c('cll','mll2'))

# Estimaçao par-a-par
system.time({
m0.1<-pair.mglmm:::glmmMultiCP(value~-1+variable+(-1+variable|ID),
                               dat$variable,dat,TRUE,cl)

JKJ.th<-rcov2(m0.1,formula=value~-1+variable+(-1+variable|ID))
})

pair.mglmm:::format1CP(m0.1,value~-1+variable+(-1+variable|ID),dat)$est
sqrt(diag(JKJ.th))[1:3]

stopCluster(cl)


# Modelo trivariado
system.time(m1<-cpglmm(value~-1+variable+(-1+variable|ID),data=dat))
m1

