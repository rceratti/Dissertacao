library(pair.mglmm)
library(numDeriv)


# Simulação de dados
beta.c1<-c(0.70,1.45,1.65,1.90)
beta.c2<-c(0.96,1.39,0.40,1.19)
beta<-matrix(c(beta.c1,beta.c2),4,2)

phi<-1; p<-1.6 

mydat<-data.sim(2,'CP',exp,beta,xi=p,phi=phi)
dat<-mydat$Data

detach('package:pair.mglmm')


m0<-cpglmm(value~-1+variable+variable:period+(-1+variable|ID),data=dat)
m1<-cpglmm(value~-1+variable+(-1+variable|ID),data=dat)

system.time(jkl0<-llik.ep(m0))
rcov(jkl0)

system.time(jkl1<-llik.ep(m1))
rcov(jkl1)

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


# P/ log-lik final
cll<-function(u,r,beta,L,phi,p){
  eta<-r$x%*%beta+r$z%*%L%*%u
  mu<-exp(eta)
  
  p1<-log(tweedie::dtweedie.series(r$y,p,mu,phi))
  p2<-dnorm(u,0,sqrt(phi),log=T)
  sum(p1)+sum(p2)
}


mll2<-function(th,r,u){
  nx<-ncol(r$x); nz<-ncol(r$z)
  
  beta<-th[1:nx]
  
  L<-matrix(0,nz,nz); nLlt<-length(L[lower.tri(L,diag=T)])
  L[lower.tri(L,diag=T)]<-th[(nx+1):(nx+nLlt)]
  
  phi<-th[nx+nLlt+1]; p<-th[nx+nLlt+2]

  H.hat<-hessian(cll,u,r=r,beta=beta,L=L,phi=phi,p=p)
  cll(u,r,beta,L,phi,p)-.5*log(abs(det(H.hat)))
}


llik.ep<-function(mod){
  beta<-fixef(mod); phi<-mod@phi; p<-mod@p
  L<-t(.Call('mer_ST_chol',mod)[[1]])
  Llt<-L[lower.tri(L,diag=T)]
  th<-c(beta,Llt,phi,p)
  
  mf.ind<-split(mod@frame,mod@flist[[1]])
  mf.r<-lapply(mf.ind,f0,mod=mod)

  u<-matrix(mod@u,ncol=2,byrow=T)

  res<-lapply(1:length(mf.r),function(i){
    g<-grad(mll2,th,r=mf.r[[i]],u=c(u[i,]))
    h<-hessian(mll2,th,r=mf.r[[i]],u=c(u[i,]))
    list(g=g,h=h)
  })

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






#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Código com RcppArmadillo 
# Gradient vector function
colvec grad(colvec u, function f, List data){
  double y1, y2;
  colvec xi()
}


# Conditional log-likelihood
double cll(colvec u, mat x, mat z, colvec y, 
           colvec beta, mat L, double phi, double p) {

  colvec eta = x*beta + z*u;
  
  int n = eta.size();
  colvec mu;
  for(int i=0; i<n; i++) mu(i) = exp(eta(i));

  colvec p1=0, p2=0;
  int m = u.size();
  for(int i=0; i<n; i++) p1 += log(dcp(y(i), mu(i), phi, p));
  for(int j=0; j<m; j++) p2 += Rf_dnorm4(u(j), 0, sqrt(phi), 1);

  return wrap(p1+p2);
}


# Marginal log-likelihood
SEXP mll(SEXP x, SEXP z, SEXP y, SEXP u,
         SEXP beta, SEXP L, SEXP phi, SEXP p) {

  mat    x    = as<mat>(x),
         z    = as<mat>(z),
         L    = as<mat>(L);

  colvec y    = as<colvec>(y),
         beta = as<colvec>(beta),
         phi  = as<colvec>(phi),
         p    = as<colvec>(p);

  mat    Hhat = hess();
}

