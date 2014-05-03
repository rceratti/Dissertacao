library(pair.mglmm)
library(numDeriv)
library(bayespack)

# Simulação de dados
beta.c1<-c(0.70,1.45,1.65,1.90)
beta.c2<-c(0.96,1.39,0.40,1.19)
beta.c3<-c(1.25,1.86,0.19,-0.39)
beta<-matrix(c(beta.c1,beta.c2,beta.c3),4,3)

phi<-1; p<-1.6 

mydat<-data.sim(2,'CP',exp,xi=p,phi=phi)
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


mll2<-function(th,r){
  nx<-ncol(r$x); nz<-ncol(r$z)
  
  beta<-th[1:nx]
  
  L<-matrix(0,nz,nz); nLlt<-length(L[lower.tri(L,diag=T)])
  L[lower.tri(L,diag=T)]<-th[(nx+1):(nx+nLlt)]
  
  phi<-th[nx+nLlt+1]; p<-th[nx+nLlt+2]
 
  q<-ncol(L)
  banint(cll,start=rep(0,q),r=r,beta=beta,L=L,phi=phi,p=p,trans='none')$nrmcon[1]
}


llik.ep<-function(mod){
  beta<-fixef(mod); phi<-mod@phi; p<-mod@p
  L<-t(.Call('mer_ST_chol',mod)[[1]])
  Llt<-L[lower.tri(L,diag=T)]
  th<-c(beta,Llt,phi,p)
  
  mf.ind<-split(mod@frame,mod@flist[[1]])
  mf.r<-lapply(mf.ind,f0,mod=mod)

  res<-lapply(mf.r,function(r){
    g<-grad(mll2,th,r=r)
    h<-hessian(mll2,th,r=r)
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


# Para lista de modelos par-a-par
rcov2<-function(lmods){
  jkl<-lapply(lmods,llik.ep)
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


## O código acima *parece* funcionar bem, mas é muito, MUITO demorado.
## Segue abaixo uma tentativa de otimizar as partes do código que demoram mais.
## Especificamente, as chamadas à função mll2() que, por sua vez, chama uma 
## integração numérica via função banint(). 

# Simplificando a função banint(), ajustando-a aos meus interesses
function (logPost, mnFns = function(x, ...) NULL, mode = NULL, 
    start = NULL, method = "Gauss-Hermite", trans = "split-t", 
    control = list(), optimctrl = list(), optimMethod = "BFGS", 
    ...) 
{
    if (is.null(start) & is.null(mode)) 
        stop("either mode or start argument needed")
    if (is.null(mode)) {
        fn <- function(x) -logPost(x, ...)
        optobj <- optim(start, fn, method = optimMethod, control = optimctrl)
        if (optobj$convergence > 0) 
            stop("failure in optimization")
        mode <- optobj$par
    }

    con <- list(relreq = 0.01, maxvls = 10000, numtrn = 10, df = 4, 
                print = FALSE, method = "Gauss-Hermite", mthd = 30)

  # Continuar com o restante dos argumentos...
}


# Código com RcppArmadillo 
double cll(vec u, mat x, mat z, vec y, 
           vec beta, mat L, double phi, double p) {

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


SEXP mll(SEXP x, SEXP z, SEXP y,
         SEXP beta, SEXP L, SEXP phi, SEXP p) {

  mat    x    = as<mat>(x),
         z    = as<mat>(z),
         L    = as<mat>(L);

  colvec y    = as<colvec>(y),
         beta = as<colvec>(beta),
         phi  = as<colvec>(phi),
         p    = as<colvec>(p);

  // resolver como incorporar método de integração
  // do pacote banint -- se é que é possível
}
