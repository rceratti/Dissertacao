library(MASS)
library(pair.mglmm)
library(Rcpp)
library(inline)


# Simulação de dados
beta.c1<-c(0.70,1.45,1.65,1.90)
beta.c2<-c(0.96,1.39,0.40,1.19)
beta.c3<-c(1.25,1.86,0.19,-0.39)
beta<-matrix(c(beta.c1,beta.c2,beta.c3),4,3)

phi<-1; p<-1.6 

mydat<-data.sim(5,'CP',exp,xi=p,phi=phi)
dat<-mydat$Data

detach('package:pair.mglmm')


# Ajuste do modelo glmm CP
m0<-cpglmm(value~-1+variable+variable:period+(-1+variable|ID),data=dat)
m1<-cpglmm(value~-1+variable+(-1+variable|ID),data=dat)

# Modelos apenas com o esqueleto do modelo
m0.1<-cpglmm(value~-1+variable+variable:period+(-1+variable|ID),data=dat,doFit=F)
m1.1<-cpglmm(value~-1+variable+(-1+variable|ID),data=dat,doFit=F)

# Aproximaçao da verossimilhança via MC - puro R
system.time({
  (llm0<-llik.fim(m0.1,fixef(m0),VarCorr(m0)[[1]],m0@phi,m0@p,1e4))
})

system.time({
  (llm1<-llik.fim(m1.1,fixef(m1),VarCorr(m1)[[1]],m1@phi,m1@p,1e4))
})

# Aproximaçao da verossimilhança via MC - R e Rcpp
system.time({
  (llm0.1<-llik.fim2(m0.1,fixef(m0),VarCorr(m0)[[1]],m0@phi,m0@p,1e4))
})

system.time({
  (llm1.1<-llik.fim2(m1.1,fixef(m1),VarCorr(m1)[[1]],m1@phi,m1@p,1e4))
})


# Comparaçao dos resultados          
llm0-llm1
llm0.1-llm1.1
logLik(m0)-logLik(m1)



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


# P/ log-lik final - integração por Monte Carlo
cll<-function(u,r,beta,phi,p){
  eta<-r$x%*%beta+r$z%*%u
  mu<-exp(eta)
  prod(tweedie::dtweedie.series(r$y,p,mu,phi))
}

mll<-function(r,beta,u,phi,p){
  B<-ncol(u)
  l.i<-apply(u,1,cll,r=r,beta=beta,phi=phi,p=p)
  log(mean(l.i))
}

llik.fim<-function(mod,beta,S,phi,p,B){
  mf.ind<-split(mod@frame,mod@flist[[1]])
  mf.r<-lapply(mf.ind,f0,mod=mod)

  q<-ncol(S)
  U<-mvrnorm(B,rep(0,q),S)
  
  ll<-sapply(mf.r,function(r) mll(r,beta,U,phi,p))
  
  sum(ll)
}


## Versão utilizando Rcpp
## Parte 1: Rcpp
# Densidade da distribuição Poisson composta
inclCP<-"
double dcppos(double y, double mu, double phi, double p){
  double a, a1, logz, drop = 37, jmax, j, cc, wmax, estlogw;
  double wm = -1.0E16, sum_ww = 0, *ww, ld;
  int k, lo_j, hi_j;
    
  a = (2-p)/(1-p);
  a1 = 1 - a ;
  logz = -a*log(y)+a*log(p-1)-a1*log(phi)-log(2-p);
  jmax = R_pow(y,2-p)/(phi*(2-p));

  jmax = Rf_fmax2(1.0,jmax);
  j = jmax;
  cc = logz+a1+a*log(-a);
  wmax = a1*jmax;
  estlogw = wmax;

  while(estlogw > (wmax - drop)){
    j += 2.0;
    estlogw = j*(cc-a1*log(j)) ;
  }

  hi_j = ceil(j);
  j = jmax;
  estlogw = wmax;

  while((estlogw > (wmax - drop)) && (j >= 2)){
    j = Rf_fmax2(1,j-2);
    estlogw = j*(cc-a1*log(j));
  }

  lo_j = Rf_imax2(1,floor(j));
  ww = Calloc(hi_j-lo_j+1, double);

  for(k=lo_j; k<hi_j+1; k++){
    ww[k-lo_j] = k*logz-lgamma(1+k)-lgamma(-a*k);
    wm = Rf_fmax2(wm,ww[k-lo_j]);
  }

  for(k=lo_j; k<hi_j+1; k++)
    sum_ww += exp(ww[k-lo_j]-wm);

  ld = -y/(phi*(p-1)*R_pow(mu, p-1))-
       (R_pow(mu, 2-p)/(phi*(2-p)))-log(y)+
        log(sum_ww)+wm;

  Free(ww);
  return exp(ld);
}

double dcp(double y, double mu, double phi, double p){
  return (y > 0) ? dcppos(y, mu, phi, p) :
                   exp((-1)*R_pow(mu, 2-p)/(phi*(2-p))); 
}
"


codeCP<-"
  arma::colvec y   = as<arma::colvec>(Y);    // Vetor de respostas: n_i x 1
  arma::colvec b   = as<arma::colvec>(est);  // Efeitos fixos: p x 1
  arma::mat    x   = as<arma::mat>(X);       // Matriz de delineamento EFs: n_1 x p
  arma::mat    z   = as<arma::mat>(Z);       // Matriz de delineamento EAs: n_i x m
  arma::mat    u   = as<arma::mat>(U);       // Matriz de vetores MVN: B x m
  double       phi = as<double>(phi_r),
               p   = as<double>(p_r);              

  arma::mat    ZUt = z*u.t();
  arma::colvec Xb  = x*b;

  int B = u.n_rows, ni = x.n_rows, i, j;
  double mu, res = 0;
  arma::colvec res_i = arma::ones(B,1);

  for(i=0; i<B; i++){
    for(j=0; j<ni; j++){
      mu = exp(Xb(j)+ZUt(j,i));
      res_i(i) *= dcp(y(j),mu,phi,p);
    }
    res += res_i(i);
  }

  res = res/B;

  return wrap(res);
"

sigCP<-signature(Y="numeric",est="numeric",phi_r="numeric",p_r="numeric",
                 X="numeric",Z="numeric",U="numeric")

mllCP<-cxxfunction(sigCP, codeCP, 'RcppArmadillo', inclCP)


# Parte 2: R
llik.fim2<-function(mod,beta,S,phi,p,B){
  mf.ind<-split(mod@frame,mod@flist[[1]])
  mf.r<-lapply(mf.ind,f0,mod=mod)

  q<-ncol(S)
  U<-mvrnorm(B,rep(0,q),S)
  
  ll<-sapply(mf.r,function(r) mllCP(r$y,beta,phi,p,r$x,r$z,U))
  
  sum(log(ll))
}






#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

# P/ log-lik final
cll<-function(u,r,beta,S,phi,p){
  eS<-eigen(S,symmetric=TRUE,EISPACK=TRUE)
  ev<-eS$values
  L<-eS$vectors%*%diag(sqrt(pmax(ev,0)))
  L<-L/phi
  
  eta<-r$x%*%beta+r$z%*%L%*%u
  mu<-exp(eta)
  
  p1<-log(tweedie::dtweedie.series(r$y,p,mu,phi))
  p2<-dnorm(u,0,sqrt(phi),log=T)
  sum(p1)+sum(p2)
}

mll<-function(r,beta,S,phi,p){
  op<-optim(par=rep(0,ncol(r$z)),fn=cll,method="BFGS",hessian=T,control=list(fnscale=-1),
            r=r,beta=beta,S=S,phi=phi,p=p)
  
  op$value-.5*log(abs(det(op$hessian)))
}

llik.fim<-function(mod,beta,S,phi,p){
  mf.ind<-split(mod@frame,mod@flist[[1]])
  mf.r<-lapply(mf.ind,f0,mod=mod)
  
  ll<-sapply(mf.r,function(r) mll(r,beta,S,phi,p))
  
  sum(ll)
}

library(numDeriv)
ep<-hessian(llik.fim,fixef(m0),mod=m0.1,S=VarCorr(m0)[[1]],phi=m0@phi,p=m0@p)

llm0<-llik.fim(m0.1,fixef(m0),VarCorr(m0)[[1]],m0@phi,m0@p)
llm1<-llik.fim(m1.1,fixef(m1),VarCorr(m1)[[1]],m1@phi,m1@p)

llm0-llm1








cll<-function(u,r,beta,phi,p){
  eta<-r$x%*%beta+r$z%*%u
  mu<-exp(eta)
  
  p1<-log(tweedie::dtweedie.series(r$y,p,mu,phi))
  p2<-dnorm(u,0,sqrt(phi),log=T)
  
  sum(p1)+sum(p2)
}

gr.cll<-function(u,r,beta,phi,p){
  eta<-r$x%*%beta+r$z%*%u
  mu<-exp(as.vector(eta))
  
  q<-(r$y-mu)/mu
  W<-diag(mu^(2-p))
  
  gr<-((t(r$z)%*%W%*%q)-u)/phi
  as.vector(gr)
}

mll<-function(r,beta,phi,p){
  op<-optim(par=rep(0,ncol(r$z)),fn=cll,method="BFGS",hessian=T,control=list(fnscale=-1),
            r=r,beta=beta,phi=phi,p=p,gr=gr.cll)
  
  op$value-.5*log(abs(det(op$hessian)))
}

llik.fim<-function(mod,beta,S,phi,p){
  r<-list()
  r$y<-mod@y; r$x<-mod@X; r$z<-t(mod@Zt)
  
  lev.re<-unique(mod@flist[[1]])
  nlev<-length(lev.re)
  
  eS<-eigen(S,symmetric=TRUE,EISPACK=TRUE)
  ev<-eS$values
  Lamb0<-eS$vectors%*%diag(sqrt(pmax(ev,0)))
  Lamb0<-Lamb0/phi
  Lamb<-kronecker(diag(rep(1,nlev)),Lamb0)
  r$z<-as.matrix(r$z)%*%Lamb
  
  mll(r,beta,phi,p)
}