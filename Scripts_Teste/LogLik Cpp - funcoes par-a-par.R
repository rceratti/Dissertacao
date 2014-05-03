


# Integraçao monte carlo da funçao de verossimilhança p/ mglmm poisson
codePois<-"
arma::colvec y=Rcpp::as<arma::colvec>(Y);      // Vetor de respostas: n_i x 1
arma::colvec b=Rcpp::as<arma::colvec>(est);    // Efeitos fixos: p x 1
arma::mat    x=Rcpp::as<arma::mat>(X);         // Matriz de delineamento EFs: n_1 x p
arma::mat    z=Rcpp::as<arma::mat>(Z);         // Matriz de delineamento EAs: n_i x m
arma::mat    u=Rcpp::as<arma::mat>(U);         // Matriz de vetores MVN: B x m

arma::mat    ZUt=z*u.t();
arma::colvec Xb=x*b;

int B=u.n_rows, ni=x.n_rows, i;
double mu, res=0;
arma::colvec res_i=arma::ones(B,1);

for(i=0; i<B; i++){
  for(int j=0; j<ni; j++){
    mu=R_pow(M_E,Xb(j)+ZUt(j,i));
    res_i(i)*=Rf_dpois(y(j),mu,0);
  }

  res+=res_i(i);
}

res=res/B;

return wrap(res);
"

codeGauss<-"
arma::colvec y=Rcpp::as<arma::colvec>(Y);      // Vetor de respostas: n_i x 1
arma::colvec b=Rcpp::as<arma::colvec>(est);    // Efeitos fixos: p x 1
arma::mat    x=Rcpp::as<arma::mat>(X);         // Matriz de delineamento EFs: n_1 x p
arma::mat    z=Rcpp::as<arma::mat>(Z);         // Matriz de delineamento EAs: n_i x m
arma::mat    u=Rcpp::as<arma::mat>(U);         // Matriz de vetores MVN: B x m

double SD=Rcpp::as<double>(sd);

arma::mat    ZUt=z*u.t();
arma::colvec Xb=x*b;

int B=u.n_rows, ni=x.n_rows, i;
double mu, res=0;
arma::colvec res_i=arma::ones(B,1);

for(i=0; i<B; i++){
  for(int j=0; j<ni; j++){
    mu=Xb(j)+ZUt(j,i);
    res_i(i)*=Rf_dnorm4(y(j),mu,SD,0);
  }

  res+=res_i(i);
}

res=res/B;

return wrap(res);
"

codeBern<-"
arma::colvec y=Rcpp::as<arma::colvec>(Y);      // Vetor de respostas: n_i x 1
arma::colvec b=Rcpp::as<arma::colvec>(est);    // Efeitos fixos: p x 1
arma::mat    x=Rcpp::as<arma::mat>(X);         // Matriz de delineamento EFs: n_1 x p
arma::mat    z=Rcpp::as<arma::mat>(Z);         // Matriz de delineamento EAs: n_i x m
arma::mat    u=Rcpp::as<arma::mat>(U);         // Matriz de vetores MVN: B x m

arma::mat    ZUt=z*u.t();
arma::colvec Xb=x*b;

int B=u.n_rows, ni=x.n_rows, i;
double mu, res=0;
arma::colvec res_i=arma::ones(B,1);

for(i=0; i<B; i++){
  for(int j=0; j<ni; j++){
    mu=1/(1+R_pow(M_E,-Xb(j)-ZUt(j,i)));
    res_i(i)*=Rf_dbinom(y(j),1,mu,0);
  }

  res+=res_i(i);
}

res=res/B;

return wrap(res);
"


sigPois<-signature(Y="numeric",est="numeric",X="numeric",Z="numeric",U="numeric")
sigGaus<-signature(Y="numeric",est="numeric",X="numeric",Z="numeric",U="numeric",sd="numeric")
sigBern<-signature(Y="numeric",est="numeric",X="numeric",Z="numeric",U="numeric")
funcs<-cxxfunction(list(mc_pois=sigPois,mc_gauss=sigGaus,mc_bern=sigBern),
                   list(codePois,codeGauss,codeBern),plugin='RcppArmadillo')


# Rebuild the log-likelihood
llik.glmm2<-function(mod,beta,S,phi=NULL,par=FALSE,cl=NULL,method){
  # List with elements to compute the likelihood for each subject
  f0<-function(mod,data){
    formula<-eval(mod$formula)
    
    mc<-match.call()
    fr<-lme4:::lmerFrames(mc,formula,NULL)
    fl<-lme4:::lmerFactorList(formula,fr,0L,0L)
    
    y<-fr$Y; x<-fr$X; z<-t(as.matrix(fl$trms[[1]]$Zt))
    
    ffun<-if(is.null(mod$glmFit)) 
      gaussian() else mod$glmFit$family
    
    fam<-ffun$family
    n<-0
    if(fam=='binomial') n<-mod$glmFit$prior.weights
    
    list(y=y,x=x,z=z,n=n,fam=fam,ffun=ffun)
  }
  
  # Log-likelihood function -- method='GH' via banint()
  f1<-function(r,beta,u,S,phi){
    if(class(u)=='data.frame') u<-as.matrix(u)
    if(is.vector(u)) u<-matrix(u,nc=ncol(S))
    u.vec<-matrix(u,nc=1)
    
    eta<-r$x%*%beta+r$z%*%u.vec
    mu<-r$ffun$linkinv(eta)
    
    l.u<-sum(dmvnorm(u,sigma=S,log=T))
    
    res<-switch(r$fam,
                poisson=sum(dpois(r$y,mu,log=T))+l.u,
                gaussian=sum(dnorm(r$y,mu,phi,log=T))+l.u,
                binomial=sum(dbinom(r$y,r$n,mu,log=T))+l.u
    )
    return(res)
  }
  
  # Log-likelihood function -- method='MC'
  mint<-function(method,...,dist=NULL){
    if(method=='GH') banint(f1,...,trans='none')$nrmcon[1]
    
    if(method=='MC'){
      l<-switch(dist,
                poisson=.Call('mc_pois',...,PACKAGE="pair.mglmm"),
                gaussian=.Call('mc_gauss'...,PACKAGE="pair.mglmm"),
                binomial=.Call('mc_bern',...,PACKAGE="pair.mglmm")
      )
      log(l)
    }
  }
  
  if(method=='MC') U<-mvrnorm(1e5,rep(0,ncol(S)),S)
  
  
  # Splitting the original data frame by subject
  mf.ind<-split(mod$fr$mf,mod$FL$fl)
  nr<-length(mf.ind)
  
  # subject-wise likelihood evaluation
  if(par && getDoParWorkers()>1){
    b<-foreach(i=1:nr,.combine=c) %dopar% {
      r<-f0(mod,mf.ind[[i]])
      b.i<-switch(method, 
                  GH=mint('GH',start=rep(0,ncol(r$z)),r=r,beta=beta,S=S,phi=phi),
                  MC=mint('MC',r$y,beta,r$x,r$z,U,dist=r$fam)
           )
    }
  }
  
  if(!par){
    b<-foreach(i=1:nr,.combine=c) %do% {
      r<-f0(mod,mf.ind[[i]])
      b.i<-switch(method, 
                  GH=mint('GH',start=rep(0,ncol(r$z)),r=r,beta=beta,S=S,phi=phi),
                  MC=mint('MC',r$y,beta,r$x,r$z,U,dist=r$fam)
      )
    }
  }
  
  sum(b)
}



# Combining the above functions inside the main MGLMM function (to be called)
mglmm2<-function(formula,id,family,data,llik=TRUE,par=FALSE,cl=NULL,method='GH'){
  # Pair-wise fitting
  m0<-glmmMulti(formula,id,family,data,par,cl)
  
  # Averaging the estimates from the pairwise models
  prm<-format1(m0,formula,data,family)
  
  # Log-likelihood
  if(llik)
    LL<-llik.glmm2(prm$m1,prm$est,prm$SigHat,prm$phi,par,cl,method)
  else
    LL<-NA
  
  phi<-ifelse(prm$phi>0,prm$phi,NA)
  
  # Output - Fixed effects and var-cov estimates and log-likelihood
  list(Estimates=list(fixef=prm$est,VarCov=prm$SigHat,phi=phi),logLik=LL)
}


setwd('~/Pacotes R/pair.mglmm_1.0-1')
RcppArmadillo.package.skeleton("pair.mglmm",'mint',example_code=FALSE)

setwd('~/Pacotes R/pair.mglmm_1.0-1/pair.mglmm/src')
cat(funcs[[1]]@code,file='pair.glmm.cpp')
