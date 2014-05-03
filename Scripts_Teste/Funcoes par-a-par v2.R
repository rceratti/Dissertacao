library(cplm)
library(MASS)
library(reshape)
library(mvtnorm)
library(bayespack)
library(doParallel)


# Function that fits all the bivariate models 
glmmMulti<-function(formula,id,family,data,par=FALSE,cl=NULL){
  # All pairwise combinations of response variables
  lev<-paste(unique(id))
  tes<-combn(lev,m=2)

  glmm.fit<-function(x,data,formula,family){
    ind<-id %in% x
    dat<-subset(data,ind)

    glmer(formula,family=family,data=dat)
  }

  # Parallel version
  if(par && getDoParWorkers()>1){
    #if(!is.null(cl)) clusterEvalQ(cl,library(lme4))
    res<-foreach(i=1:ncol(tes)) %dopar% glmm.fit(tes[,i],data,formula,family)
  }

  # Serial version
  if(!par) res<-foreach(i=1:ncol(tes)) %do% glmm.fit(tes[,i],data,formula,family)

  return(res)
}



# Re-arranging the estimates (fixed effects, residual std. deviation and var-cov)
# from an object fitted with glmer() (class 'mer')
format0<-function(mod){
  betas<-fixef(mod)
  df.b<-melt(betas)
  df.b$Parametro<-rownames(df.b)
  rownames(df.b)<-1:nrow(df.b)

  df.b<-df.b[,c(2,1)]

  S<-VarCorr(mod)
  uS<-unlist(S)
  triS<-unique(uS)

  nameS<-unlist(attr(S[[1]],"dimnames"))
  nameS<-combn(sort(nameS),2)
  nameS<-unique(t(nameS))
  nameS<-apply(nameS,1,paste,collapse=":")

  df.s<-data.frame(Parametro=nameS,value=triS)

  phi<-as.numeric(attr(S,'sc'))
  phi<-data.frame(Parametro='phi',value=ifelse(is.na(phi),0,phi))

  rbind(df.b,df.s,phi)
}



# Function to average the estimates for a MGLMM. Uses the results from glmmMulti()
format1<-function(mod.list,formula,data,family){
  # Calling format0 and averaging the estimates
  tes<-lapply(mod.list,format0)
  tes.1<-do.call(rbind,tes)
  df.m<-aggregate(value~Parametro,tes.1,mean)

  # 'mer' structure for the m-dimensional model. Also used in llik.glmm().
  m1<-glmer(formula,data,family,doFit=F)
  m1$formula<-formula

  # Ordered fixed effects, residual std. dev. and var-cov matrix
  matn<-attr(m1$fr$fixef,"names")
  ind<-match(matn,df.m$Parametro)
  est<-as.vector(df.m[ind,2]); names(est)<-df.m[ind,1]

  ind.p<-match('phi',df.m$Parametro)
  phi<-df.m[ind.p,2]

  SigHat<-m1$FL$trms[[1]]$ST
  indS<-!seq(1,nrow(df.m)) %in% c(ind,ind.p)
  SigHat[lower.tri(SigHat,diag=T)]<-df.m[indS,2]
  SigHat<-SigHat+t(lower.tri(SigHat)*SigHat)

  list(m1=m1,est=est,phi=phi,SigHat=SigHat)
}



# Rebuild the log-likelihood
llik.glmm<-function(mod,beta,S,phi=NULL,par=FALSE,cl=NULL){
  # List with elements to compute the likelihood for each subject
  f0<-function(mod,data){
    formula<-eval(mod$formula)

    fix.form<-lme4:::nobars(mod$formula)
    fr<-model.frame(fix.form,data=data)
    Fr<-list()
    Fr$Y<-model.response(fr); Fr$X<-model.matrix(attr(fr,"terms"),data=fr)
    Fr$wts<-mod$fr$wts; Fr$off<-mod$fr$off
    Fr$mf<-data; Fr$fixef<-mod$fr$fixef

    fl<-lme4:::lmerFactorList(formula,Fr,0L,0L)

    x<-Fr$X; y<-Fr$Y; z<-t(as.matrix(fl$trms[[1]]$Zt))

    ffun<-if(is.null(mod$glmFit)) 
            gaussian() else mod$glmFit$family

    fam<-ffun$family

    list(y=y,x=x,z=z,fam=fam,ffun=ffun)#
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
  
  
  # Splitting the original data frame by subject
  mf.ind<-split(mod$fr$mf,mod$FL$fl)
  nr<-length(mf.ind)

  fam<-mod$glmFit$family$family
  if(fam=='binomial') ndat<-split(mod$glmFit$prior.weights,mod$FL$fl) 

  # Log-likelihood evaluation
  if(par && getDoParWorkers()>1){
    b<-foreach(i=1:nr,.combine=c) %dopar% {
      r<-f0(mod,mf.ind[[i]])
      if(fam=='binomial') r$n<-ndat[[i]]
      b.i<-banint(f1,start=runif(ncol(r$z),-5,5),r=r,beta=beta,
                  S=S,phi=phi,trans='none')
      b.i$nrmcon[1]
    }
  }

  if(!par){
    b<-foreach(i=1:nr,.combine=c) %do% {
      r<-f0(mod,mf.ind[[i]])
      if(fam=='binomial') r$n<-ndat[[i]]
      b.i<-banint(f1,start=runif(ncol(r$z),-5,5),r=r,beta=beta,
                  S=S,phi=phi,trans='none')
      b.i$nrmcon[1]
    }
  }

  sum(b)
}



# Combining the above functions inside the main MGLMM function (to be called)
mglmm<-function(formula,id,family,data,llik=TRUE,par=FALSE,cl=NULL){
  # Pair-wise fitting
  m0<-glmmMulti(formula,id,family,data,par,cl)

  # Averaging the estimates from the pairwise models
  prm<-format1(m0,formula,data,family)

  # Log-likelihood
  if(llik)
    LL<-llik.glmm(prm$m1,prm$est,prm$SigHat,prm$phi,par,cl)
  else
    LL<-NA

  SD<-ifelse(prm$phi>0,prm$phi,NA)

  # Output - Fixed effects and var-cov estimates and log-likelihood
  list(Estimates=list(fixef=prm$est,VarCov=prm$SigHat,SD=SD),logLik=LL)
}