library(cplm)
library(MASS)
library(reshape)
library(mvtnorm)
library(bayespack)
library(doParallel)



# Simulação de dados
beta.c1<-c(0.70,1.45,1.65,1.90)
beta.c2<-c(0.96,1.39,0.40,1.19)
beta.c3<-c(1.25,1.86,0.19,-0.39)
beta<-matrix(c(beta.c1,beta.c2,beta.c3),4,3)

phi<-1; p<-1.6 

mydat<-data.sim(3,'CP',exp,beta,xi=p,phi=phi)
dat<-mydat$Data



# 
gm0.0<-cpglmm(value~-1+period:variable+(-1+variable|ID),data=dat)
gm0.1<-cpglmm(value~-1+variable+(-1+variable|ID),data=dat)


#
cl<-makeCluster(4)
registerDoParallel(cl)
clusterEvalQ(cl,c(library(cplm),library(bayespack),library(mvtnorm)))

gm1.0<-mglmmCP(value~-1+period:variable+(-1+variable|ID),dat$variable,dat,T,T,cl)
gm1.1<-mglmmCP(value~-1+variable+(-1+variable|ID),dat$variable,dat,T,T,cl)


# Comparação das diferenças nas verossimilhanças
logLik(gm0.1)-logLik(gm0.0)
gm1.1$logLik-gm1.0$logLik



# Função para ajuste par-a-par do modelo multivariado
glmmMultiCP<-function(formula,id,data,par=FALSE,cl=NULL){
  # Combinação dos níveis de composto
  lev<-paste(unique(id))
  tes<-combn(lev,m=2)

  glmm.fit<-function(x,data,formula){
    ind<-id %in% x
    dat0<-subset(data,ind)

    cpglmm(formula,data=dat0)
  }

  # Em paralelo
  if(par==TRUE){
    #clusterEvalQ(cl,library(cplm))
    res<-foreach(i=1:ncol(tes)) %dopar% glmm.fit(tes[,i],data,formula)
  }

  # Em série
  if(par==FALSE) res<-foreach(i=1:ncol(tes)) %do% glmm.fit(tes[,i],data,formula)

  return(res)
}



# Função de formatação das estimativas (parte 1)
format0CP<-function(mod){
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

  phi<-as.numeric(attr(S,'sc'))^2
  phi<-data.frame(Parametro='phi',value=ifelse(is.na(phi),0,phi))

  p<-data.frame(Parametro='p',value=mod@p)

  rbind(df.b,df.s,phi,p)
}



# Função de formatação das estimativas (parte 2)
format1CP<-function(mod.list,formula,data){
  # Estimativas dos parâmetros e formatação (parte 1)
  tes<-lapply(mod.list,format0CP)
  tes.1<-do.call(rbind,tes)
  df.m<-aggregate(value~Parametro,tes.1,mean)

  # Formatação (parte 2):
  # Efeitos fixos (ordenados conforme matriz de delineamento)
  m1<-cpglmm(formula,data=data,doFit=F)

  matn<-names(fixef(m1))
  ind<-match(matn,df.m$Parametro)
  est<-as.vector(df.m[ind,2]); names(est)<-df.m[ind,1]

  ind.phi<-match('phi',df.m$Parametro)
  phi<-df.m[ind.phi,2]

  ind.p<-match('p',df.m$Parametro)
  p<-df.m[ind.p,2]

  SigHat<-m1@ST[[1]]
  indS<-!seq(1,nrow(df.m)) %in% c(ind,ind.phi,ind.p)
  SigHat[lower.tri(SigHat,diag=T)]<-df.m[indS,2]
  SigHat<-SigHat+t(lower.tri(SigHat)*SigHat)

  list(m1=m1,est=est,phi=phi,p=p,SigHat=SigHat)
}



# Log-Verossimilhança reconstruída
llik.glmmCP<-function(mod,beta,S,phi,p,formula,par=FALSE,cl=NULL){

  f0<-function(mod,data){
    formula<-eval(mod@call$formula)

    mc<-match.call()
    fr<-lme4:::lmerFrames(mc,formula,NULL)
    fl<-lme4:::lmerFactorList(formula,fr,0L,0L)

    y<-fr$Y; x<-fr$X; z<-t(as.matrix(fl$trms[[1]]$Zt))

    list(y=y,x=x,z=z)
  }

  f1<-function(r,beta,u,S,phi,p){
    if(class(u)=='data.frame') u<-as.matrix(u)
    if(is.vector(u)) u<-matrix(u,nc=ncol(S))
    u.vec<-matrix(u,nc=1)

    eta<-r$x%*%beta+r$z%*%u.vec
    mu<-exp(eta)

    l.u<-sum(dmvnorm(u,sigma=S,log=T))

    res<-sum(log(tweedie::dtweedie(r$y,p,mu,phi)))+l.u

    return(res)
  }

  mf.ind<-split(mod@frame,mod@flist[[1]])
  nr<-length(mf.ind)

  if(par && getDoParWorkers()>1){
    b<-foreach(i=1:nr,.combine=c) %dopar% {
      r<-f0(mod,mf.ind[[i]])
      ini<-rnorm(ncol(r$z))
      b.i<-banint(f1,start=ini,r=r,beta=beta,S=S,phi=phi,p=p,trans='none')
      b.i$nrmcon[1]
    }
  }
  
  if(!par){
    b<-foreach(i=1:nr,.combine=c) %do% {
      r<-f0(mod,mf.ind[[i]])
      ini<-rnorm(ncol(r$z))
      b.i<-banint(f1,start=ini,r=r,beta=beta,S=S,phi=phi,p=p,trans='none')
      b.i$nrmcon[1]
    }
  }

  sum(b)
}



# Juntando as funções de estimação par-a-par, formatação e log-verossimilhança
mglmmCP<-function(formula,id,data,llik=TRUE,par=FALSE,cl=NULL){
  # Ajuste dos modelos par-a-par
  m0<-glmmMultiCP(formula,id,data,par,cl)

  # Estimativas dos parâmetros e formatação
  prm<-format1CP(m0,formula,data)

  if(llik)
    LL<-llik.glmmCP(prm$m1,prm$est,prm$SigHat,prm$phi,prm$p,formula)
  else
    LL<-NA

  # Resultado final
  list(Estimates=list(fixef=prm$est,VarCov=prm$SigHat,phi=prm$phi,p=prm$p),
       logLik=LL)
}