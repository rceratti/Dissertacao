library(lme4)
library(MASS)
library(Matrix)
library(reshape)
library(mvtnorm)
library(bayespack)
library(doParallel)


# Simulação de dados
beta.c1<-c(0.70,1.45,1.65,1.90)
beta.c2<-c(0.96,1.39,0.40,1.19)
beta.c3<-c(1.25,1.86,0.19,-0.39)
beta<-matrix(c(beta.c1,beta.c2,beta.c3),4,3)

mydat<-data.sim(3,'poisson',exp,beta)
dat<-mydat$Data


## Ajuste dos modelos concorrentes com 3 dimensões
system.time({
gm0.0<-glmer(value~-1+period:variable+(-1+variable|ID),family=poisson,data=dat)
})

system.time({
gm0.1<-glmer(value~-1+variable+(-1+variable|ID),family=poisson,data=dat)
})


## Ajuste par-a-par dos modelos concorrentes -- rodar as funções abaixo antes!
cl<-makeCluster(4)
registerDoParallel(cl)
clusterEvalQ(cl,c(library(lme4),library(bayespack),library(mvtnorm)))

system.time({
gm1.0<-mglmm(value~-1+period:variable+(-1+variable|ID),dat$variable,
             poisson,dat,par=T,cl=cl)
})

system.time({
gm1.1<-mglmm(value~-1+variable+(-1+variable|ID),dat$variable,
             poisson,dat,par=T,cl=cl)
})

stopCluster(cl)

# Comparação das diferenças nas verossimilhanças
logLik(gm0.1)-logLik(gm0.0)
gm1.1$logLik-gm1.0$logLik



# Function to simulate new data
data.sim<-function(m,distr,link.inv,beta=NULL,...){
  if(m<=0) 
    stop("\'m\' must be >= 1")

  if(!is.null(beta) & !all(dim(beta)==c(4,m))) 
    stop("dim(beta) must be 4 by m")

  if(!distr %in% c('poisson','binomial','gaussian','CP'))
    stop("Please specify either 'poisson', 'binomial', 'gaussian' or 'CP'")

  # Fixed and random effects factors
  period<-factor(1:4)
  ID<-factor(1:30)

  dat<-expand.grid(period=period,ID=ID)

  if(is.null(beta))
    beta<-matrix(runif(4*m,-2,2),4,m)
  X<-model.matrix(~-1+period,dat)

  # Function to generate a PD matrix for the RE's
  f.S<-function(n){
    S<-matrix(0,n,n)
    diag(S)<-rep(1,n)

    corr<-runif(n*(n-1)/2,-.9,.9)
    S[lower.tri(S)]<-S[upper.tri(S)]<-corr

    nearPD(S,corr=T)$mat
  }

  S<-f.S(m)
  u<-mvrnorm(length(levels(dat$ID)),rep(0,m),S)
  Z<-model.matrix(~-1+ID,dat)


  Eta<-X%*%beta+Z%*%u
  Mu<-link.inv(Eta)

  Mu<-data.frame(Mu)
  names(Mu)<-paste("C",1:m,sep="")

  dat<-cbind(dat,Mu)
  dat<-melt(dat,id=c('period','ID'))
  dat$value<-switch(distr,
                    poisson=rpois(nrow(dat),dat$value),
                    binomial=rbinom(nrow(dat),prob=dat$value,...),
                    gaussian=rnorm(nrow(dat),mean=dat$value,...),
                    CP=tweedie::rtweedie(nrow(dat),mu=dat$value,...)
              )

  list(Data=dat,beta=beta,S=S)
}



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
  
  
  # Splitting the original data frame by subject
  mf.ind<-split(mod$fr$mf,mod$FL$fl)
  nr<-length(mf.ind)

  # Log-likelihood evaluation
  if(par && getDoParWorkers()>1){
    b<-foreach(i=1:nr,.combine=c) %dopar% {
      r<-f0(mod,mf.ind[[i]])
      b.i<-banint(f1,start=runif(ncol(r$z),-5,5),r=r,beta=beta,
                  S=S,phi=phi,trans='none')
      b.i$nrmcon[1]
    }
  }

  if(!par){
    b<-foreach(i=1:nr,.combine=c) %do% {
      r<-f0(mod,mf.ind[[i]])
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






package.skeleton("pair.mglmm",c('data.sim','format0','format1','llik.glmm','mglmm',
                                'format0CP','format1CP','llik.glmmCP','mglmmCP',
                                'glmmMulti','glmmMultiCP'))



## Testes com especificação do modelo -- Usuário especifica como um modelo 
## univariado e a função muda a fórmula para formato multivariado adequado.
# Edit 1: Talvez seja melhor pensar em uma forma de especificar em três partes:
#         1. Parâmetros que variam conforme a variável resposta;
#         2. Parâmetros com valor comum em todas as variáveis resposta;
#         3. Variável identificadora das variáveis resposta.
form<-resp~trt+tempo*trt+(1|id1)+(tes|id2)
id<-'composto'

form.ch<-paste(form)

pred.ch<-form.ch[grep('\\+',form.ch)]
pred.ch.1<-strsplit(pred.ch,"\\+")[[1]]

re.ind<-grep("\\|",pred.ch.1)
fix.ch<-pred.ch.1[-re.ind]
re.ch<-pred.ch.1[re.ind]

paste(pred.ch.1,":",id)

til<-form.ch[grep('\\~',form.ch)]

resp~-1+trt:composto+tempo:composto+(-1+composto|id)






## Código para criação da matriz A. "A ver" se é necessário mesmo ou se 
## basta manusear as matrizes J e K tal que aggregate() se aplica.
## De qquer forma, fica o código.
gm1.2<-glmmMulti(value~-1+variable+(-1+variable|herd),dat$variable,
                 poisson,dat,par=F)

tes<-lapply(gm1.2,format0)
tes.1<-do.call(rbind,tes)

# Criação da matriz A
df.l<-sort(unique(tes.1$Parametro))

mat<-lapply(df.l,'==',tes.1$Parametro)
mat<-do.call(rbind,mat)

mat<-sweep(mat,1,rowSums(mat),"/") # matriz A(!)

data.frame(Parametro=df.l,value=mat%*%tes.1$value) # ok!
aggregate(value~Parametro,tes.1,mean)






