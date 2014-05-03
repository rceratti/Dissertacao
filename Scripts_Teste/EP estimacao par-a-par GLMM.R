library(lme4)
library(doParallel)
library(mvtnorm)

## Dados de questionário - Fieuws et al.
base<-'~/Rubem_Ceratti/Outros/pairwiseMGLMM/'  # IPEA
#base<-'C:/Users/Rubem/Dropbox/Dissertação/'    # Casa

quest<-read.table(paste(base,'questiondata.txt',sep=""),h=F,na.strings='.')

names(quest)<-c('id','DC','resp','physical','selfesteem','barriers','psychic',
                'selfperception','selfeffic','motivation')


varLevels<-sapply(1:nrow(quest),function(i,x){
  x<-x[i,]
  names(x)[which(x[,-(1:3)]==1)+3]
},x=quest)

quest1<-quest[,1:3]
quest1$item<-as.factor(varLevels)

quest2<-quest1[complete.cases(quest1),]
quest2$id<-as.factor(quest2$id)
quest2$DC<-as.factor(quest2$DC)


# Estimativas dos parâmetros do modelo com todas as observações
cl<-makeCluster(4)
registerDoParallel(cl)
clusterEvalQ(cl,c(library(lme4),library(numDeriv),library(mvtnorm)))

gm1<-glmmMulti(resp~-1+item+DC:item+(-1+item|id),quest2$item,binomial,quest2,T,cl)

stopCluster(cl)

# Lista apenas com os modelos
mods<-lapply(gm1,'[[','mod')

# Lista com os gradientes e hessianas
jkl<-lapply(gm1,'[[','jkl')
l.jkl<-length(gm1[[1]][[2]])

# Matriz K
gs0<-lapply(1:l.jkl,function(i,x){
  tes0<-lapply(x,'[[',i)
  tes1<-lapply(tes0,'[[','g')
  tes2<-lapply(1:length(tes1),function(i) lapply(tes1,outer,tes1[[i]]))
  tes3<-lapply(tes2,do.call,what=cbind)
  tes4<-do.call(rbind,tes3)
  as.matrix(tes4)
},x=jkl)

K<-(1/length(gs0))*Reduce('+',gs0)


# Matriz J
hs0<-lapply(jkl,function(x){
  tes0<-lapply(x,'[[','h')
  (-1/length(x))*Reduce('+',tes0)
})

J<-bdiag(hs0)
Jinv<-solve(J)

JKJ<-(Jinv%*%K%*%Jinv)/l.jkl


# Criação da matriz A
tes<-lapply(mods,pair.mglmm:::format0)
tes.1<-do.call(rbind,tes)
tes.1<-subset(tes.1,tes.1$Parametro!='phi')

df.l<-sort(unique(tes.1$Parametro))

mat<-lapply(df.l,'==',tes.1$Parametro)
mat<-do.call(rbind,mat)

mat<-sweep(mat,1,rowSums(mat),"/") # matriz A(!)

est<-data.frame(Parametro=df.l,value=mat%*%tes.1$value) # ok!

mat<-mat[c(4:6,8:9,11,3,7,10,2,1),]
EPmat<-mat%*%JKJ%*%t(mat)

ep<-sqrt(diag(EPmat))[c(11,10,7,1:3,8,4:5,9,6)]
est$ep<-ep


##
# Function that fits all the bivariate models 
glmmMulti<-function(formula,id,family,data,par=FALSE,cl=NULL){
  # All pairwise combinations of response variables
  lev<-paste(unique(id))
  tes<-combn(lev,m=2)

  # log-likelihood function 
  ll<-function(th,fr){
    y<-fr$y; X<-fr$X; z<-fr$z
    beta<-th[1:ncol(X)]
    S<-matrix(0,ncol(z),ncol(z))
    S[lower.tri(S,diag=T)]<-th[(ncol(X)+1):length(th)]
    S<-S+t(lower.tri(S)*S)
  
    Q.u<-function(u){
      mu<-plogis(X%*%beta+z%*%u)
      f.y<-dbinom(y,1,as.vector(mu),log=T)
      l.yu<-sum(f.y)+sum(dmvnorm(u,sigma=S,log=T))
      return(-l.yu)
    }
  
    uhat<-nlminb(rep(0,ncol(z)),Q.u)
    hess.u<-hessian(Q.u,uhat$par)
  
    l.y<--uhat$objective-.5*log(abs(det(hess.u)))
  
    return(-l.y)
  }

  # Frames
  lmerFrames2 <- function (mc, formula, contrasts, vnms = character(0)) {
    mf <- mc
    m <- match(c("data", "subset", "weights", "na.action", "offset"), 
               names(mf), 0)
    mf <- mf[c(1, m)]
    frame.form <- lme4:::subbars(formula)
    if (length(vnms) > 0) 
      frame.form[[3]] <- substitute(foo + bar, 
                                    list(foo = parse(text = paste(vnms,
                                                     collapse = " + "))[[1]], 
                                                     bar = frame.form[[3]]))
    fixed.form <- lme4:::nobars(formula)
    if (!inherits(fixed.form, "formula")) 
      fixed.form <- as.formula(substitute(foo ~ 1, list(foo = fixed.form)))
    environment(fixed.form) <- environment(frame.form) <- environment(formula)
    mf$formula <- frame.form
    mf$drop.unused.levels <- FALSE
    mf[[1]] <- as.name("model.frame")
    fe <- mf
    mf <- eval(mf, parent.frame(2))
    fe$formula <- fixed.form
    fe <- eval(fe, parent.frame(2))
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1) {
      nm <- rownames(Y)
      dim(Y) <- NULL
      if (!is.null(nm)) 
        names(Y) <- nm
    }
    mt <- attr(fe, "terms")
    X <- if (!is.empty.model(mt)) 
           model.matrix(mt, mf, contrasts)
         else matrix(, NROW(Y), 0)
    storage.mode(X) <- "double"
    fixef <- numeric(ncol(X))
    names(fixef) <- colnames(X)
    dimnames(X) <- NULL
    wts <- model.weights(mf)
    if (is.null(wts)) 
      wts <- numeric(0)
    off <- model.offset(mf)
    if (is.null(off)) 
      off <- numeric(0)
    if (any(wts <= 0)) 
      stop(gettextf("negative weights or weights of zero are not allowed"))
    if (length(off) && length(off) != NROW(Y)) 
      stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                    length(off), NROW(Y)))
    attr(mf, "terms") <- mt
    list(Y = Y, X = X, wts = as.double(wts), off = as.double(off), 
         mf = mf, fixef = fixef)
  }

  frames<-function(formula,data){
    contrasts<-NULL
    mc<-match.call()
    fr<-lmerFrames2(mc,formula,contrasts)
    fl<-lme4:::lmerFactorList(formula,fr,0L,0L)
  
    y<-fr$Y; x<-fr$X; z<-t(as.matrix(fl$trms[[1]]$Zt))
  
    list(y=y,X=x,z=z)
  }

  # Função para retornar o gradiente e hessiana
  swGHfun<-function(mod,data,formula){
    dat.s<-split(data,mod@flist[[1]])
    theta<-c(fixef(mod),unique(unlist(VarCorr(mod))))

    res<-lapply(dat.s,function(x){
      fr<-frames(formula,x)
      g<-grad(ll,theta,fr=fr)
      h<-hessian(ll,theta,fr=fr)
      list(g=g,h=h)
    })

    return(res)
  }

  glmm.fit<-function(x,data,formula,family){
    ind<-id %in% x
    dat<-subset(data,ind)

    mod<-glmer(formula,family=family,data=dat)
    jkl<-swGHfun(mod,mod@frame,formula)

    list(mod=mod,jkl=jkl)
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