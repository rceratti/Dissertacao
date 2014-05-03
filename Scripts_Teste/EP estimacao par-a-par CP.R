library(pair.mglmm)
library(numDeriv)


# Simulação de dados
beta.c1<-c(0.70,1.45,1.65,1.90)
beta.c2<-c(0.96,1.39,0.40,1.19)
beta.c3<-c(1.25,1.86,0.19,-0.39)
beta<-matrix(c(beta.c1,beta.c2,beta.c3),4,3)

phi<-1; p<-1.6 

mydat<-data.sim(3,'CP',exp,beta,xi=p,phi=phi)
dat<-mydat$Data

detach('package:pair.mglmm')


# Estimativas dos parâmetros do modelo com todas as observações
cl<-makeCluster(4)
registerDoParallel(cl)
clusterEvalQ(cl,c(library(cplm),library(numDeriv),library(foreach)))

gm1<-glmmMultiCP(value~-1+variable+(-1+variable|ID),dat$variable,dat,T,cl)

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
tes<-lapply(mods,pair.mglmm:::format0CP)
tes.1<-do.call(rbind,tes)

df.l<-sort(unique(tes.1$Parametro))

mat<-lapply(df.l,'==',tes.1$Parametro)
mat<-do.call(rbind,mat)

mat<-sweep(mat,1,rowSums(mat),"/") # matriz A(!)

est<-data.frame(Parametro=df.l,value=mat%*%tes.1$value) # ok!

mat<-mat[c(4:6,8:9,11,3,7,10,2,1),]
EPmat<-mat%*%JKJ%*%t(mat)

ep<-sqrt(diag(EPmat))[c(11,10,7,1:3,8,4:5,9,6)]
est$ep<-ep

m0<-cpglmm(value~-1+variable+(-1+variable|ID),data=dat)


##
# Função para ajuste par-a-par do modelo multivariado
glmmMultiCP<-function(formula,id,data,par=FALSE,cl=NULL){
  # Combinação das variáveis multivariadas
  lev<-paste(unique(id))
  tes<-combn(lev,m=2)

  # Função para retornar o gradiente e hessiana
  swGHfun<-function(mod,data,formula){
    # Estimates
    th<-t(.Call(lme4:::mer_ST_chol,mod)[[1]])
    th<-th[lower.tri(th,diag=T)]
    beta<-fixef(mod); lphi<-log(mod@phi); p<-mod@p

    # Splitting data by subject
    gm2<-cpglmm(formula,data=data,doFit=F)
    dat.s<-split(gm2@frame,gm2@flist[[1]])

    # log-likelihood function 
    llfun<-function(mod.up,pars) -.5*.Call('cpglmm_update_dev',mod.up,pars)

    res<-foreach(i=1:length(dat.s)) %do% {
      gm3<-cpglmm(formula,data=dat.s[[i]],doFit=F)
      g<-grad(llfun,c(th,beta,lphi,p),mod.up=gm3)
      h<-hessian(llfun,c(th,beta,lphi,p),mod.up=gm3)
      list(g=g,h=h)
    }

    return(res)
  }

  glmm.fit<-function(x,data,formula){
    ind<-id %in% x
    dat0<-subset(data,ind)

    mod<-cpglmm(formula,data=dat0)
    jkl<-swGHfun(mod,mod@frame,mod@formula)

    list(mod=mod,jkl=jkl)
  }

  # Em paralelo
  if(par){
    res<-foreach(i=1:ncol(tes)) %dopar% glmm.fit(tes[,i],data,formula)
  }

  # Em série
  if(!par) res<-foreach(i=1:ncol(tes)) %do% glmm.fit(tes[,i],data,formula)

  return(res)
}







