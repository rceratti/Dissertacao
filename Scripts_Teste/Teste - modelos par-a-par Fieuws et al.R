library(pair.mglmm)


## Dados de questionário - Fieuws et al.
base<-'~/Rubem_Ceratti/Outros/pairwiseMGLMM/'  # IPEA
base<-'C:/Users/Rubem/Dropbox/Dissertação/'    # Casa

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


#
cl<-makeCluster(detectCores())
registerDoParallel(cl)
clusterEvalQ(cl,c(library(cplm),library(bayespack),library(mvtnorm)))


# Modelo 'A' -- a_k + b_k*DC_i + u_ik
gm1.0<-mglmm(resp~-1+item+DC:item+(-1+item|id),quest2$item,binomial,quest2,TRUE,TRUE,cl)

# Modelo 'B' -- a_k + b*DC_i + u_ik (efeito de DC comum à todos os itens)
gm1.1<-mglmm(resp~-1+item+DC+(-1+item|id),quest2$item,binomial,quest2,TRUE,TRUE,cl)

# Comparaçao da log-verossimilhança -- Artigo: Chi^2_obs = 12.33
2*(gm1.0$logLik-gm1.1$logLik)  # 18.7521

# Modelos A e B via glmer()
m0<-glmer(resp~-1+item+DC:item+(-1+item|id),quest2,binomial)
m1<-glmer(resp~-1+item+DC+(-1+item|id),quest2,binomial)

2*(logLik(m0)-logLik(m1)) # 11.014 -- uh-oh...



# Modelo 'C' -- a_k + b_k*DC_i + u_i (efeito aleatório univariado)
gm1.2<-mglmm(resp~-1+item+DC:item+(1|id),quest2$item,binomial,quest2,F,TRUE,cl)

stopCluster(cl)


# PCA gm1.0
pr<-eigen(cov2cor(gm1.0$Estimates$VarCov),only.values=T)
pr$values/sum(pr$values)
