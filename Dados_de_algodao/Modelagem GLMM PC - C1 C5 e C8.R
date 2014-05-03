library(reshape)
library(lattice)
library(car)
library(snow)
library(MASS)
library(boot)
library(statmod)
library(cplm)
library(ggplot2)
library(xtable)
library(data.table)


base<-"~/Rubem_Ceratti/Outros/Dados - algodão"

setwd(base)

# Leitura da base
a <- read.csv('Reprodutivo_Final.csv', h = T, stringsAsFactors = T)

a <- subset(a, !a$composto %in% 'C26') # Remoção do C26 - aparece apenas para alguns



# Seleção de alguns compostos
ind <- paste("C", c(1, 5, 8), sep = "")
a.1 <- subset(a, a$composto %in% ind)
a.1$resp <- a.1$resp*1e3
a.1$tmp <- as.integer(substr(a.1$tempo, 1, 2))



# GLMM Poisson Composta - Especificação multivariada - TEMPO VARIÁVEL CONTÍNUA
# Com interaçao entre tempo e tratamento
modPCM <- cpglmm(resp ~ -1+composto:trt+composto:tmp+composto:trt:tmp+
                 (-1+composto|id), data = a.1)
summary(modPCM)  # AIC=2392, BIC=2672


# Apenas os efeitos principais
modPCM1 <- cpglmm(resp ~ -1+composto:trt+tmp:composto+(-1+composto|id), data = a.1)
summary(modPCM1)


# Apenas efeito de tratamento
modPCM2 <- cpglmm(resp ~ -1+composto:trt+(-1+composto|id), data = a.1)
summary(modPCM2)


# Teste da razão de verossimilhança
anova(modPCM, modPCM1, modPCM2)  

#         Df    AIC    BIC  logLik   Chisq  ChiDf  Pr(>Chisq)
# modPCM2 22 2406.9 2498.7 -1181.4
# modPCM1 25 2373.3 2477.6 -1161.6 39.5693      3   1.315e-08 ***                         
# modPCM  37 2387.6 2542.0 -1156.8  9.7185     12      0.6406


# Gráfico de dispersão dos efeitos aleatórios
plotmatrix(ranef(modPCM1)$id, colour = 'darkblue')
attr(VarCorr(modPCM1)$id, "correlation")

modPCM <- modPCM1



# Organizaçao da matriz de coeficientes do modPCM
betas<-fixef(modPCM)
df.b<-melt(betas)
df.b$Parametro<-rownames(df.b)
rownames(df.b)<-1:nrow(df.b)
df.b<-df.b[,c(2,1)]

se<-sqrt(diag(modPCM@vcov))
df.b$se<-se

df.b$composto<-substr(df.b$Parametro,1,10)

n<-nchar(df.b$Parametro)
df.b$Parametro<-substr(df.b$Parametro,12,n)

est<-recast(df.b,Parametro+variable~composto,id.var=c(1,4))
est[1:12,3:5]<-round(est[1:12,3:5],3)

est.1<-matrix(NA,6,3)
for(i in 1:6){
  est[2*i,3:5]<-paste("(",est[2*i,3:5],")",sep="")
  est.1[i,]<-paste(est[2*i-1,3:5],est[2*i,3:5])
}

est.2<-data.frame(Parametro=unique(est$Parametro),est.1)
names(est.2)[2:4]<-names(est)[3:5]

print(xtable(est.2),include.rownames=FALSE)

modPCM@p
modPCM@phi

sqrt(diag(attr(modPCM@vcov,"phi_p")))


a.2<-a.1[,c(1,2,4)]
a.2$pred<-exp(model.matrix(modPCM)%*%fixef(modPCM))
a.2<-unique(a.2)


p<-ggplot(a.2,aes(tempo,pred,group=trt,colour=trt))
p<-p+geom_line()
(p<-p+facet_grid(.~composto)+xlab("Tempo")+ylab("Massa"))

pdf("pred_PCmult.pdf",w=11)
print(p)
dev.off()


p<-ggplot(a.1,aes(tempo,resp,group=trt,colour=trt))
p<-p+geom_point()+geom_line(data=a.2,aes(tempo,pred,group=trt,colour=trt))
p+facet_grid(.~composto)+xlab("Tratamento")+ylab("Massa")


# Error plot dos coeficientes - gráfico ruim...
p<-ggplot(df.b,aes(Parametro,value))
p+geom_errorbar(aes(x=Parametro,y=NULL,ymin=value-2*se,ymax=value+2*se),width=.2)+
  geom_point(aes(Parametro,value))+facet_grid(.~composto)



#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

# Gráfico meio normal
library(doSNOW)
library(grid)


## Gráfico meio-normal com envelope simulado para CP-GLMM (cplm)
hnplot<-function(glmfit,nsim=20,par=FALSE,type,Var){
  
  myres<-function(glmfit,type){   
    y<-glmfit@y
          
    if(type=='marg'){
      X<-model.matrix(glmfit)
      eta<-as.vector(X%*%fixef(glmfit))
      mu<-exp(eta)
    }
    if(type=='cond'){
      mu<-fitted(glmfit)
    }
    
    V<-mu^glmfit@p  
    (y-mu)/sqrt(V)
  }
    
  res<-myres(glmfit,type)
  res<-split(res,Var,drop=T)
  res.1<-lapply(res,function(x) sort(abs(x)))


  sim.hnp<-function(i,glmfit,type,Var){
    dat.1<-glmfit@frame
    dat.1[,1]<-rtweedie(nrow(dat.1),glmfit$p,fitted(glmfit),glmfit$phi)
  
    mfun<-update(glmfit,data=dat.1)

    rp<-myres(mfun,type)
    rp<-split(rp,Var,drop=T)
    lapply(rp,function(x) sort(abs(x)))
  }


  if(par){
    tes<-foreach(i=1:nsim) %dopar% sim.hnp(i,glmfit,type,Var)
    tes<-lapply(1:length(res.1),function(i) t(sapply(tes,"[[",i)))
  } else {
    tes<-lapply(1:nsim,sim.hnp,glmfit=glmfit)
    tes<-do.call(rbind,tes)
  }

  linf<-lapply(tes,function(x) apply(x,2,quantile,.01))
  lsup<-lapply(tes,function(x) apply(x,2,quantile,.99))
  meio<-lapply(tes,function(x) apply(x,2,mean))

  nobs<-ncol(tes[[1]]); t.i<-1:nobs
  ind<-qnorm((t.i+nobs-1/8)/(2*nobs+1/2))
  
  tes.1<-lapply(1:length(res.1),function(i){
    data.frame(ind,linf=linf[[i]],meio=meio[[i]],lsup=lsup[[i]],res.1=res.1[[i]])
  })

  return(tes.1)
}


## Inicialização do 'cluster' e geração do gráfico meio-normal
cl<-makeCluster(4)            # library(doMC); registerDoMC()
registerDoSNOW(cl)
clusterEvalQ(cl,{library(tweedie);library(cplm)})

graf0<-hnplot(modPCM1,99,T,'cond',a.1$composto)
graf1<-hnplot(modPCM1,99,T,'marg',a.1$composto)

stopCluster(cl)


pdf("hnp_cpmult_cond.pdf",w=11)
pushViewport(viewport(layout=grid.layout(1,3))) 

for(i in 1:3){
  p<-ggplot(graf0[[i]],mapping=aes(ind,meio))+geom_line()
  p<-p+geom_line(mapping=aes(ind,linf),linetype="dashed")
  p<-p+geom_line(mapping=aes(ind,lsup),linetype="dashed")
  p<-p+geom_point(mapping=aes(ind,res.1),colour="red")
  p<-p+xlab("Valor esperado do quantil \nmeio normal")
  p<-p+ylab("Resíduos absolutos ordenados")
  print(p,vp=viewport(layout.pos.row=1,layout.pos.col=i))
}

dev.off()



pdf("hnp_cpmult_marg.pdf",w=11)
pushViewport(viewport(layout=grid.layout(1,3))) 

for(i in 1:3){
  p<-ggplot(graf1[[i]],mapping=aes(ind,meio))+geom_line()
  p<-p+geom_line(mapping=aes(ind,linf),linetype="dashed")
  p<-p+geom_line(mapping=aes(ind,lsup),linetype="dashed")
  p<-p+geom_point(mapping=aes(ind,res.1),colour="red")
  p<-p+xlab("Valor esperado do quantil \nmeio normal")
  p<-p+ylab("Resíduos absolutos ordenados")
  print(p,vp=viewport(layout.pos.row=1,layout.pos.col=i))
}

dev.off()


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

# Comparações múltiplas
library(contrast)
library(multcomp)


betas<-fixef(modPCM)
df.b<-melt(betas)
df.b$Parametro<-rownames(df.b)
rownames(df.b)<-1:nrow(df.b)
df.b<-df.b[,c(2,1)]

df.b$composto<-substr(df.b$Parametro,9,10)

n<-nchar(df.b$Parametro)
df.b$Parametro<-ifelse(n==14,substr(df.b$Parametro,12,n),
                             substr(df.b$Parametro,15,n))

df.b<-df.b[,-2]


lhtmer<-function(df,fit,a,b,...){
  a1.1<-expand.grid(a)
  b1.1<-expand.grid(b)

  for(i in 1:nrow(a1.1)){
    x<-a1.1[i,]
    ind.1<-which(x$composto==df$composto & x$trt==df$Parametro)

    for(j in 1:nrow(b1.1)){ 
      y<-b1.1[j,] 
      ind.2<-which(y$composto==df$composto & y$trt==df$Parametro)
      ind.3<-!seq(1,nrow(df)) %in% c(ind.1,ind.2)

      df.b[ind.1,(j+2)]<- 1
      df.b[ind.2,(j+2)]<--1
      df.b[ind.3,(j+2)]<- 0
    }
  }

  tes<-as.matrix(df.b[,-(1:2)])
  tes<-t(tes)
  rownames(tes)<-colnames(tes)<-NULL

  glht(fit,linfct=tes,...)
}

a1<-list(composto="C1",trt="grandis")
b1<-list(composto="C1",trt=c("Controle","frugiperda","heros","mecânico"))

a2<-list(composto="C5",trt="grandis")
b2<-list(composto="C5",trt=c("Controle","frugiperda","heros","mecânico"))

a3<-list(composto="C8",trt="grandis")
b3<-list(composto="C8",trt=c("Controle","frugiperda","heros","mecânico"))

con.list<-list(list(a=a1,b=b1),list(a=a2,b=b2),list(a=a3,b=b3))



tes<-lapply(1:3,function(i,x,y,cx,df){
  y<-y[[i]]
  z.0<-lhtmer(df,x,y$a,y$b,alternative="two.sided")
  z.1<-summary(z.0,test=adjusted("holm"))
  
  pval<-z.1$test$pvalues
  Est<-round(z.1$test$coefficients,3)
  pval<-ifelse(pval<1e-4,paste(Est,"(<1e-4)"),
               paste(Est,paste("(",round(pval,4),")",sep="")))
  
  cont<-c('trtG - trtC','trtG - trtF',
          'trtG - trtH','trtG - trtM')
  
  z<-data.frame(Contraste=cont,p.valor=pval)
  z$Composto<-cx[[i]]
  return(z)
},x=modPCM1,y=con.list,cx=c('C1','C5','C8'),df=df.b)

tes<-do.call(rbind,tes)
tes1<-cast(tes,Composto~Contraste,value='p.valor')
print(xtable(tes1),include.rownames=FALSE)
