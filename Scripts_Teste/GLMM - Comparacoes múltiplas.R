library(multcomp)
library(lme4)
library(contrast)
 
warpbreaks$id<-rep(1:27,each=2)

mod<-glmer(breaks~wool+tension+(1|id),warpbreaks,poisson)
# mod<-glmer(cbind(breaks,70-breaks)~wool+tension+(1|id),warpbreaks,binomial)

a<-list(wool="A",tension="M")
b<-list(wool=c("A","B"),tension=c("L","H"))

contrMer(mod,a,b,'less',tes=adjusted(type='fdr'))
contrMer(mod,a,b,tes=adjusted(type='fdr'))



contrMer<-function(fit,a,b,Ha='two.sided',...){
  fr<-fit@frame; fam<-paste(fit@call$family)
  form<-attr(fr,"terms"); D<-dim(fr[,1])

  if(!is.null(D)){
    y<-model.response(fr); fr<-fr[,-1]
    fr<-data.frame(y[,1],fr)
    names(fr)[1]<-names(as.data.frame(y))[1]
  }

  tes<-contrast(glm(form,fam,fr),a=a,b=b)


  rename<-function(x){
    nomes<-names(x)  

    b.1<-lapply(nomes,function(x,lista){
      p1<-substr(x,1,1); p2<-lista[[x]]
      p3<-paste(p1,p2,sep="")
    },lista=x)

    b.2<-expand.grid(b.1)

    new<-b.2[,1]; n<-ncol(b.2)
    for(i in 2:n){
      old<-new
      new<-paste(old,b.2[,i],sep="_")
    }
    return(new)
  }

  n.a<-rename(a); n.b<-rename(b)
  new.names<-outer(n.a,n.b,function(x,y) paste(x,y,sep=" - "))
  rownames(tes$X)<-c(new.names)

  summary(glht(mod,linfct=tes$X,alternative=Ha),...)
}