y<-rtweedie(10,1.2,.8,1)

dcp<-function(y,p,mu,phi){
  if(any(p<=1)|any(p>=2))  stop("p deve estar entre 1 e 2")
  if(any(phi<=0))          stop("phi deve ser maior que 0")
  if(any(mu <=0))          stop("mu deve ser maior que 0")
  if(any(y  < 0))          stop("y deve ser não negativo")

  if(length(mu) ==1) mu <-rep(mu, length(y))
  if(length(phi)==1) phi<-rep(phi,length(y))

  dev<-(y^(2-p))/((1-p)*(2-p))-(y*(mu^(1-p)))/(1-p)+(mu^(2-p))/(2-p)
  dev<-2*dev

  eps<-1/10
  y.eps<-y+eps

  dens<-y
  dens<-(2*pi*phi*y.eps^p)^(-1/2)*exp(-dev/(2*phi))

  if(any(y==0)){
    y0<-(y==0)
    lambda<-mu[y0]^(2-p)/(phi[y0]*(2-p))
    dens[y0]<-exp(-lambda)
  }

  return(dens)
}

dtweedie(y,1.2,.8,1)
dtweedie.saddle(y,1.2,.8,1)
dcp(y,1.2,.8,1)


x<-seq(0,10,l=200)
y<-dtweedie(x,1.2,.8,1)
plot(y~x,ty='l')

x<-rtweedie(1e6,1.2,.8,1)
hist(x,fr=F,br=40)