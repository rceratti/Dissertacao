library(lme4) 
library(numDeriv)

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


ll<-function(th,fr){
  y<-fr$y; X<-fr$X; z<-fr$z
  beta<-th[1:ncol(X)]
  sig<-th[(ncol(X)+1)]
  
  Q.u<-function(u){
    eta<-exp(X%*%beta+z%*%u)
    f.y<-dpois(y,as.vector(eta),log=T)
    l.yu<-sum(f.y)+sum(dnorm(u,sd=sig,log=T))
    return(-l.yu)
  }
  
  uhat<-nlminb(rep(0,ncol(z)),Q.u)
  hess.u<-hessian(Q.u,uhat$par)
  
  l.y<--uhat$objective-.5*log(abs(det(hess.u)))
  
  return(-l.y)
}


main <- size ~ period+(1|herd)
m1 <- glmer(main, cbpp, poisson)
theta <- c(fixef(m1),0.1947178)

cbpp.list <- split(m1@frame, m1@flist)

m1gh <- lapply(cbpp.list, function(x) {
  fr <- frames(main, x)
  g <- grad(ll, theta, fr=fr)
  h <- hessian(ll, theta, fr=fr)
  list(g=g, h=h)
})



# Lista com os gradientes e hessianas
jkl<-m1gh
l.jkl<-length(jkl)

# Matriz K
gs0<-lapply(jkl,function(x){
  tes1<-x[['g']]
  tes2<-tcrossprod(tes1)
  as.matrix(tes2)
})

K<-(1/length(gs0))*Reduce('+',gs0)


# Matriz J
hs0<-(-1/length(jkl))*Reduce('+',lapply(jkl,"[[",'h'))

J<-hs0
Jinv<-solve(J)

JKJ<-(Jinv%*%K%*%Jinv)/length(jkl)
sqrt(diag(JKJ))
m1
