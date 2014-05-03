library(lme4)  # Vesao de desenvolvimento
library(numDeriv)

(m1 <- glmer(size ~ period+(1|herd), cbpp, poisson))

m1dev <- update(m1, devFunOnly=T)

m1th <- getME(m1, "theta")
grad(m1dev, m1th)
hessian(m1dev, m1th)




cbpp.list <- split(m1@frame, m1@flist)

m1.list <- lapply(cbpp.list, function(x) update(m1, data = x))


th <- t(.Call(lme4:::mer_ST_chol, m1)[[1]])
th <- th[lower.tri(th, diag=T)]

beta<-fixef(m1)


g<-grad(llfun,c(th,beta),mod.up=gm3)
 


