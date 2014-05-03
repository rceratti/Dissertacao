library(cplm)
library(statmod)
library(tweedie)
library(MASS)


FineRoot$Zone <- factor(FineRoot$Zone)
FineRoot$Stock <- factor(FineRoot$Stock)

fit1 <- cpglm(RLD ~ Zone * Spacing, data = FineRoot)
fit2 <- cpglm(RLD ~ Zone + Spacing, data = FineRoot)
fit3 <- cpglm(RLD ~ Zone, data = FineRoot)


devBootPvalue <- function(fit1, fit2, B) {
  #X2.obs <- deviance(fit2) - deviance(fit1)
  d1 <- fit1@model.frame; d2 <- fit2@model.frame
  muhat.fit2 <- fitted(fit2); n <- length(muhat.fit2)

  modboot <- sapply(1:B, function(i) {
    d1[, 1] <- d2[, 1] <- rtweedie(n, fit2@p, muhat.fit2, fit2@phi)
    fit1.1 <- update(fit1, data = d1)
    fit2.1 <- update(fit2, data = d2)
    deviance(fit2.1) - deviance(fit1.1)
  })
  modboot
}


tes <- devBootPvalue(fit1, fit2, 1e2)
tes1 <- devBootPvalue(fit2, fit3, 1e2)

hist(tes, fr = F)
abline(v = deviance(fit2) - deviance(fit1), lty = 2)

x <- seq(min(tes), max(tes), l = 1e2)
lines(x, dchisq(x, fit2$df.residual-fit1$df.residual))

m <- fitdistr(tes, "chi-squared", start = list(df = 2), method = "L-BFGS-B")
lines(x, dchisq(x, m$estimate), lty = 2)
