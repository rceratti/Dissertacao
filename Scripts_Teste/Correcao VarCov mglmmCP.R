library(pair.mglmm)

# Simulated data
phi <- 1; p <- 1.6 

mydat <- data.sim(3, 'CP', exp, xi=p, phi=phi)
dat <- mydat$Data

dat$x <- rnorm(nrow(dat))

mainForm0 <- value ~ -1 + variable + variable:period + (-1+variable|ID)
mainForm1 <- value ~ -1 + variable + (-1+variable|ID)


# Pairwise models
cl<-makeCluster(4)
registerDoParallel(cl)
clusterEvalQ(cl, library(pair.mglmm))

system.time(m0.1 <- mglmmCP(mainForm0, dat$variable, dat))
system.time(m1.1 <- mglmmCP(mainForm1, dat$variable, dat))

stopCluster(cl)


# Multivariate models
system.time(m0 <- cpglmm(mainForm0, data=dat))
system.time(m1 <- cpglmm(mainForm1, data=dat))


# Model comparison
m0
summary.cp(m0.1)

m1
summary.cp(m1.1)


mtes <- lmer(mainForm0, data = dat, doFit = F)
sqrt(diag(rcov2(mtes, mainForm0, VarCorr(m0)[[1]], m0@phi, m0@p, fitted(m0))))
sqrt(diag(vcov(m0)))

sqrt(diag(rcov2(mtes, mainForm0, m0.1$VarCov, m0.1$phi, m0.1$p, fitted.cp(m0.1))))


mtes <- lmer(mainForm1, data = dat, doFit = F)
sqrt(diag(rcov2(mtes, mainForm1, m1.1$VarCov, m1.1$phi, m1.1$p, fitted.cp(m1.1))))
sqrt(diag(rcov2(mtes, mainForm1, VarCorr(m1)[[1]], m1@phi, m1@p, fitted(m1))))
sqrt(diag(vcov(m1)))


#::::::::::::::::#
rcov2 <-
function (mod, formula, S, phi, p, fit) 
{
    q <- length(unique(mod$FL$fl[[1]]))
    Lambdat <- chol(S/phi)
    Lambdat <- lapply(1:q, function(i) Lambdat)
    Lambdat <- bdiag(Lambdat)
    fzt <- lme4:::findbars(formula[[3]])[[1]]
    formZt <- paste("~", deparse(fzt[[2]]), ":", fzt[[3]], collapse = "")
    Zt <- t(model.matrix(as.formula(formZt), mod$fr$mf))
    W <- diag(fit^(2 - p))
    LtZt <- Lambdat %*% Zt
    LtZtWZL <- LtZt %*% W %*% t(LtZt)
    q <- ncol(LtZtWZL)
    L_z <- t(chol(LtZtWZL + diag(1, q)))
    X <- mod$fr$X
    L_xz <- t(solve(L_z) %*% LtZt %*% W %*% X)
    L_x <- chol(t(X) %*% W %*% X - L_xz %*% t(L_xz))
    phi * chol2inv(L_x)
}