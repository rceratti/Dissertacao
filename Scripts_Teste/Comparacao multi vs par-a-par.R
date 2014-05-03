library(pair.mglmm)

# Simulação de dados
phi <- 1; p <- 1.6 

mydat <- data.sim(3, 'CP', exp, xi=p, phi=phi)
dat <- mydat$Data


cl<-makeCluster(4)
registerDoParallel(cl)
clusterEvalQ(cl, library(pair.mglmm))

mainForm0 <- value~-1+variable:period+(-1+variable|ID)
mainForm1 <- value~-1+variable+(-1+variable|ID)

system.time(m0.1 <- mglmmCP(mainForm0, dat$variable, dat))
system.time(m1.1 <- mglmmCP(mainForm1, dat$variable, dat))

stopCluster(cl)


system.time(m0 <- cpglmm(mainForm0, data=dat))
system.time(m1 <- cpglmm(mainForm1, data=dat))

m0
m0.1

m1
m1.1

logLik(m0)-logLik(m1)
m0.1$logLik-m1.1$logLik