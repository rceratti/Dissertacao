library(doParallel)
library(MASS)
library(statmod)
library(cplm)



base <- "~/Rubem_Ceratti/Outros/Dados - algodão"
setwd(base)

# Leitura da base
a <- read.csv('Reprodutivo_Final.csv', h = T, stringsAsFactors = T)
a <- subset(a, !a$composto %in% 'C26') # Remoção do C26 - aparece apenas para alguns



# Seleção de alguns compostos
ind <- paste("C", c(1, 5, 8), sep = "")
a.1 <- subset(a, a$composto %in% ind)
a.1$resp <- a.1$resp*1e3
a.1$tmp <- as.integer(substr(a.1$tempo, 1, 2))


devBootPvalue <- function(fit1, fit2, B) {
  X.obs <- deviance(fit2) - deviance(fit1)  
  muhat.fit2 <- fitted(fit2); n <- length(muhat.fit2)

  bootDev <- vector('numeric', B)

  for(i in 1:B) {
    d1 <- fit1@model.frame; d2 <- fit2@model.frame
    d1[, 1] <- d2[, 1] <- rtweedie(n, fit2@p, muhat.fit2, fit2@phi)
    assign("d1", d1, envir = .GlobalEnv)
    assign("d2", d2, envir = .GlobalEnv)
    fit1.1 <- tryCatch(cpglm(as.formula(fit1@call$formula), data = d1), 
                       error = function(e) NA)
    fit2.1 <- tryCatch(cpglm(as.formula(fit2@call$formula), data = d2),
                       error = function(e) NA)
    if(is.na(fit1.1) | is.na(fit2.1))
      bootDev[i] <- NA
    else
      bootDev[i] <- deviance(fit2.1) - deviance(fit1.1)
    rm(d1, envir = .GlobalEnv)
    rm(d2, envir = .GlobalEnv)
  }
  mean(na.omit(bootDev) > X.obs)
}


cl <- makeCluster(3)
registerDoParallel(cl)
clusterExport(cl, c('a.1.2', 'devBootPvalue'))
clusterEvalQ(cl, {library(cplm); library(tweedie)})
snow::clusterSetupRNGstream(cl, seed = 2645751)


modPC <- foreach(i = 1:3) %dopar% {
  x <- a.1.2[[i]]
  
  modPC0 <- cpglm(resp~trt*tmp, data = x)
  modPC1 <- cpglm(resp~trt+tmp, data = x)
  modPC2 <- cpglm(resp~trt, data = x)

  tes.01 <- devBootPvalue(modPC0, modPC1, 5e2)
  tes.12 <- devBootPvalue(modPC1, modPC2, 5e2)

  d0 <- deviance(modPC0)
  d1 <- deviance(modPC1)
  d2 <- deviance(modPC2)

  data.frame(Modelo = c('mPC0', 'mPC1', 'mPC2'),
             GL = c(modPC0$df.residual, modPC1$df.residual, modPC2$df.residual),
             Deviance = c(d0, d1, d2),
             Di_Dj = c(NA, d1 - d0, d2 - d1),
             P_valor = c(NA, tes.01, tes.12))
}

stopCluster(cl)

modPC