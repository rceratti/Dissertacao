library(pair.mglmm)

base <- "~/Rubem_Ceratti/Outros/Dados - algodão"

setwd(base)

# Leitura da base
a <- read.csv('Reprodutivo_Final.csv', h=T, stringsAsFactors=T)
a <- subset(a, !a$composto %in% 'C26') # Remoção do C26 


# Seleção de alguns compostos
a$resp <- a$resp*1e3
a$tmp <- as.integer(substr(a$tempo, 1, 2))

mainForm0 <- resp ~ -1 + composto:trt + (-1 + composto|id)
m <- c(3, 5, 10, 15)

tempo <- lapply(m[3], function(i) {
  ind <- paste("C", 1:i, sep="")
  a.1 <- subset(a, a$composto %in% ind)

  t1 <- system.time(m0 <- cpglmm(mainForm0, data = a.1))

  cl<-makeCluster(4)
  registerDoParallel(cl)
  clusterEvalQ(cl, library(pair.mglmm))

  t2 <- system.time(m0.1 <- mglmmCP(mainForm0, a.1$composto, a.1))

  stopCluster(cl)

  data.frame(m = i, Multi = t1[3], Par = t2[3])
})


tempo1 <- do.call(rbind, tempo)
rownames(tempo1) <- NULL

tempo1$Razao <- with(tempo1, Par/Multi)

#   m  Multi    Par    Razao
# 1 3  29.97 239.06 7.976643
# 2 5 148.35 858.54 5.787260


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#


library(pair.mglmm)

# Simulação de dados
phi <- 1; p <- 1.6 

mainForm0 <- value~-1+variable:period+(-1+variable|ID)
m <- c(3, 7, 10, 15)

tempo <- lapply(m, function(i) {
  a.1 <- data.sim(i, 'CP', exp, xi=p, phi=phi)$Data

  t1 <- system.time(m0 <- cpglmm(mainForm0, data = a.1))

  cl<-makeCluster(4)
  registerDoParallel(cl)
  clusterEvalQ(cl, library(pair.mglmm))

  t2 <- system.time(m0.1 <- mglmmCP(mainForm0, a.1$variable, a.1))

  stopCluster(cl)

  data.frame(m = i, Multi = t1[3], Par = t2[3])
})


tempo1 <- do.call(rbind, tempo)
rownames(tempo1) <- NULL

tempo1$Razao <- with(tempo1, Par/Multi)



# CPU: Intel Core i5 750 (IPEA)
# 30 níveis do fator de efeito aleatório
#    m   Multi     Par     Razao
# 1  3   18.13  113.82 6.2779923
# 2  5  104.55  414.09 3.9606887
# 3  7  294.98  814.11 2.7598820
# 4 10 1740.00 1776.57 1.0210172
# 5 13 4684.29 3174.71 0.6777356
# 6 15 8509.26 3960.32 0.4654130
# 7 20      NA 7354.24        NA

tempo1 <- structure(list(m = c(3, 5, 7, 10, 13, 15, 20), Multi = c(18.1300000000001, 
                    104.55, 294.98, 1740, 4684.29, 8509.26, NA), 
                    Par = c(113.82, 414.09, 814.11, 1776.57, 3174.71, 3960.32, 
                    7354.24), Razao = c(6.27799227799224, 3.96068866571019, 
                    2.75988202590006, 1.02101724137931, 0.677735579991845, 
                    0.46541297363108, NA)), .Names = c("m", "Multi", "Par", "Razao"
                    ), row.names = c(NA, -7L), class = "data.frame")


plot(tempo1$m, tempo1$Multi, pch = 20)
points(tempo1$m, tempo1$Par, pch = 20, col = 4)

plot(tempo1$m, tempo1$Razao, pch = 20)

plot(choose(tempo1$m, 2), tempo1$Par, pch = 20, col = 4)
summary(m0 <- lm(tempo1$Par~0+choose(tempo1$m, 2)))
lines(choose(tempo1$m, 2), fitted(m0), col = 4)



# CPU: Intel Core i5 2430M (casa)
# 30 níveis do fator de efeito aleatório
#    m   Multi     Par     Razao
# 1  3   21.15  167.57 7.9229314
# 2  7  336.27 1077.93 3.2055491
# 3 10  973.26 2370.23 2.4353513
# 4 15 9124.22 5615.68 0.6154696

tempo1 <- structure(list(m = c(3, 7, 10, 15), Multi = c(21.15, 336.27, 
          973.26, 9124.22), Par = c(167.57, 1077.93, 2370.23, 5615.68), 
          Razao = c(7.92293144208038, 3.20554911232046, 2.43535129359061, 
          0.615469596305218)), .Names = c("m", "Multi", "Par", "Razao"
          ), row.names = c(NA, -4L), class = "data.frame")
