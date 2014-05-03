library(pair.mglmm)  # pair.mglmm_1.0-3

# Simulação de dados
phi <- 1; p <- 1.6 

mainForm0 <- value~-1+variable:period+(-1+variable|ID)
m <- c(3, 5, 7, 10, 13, 15, 20)

cl<-makeCluster(4)
registerDoParallel(cl)
clusterEvalQ(cl, library(pair.mglmm))

tempo <- lapply(m, function(i) {
  a.1 <- data.sim(i, 'CP', exp, xi=p, phi=phi)$Data
  t2 <- system.time(m0.1 <- mglmmCP(mainForm0, a.1$variable, a.1))
  data.frame(m = i, Par = t2[3])
})

stopCluster(cl)


tempo1 <- do.call(rbind, tempo)
rownames(tempo1) <- NULL

tempo2 <- tempo1
plot(tempo2)
plot(choose(tempo2$m, 2), tempo2$Par)

m0 <- lm(Par ~ -1 + choose(m, 2), tempo2)
summary(m0) # ~ 3.19 (seg/par)

lines(choose(tempo2$m, 2), fitted(m0))


tempo1 <- structure(list(m = c(3, 5, 7, 10, 13, 15, 20), Multi = c(18.1300000000001, 
                    104.55, 294.98, 1740, 4684.29, 8509.26, NA), 
                    Par = c(113.82, 414.09, 814.11, 1776.57, 3174.71, 3960.32, 
                    7354.24), Razao = c(6.27799227799224, 3.96068866571019, 
                    2.75988202590006, 1.02101724137931, 0.677735579991845, 
                    0.46541297363108, NA)), .Names = c("m", "Multi", "Par1", "Razao"
                    ), row.names = c(NA, -7L), class = "data.frame")


tempo3 <- cbind(tempo1[,1:3], Par2 = tempo2[,2], Razao1 = tempo1[,4])


tempo3 <- structure(list(m = c(3, 5, 7, 10, 13, 15, 20), Multi = c(18.1300000000001, 
          104.55, 294.98, 1740, 4684.29, 8509.26, NA), Par1 = c(113.82, 
          414.09, 814.11, 1776.57, 3174.71, 3960.32, 7354.24), Par2 = c(18.66, 
          43.08, 79.15, 150.72, 260.01, 331.05, 600.26), Razao1 = c(6.27799227799224, 
          3.96068866571019, 2.75988202590006, 1.02101724137931, 0.677735579991845, 
          0.46541297363108, NA), Razao2 = c(1.02923331494759, 0.412051649928264, 
          0.268323276154315, 0.0866206896551724, 0.0555068110642168, 
          0.0389046756122154, 
          NA)), .Names = c("m", "Multi", "Par1", "Par2", "Razao1", "Razao2"
          ), row.names = c(NA, -7L), class = "data.frame")

tempo3$Razao2 <- tempo3$Par2/tempo3$Multi

plot(tempo3$m, tempo3$Multi, pch = 20)
points(tempo3$m, tempo3$Par1, pch = 20, col = 2)
points(tempo3$m, tempo3$Par2, pch = 20, col = 4)


library(reshape)
library(ggplot2)

tempo4 <- tempo3[, 1:4]
tempo4 <- melt(tempo4, id = 'm')

p <- ggplot(tempo4, aes(m, value, colour = variable, group = variable))
pdf('Multivariado vs par-a-par.pdf', w = 11)
print(p + geom_point())
dev.off()