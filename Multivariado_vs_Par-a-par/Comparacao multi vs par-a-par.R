library(pair.mglmm)

# Simulação de dados
phi <- 1; p <- 1.6 

mainForm0 <- value ~ -1 + variable:period + (-1 + variable|ID)
m <- c(3, 5, 7, 10, 13, 15, 20)

# Ajuste dos modelos multivariado e par-a-par
tempo <- lapply(m, function(i) {
  a.1 <- data.sim(i, 'CP', exp, xi=p, phi=phi)$Data

  t1 <- system.time(m0 <- cpglmm(mainForm0, data = a.1))

  cl <- makeCluster(4)
  registerDoParallel(cl)
  clusterEvalQ(cl, library(pair.mglmm))

  t2 <- system.time(m0.1 <- mglmmCP(mainForm0, a.1$variable, a.1))

  stopCluster(cl)

  data.frame(m = i, Multi = t1[3], Par = t2[3])
})


tempo1 <- do.call(rbind, tempo)
rownames(tempo1) <- NULL

# Tempo de execução do modelo par-a-par em função do número de pares de variáveis
m0 <- lm(Par ~ -1 + choose(m, 2), tempo1)
summary(m0) # ~ 3.19 (seg/par)

p <- ggplot(tempo2, aes(choose(m, 2), Par2))
p <- p + geom_point() + geom_smooth(aes(choose(m, 2), fitted(m0)), stat="identity")
p <- p + xlab("Número de pares") + ylab("Tempo (s)")
lab <- as.character(expression(Tempo %~~% 3.19~m(m-1)/2))
p <- p + geom_text(aes(50, 500, label = lab), parse = T)
(p <- p + opts(legend.position = "none"))

pdf('Tempo_m_Par.pdf', w = 11)
print(p)
dev.off()


names(tempo1)[3] <- 'Par-a-par'

tempo2 <- melt(tempo1, id = 'm')
names(tempo2)[2] <- 'Estimação'

# Comparação entre os tempos de execução dos modelos multivariado e par-a-par
p <- ggplot(tempo2, aes(m, value, colour = Estimação, group = Estimação))
p <- p + geom_point()
(p <- p + xlab("Dimensão") + ylab("Tempo (s)"))


pdf('Multi_vs_Par.pdf', w = 11)
print(p)
dev.off()
