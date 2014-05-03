library(pair.mglmm)
library(reshape)
library(xtable)
library(ggplot2)
library(gdata)
library(calibrate)
library(corrplot)


base <- "~/Rubem_Ceratti/Outros/Dados - algodão"
setwd(base)

# Leitura da base e correspondência dos compostos (código vs nome)
a <- read.csv('Reprodutivo_Final.csv', h=T, stringsAsFactors=T)
a <- subset(a, !a$composto %in% 'C26') # Remoção do C26

codcomp <- read.csv('Codigos_compostos.csv')


# Seleção de alguns compostos
ind0 <- (1:25)[-c(13, 17)]
ind <- paste("C", ind0, sep="") 
a.1 <- subset(a, a$composto %in% ind)
a.1$resp <- a.1$resp * 1e3
a.1$tmp <- as.integer(substr(a.1$tempo, 1, 2))


# Fórmulas
formPCM0 <- resp ~ 0 + composto:(trt * tmp) + (0 + composto|id)
formPCM1 <- resp ~ 0 + composto:(trt + tmp) + (0 + composto|id)
formPCM2 <- resp ~ 0 + composto:trt + (0 + composto|id)

# Ajuste dos modelos par-a-par
cl <- makeCluster(3)
registerDoParallel(cl)
clusterEvalQ(cl, library(pair.mglmm))

system.time(m0.1 <- mglmmCP(formPCM0, a.1$composto, a.1)) # 4855.66 seg
system.time(m1.1 <- mglmmCP(formPCM1, a.1$composto, a.1)) # 1554.31 seg (Aviso Hess)
system.time(m2.1 <- mglmmCP(formPCM2, a.1$composto, a.1)) # 1215.37 seg

stopCluster(cl)

save(m0.1, file = 'm0.1_Ex_13_17.RData')
save(m1.1, file = 'm1.1_Ex_13_17.RData')
save(m2.1, file = 'm2.1_Ex_13_17.RData')

load('m0.1_Ex_13_17.RData')
load('m1.1_Ex_13_17.RData')
load('m2.1_Ex_13_17.RData')


## Seleção do modelo
# Teste RV
(ll0 <- m0.1$logLik - m1.1$logLik)
(ll1 <- m1.1$logLik - m2.1$logLik)

pv1 <- 1 - pchisq(2*ll0, m0.1$df-m1.1$df)
pv2 <- 1 - pchisq(2*ll1, m1.1$df-m2.1$df)

# AIC e BIC
aic1 <- -2 * (m0.1$logLik - m0.1$df)
bic1 <- -2 * m0.1$logLik + log(nrow(a.1)) * m0.1$df
aic2 <- -2 * (m1.1$logLik - m1.1$df)
bic2 <- -2 * m1.1$logLik + log(nrow(a.1)) * m1.1$df
aic3 <- -2 * (m2.1$logLik - m2.1$df)
bic3 <- -2 * m2.1$logLik + log(nrow(a.1)) * m2.1$df

# Tabela dos testes de razão de verossimlhanças
lrt.tab <- data.frame(Modelo = c('trt*tmp', 'trt+tmp', 'trt'),
                      GL = c(m0.1$df, m1.1$df, m2.1$df),
                      AIC = c(aic1, aic2, aic3),
                      BIC = c(bic1, bic2, bic3),
                      logLik = c(m0.1$logLik, m1.1$logLik, m2.1$logLik),
                      X2_obs = round(c(NA, 2*ll0, 2*ll1), 3),
                      GL_X2 = c(NA, m0.1$df-m1.1$df, m1.1$df-m2.1$df),
                      P_valor = c(NA, pv1, pv2))

lrt.tab
print(xtable(lrt.tab, dig = 3), include.rownames = F)


# Modelo selecionado: m0.1 (Interação entre tratamento e tempo)
# Tabela com as estimativas e erro padrão do modelo selecionado
est.mod <- m0.1$fixef
tes0 <- strsplit(est.mod$Parameter, ":")
tes0 <- lapply(tes0, function(x){
  if(length(x) == 3) x[2] <- paste(x[2], x[3], sep = ":")
  x[1:2]
})
tes <- do.call(rbind, tes0)
tes <- as.data.frame(tes)

beta.mod <- round(est.mod$Estimate, 3)
bmc <- as.character(beta.mod)
ep.mod <- round(est.mod$stdErr, 3)
emc <- as.character(ep.mod)
tes$value <- paste(bmc, paste('(', emc, ')', sep = ''), sep = " ")

tes.1 <- cast(tes, V1 ~ V2)
names(tes.1)[1] <- 'Composto'

cc <- as.character(tes.1$Composto)
nc <- nchar(cc)
linf <- ifelse(nc == 10, nc, nc - 1)
ind.1 <- order(as.integer(substr(cc, linf, nc)))
tes.1 <- tes.1[ind.1,]
tes.1$Composto <- substr(tes.1$Composto, 9, 11)

print(xtable(tes.1), include.rownames = F)

vc <- m0.1$VarCov[ind.1, ind.1]
m0.1$phi; m0.1$p


vc.df <- melt(diag(vc))
vc.df$Composto <- substr(rownames(vc.df), 9, 11)
vc.df$Composto <- factor(vc.df$Composto, rev(ind))
p0 <- ggplot(vc.df, aes(Composto, value))
p0 <- p0 + geom_bar() + coord_flip() + ylab(expression(delta[j]^2))
pdf("VCdelta2.pdf", w = 11)
print(p0)
dev.off()


corr <- cov2cor(vc)
rownames(corr) <- colnames(corr) <- substr(colnames(corr), 9, 11)
pdf("CorrVC.pdf")
corrplot(corr, "ellipse", "lower", tl.col = 1, tl.srt = 0, 
         order = "hclust", hclust.method = "average")
dev.off()


# Gráfico de médias
a.2 <- aggregate(resp ~ trt + composto + tempo, a.1, mean)
a.2$composto <- factor(a.2$composto, levels = ind)
## Com os nomes dos compostos:
# ind.comp <- match(a.2$composto, codcomp$composto)
# a.2$composto <- factor(codcomp$descricao[ind.comp], levels = codcomp$descricao[1:25])

p <- ggplot(a.2, aes(tempo, resp, group = trt, colour = trt))
p <- p + geom_line() + geom_point()
p <- p + facet_wrap(~composto, ncol = 5) + xlab("Tempo") + ylab("Massa")

pdf('perfil_C1_C25_Ex13_17.pdf', w = 11)
print(p)
dev.off()


# Gráfico de valores preditos
a.3 <- a.1[, c(1, 2, 4)]
a.3$pred <- exp(m0.1$frames$fr$X %*% m0.1$fixef$Estimate)
a.3 <- unique(a.3)
a.3$composto <- factor(a.3$composto, levels = ind)
## Com os nomes dos compostos:
# ind.comp <- match(a.3$composto, codcomp$composto)
# a.3$Composto <- factor(codcomp$descricao[ind.comp], levels = codcomp$descricao[1:25])

p <- ggplot(a.3, aes(tempo, pred, group = trt, colour = trt))
p <- p + geom_line()
(p <- p + facet_wrap(~composto, ncol = 5) + xlab("Tempo") + ylab("Predito"))

pdf("pred_C1_C25_Ex13_17.pdf", w = 11)
print(p)
dev.off()


# Análise multivariada
acp <- eigen(vc)
va <- acp$values
(irel <- va/sum(va))
cumsum(irel)

ve <- acp$vector[, 1:2]
rownames(ve) <- rownames(vc)
colnames(ve) <- c('PCA1', 'PCA2')

# Apenas a projeção das componentes principais
ve.df <- as.data.frame(ve)
ve.df$cmp <- substr(rownames(ve.df), 9, 11)
(p2 <- ggplot(ve.df, aes(PCA1, PCA2, label = cmp)) +
       geom_segment(aes(x = 0, y = 0, xend = PCA1, yend = PCA2)) +
       geom_text(size = 3.5) +  xlab('CP 1 (56.8%)') +
       ylab('CP 2 (18.2%)'))

pdf("ACP_C1C25.pdf")
print(p2)
dev.off()

# Versão alternativa (com os nomes dos compostos)
ve.df.1 <- merge(ve.df, codcomp, by.x = 'cmp', by.y = 'composto')

(p2.1 <- ggplot(ve.df.1, aes(PCA1, PCA2, label = descricao)) +
         geom_segment(aes(x = 0, y = 0, xend = PCA1, yend = PCA2)) +
         geom_text(size = 3.5) +  xlab('CP 1 (56.8%)') +
         ylab('CP 2 (18.2%)'))

pdf("ACP_C1C25_v2.pdf")
print(p2.1)
dev.off()
      
# Componentes principais e pontos (compostos quase indistinguíveis)
ve.obs <- m0.1$ranef %*% ve
cod <- rep(1:5, each = 8)
plot(ve.obs[,1], ve.obs[,2], pch = 20, col = cod)
textxy(ve.obs[,1], ve.obs[,2], rownames(ve.obs), cx = .6)
segments(0, 0, ve[,1], ve[,2])
textxy(ve[,1], ve[,2], substr(rownames(ve), 9, 11), cx = .6)

# Apenas pontos, visualização com ggplot2
ve.obs.1 <- as.data.frame(ve.obs)
ve.obs.1$trt <- tolower(substr(rownames(ve.obs), 1, 2))
g <- ggplot(ve.obs.1, aes(PCA1, PCA2, colour = trt))
g + geom_point()


# Análise de agrupamentos dos compostos (efeitos aleatórios)
t.ranef <- t(m0.1$ranef)
rownames(t.ranef) <- substr(colnames(m0.1$VarCov), 9, 11)
t.ranef <- t.ranef[ind.1, ]

# Método hieráriquico com distâncias médias
ha.cl <- hclust(dist(t.ranef), 'average')

pdf("clusters_C1C25.pdf", w = 10)
par(mar = c(2, 5, 3, 2))
plot(ha.cl, main = "", ylab = "Distância", xlab = "Compostos")
abline(h = 5, lty = 2)
dev.off()


# Versão alternativa (com os nomes dos compostos)
rownames(t.ranef) <- codcomp$descricao[ind0]

# Método hieráriquico com distâncias médias
ha.cl <- hclust(dist(t.ranef), 'average')

pdf("clusters_C1C25_v2.pdf", w = 10)
par(mar = c(2, 5, 3, 2))
plot(ha.cl, main = "", ylab = "Distância", xlab = "Compostos", cex = .9)
abline(h = 5, lty = 2)
dev.off()


# Comparações múltiplas - Para cada composto, tem-se comparações entre os 
# tratamentos fixando-se um tempo de cada vez. Da mesma forma, tem-se comparações 
# entre tempos para um tratamento fixado. Assim, tem-se uma tabela com os trt's 
# nas linhas e tempos nas colunas, com 40 comparações para cada composto.

cntr.fun <- function(ac, X.mod, id.level, mod) {
  ind.c1 <- ac$composto == id.level
  ac.c1 <- ac[ind.c1,]
  ac.c1$trt.tmp <- with(ac.c1, paste(substr(trt, 1, 3), tmp, sep = ":"))
  X.c1 <- X.mod[ind.c1,]
  indice <- seq_len(nrow(ac.c1))

  trt.l <- split(indice, ac.c1$trt)
  trt.comb <- lapply(trt.l, combn, 2)
  trt.bind <- do.call(cbind, trt.comb)

  tmp.l <- split(indice, ac.c1$tmp)
  tmp.comb <- lapply(tmp.l, combn, 2)
  tmp.bind <- do.call(cbind, tmp.comb)

  p.c1 <- cbind(trt.bind, tmp.bind)

  r.c1 <- lapply(1:ncol(p.c1), function(j) {
    p.c1.j <- p.c1[, j]
    x.j <- matrix(X.c1[p.c1.j[1],] - X.c1[p.c1.j[2],], nrow = 1)

    ctr1 <- ac.c1$trt.tmp[p.c1.j[1]]
    ctr2 <- ac.c1$trt.tmp[p.c1.j[2]]
    rownames(x.j) <- paste(ctr1, ctr2, sep = " - ")

    x.j
  })

  r.c1 <- do.call(rbind, r.c1)

  est <- r.c1 %*% as.vector(mod$fixef[,2])
  vc <- mod$VC.fixef
  se <- sqrt(apply(r.c1, 1, function(x) t(x) %*% as.matrix(vc) %*% x))

  p.values <- 2 * pnorm(abs(est/se), lower = F)
  p.val.adj <- p.adjust(p.values, 'holm')

  cbind(Estimate = c(est), stdError = se, pValueAdj = p.val.adj)  
}


ac <- unique(a.1[, c('composto', 'trt', 'tmp')])
rownames(ac) <- NULL
ac <- drop.levels(ac)

X.mod <- model.matrix(~ 0 + composto:(trt * tmp), ac)

comps <- ind
cntr.0 <- lapply(comps, cntr.fun, ac = ac, X.mod = X.mod, mod = m0.1)

cntr.1 <- lapply(ind.1, function(i) {
  x <- as.data.frame(cntr.0[[i]])
  #x[order(x$pValueAdj),]
})


## IMPORTANTE!!
## Organização dos contrastes em algo mais palatável para colocar na dissertação

cntr.fun.2 <- function(x) {
  x.1 <- strsplit(rownames(x), ' - ')
  x.1 <- do.call(rbind, x.1)
  x.1 <- as.data.frame(x.1)
  x.1$pvalor <- x$pValueAdj

  split.cntr <- with(x.1, substr(V1, 1, 3) == substr(V2, 1, 3))
  x.l <- split(x.1, split.cntr)

  # Tempo fixado
  x.tmp <- x.l$'FALSE'

  x.tmp$tmp <- substr(x.tmp$V1, 5, 6)
  x.tmp$V1 <- substr(x.tmp$V1, 1, 3)
  x.tmp$V2 <- substr(x.tmp$V2, 1, 3)

  ord.trt <- apply(x.tmp[, c('V1', 'V2')], 1, sort)
  x.tmp$V1 <- ord.trt[1,]
  x.tmp$V2 <- ord.trt[2,]

  x.tmp.1 <- cast(x.tmp, V1 + V2 ~ tmp, value = 'pvalor')


  # tratamento fixado
  x.trt <- x.l$'TRUE'

  x.trt$tmp1 <- substr(x.trt$V1, 5, 6)
  x.trt$tmp2 <- substr(x.trt$V2, 5, 6)
  x.trt$V1 <- substr(x.trt$V1, 1, 3)

  ord.tmp <- apply(x.trt[, c('tmp1', 'tmp2')], 1, sort)
  x.trt$tmp1 <- ord.tmp[1,]
  x.trt$tmp2 <- ord.tmp[2,]

  x.trt.1 <- cast(x.trt, tmp1 + tmp2 ~ V1, value = 'pvalor')

  # lista com as tabelas por tempo e tratamento
  list(tmp = x.tmp.1, trt = x.trt.1)
}

cntr.2 <- lapply(cntr.1, cntr.fun.2)

# Tabela cruzada de calores preditos para cada composto
a.3.1 <- a.3
a.3.1$pred <- round(a.3.1$pred, 2)
a.3.1l <- split(a.3.1, a.3.1$composto)


cntr.3 <- lapply(cntr.2, function(x) {
  y <- x$tmp; z <- x$trt
  
  y.1 <- y[, -(1:2)]
  z.1 <- z[, -(1:2)]

  y.lab <- paste(y$V1, y$V2, sep = " - ")
  z.lab <- paste(z$tmp1, z$tmp2, sep = " - ")

  y.2 <- apply(y.1, 2, function(x){
    x.1 <- x < 0.05
    names(x.1) <- y.lab
    multcomp:::insert_absorb(x.1, decreasing = FALSE)$Letters
  })

  z.2 <- apply(z.1, 2, function(x){
    x.1 <- x < 0.05
    names(x.1) <- z.lab
    multcomp:::insert_absorb(x.1, decreasing = FALSE)$Letters
  })

  y.2 <- as.data.frame(y.2)
  y.2 <- data.frame(trt = rownames(y.2), y.2)
  y.2 <- melt(y.2, id = 'trt', variable_name = "tempo")
  levels(y.2$tempo) <- paste(names(y.1), "h", sep = "")

  z.2 <- as.data.frame(z.2)
  names(z.2) <- names(z.1)
  z.2 <- data.frame(tempo = paste(rownames(z.2), "h", sep = ""), z.2)
  z.2 <- melt(z.2, id = 'tempo', variable_name = 'trt')

  x.2 <- merge(z.2, y.2, by = c('trt', 'tempo'))
  x.2$cod <- paste(x.2$value.x, x.2$value.y, sep = ", ")
  levels(x.2$trt) <- levels(a.3.1$trt)
  x.2[, c('trt', 'tempo', 'cod')]
})


cntr.4 <- lapply(1:23, function(i){
  x <- merge(a.3.1l[[i]], cntr.3[[i]], by = c('trt', 'tempo'))
  x$value <- with(x, paste("$", pred, "^{", cod, "}$", sep = ""))
  cast(x, composto + trt ~ tempo, value = 'value')
})


cntr.5 <- do.call(rbind, cntr.4)

print(xtable(cntr.5), include.rownames = F, 
      sanitize.text.function = function(x){x}, 
      tabular.environment = 'longtable', floating = F)