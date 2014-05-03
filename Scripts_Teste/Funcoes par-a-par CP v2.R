glmmMultiCP <- function(formula, id, data, cl = NULL) {
  lev <- paste(unique(id))
  tes <- combn(lev, m = 2)

  glmm.fit <- function(x, data, formula) {
    ind <- id %in% x
    dat0 <- subset(data, ind)
    cpglmm(formula, data = dat0)
  }

  res <- foreach(i = 1:ncol(tes)) %dopar% {
           glmm.fit(tes[,i], data, formula)
         }

  return(res)
}



format0CP <- function(mod) {
  betas <- fixef(mod)
  df.b <- melt(betas)
  df.b$Parametro <- rownames(df.b)
  rownames(df.b) <- 1:nrow(df.b)

  df.b <- df.b[,c(2,1)]

  S <- VarCorr(mod)
  uS <- unlist(S)
  triS <- unique(uS)

  nameS <- unlist(attr(S[[1]], "dimnames"))
  nameS <- combn(sort(nameS), 2)
  nameS <- unique(t(nameS))
  nameS <- apply(nameS, 1, paste, collapse = ":")

  df.s <- data.frame(Parametro = nameS, value = triS)

  phi <- as.numeric(attr(S, 'sc'))^2
  phi <- data.frame(Parametro = 'phi', 
                    value = ifelse(is.na(phi), 1, phi))

  p <- data.frame(Parametro = 'p', value = mod@p)

  rbind(df.b, df.s, phi, p)
}



estAveFun <- function(mod.list, jkj) {
  tes <- lapply(mod.list, format0CP)
  tes.1 <- do.call(rbind, tes)

  df.l <- sort(unique(tes.1$Parametro))

  mat <- lapply(df.l, '==', tes.1$Parametro)
  mat <- do.call(rbind, mat)

  mat <- sweep(mat, 1, rowSums(mat), "/")

  est <- mat %*% tes.1$value
  ep0 <- mat %*% jkj %*% t(mat)
  ep1 <- sqrt(diag(ep0))

  data.frame(Parameter = df.l,
             Estimate = est,
             stdError = ep1)
}



format1CP <- function(df.m, formula, data) {
  m1 <- lmer(formula, data = data, doFit = F)

  matn <- names(m1$fr$fixef)
  ind <- match(matn, df.m$Parameter)
  est <- df.m[ind,]

  ind.phi <- match('phi', df.m$Parameter)
  phi <- df.m[ind.phi,-1]

  ind.p <- match('p', df.m$Parameter)
  p <- df.m[ind.p,-1]

  SigHat <- m1$FL$trms[[1]]$ST
  indS <- !seq(1, nrow(df.m)) %in% c(ind, ind.phi, ind.p)
  SigHat[lower.tri(SigHat, diag = T)] <- df.m[indS,2]
  SigHat <- SigHat + t(lower.tri(SigHat)*SigHat)
  SigHat <- nearPD(SigHat)$mat

  list(m1 = m1, est = est, phi = phi,
       p = p, SigHat = SigHat)
}



f0.ep <- function (data, formula) {
    formula0 <- eval(formula)
    fix.form <- lme4:::nobars(formula)
    fr <- model.frame(fix.form, data = data)

    Fr <- list()
    Fr$Y <- model.response(fr)
    Fr$X <- model.matrix(attr(fr, "terms"), data = fr)
    Fr$wts <- rep(1, nrow(Fr$X))
    Fr$off <- rep(0, nrow(Fr$X))
    Fr$mf <- data
    Fr$fixef <- rep(0, ncol(Fr$X))

    fl <- lme4:::lmerFactorList(formula0, Fr, 0L, 0L)

    x <- Fr$X; y <- Fr$Y
    z <- t(as.matrix(fl$trms[[1]]$Zt))

    list(y = y, x = x, z = z)
}



mll2 <- function (th, r) {
    nx <- ncol(r$x)
    nz <- ncol(r$z)

    beta <- th[1:nx]
    L <- matrix(0, nz, nz)
    nLlt <- length(L[lower.tri(L, diag = T)])
    L[lower.tri(L, diag = T)] <- th[(nx + 1):(nx + nLlt)]
    phi <- th[nx + nLlt + 1]
    p <- th[nx + nLlt + 2]

    pr.grid <- createProductRuleGrid("GQN", 2, 11, sym = FALSE)

    res <- .Call("mll", r, beta, L, phi, p, pr.grid$nodes, 
                 pr.grid$weights, PACKAGE = "pair.mglmm")

    as.vector(res)
}



llik.ep <- function (mod, formula) {
    beta <- fixef(mod)
    phi <- mod@phi
    p <- mod@p
    L <- t(.Call("mer_ST_chol", mod)[[1]])
    Llt <- L[lower.tri(L, diag = T)]
    th <- c(beta, Llt, phi, p)

    mf.ind <- split(mod@frame, mod@flist[[1]])
    mf.r <- lapply(mf.ind, f0.ep, formula = formula)

    res <- foreach(i = 1:length(mf.r)) %dopar% {
        g <- grad(mll2, th, r = mf.r[[i]])
        h <- hessian(mll2, th, r = mf.r[[i]])
        list(g = g, h = h)
    }

    return(res)
}



rcov2 <- function (lmods, formula) {
    jkl <- lapply(lmods, llik.ep, formula = formula)
    l.jkl <- length(jkl[[1]])

    gs0 <- lapply(1:l.jkl, function(i, x) {
        p0 <- lapply(x, "[[", i)
        p1 <- lapply(p0, "[[", "g")
        p2 <- lapply(1:length(p1), 
                     function(i) lapply(p1, outer, p1[[i]]))
        p3 <- lapply(p2, do.call, what = cbind)
        do.call(rbind, p3)
    }, x = jkl)

    K <- (1/length(gs0)) * Reduce("+", gs0)

    hs0 <- lapply(jkl, function(x) {
        p0 <- lapply(x, "[[", "h")
        (-1/length(x)) * Reduce("+", p0)
    })

    J <- bdiag(hs0)
    Jinv <- solve(J)
    (Jinv %*% K %*% Jinv)/l.jkl
}



llik.fim <- function (mod, formula, beta, S, phi, p, B = 10000, cl = NULL) {
    mf.ind <- split(mod$fr$mf, mod$FL$fl[[1]])
    mf.r <- lapply(mf.ind, f0.ep, formula = formula)

    q <- ncol(S)
    U <- mvrnorm(B, rep(0, q), S)

    ll <- foreach(i = 1:length(mf.r), .combine = c) %dopar% {
        r <- mf.r[[i]]
        .Call("mllCPMC", r$y, beta, r$x, r$z, U, phi, p, PACKAGE = "pair.mglmm")
    }

    sum(log(ll))
}



mglmmCP <- function(formula, id, data, cl = NULL) {
    # Ajuste dos modelos par-a-par
    m0 <- glmmMultiCP(formula, id, data)
    EP <- rcov2(m0, formula)

    # Estimativas dos parâmetros e formatação
    prm0 <- estAveFun(m0, EP)
    prm1 <- format1CP(prm0, formula, data)

    # Log-verossimilhança
    LL <- llik.fim(prm1$m1, formula, prm1$est[,2],
                   prm1$SigHat, prm1$phi[,2], prm1$p[,2])

    # Resultado final
    list(fixef = prm$est, VarCov = prm1$SigHat, phi = prm1$phi,
         p = prm1$p, logLik = LL)
}