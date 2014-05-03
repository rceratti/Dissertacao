(m1 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject),
             data = sleepstudy))

##
Lambdat <- t(with(expand(m1), S %*% T))
Zt <- m1@Zt

P <- expand(m1)$P

LtZt <- Lambdat %*% Zt

LtZtZL <- LtZt %*% t(LtZt)
q <- ncol(LtZtZL)

L_z <- t(chol(P %*% (LtZtZL + diag(1, q)) %*% t(P)))

X <- m1@X
L_xz <- t(solve(L_z) %*% P %*% LtZt %*% X)

L_x <- chol(t(X) %*% X - L_xz %*% t(L_xz))

sqrt(diag(sigma(m1)^2 * chol2inv(L_x)))
vcov(m1)



(m1 <- glmer(size ~ period + (-1 + period | herd), cbpp, poisson))

##
q <- length(unique(unlist(m1@flist)))

Lambdat <- .Call("mer_ST_chol", m1)[[1]]
Lambdat <- lapply(1:q, function(i) Lambdat)
Lambdat <- bdiag(Lambdat)

Zt <- t(model.matrix(~-1 + period:herd, cbpp))

W <- diag(fitted(m1))  # Poisson Composta: W <- diag(fitted(m1)^(2-m1@p))

LtZt <- Lambdat %*% Zt

LtZtWZL <- LtZt %*% W %*% t(LtZt)
q <- ncol(LtZtWZL)
 
L_z <- t(chol(LtZtWZL + diag(1, q)))

X <- m1@X
L_xz <- t(solve(L_z) %*% LtZt %*% W %*% X)

L_x <- chol(t(X) %*% W %*% X - L_xz %*% t(L_xz))

sigma(m1)^2 * chol2inv(L_x)
vcov(m1)
