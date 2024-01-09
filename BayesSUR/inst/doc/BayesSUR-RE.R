## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = FALSE)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----eval=TRUE----------------------------------------------------------------
library("BayesSUR")

## -----------------------------------------------------------------------------
#  library("gRbase")
#  sim.ssur <- function(n, s, p, t0 = 0, seed = 123, mv = TRUE,
#                       t.df = Inf, random.intercept = 0, intercept = TRUE) {
#    # set seed to fix coefficients
#    set.seed(7193)
#    sd_b <- 1
#    mu_b <- 1
#    b <- matrix(rnorm((p + ifelse(t0 == 0, 1, 0)) * s, mu_b, sd_b), p + ifelse(t0 == 0, 1, 0), s)
#  
#    # design groups and pathways of Gamma matrix
#    gamma <- matrix(FALSE, p + ifelse(t0 == 0, 1, 0), s)
#    if (t0 == 0) gamma[1, ] <- TRUE
#    gamma[2:6 - ifelse(t0 == 0, 0, 1), 1:5] <- TRUE
#    gamma[11:21 - ifelse(t0 == 0, 0, 1), 6:12] <- TRUE
#    gamma[31:51 - ifelse(t0 == 0, 0, 1), 1:5] <- TRUE
#    gamma[31:51 - ifelse(t0 == 0, 0, 1), 13:15] <- TRUE
#    gamma[52:61 - ifelse(t0 == 0, 0, 1), 1:12] <- TRUE
#    gamma[71:91 - ifelse(t0 == 0, 0, 1), 6:15] <- TRUE
#    gamma[111:121 - ifelse(t0 == 0, 0, 1), 1:15] <- TRUE
#    gamma[122 - ifelse(t0 == 0, 0, 1), 16:18] <- TRUE
#    gamma[123 - ifelse(t0 == 0, 0, 1), 19] <- TRUE
#    gamma[124 - ifelse(t0 == 0, 0, 1), 20] <- TRUE
#  
#    G_kron <- matrix(0, s * p, s * p)
#    G_m <- bdiag(matrix(1, ncol = 5, nrow = 5),
#                 matrix(1, ncol = 7, nrow = 7),
#                 matrix(1, ncol = 8, nrow = 8))
#    G_p <- bdiag(matrix(1, ncol = 5, nrow = 5), diag(3),
#                 matrix(1, ncol = 11, nrow = 11), diag(9),
#                 matrix(1, ncol = 21, nrow = 21),
#                 matrix(1, ncol = 10, nrow = 10), diag(9),
#                 matrix(1, ncol = 21, nrow = 21), diag(19),
#                 matrix(1, ncol = 11, nrow = 11), diag(181))
#    G_kron <- kronecker(G_m, G_p)
#  
#    combn11 <- combn(rep((1:5 - 1) * p, each = length(1:5)) +
#                       rep(1:5, times = length(1:5)), 2)
#    combn12 <- combn(rep((1:5 - 1) * p, each = length(30:60)) +
#                       rep(30:60, times = length(1:5)), 2)
#    combn13 <- combn(rep((1:5 - 1) * p, each = length(110:120)) +
#                       rep(110:120, times = length(1:5)), 2)
#    combn21 <- combn(rep((6:12 - 1) * p, each = length(10:20)) +
#                       rep(10:20, times = length(6:12)), 2)
#    combn22 <- combn(rep((6:12 - 1) * p, each = length(51:60)) +
#                       rep(51:60, times = length(6:12)), 2)
#    combn23 <- combn(rep((6:12 - 1) * p, each = length(70:90)) +
#                       rep(70:90, times = length(6:12)), 2)
#    combn24 <- combn(rep((6:12 - 1) * p, each = length(110:120)) +
#                       rep(110:120, times = length(6:12)), 2)
#    combn31 <- combn(rep((13:15 - 1) * p, each = length(30:50)) +
#                       rep(30:50, times = length(13:15)), 2)
#    combn32 <- combn(rep((13:15 - 1) * p, each = length(70:90)) +
#                       rep(70:90, times = length(13:15)), 2)
#    combn33 <- combn(rep((13:15 - 1) * p, each = length(110:120)) +
#                       rep(110:120, times = length(13:15)), 2)
#    combn4 <- combn(rep((16:18 - 1) * p, each = length(121)) +
#                      rep(121, times = length(16:18)), 2)
#    combn5 <- matrix(rep((19 - 1) * p, each = length(122)) +
#                       rep(122, times = length(19)), nrow = 1, ncol = 2)
#    combn6 <- matrix(rep((20 - 1) * p, each = length(123)) +
#                       rep(123, times = length(20)), nrow = 1, ncol = 2)
#  
#    combnAll <- rbind(t(combn11), t(combn12), t(combn13),
#                      t(combn21), t(combn22), t(combn23), t(combn24),
#                      t(combn31), t(combn32), t(combn33),
#                      t(combn4), combn5, combn6)
#  
#    set.seed(seed + 7284)
#    sd_x <- 1
#    x <- matrix(rnorm(n * p, 0, sd_x), n, p)
#  
#    if (t0 == 0 & intercept) x <- cbind(rep(1, n), x)
#    if (!intercept) {
#      gamma <- gamma[-1, ]
#      b <- b[-1, ]
#    }
#    xb <- matrix(NA, n, s)
#    if (mv) {
#      for (i in 1:s) {
#        if (sum(gamma[, i]) >= 1) {
#          if (sum(gamma[, i]) == 1) {
#            xb[, i] <- x[, gamma[, i]] * b[gamma[, i], i]
#          } else {
#            xb[, i] <- x[, gamma[, i]] %*% b[gamma[, i], i]
#          }
#        } else {
#          xb[, i] <- sapply(1:s, function(i) rep(1, n) * b[1, i])
#        }
#      }
#    } else {
#      if (sum(gamma) >= 1) {
#        xb <- x[, gamma] %*% b[gamma, ]
#      } else {
#        xb <- sapply(1:s, function(i) rep(1, n) * b[1, i])
#      }
#    }
#  
#    corr_param <- 0.9
#    M <- matrix(corr_param, s, s)
#    diag(M) <- rep(1, s)
#  
#    ## wanna make it decomposable
#    Prime <- list(c(1:(s * .4), (s * .8):s),
#                  c((s * .4):(s * .6)),
#                  c((s * .65):(s * .75)),
#                  c((s * .8):s))
#    G <- matrix(0, s, s)
#    for (i in 1:length(Prime)) {
#      G[Prime[[i]], Prime[[i]]] <- 1
#    }
#  
#    # check
#    dimnames(G) <- list(1:s, 1:s)
#    length(gRbase::mcsMAT(G - diag(s))) > 0
#  
#    var <- solve(BDgraph::rgwish(n = 1, adj = G, b = 3, D = M))
#  
#    # change seeds to add randomness on error
#    set.seed(seed + 8493)
#    sd_err <- 0.5
#    if (is.infinite(t.df)) {
#      err <- matrix(rnorm(n * s, 0, sd_err), n, s) %*% chol(as.matrix(var))
#    } else {
#      err <- matrix(rt(n * s, t.df), n, s) %*% chol(as.matrix(var))
#    }
#  
#    if (t0 == 0) {
#      b.re <- NA
#      z <- NA
#      y <- xb + err
#      if (random.intercept != 0) {
#        y <- y + matrix(rnorm(n * s, 0, sqrt(random.intercept)), n, s)
#      }
#  
#      z <- sample(1:4, n, replace = T, prob = rep(1 / 4, 4))
#  
#      return(list(y = y, x = x, b = b, gamma = gamma, z = model.matrix(~ factor(z) + 0)[, ],
#                  b.re = b.re, Gy = G, mrfG = combnAll))
#    } else {
#      # add random effects
#      z <- t(rmultinom(n, size = 1, prob = c(.1, .2, .3, .4)))
#      z <- sample(1:t0, n, replace = T, prob = rep(1 / t0, t0))
#      set.seed(1683)
#      b.re <- rnorm(t0, 0, 2)
#      y <- matrix(b.re[z], nrow = n, ncol = s) + xb + err
#  
#      return(list(
#        y = y, x = x, b = b, gamma = gamma, z = model.matrix(~ factor(z) + 0)[, ],
#        b.re = b.re, Gy = G, mrfG = combnAll
#      ))
#    }
#  }

## -----------------------------------------------------------------------------
#  library("Matrix")
#  n <- 250
#  s <- 20
#  p <- 300
#  sim1 <- sim.ssur(n, s, p, seed = 1)

## -----------------------------------------------------------------------------
#  t0 <- 4
#  sim2 <- sim.ssur(n, s, p, t0, seed = 1) # learning data
#  sim2.val <- sim.ssur(n, s, p, t0, seed=101) # validation data

## -----------------------------------------------------------------------------
#  hyperpar <- list(mrf_d = -2, mrf_e = 1.6, a_w0 = 100, b_w0 = 500, a_w = 15, b_w = 60)
#  set.seed(1038)
#  fit2 <- BayesSUR(
#    data = cbind(sim2$y, sim2$z, sim2$x),
#    Y = 1:s,
#    X_0 = s + 1:t0,
#    X = s + t0 + 1:p,
#    outFilePath = "sim2_mrf_re",
#    hyperpar = hyperpar,
#    gammaInit = "0",
#    betaPrior = "reGroup",
#    nIter = 300, burnin = 100,
#    covariancePrior = "HIW",
#    standardize = F,
#    standardize.response = F,
#    gammaPrior = "MRF",
#    mrfG = sim2$mrfG,
#    output_CPO = T
#  )

## -----------------------------------------------------------------------------
#  summary(fit2)

## -----------------------------------------------------------------------------
#  # compute accuracy, sensitivity, specificity of variable selection
#  gamma <- getEstimator(fit2)
#  (accuracy <- sum(data.matrix(gamma > 0.5) == sim2$gamma) / prod(dim(gamma)))

## -----------------------------------------------------------------------------
#  (sensitivity <- sum((data.matrix(gamma > 0.5) == 1) & (sim2$gamma == 1)) / sum(sim2$gamma == 1))

## -----------------------------------------------------------------------------
#  (specificity <- sum((data.matrix(gamma > 0.5) == 0) & (sim2$gamma == 0)) / sum(sim2$gamma == 0))

## -----------------------------------------------------------------------------
#  # compute RMSE and RMSPE for prediction performance
#  beta <- getEstimator(fit2, estimator = "beta", Pmax = .5, beta.type = "conditional")
#  (RMSE <- sqrt(sum((sim2$y - cbind(sim2$z, sim2$x) %*% beta)^2) / prod(dim(sim2$y))))

## -----------------------------------------------------------------------------
#  (RMSPE <- sqrt(sum((sim2.val$y - cbind(sim2.val$z, sim2.val$x) %*% beta)^2) / prod(dim(sim2.val$y))))

## -----------------------------------------------------------------------------
#  # compute bias of beta estimates
#  b <- sim2$b
#  b[sim2$gamma == 0] <- 0
#  (beta.l2 <- sqrt(sum((beta[-c(1:4), ] - b)^2) / prod(dim(b))))

## -----------------------------------------------------------------------------
#  g.re <- getEstimator(fit2, estimator = "Gy")
#  (g.accuracy <- sum((g.re > 0.5) == sim2$Gy) / prod(dim(g.re)))

## -----------------------------------------------------------------------------
#  (g.sensitivity <- sum(((g.re > 0.5) == sim2$Gy)[sim2$Gy == 1]) / sum(sim2$Gy == 1))

## -----------------------------------------------------------------------------
#  (g.specificity <- sum(((g.re > 0.5) == sim2$Gy)[sim2$Gy == 0]) / sum(sim2$Gy == 0))

