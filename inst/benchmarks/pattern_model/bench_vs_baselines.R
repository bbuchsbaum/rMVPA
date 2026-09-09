# Head-to-head predictive benchmark for pattern_model (the comparison the
# Phase 3 plan left unrun).
#
# Question: on identical outer splits, how does the whole-brain pattern model
# decode against the estimators a practitioner would actually reach for?
#
#   Classification (Part A):
#     * pattern_model, rank tuned, no penalty
#     * pattern_model, rank + sparse penalty tuned
#     * pattern_model, rank + sparse + signed_smooth tuned (grid graph)
#     * spacenet_tvl1 (TV-L1 logistic, alpha tuned by its own inner CV)
#     * whole-brain shrinkage LDA
#     * searchlight shrinkage LDA, honest: the sphere is chosen by inner
#       blocked CV on the training rows only, then evaluated on the outer fold
#     * searchlight oracle: the sphere with the best OUTER-fold accuracy.
#       This peeks at the test data; it is an upper bound on any sphere
#       selection rule, not an estimator.
#
#   Continuous decoding (Part B):
#     * pattern_model, rank tuned, decode direction
#     * PLS regression brain -> features with ncomp chosen by inner blocked CV
#       (the same kernel feature_rsa_model(method = "pls") uses)
#
# The simulation plants four kinds of structure on a 3D grid: a smooth
# informative blob (component 1), a fine-grained sign-flipping informative
# region (component 2), a region carrying only localized low-rank nuisance
# (which the global diag+lowrank noise model cannot fully represent -- the
# planned failure mode), and a redundant copy of component 1. SNR is set low
# enough that no method saturates.
#
# Run:  Rscript inst/benchmarks/pattern_model/bench_vs_baselines.R [quick]
# Record results in adocs/pattern-model-benchmarks.md.

if (!isNamespaceLoaded("rMVPA")) suppressPackageStartupMessages(library(rMVPA))
suppressPackageStartupMessages(library(Matrix))
set.seed(20260909)
quick <- "quick" %in% commandArgs(trailingOnly = TRUE)

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---------------------------------------------------------------------------
# simulation on a 3D grid
# ---------------------------------------------------------------------------

make_sim <- function(n, dims, K = 3, snr = 0.20, q = 5, seed = 1) {
  set.seed(seed)
  p <- prod(dims)
  coords <- arrayInd(seq_len(p), dims)

  half_x <- floor(dims[1] / 2); half_y <- floor(dims[2] / 2)
  quadrant <- function(xr, yr) {
    which(coords[, 1] >= xr[1] & coords[, 1] <= xr[2] &
          coords[, 2] >= yr[1] & coords[, 2] <= yr[2])
  }
  reg <- list(blob  = quadrant(c(1, half_x), c(1, half_y)),
              flip  = quadrant(c(half_x + 1, dims[1]), c(1, half_y)),
              nuis  = quadrant(c(1, half_x), c(half_y + 1, dims[2])),
              redun = quadrant(c(half_x + 1, dims[1]), c(half_y + 1, dims[2])))

  A <- matrix(0, p, 2)
  centre1 <- colMeans(coords[reg$blob, , drop = FALSE])
  d2 <- rowSums(sweep(coords[reg$blob, , drop = FALSE], 2L, centre1, "-")^2)
  A[reg$blob, 1] <- exp(-d2 / (2 * (max(dims) / 5)^2))          # smooth blob
  cf <- coords[reg$flip, , drop = FALSE]
  A[reg$flip, 2] <- ifelse((cf[, 1] + cf[, 2] + cf[, 3]) %% 2 == 0, 1, -1) *
    (1 + abs(rnorm(length(reg$flip), sd = 0.2)))                # sign-flipping code
  A[reg$redun, 1] <- rnorm(length(reg$redun), sd = 0.5)         # redundant copy of comp 1
  A <- A * snr

  U <- matrix(0, p, 3)                                          # localized nuisance
  for (j in 1:3) U[reg$nuis, j] <- rnorm(length(reg$nuis)) * 1.2

  # class targets (Part A) and continuous feature targets (Part B) share A
  y <- factor(rep(letters[seq_len(K)], length.out = n))
  Yc <- scale(stats::model.matrix(~ y - 1), scale = FALSE)
  S <- crossprod(Yc) / n
  eg <- eigen(S, symmetric = TRUE)
  keep <- eg$values > 1e-8
  Wy <- eg$vectors[, keep, drop = FALSE] %*% diag(1 / sqrt(eg$values[keep]), nrow = sum(keep))
  Cc <- qr.Q(qr(matrix(rnorm(sum(keep) * 2), sum(keep), 2)))
  T_class <- (Yc %*% Wy) %*% Cc

  Yq <- matrix(rnorm(n * q), n, q, dimnames = list(NULL, paste0("f", seq_len(q))))
  Cq <- qr.Q(qr(matrix(rnorm(q * 2), q, 2)))
  T_cont <- scale(Yq, scale = FALSE) %*% Cq

  noise <- function(m) matrix(rnorm(m * p), m, p) +
    matrix(rnorm(m * ncol(U)), m, ncol(U)) %*% t(U)
  X_class <- T_class %*% t(A) + noise(n)
  X_cont  <- T_cont %*% t(A) + noise(n)

  blocks <- as.integer(cut(seq_len(n), 3, labels = 1:3))
  mask <- neuroim2::NeuroVol(array(1, dims), neuroim2::NeuroSpace(dims, c(1, 1, 1)))
  list(X_class = X_class, y = y, X_cont = X_cont, Yq = Yq,
       blocks = blocks, coords = coords, dims = dims, mask = mask, p = p, reg = reg)
}

# ---------------------------------------------------------------------------
# baselines
# ---------------------------------------------------------------------------

shrinkage_lda_fit <- function(Xtr, ytr, lambda = 0.1) {
  lev <- levels(droplevels(ytr))
  mu <- t(vapply(lev, function(l) colMeans(Xtr[ytr == l, , drop = FALSE]), numeric(ncol(Xtr))))
  Xc <- Xtr - mu[as.integer(droplevels(ytr)), , drop = FALSE]
  S <- crossprod(Xc) / max(nrow(Xtr) - length(lev), 1)
  S <- (1 - lambda) * S + lambda * mean(diag(S)) * diag(ncol(Xtr))
  W <- solve(S, t(mu))
  list(W = W, b = 0.5 * colSums(t(mu) * W), lev = lev)
}
shrinkage_lda_predict <- function(fit, Xte) {
  sc <- Xte %*% fit$W - matrix(fit$b, nrow(Xte), length(fit$lev), byrow = TRUE)
  factor(fit$lev[max.col(sc, ties.method = "first")], levels = fit$lev)
}

# searchlight spheres on the grid (radius in voxels)
make_spheres <- function(coords, radius = 2) {
  offs <- as.matrix(expand.grid(dx = -radius:radius, dy = -radius:radius, dz = -radius:radius))
  offs <- offs[rowSums(offs^2) <= radius^2, , drop = FALSE]
  dims <- apply(coords, 2L, max)
  lin <- function(cc) {
    ok <- cc[, 1] >= 1 & cc[, 1] <= dims[1] & cc[, 2] >= 1 & cc[, 2] <= dims[2] &
          cc[, 3] >= 1 & cc[, 3] <= dims[3]
    (cc[ok, 3] - 1) * dims[1] * dims[2] + (cc[ok, 2] - 1) * dims[1] + cc[ok, 1]
  }
  lapply(seq_len(nrow(coords)), function(v) {
    cc <- sweep(offs, 2L, as.numeric(coords[v, ]), "+")
    lin(cc)
  })
}

sphere_accs <- function(X, ytr, tr, te, spheres) {
  vapply(spheres, function(vox) {
    f <- shrinkage_lda_fit(X[tr, vox, drop = FALSE], ytr[tr])
    mean(shrinkage_lda_predict(f, X[te, vox, drop = FALSE]) == ytr[te])
  }, numeric(1))
}

# ---------------------------------------------------------------------------
# Part A: classification on identical outer splits
# ---------------------------------------------------------------------------

n <- if (quick) 120 else 240
dims <- if (quick) c(10, 10, 4) else c(14, 14, 8)
snr_grid <- if (quick) 0.35 else c(0.20, 0.10)

setup_sim <- function(snr) {
  sim <<- make_sim(n, dims, snr = snr, seed = 42)
  folds <<- lapply(1:3, function(b) list(train = which(sim$blocks != b),
                                         test = which(sim$blocks == b)))
  graph <<- spatial_graph(list(A = {
    # 6-neighbour grid adjacency
    i <- integer(0); j <- integer(0)
    for (ax in 1:3) {
      step <- c(1, sim$dims[1], sim$dims[1] * sim$dims[2])[ax]
      ok <- sim$coords[, ax] < sim$dims[ax]
      i <- c(i, which(ok)); j <- c(j, which(ok) + step)
    }
    Matrix::sparseMatrix(i = c(i, j), j = c(j, i), x = 1, dims = c(sim$p, sim$p))
  }))
}
setup_sim(snr_grid[1])
cat(sprintf("## Head-to-head classification (p = %d, n = %d, K = 3, 3 blocked outer folds)\n",
            sim$p, n))

pattern_run <- function(penalty = NULL, use_graph = FALSE) {
  ctrl <- pattern_control(max_rank = 3)
  pred <- factor(rep(NA_character_, n), levels = levels(sim$y))
  secs <- system.time(for (f in folds) {
    tr <- f$train
    sel <- rMVPA:::.pattern_select_config(sim$X_class[tr, ], sim$y[tr], sim$blocks[tr],
                                          ctrl, penalty = penalty,
                                          graph = if (use_graph) graph else NULL)
    fit <- rMVPA:::.pattern_fit(sim$X_class[tr, ], sim$y[tr], rank = sel$rank, control = ctrl,
                                penalty = sel$penalty, graph = if (use_graph) graph else NULL,
                                cap_rank = TRUE)
    pred[f$test] <- predict(fit, sim$X_class[f$test, ], type = "class")
  })[["elapsed"]]
  list(acc = mean(pred == sim$y), secs = secs)
}

spacenet_run <- function() {
  m <- load_model("spacenet_tvl1")
  param <- m$grid(sim$X_class, sim$y, 1)
  pred <- factor(rep(NA_character_, n), levels = levels(sim$y))
  secs <- system.time(for (f in folds) {
    fit <- m$fit(sim$X_class[f$train, ], sim$y[f$train], wts = NULL, param = param,
                 lev = levels(sim$y), last = TRUE, weights = NULL, classProbs = TRUE,
                 feature_ids = seq_len(sim$p), spatial_mask = sim$mask)
    pred[f$test] <- m$predict(fit, sim$X_class[f$test, ])
  })[["elapsed"]]
  list(acc = mean(pred == sim$y), secs = secs)
}

lda_run <- function() {
  pred <- factor(rep(NA_character_, n), levels = levels(sim$y))
  secs <- system.time(for (f in folds) {
    fit <- shrinkage_lda_fit(sim$X_class[f$train, ], sim$y[f$train])
    pred[f$test] <- shrinkage_lda_predict(fit, sim$X_class[f$test, ])
  })[["elapsed"]]
  list(acc = mean(pred == sim$y), secs = secs)
}

searchlight_run <- function() {
  spheres <- make_spheres(sim$coords, radius = 2)
  pred <- factor(rep(NA_character_, n), levels = levels(sim$y))
  oracle_hits <- 0L
  secs <- system.time(for (f in folds) {
    tr <- f$train
    # honest selection: best mean accuracy over inner leave-one-block-out
    # folds of the training rows; the outer test rows are never touched
    tr_blocks <- unique(sim$blocks[tr])
    inner <- lapply(tr_blocks, function(b) list(train = tr[sim$blocks[tr] != b],
                                                test = tr[sim$blocks[tr] == b]))
    inner_acc <- rowMeans(vapply(inner, function(g)
      sphere_accs(sim$X_class, sim$y, g$train, g$test, spheres), numeric(sim$p)))
    best <- which.max(inner_acc)
    vox <- spheres[[best]]
    fit <- shrinkage_lda_fit(sim$X_class[tr, vox, drop = FALSE], sim$y[tr])
    pred[f$test] <- shrinkage_lda_predict(fit, sim$X_class[f$test, vox, drop = FALSE])
    # oracle: the sphere that happens to do best on this outer fold's test rows
    outer_acc <- sphere_accs(sim$X_class, sim$y, tr, f$test, spheres)
    oracle_hits <- oracle_hits + max(outer_acc) * length(f$test)
  })[["elapsed"]]
  list(acc = mean(pred == sim$y), oracle = oracle_hits / n, secs = secs)
}

for (snr_a in snr_grid) {
  setup_sim(snr_a)
  cat(sprintf("\n### snr = %.2f\n\n", snr_a))
  cat("| method | accuracy | seconds |\n|---|---|---|\n")
  r <- pattern_run()
  cat(sprintf("| pattern_model (rank auto) | %.3f | %.1f |\n", r$acc, r$secs))
  r <- pattern_run(penalty = list(sparse = "auto"))
  cat(sprintf("| pattern_model (sparse auto) | %.3f | %.1f |\n", r$acc, r$secs))
  r <- pattern_run(penalty = list(sparse = "auto", signed_smooth = "auto"), use_graph = TRUE)
  cat(sprintf("| pattern_model (sparse + smooth auto) | %.3f | %.1f |\n", r$acc, r$secs))
  r <- spacenet_run()
  cat(sprintf("| spacenet_tvl1 | %.3f | %.1f |\n", r$acc, r$secs))
  r <- lda_run()
  cat(sprintf("| whole-brain shrinkage LDA | %.3f | %.1f |\n", r$acc, r$secs))
  r <- searchlight_run()
  cat(sprintf("| searchlight LDA (honest best sphere) | %.3f | %.1f |\n", r$acc, r$secs))
  cat(sprintf("| searchlight LDA (oracle sphere; peeks at test) | %.3f | - |\n", r$oracle))
}
cat(sprintf("\nChance is %.3f.\n", 1 / 3))
setup_sim(snr_grid[1])   # Part B uses the first condition

# ---------------------------------------------------------------------------
# Part B: continuous decoding (brain -> feature vector) on the same splits
# ---------------------------------------------------------------------------

cat("\n## Head-to-head continuous decoding (q = 5 responses, mean predictive R^2)\n\n")

r2_of <- function(pred_list) {
  # pooled R^2 per response against each fold's training-mean baseline
  sse <- 0; sst <- 0
  for (f in seq_along(folds)) {
    te <- folds[[f]]$test; tr <- folds[[f]]$train
    Y <- sim$Yq[te, , drop = FALSE]
    mu <- colMeans(sim$Yq[tr, , drop = FALSE])
    sse <- sse + colSums((Y - pred_list[[f]])^2)
    sst <- sst + colSums(sweep(Y, 2L, mu, "-")^2)
  }
  mean(1 - sse / sst)
}

pattern_decode_run <- function() {
  ctrl <- pattern_control(max_rank = 3)
  preds <- vector("list", length(folds))
  secs <- system.time(for (f in seq_along(folds)) {
    tr <- folds[[f]]$train
    sel <- rMVPA:::.pattern_select_config(sim$X_cont[tr, ], sim$Yq[tr, ], sim$blocks[tr], ctrl)
    fit <- rMVPA:::.pattern_fit(sim$X_cont[tr, ], sim$Yq[tr, ], rank = sel$rank,
                                control = ctrl, cap_rank = TRUE)
    preds[[f]] <- predict(fit, sim$X_cont[folds[[f]]$test, ], type = "decode")
  })[["elapsed"]]
  list(r2 = r2_of(preds), secs = secs)
}

pls_run <- function(max_comp = 8) {
  preds <- vector("list", length(folds))
  secs <- system.time(for (f in seq_along(folds)) {
    tr <- folds[[f]]$train
    # ncomp by inner leave-one-block-out CV on the training rows
    tr_blocks <- unique(sim$blocks[tr])
    sse <- numeric(max_comp)
    for (b in tr_blocks) {
      itr <- tr[sim$blocks[tr] != b]; ite <- tr[sim$blocks[tr] == b]
      pf <- pls::plsr(Y ~ X, ncomp = max_comp,
                      data = list(Y = sim$Yq[itr, , drop = FALSE],
                                  X = sim$X_cont[itr, , drop = FALSE]))
      ph <- predict(pf, newdata = list(X = sim$X_cont[ite, , drop = FALSE]))
      for (k in seq_len(max_comp)) sse[k] <- sse[k] + sum((sim$Yq[ite, ] - ph[, , k])^2)
    }
    k_best <- which.min(sse)
    pf <- pls::plsr(Y ~ X, ncomp = k_best,
                    data = list(Y = sim$Yq[tr, , drop = FALSE],
                                X = sim$X_cont[tr, , drop = FALSE]))
    ph <- predict(pf, newdata = list(X = sim$X_cont[folds[[f]]$test, , drop = FALSE]))
    preds[[f]] <- ph[, , k_best]
  })[["elapsed"]]
  list(r2 = r2_of(preds), secs = secs)
}

cat("| method | mean R^2 | seconds |\n|---|---|---|\n")
r <- pattern_decode_run()
cat(sprintf("| pattern_model (rank auto, decode) | %.3f | %.1f |\n", r$r2, r$secs))
r <- pls_run()
cat(sprintf("| PLS (ncomp by inner CV) | %.3f | %.1f |\n", r$r2, r$secs))
cat("\nBaseline (training-mean prediction) has R^2 = 0 by construction.\n")
