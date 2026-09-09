# Benchmarks for pattern_model (Phase 3).
#
# Not run by the test suite. Two questions:
#
#   1. Cost. How long does the *whole* workflow take as the feature count
#      grows -- tuning, covariance estimation, fitting, and prediction, never
#      a single favourable solver iteration.
#   2. Where accuracy is lost. When the model decodes worse than an
#      independently fitted local model, is that a poor residual covariance or
#      a missed signal subspace? The regional comparison below separates them
#      by holding the learned patterns fixed and swapping only the prediction
#      head.
#
# Run:  Rscript inst/benchmarks/pattern_model/bench_pattern_model.R [quick]
# Writes a markdown table to stdout; record results in
# .planning/pattern-model-benchmarks.md.

suppressPackageStartupMessages(library(rMVPA))
set.seed(11)
quick <- "quick" %in% commandArgs(trailingOnly = TRUE)

`%||%` <- function(x, y) if (is.null(x)) y else x

# time and peak R heap in one pass
timeit <- function(expr) {
  gc(reset = TRUE, full = TRUE)
  t0 <- proc.time()[["elapsed"]]
  val <- force(expr)
  secs <- proc.time()[["elapsed"]] - t0
  g <- gc(full = TRUE)
  list(value = val, seconds = secs, peak_mb = sum(g[, "max used"] * c(8, 8)) / 1024^2)
}

# ---------------------------------------------------------------------------
# simulation: known truth, with localized nuisance and redundant signal
# ---------------------------------------------------------------------------

make_sim <- function(n, p, K = 3, r = 2, snr = 0.7, n_local_nuisance = 3, seed = 1) {
  set.seed(seed)
  y <- factor(rep(letters[seq_len(K)], length.out = n))
  Y <- model.matrix(~ y - 1)
  Yc <- scale(Y, scale = FALSE)
  S <- crossprod(Yc) / n
  eg <- eigen(S, symmetric = TRUE)
  keep <- eg$values > 1e-8 * eg$values[1]
  Wy <- eg$vectors[, keep, drop = FALSE] %*% diag(1 / sqrt(eg$values[keep]), nrow = sum(keep))
  C <- qr.Q(qr(matrix(rnorm(ncol(Wy) * r), ncol(Wy), r)))[, seq_len(r), drop = FALSE]

  region <- rep(seq_len(4), each = ceiling(p / 4))[seq_len(p)]   # four contiguous parcels
  A <- matrix(0, p, r)
  A[region == 1, 1] <- rnorm(sum(region == 1)) * snr
  A[region == 2, 2] <- rnorm(sum(region == 2)) * snr
  A[region == 4, 1] <- rnorm(sum(region == 4)) * snr  # redundant with region 1

  # nuisance confined to region 3: local covariance the global model cannot see
  U <- matrix(0, p, n_local_nuisance)
  idx3 <- which(region == 3)
  for (j in seq_len(n_local_nuisance)) U[idx3, j] <- rnorm(length(idx3)) * 1.2

  make_X <- function(m, yy) {
    Ym <- model.matrix(~ yy - 1)
    Tm <- sweep(Ym, 2L, colMeans(Y), "-") %*% Wy %*% C
    E <- matrix(rnorm(m * p), m, p) + matrix(rnorm(m * ncol(U)), m, ncol(U)) %*% t(U)
    Tm %*% t(A) + E
  }
  X <- make_X(n, y)
  blocks <- as.integer(cut(seq_len(n), 3, labels = 1:3))
  list(X = X, y = y, A = A, region = region, blocks = blocks, p = p, n = n)
}

# ---------------------------------------------------------------------------
# 1. cost of the whole workflow
# ---------------------------------------------------------------------------

cat("## Whole-workflow cost\n\n")
cat("Rank is tuned in every row; the `penalty` column says whether the penalty\n")
cat("strength is also tuned.\n\n")
cat("| p | n | penalty | seconds | peak MB |\n|---|---|---|---|---|\n")

sizes <- if (quick) list(c(2000, 200)) else list(c(5000, 400), c(30000, 400), c(100000, 400))
for (sz in sizes) {
  p <- sz[1]; n <- sz[2]
  sim <- make_sim(n, p, seed = 2)
  graph <- NULL
  for (pen in list(NULL, list(sparse = 0.2), list(sparse = "auto"),
                   list(sparse = 0.2, signed_smooth = 4))) {
    label <- if (is.null(pen)) "none" else paste(names(pen), unlist(pen), sep = "=", collapse = ",")
    if (!is.null(pen$signed_smooth) && is.null(graph)) {
      # a chain graph over the feature order: enough to exercise the Laplacian
      graph <- spatial_graph(list(A = Matrix::bandSparse(
        p, k = 1, diagonals = list(rep(1, p - 1)))))
    }
    tm <- timeit({
      folds <- split(seq_len(n), sim$blocks)
      for (k in seq_along(folds)) {
        tr <- setdiff(seq_len(n), folds[[k]])
        sel <- rMVPA:::.pattern_select_config(sim$X[tr, , drop = FALSE], sim$y[tr],
                                              sim$blocks[tr], pattern_control(max_rank = 3),
                                              penalty = pen, graph = graph)
        fit <- rMVPA:::.pattern_fit(sim$X[tr, , drop = FALSE], sim$y[tr], rank = sel$rank,
                                    control = pattern_control(max_rank = 3),
                                    penalty = sel$penalty, graph = graph, cap_rank = TRUE)
        invisible(predict(fit, sim$X[folds[[k]], , drop = FALSE], type = "prob"))
      }
    })
    cat(sprintf("| %d | %d | %s | %.1f | %.0f |\n", p, n, label, tm$seconds, tm$peak_mb))
  }
}

# ---------------------------------------------------------------------------
# 2. regional accuracy: covariance model vs signal subspace
# ---------------------------------------------------------------------------
#
# All rows use the SAME assessment rows and the SAME region. Rows 1-3 keep the
# whole-brain patterns A_R fixed and change only the prediction head, so a gap
# among them is a covariance-model effect. Row 4 fits an independent local
# model, so a gap between the best of 1-3 and row 4 is a missed subspace.

cat("\n## Regional accuracy: where is the loss?\n\n")
cat("| region | restricted (diag) | restricted (diag+lowrank) | adapted (local shrinkage) | independent shrinkage LDA |\n")
cat("|---|---|---|---|---|\n")

shrinkage_lda <- function(Xtr, ytr, Xte, lambda = 0.1) {
  lev <- levels(ytr)
  mu <- t(vapply(lev, function(l) colMeans(Xtr[ytr == l, , drop = FALSE]), numeric(ncol(Xtr))))
  Xc <- Xtr - mu[as.integer(ytr), , drop = FALSE]
  S <- crossprod(Xc) / max(nrow(Xtr) - length(lev), 1)
  S <- (1 - lambda) * S + lambda * mean(diag(S)) * diag(ncol(Xtr))
  W <- solve(S, t(mu))
  sc <- Xte %*% W - matrix(0.5 * colSums(t(mu) * W), nrow(Xte), length(lev), byrow = TRUE)
  factor(lev[max.col(sc, ties.method = "first")], levels = lev)
}

local_head <- function(fit, Xtr, ytr, Xte, keep, mode) {
  A_R <- fit$A[keep, , drop = FALSE]
  noise_R <- switch(mode,
    diag = rMVPA:::new_pattern_noise(fit$noise$D[keep], NULL),
    lowrank = rMVPA:::.noise_restrict(fit$noise, keep),
    adapted = {
      # re-estimate the residual covariance inside the region, on training rows
      Tm <- rMVPA:::.pattern_targets_apply(fit$y_transform, ytr) %*% fit$C
      Xc <- sweep(Xtr[, keep, drop = FALSE], 2L, fit$x_transform$mu[keep], "-")
      E <- Xc - Tm %*% t(A_R)
      rMVPA:::.estimate_pattern_noise(E, type = "diag_lowrank", rank = "auto", max_rank = 5)
    })
  sub <- fit
  sub$A <- A_R
  sub$noise <- noise_R
  sub$precision_A <- rMVPA:::.noise_apply_precision(noise_R, A_R)
  sub$G <- crossprod(A_R, sub$precision_A)
  sub$x_transform <- list(mu = fit$x_transform$mu[keep], sd = fit$x_transform$sd[keep],
                          scale = fit$x_transform$scale)
  sub$feature_index <- seq_along(keep)
  sub$p_input <- length(keep)
  predict(sub, Xte[, keep, drop = FALSE], type = "class")
}

# snr chosen so accuracies land away from both chance and ceiling; at higher
# snr every head scores 1.000 and the comparison says nothing.
sim <- make_sim(n = 240, p = if (quick) 400 else 2000, snr = 0.08, seed = 3)
tr <- which(sim$blocks != 3); te <- which(sim$blocks == 3)
fit <- rMVPA:::.pattern_fit(sim$X[tr, ], sim$y[tr], rank = 2, cap_rank = TRUE)
for (rg in 1:4) {
  keep <- which(sim$region == rg)
  acc <- function(p) mean(p == sim$y[te])
  a1 <- acc(local_head(fit, sim$X[tr, ], sim$y[tr], sim$X[te, ], keep, "diag"))
  a2 <- acc(local_head(fit, sim$X[tr, ], sim$y[tr], sim$X[te, ], keep, "lowrank"))
  a3 <- acc(local_head(fit, sim$X[tr, ], sim$y[tr], sim$X[te, ], keep, "adapted"))
  a4 <- acc(shrinkage_lda(sim$X[tr, keep, drop = FALSE], sim$y[tr],
                          sim$X[te, keep, drop = FALSE]))
  cat(sprintf("| %d%s | %.3f | %.3f | %.3f | %.3f |\n", rg,
              if (rg == 3) " (local nuisance)" else if (rg == 4) " (redundant)" else "",
              a1, a2, a3, a4))
}

cat("\nRows 1-3 share the whole-brain patterns and differ only in the residual\n",
    "covariance used by the prediction head, so differences among them are a\n",
    "covariance-model effect. A gap between the best of those and the last\n",
    "column is a missed signal subspace, not a covariance problem.\n", sep = "")

# ---------------------------------------------------------------------------
# 3. where does the time go inside a penalized fit?
# ---------------------------------------------------------------------------
#
# The plan gates a compiled kernel on whether the A-step dominates at
# whole-brain feature counts. Measure the share rather than guess.

cat("\n## A-step share of a penalized fit\n\n")
p <- if (quick) 5000 else 100000
n <- 400
sim <- make_sim(n, p, snr = 0.3, seed = 4)
ctrl <- pattern_control(max_rank = 2)
tm <- timeit(rMVPA:::.pattern_fit(sim$X, sim$y, rank = 2, control = ctrl,
                                  penalty = list(sparse = 0.2)))
fit <- tm$value
Xc <- scale(sim$X, scale = FALSE)
Yw <- rMVPA:::.pattern_targets_apply(fit$y_transform, sim$y)
const <- rMVPA:::.noise_quadform(fit$noise, Xc)
prob <- rMVPA:::.pattern_astep_problem(Xc, Yw %*% fit$C, fit$noise, NULL,
                                       fit$penalty$lambda_s, 0, const)
reps <- 20
t_iter <- (function() {
  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(reps)) {
    invisible(rMVPA:::.pattern_astep_gradient(fit$A, prob))
    invisible(rMVPA:::.pattern_astep_objective(fit$A, prob))
  }
  (proc.time()[["elapsed"]] - t0) / reps
})()
n_inner <- sum(vapply(fit$diagnostics$inner, `[[`, numeric(1), "iterations"))
cat(sprintf("| p | fit seconds | inner iterations | seconds per iteration | A-step share |\n"))
cat("|---|---|---|---|---|\n")
cat(sprintf("| %d | %.1f | %d | %.4f | %.0f%% |\n", p, tm$seconds, n_inner, t_iter,
            100 * n_inner * t_iter / tm$seconds))
cat("\nIf the A-step share is small, a compiled proximal kernel would not change\n",
    "the workflow's cost and is not worth the build complexity.\n", sep = "")
