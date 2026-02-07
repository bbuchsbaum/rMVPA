## REMAP-RRR diagnostics battery (standalone harness)
## Requires: rrpack, corpcor

suppressPackageStartupMessages({
  require(rrpack)
  require(corpcor)
})

row_center <- function(M) sweep(M, 2, colMeans(M), "-")

rand_orth <- function(p, r) {
  Q <- qr.Q(qr(matrix(rnorm(p * r), p, r)))
  Q[, seq_len(r), drop = FALSE]
}

whiten_shrink <- function(X, eps = 1e-7) {
  Xc <- row_center(X)
  S  <- cov.shrink(Xc, verbose = FALSE)
  ee <- eigen(S, symmetric = TRUE)
  Winvhalf <- ee$vectors %*% diag(1 / sqrt(pmax(ee$values, eps))) %*% t(ee$vectors)
  list(Xw = Xc %*% Winvhalf, W = Winvhalf, mu = colMeans(X))
}

row_cor <- function(A, v) {
  A0 <- row_center(A)
  v0 <- v - mean(v)
  den <- sqrt(rowSums(A0^2)) * sqrt(sum(v0^2))
  as.numeric((A0 %*% v0) / pmax(den, .Machine$double.eps))
}

simulate_pm <- function(K = 60, p = 500, r_lat = 10, r_map = 5,
                        n_rep_p = 6, snr_p = 2.0, snr_m = 1.5,
                        map_strength = 1, cov_shift = 0.0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Lp <- rand_orth(p, r_lat)
  S  <- matrix(rnorm(K * r_lat), K, r_lat)
  X_signal <- S %*% t(Lp)
  U <- rand_orth(p, r_map); V <- rand_orth(p, r_map)
  sv <- map_strength * seq(r_map, 1, length.out = r_map) / r_map
  B_true <- V %*% diag(sv) %*% t(U)
  Y_signal <- X_signal %*% B_true
  sigvar_p <- mean(apply(X_signal, 2, var))
  noise_sd_p <- sqrt(sigvar_p / pmax(snr_p, 1e-6))
  X_trials <- array(NA_real_, dim = c(K, n_rep_p, p))
  for (i in 1:K) for (j in 1:n_rep_p) X_trials[i, j, ] <- X_signal[i, ] + rnorm(p, 0, noise_sd_p)
  X_proto <- apply(X_trials, c(1,3), mean)
  sigvar_m <- mean(apply(Y_signal, 2, var))
  base_sd_m <- sqrt(sigvar_m / pmax(snr_m, 1e-6))
  if (cov_shift > 0) {
    Qm <- rand_orth(p, r_lat)
    scales <- 1 + cov_shift * seq(1, 2, length.out = r_lat)
    Sigma_m <- Qm %*% diag(scales^2) %*% t(Qm) + diag(p)
    Lm <- eigen(Sigma_m, symmetric = TRUE)
    A  <- Lm$vectors %*% diag(sqrt(pmax(Lm$values, 1e-8))) %*% t(Lm$vectors)
    m_noise <- function() as.numeric((base_sd_m) * (rnorm(p) %*% A))
  } else {
    m_noise <- function() rnorm(p, 0, base_sd_m)
  }
  Y_trial <- Y_signal + t(replicate(K, m_noise()))
  list(X_proto = X_proto, X_trials = X_trials, Y_trial = Y_trial, B_true = B_true)
}

decode_naive <- function(X_proto, Y_trial) {
  K <- nrow(X_proto)
  pred <- integer(K)
  for (i in 1:K) {
    scores <- row_cor(X_proto, Y_trial[i, ])
    pred[i] <- which.max(scores)
  }
  mean(pred == seq_len(K))
}

decode_remap_rrr_loio <- function(X_proto, Y_trial, rank = "auto") {
  K <- nrow(X_proto); p <- ncol(X_proto)
  correct <- logical(K)
  for (i in 1:K) {
    keep <- (seq_len(K) != i)
    Xp_fit <- X_proto[keep, , drop = FALSE]
    Ym_fit <- Y_trial[keep, , drop = FALSE]
    Wp <- whiten_shrink(Xp_fit); Xp_w <- Wp$Xw
    Wm <- whiten_shrink(Ym_fit); Ym_w <- Wm$Xw
    if (identical(rank, "auto")) {
      nfold <- max(3, min(5, floor(nrow(Xp_w)/3)))
      cv    <- rrpack::cv.rrr(Y = Ym_w, X = Xp_w, nfold = nfold)
      rsel  <- max(1, min(cv$rank, min(ncol(Xp_w), ncol(Ym_w), nrow(Xp_w)-1)))
      fit   <- rrpack::rrr.fit(Y = Ym_w, X = Xp_w, nrank = rsel)
    } else {
      rsel  <- as.integer(rank)
      fit   <- rrpack::rrr.fit(Y = Ym_w, X = Xp_w, nrank = rsel)
    }
    Cw <- fit$coef
    Xp_all_w <- (row_center(X_proto) %*% Wp$W)
    Xp_pred_w <- Xp_all_w %*% Cw
    y_i_w <- ((Y_trial[i, , drop = FALSE] - matrix(Wm$mu, 1, p)) %*% Wm$W)[1, ]
    scores <- row_cor(Xp_pred_w, y_i_w)
    pred_i <- which.max(scores)
    correct[i] <- (pred_i == i)
  }
  mean(correct)
}

run_battery <- function(n_sims = 30,
                        K = 48, p = 400, r_lat = 8, r_map = 4,
                        n_rep_p = 6, snr_p = 2.0, snr_m = 1.5,
                        map_strength = 1.0, cov_shift = 0.0,
                        rank = "auto", seed = 1) {
  set.seed(seed)
  res <- data.frame(sim = integer(), method = character(), acc = numeric())
  for (s in 1:n_sims) {
    sim <- simulate_pm(K, p, r_lat, r_map, n_rep_p, snr_p, snr_m, map_strength, cov_shift)
    acc_naive <- decode_naive(sim$X_proto, sim$Y_trial)
    acc_remap <- decode_remap_rrr_loio(sim$X_proto, sim$Y_trial, rank = rank)
    res <- rbind(res,
                 data.frame(sim = s, method = "naive", acc = acc_naive),
                 data.frame(sim = s, method = "remap_rrr", acc = acc_remap))
  }
  res
}

if (interactive()) {
  res_A <- run_battery(n_sims = 20, r_map = 0, map_strength = 0, cov_shift = 0.0, snr_p = 2.5, snr_m = 2.5, seed = 42)
  res_B <- run_battery(n_sims = 20, r_map = 4, map_strength = 1.0, cov_shift = 0.0, seed = 43)
  res_C <- run_battery(n_sims = 20, r_map = 4, map_strength = 1.0, cov_shift = 0.8, snr_m = 1.2, seed = 44)
  summarize <- function(df, label) {
    na <- subset(df, method=="naive")$acc
    rr <- subset(df, method=="remap_rrr")$acc
    cat("\n", label, "\n",
        sprintf("naive:  mean=%.3f  sd=%.3f", mean(na), sd(na)), "\n",
        sprintf("remap:  mean=%.3f  sd=%.3f", mean(rr), sd(rr)), "\n",
        sprintf("improv: mean=%.3f (paired t p=%.3g)",
                mean(rr - na), t.test(rr, na, paired=TRUE)$p.value), "\n", sep="")
  }
  summarize(res_A, "Scenario A: No domain shift (sanity)")
  summarize(res_B, "Scenario B: Low-rank mapping")
  summarize(res_C, "Scenario C: Low-rank + covariance shift")
}

