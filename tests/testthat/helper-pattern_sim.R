# Known-truth simulator for pattern_model tests.
#
# Generates X = T A' + E on a small volumetric grid with an mvpa_dataset /
# mvpa_design pair plus the ground truth used to build it.
#
# scenarios
#   basic      : r disjoint compact informative regions, diagonal noise
#   suppressor : one informative region plus a region that carries a shared
#                noise component but no task signal (a noise canceller)
#   redundant  : the same task dimension expressed in two separate regions
#   diffuse    : weak signal spread over most of the volume
#   signflip   : one contiguous informative region with alternating-sign
#                loadings (fine-scale code inside a coherent territory)
sim_pattern_data <- function(n = 90, dims = c(8, 8, 4), K = 3, q = NULL, r = NULL,
                             snr = 2, scenario = c("basic", "suppressor", "redundant",
                                                   "diffuse", "signflip"),
                             blocks = 3, noise_sd = 1, external_test = FALSE,
                             n_test = n, seed = NULL) {
  scenario <- match.arg(scenario)
  if (!is.null(seed)) set.seed(seed)
  p <- prod(dims)
  categorical <- is.null(q)
  if (is.null(r)) r <- if (categorical) min(K - 1L, 2L) else min(q, 2L)

  coords <- arrayInd(seq_len(p), dims)
  region_box <- function(xr, yr) {
    which(coords[, 1] >= xr[1] & coords[, 1] <= xr[2] &
          coords[, 2] >= yr[1] & coords[, 2] <= yr[2])
  }
  half_x <- floor(dims[1] / 2); half_y <- floor(dims[2] / 2)
  boxes <- list(region_box(c(1, half_x), c(1, half_y)),
                region_box(c(half_x + 1, dims[1]), c(1, half_y)),
                region_box(c(1, half_x), c(half_y + 1, dims[2])),
                region_box(c(half_x + 1, dims[1]), c(half_y + 1, dims[2])))

  A <- matrix(0, p, r)
  U <- NULL
  for (k in seq_len(r)) A[boxes[[k]], k] <- rnorm(length(boxes[[k]]))
  if (scenario == "suppressor") {
    A[] <- 0; A[boxes[[1]], 1] <- abs(rnorm(length(boxes[[1]]))) + 0.5
    U <- matrix(0, p, 1); U[c(boxes[[1]], boxes[[2]]), 1] <- 2
    r <- 1L; A <- A[, 1, drop = FALSE]
  } else if (scenario == "redundant") {
    A[] <- 0; A[boxes[[1]], 1] <- rnorm(length(boxes[[1]])); A[boxes[[4]], 1] <- rnorm(length(boxes[[4]]))
    r <- 1L; A <- A[, 1, drop = FALSE]
  } else if (scenario == "diffuse") {
    A[] <- 0; A[, 1] <- rnorm(p, sd = 0.3)
    r <- 1L; A <- A[, 1, drop = FALSE]
  } else if (scenario == "signflip") {
    A[] <- 0
    idx <- boxes[[1]]
    A[idx, 1] <- ifelse((coords[idx, 1] + coords[idx, 2]) %% 2 == 0, 1, -1) * (1 + abs(rnorm(length(idx), sd = 0.2)))
    r <- 1L; A <- A[, 1, drop = FALSE]
  }
  A <- A * snr

  # Raw targets. Centred one-hot codes are rank K - 1, so a direction drawn in
  # raw R^K would put part of the signal in the null space where it cancels.
  # We therefore whiten the targets on the training rows first and draw C in
  # the whitened space: the score columns then have unit variance and `snr` is
  # the per-voxel signal standard deviation relative to `noise_sd`.
  make_targets <- function(m) {
    if (categorical) {
      y <- factor(rep(letters[seq_len(K)], length.out = m))
      list(y = y, Y = model.matrix(~ y - 1))
    } else {
      Y <- matrix(rnorm(m * q), m, q, dimnames = list(NULL, paste0("f", seq_len(q))))
      list(y = Y, Y = Y)
    }
  }
  tg <- make_targets(n)
  y_mu <- colMeans(tg$Y)
  S <- crossprod(sweep(tg$Y, 2L, y_mu, "-")) / n
  eg <- eigen(S, symmetric = TRUE)
  keep_e <- eg$values > 1e-8 * max(eg$values[1], .Machine$double.eps)
  Wy <- eg$vectors[, keep_e, drop = FALSE] %*% diag(1 / sqrt(eg$values[keep_e]), nrow = sum(keep_e))
  q_eff <- ncol(Wy)
  r <- min(r, q_eff)
  A <- A[, seq_len(r), drop = FALSE]
  C <- qr.Q(qr(matrix(rnorm(q_eff * r), q_eff, r)))[, seq_len(r), drop = FALSE]

  whiten <- function(Y) sweep(Y, 2L, y_mu, "-") %*% Wy
  make_X <- function(m, tg) {
    Tm <- whiten(tg$Y) %*% C
    E <- matrix(rnorm(m * p, sd = noise_sd), m, p)
    if (!is.null(U)) E <- E + matrix(rnorm(m * ncol(U)), m, ncol(U)) %*% t(U)
    Tm %*% t(A) + E
  }
  to_vec <- function(X) {
    m <- nrow(X)
    arr <- array(t(X), c(dims, m))
    neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(c(dims, m), c(1, 1, 1)))
  }

  X <- make_X(n, tg)
  mask <- as.logical(neuroim2::NeuroVol(array(1, dims), neuroim2::NeuroSpace(dims, c(1, 1, 1))))
  block_var <- as.integer(as.character(cut(seq_len(n), blocks, labels = seq_len(blocks))))

  if (external_test) {
    tg_test <- make_targets(n_test)
    X_test <- make_X(n_test, tg_test)
    dataset <- mvpa_dataset(train_data = to_vec(X), test_data = to_vec(X_test), mask = mask)
    if (categorical) {
      design <- mvpa_design(data.frame(y = tg$y, block = block_var),
                            data.frame(y = tg_test$y),
                            y_train = ~ y, y_test = ~ y, block_var = ~ block)
    } else {
      design <- mvpa_design(data.frame(id = seq_len(n), block = block_var),
                            data.frame(id = seq_len(n_test)),
                            cv_labels = seq_len(n), targets = tg$y,
                            y_test = seq_len(n_test), targets_test = tg_test$y,
                            block_var = ~ block)
    }
  } else {
    X_test <- NULL; tg_test <- NULL
    dataset <- mvpa_dataset(train_data = to_vec(X), mask = mask)
    if (categorical) {
      design <- mvpa_design(data.frame(y = tg$y, block = block_var), y_train = ~ y, block_var = ~ block)
    } else {
      design <- mvpa_design(data.frame(id = seq_len(n), block = block_var),
                            cv_labels = seq_len(n), targets = tg$y, block_var = ~ block)
    }
  }

  list(dataset = dataset, design = design, X = X, X_test = X_test,
       targets = tg$y, targets_test = if (!is.null(tg_test)) tg_test$y else NULL,
       A = A, C = C, U = U, boxes = boxes, r = r, dims = dims, scenario = scenario,
       informative = which(rowSums(A != 0) > 0))
}
