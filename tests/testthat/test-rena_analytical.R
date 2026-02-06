context("ReNA analytical verification: exact mathematical correctness")

## These tests verify that ReNA model outputs match hand-computed expected values
## using simple, deterministic synthetic data where the correct answer is known
## analytically. This provides stronger guarantees than the heuristic-threshold
## tests in the model-specific test files.

library(testthat)
library(rMVPA)

## The 4x4x4 mask = 64 voxels. All item-pattern matrices use V=64.
NVOX <- 64L

## ---------------------------------------------------------------------------
## Helper: build a minimal ROI + model spec from explicit item-level patterns
## ---------------------------------------------------------------------------
build_rena_roi <- function(E_items, items, reps_per_item = 2L) {
  K <- length(items)
  stopifnot(ncol(E_items) == NVOX)
  n_obs <- K * reps_per_item

  ds <- gen_sample_dataset(D = c(4, 4, 4),
                           nobs = n_obs,
                           blocks = reps_per_item,
                           nlevels = 2)
  mvdes <- mvpa_design(
    train_design = ds$design$train_design,
    y_train      = ds$design$y_train,
    block_var    = ds$design$block_var
  )
  mvdes$train_design$.rownum <- rep(items, each = reps_per_item)

  vox <- which(ds$dataset$mask > 0)
  n_vox <- length(vox)
  stopifnot(n_vox == NVOX)
  samp <- data_sample(ds$dataset, vox)
  roi  <- as_roi(samp, ds$dataset)

  X_rep <- E_items[as.character(mvdes$train_design$.rownum), , drop = FALSE]
  roi$train_roi@.Data <- X_rep

  list(ds = ds, mvdes = mvdes, roi = roi, n_vox = n_vox)
}

## ===========================================================================
## SECTION 1: repnet_model analytical tests
## ===========================================================================

test_that("repnet: identical ROI and seed RDMs yield conn_raw = 1.0", {
  set.seed(42)
  K <- 8
  items <- paste0("it", seq_len(K))

  E <- matrix(rnorm(K * NVOX), nrow = K, ncol = NVOX)
  rownames(E) <- items

  ## Compute the ROI RDM that process_roi will compute: 1 - cor(t(E))
  D_roi <- 1 - cor(t(E))
  rownames(D_roi) <- colnames(D_roi) <- items

  seed_rdm <- D_roi

  env <- build_rena_roi(E, items, reps_per_item = 2L)

  rn_des <- repnet_design(
    design    = env$mvdes,
    key_var   = ~ .rownum,
    seed_rdm  = seed_rdm
  )
  mspec <- repnet_model(
    dataset    = env$ds$dataset,
    design     = env$mvdes,
    repnet_des = rn_des,
    distfun    = cordist("pearson"),
    simfun     = "pearson"
  )

  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  perf <- res$performance[[1]]

  expect_false(res$error)
  expect_equal(unname(perf["conn_raw"]), 1.0, tolerance = 1e-10)
})


test_that("repnet: conn_raw matches hand-computed correlation of lower triangles", {
  set.seed(43)
  K <- 10
  items <- paste0("it", seq_len(K))

  E <- matrix(rnorm(K * NVOX), nrow = K, ncol = NVOX)
  rownames(E) <- items

  ## Hand-compute ROI RDM
  D_roi <- 1 - cor(t(E))
  d_roi_vec <- D_roi[lower.tri(D_roi)]

  ## Seed is a different RDM
  seed_features <- matrix(rnorm(K * 8), nrow = K, ncol = 8)
  D_seed <- 1 - cor(t(seed_features))
  rownames(D_seed) <- colnames(D_seed) <- items
  d_seed_vec <- D_seed[lower.tri(D_seed)]

  ## Expected conn_raw
  expected_conn <- cor(d_roi_vec, d_seed_vec)

  env <- build_rena_roi(E, items, reps_per_item = 2L)

  rn_des <- repnet_design(
    design    = env$mvdes,
    key_var   = ~ .rownum,
    seed_rdm  = D_seed
  )
  mspec <- repnet_model(
    dataset    = env$ds$dataset,
    design     = env$mvdes,
    repnet_des = rn_des,
    distfun    = cordist("pearson"),
    simfun     = "pearson"
  )

  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  perf <- res$performance[[1]]

  expect_false(res$error)
  expect_equal(unname(perf["conn_raw"]), expected_conn, tolerance = 1e-10)
})


test_that("repnet: beta_seed matches hand-computed partial correlation", {
  set.seed(44)
  K <- 12
  items <- paste0("it", seq_len(K))

  E <- matrix(rnorm(K * NVOX), nrow = K, ncol = NVOX)
  rownames(E) <- items

  D_roi <- 1 - cor(t(E))
  d_roi_vec <- D_roi[lower.tri(D_roi)]

  ## Main seed + one confound
  F_seed <- matrix(rnorm(K * 6), nrow = K, ncol = 6)
  D_seed <- 1 - cor(t(F_seed))
  rownames(D_seed) <- colnames(D_seed) <- items
  d_seed_vec <- D_seed[lower.tri(D_seed)]

  F_conf <- matrix(rnorm(K * 5), nrow = K, ncol = 5)
  D_conf <- 1 - cor(t(F_conf))
  rownames(D_conf) <- colnames(D_conf) <- items
  d_conf_vec <- D_conf[lower.tri(D_conf)]

  ## Expected: the partial correlation of seed residualized against confound
  ## (this is what repnet_model computes as beta_seed when confounds present)
  resid_seed <- lm(d_seed_vec ~ d_conf_vec)$residuals
  expected_partial <- cor(resid_seed, d_roi_vec)

  env <- build_rena_roi(E, items, reps_per_item = 2L)

  rn_des <- repnet_design(
    design        = env$mvdes,
    key_var       = ~ .rownum,
    seed_rdm      = D_seed,
    confound_rdms = list(conf1 = D_conf)
  )
  mspec <- repnet_model(
    dataset    = env$ds$dataset,
    design     = env$mvdes,
    repnet_des = rn_des,
    distfun    = cordist("pearson"),
    simfun     = "pearson"
  )

  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  perf <- res$performance[[1]]

  expect_false(res$error)
  expect_equal(unname(perf["beta_seed"]), expected_partial, tolerance = 1e-10)
})


test_that("repnet: Spearman simfun matches hand-computed Spearman correlation", {
  set.seed(45)
  K <- 10
  items <- paste0("it", seq_len(K))

  E <- matrix(rnorm(K * NVOX), nrow = K, ncol = NVOX)
  rownames(E) <- items

  D_roi <- 1 - cor(t(E))
  d_roi_vec <- D_roi[lower.tri(D_roi)]

  seed_features <- matrix(rnorm(K * 8), nrow = K, ncol = 8)
  D_seed <- 1 - cor(t(seed_features))
  rownames(D_seed) <- colnames(D_seed) <- items
  d_seed_vec <- D_seed[lower.tri(D_seed)]

  expected_conn <- cor(d_roi_vec, d_seed_vec, method = "spearman")

  env <- build_rena_roi(E, items, reps_per_item = 2L)

  rn_des <- repnet_design(
    design    = env$mvdes,
    key_var   = ~ .rownum,
    seed_rdm  = D_seed
  )
  mspec <- repnet_model(
    dataset    = env$ds$dataset,
    design     = env$mvdes,
    repnet_des = rn_des,
    distfun    = cordist("pearson"),
    simfun     = "spearman"
  )

  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  perf <- res$performance[[1]]

  expect_false(res$error)
  expect_equal(unname(perf["conn_raw"]), expected_conn, tolerance = 1e-10)
})


test_that("repnet: negated RDM yields conn_raw = -1.0", {
  set.seed(46)
  K <- 8
  items <- paste0("it", seq_len(K))

  E <- matrix(rnorm(K * NVOX), nrow = K, ncol = NVOX)
  rownames(E) <- items

  D_roi <- 1 - cor(t(E))

  ## Negate the RDM: conn should be -1
  seed_rdm <- -D_roi
  rownames(seed_rdm) <- colnames(seed_rdm) <- items

  env <- build_rena_roi(E, items, reps_per_item = 2L)

  rn_des <- repnet_design(
    design    = env$mvdes,
    key_var   = ~ .rownum,
    seed_rdm  = seed_rdm
  )
  mspec <- repnet_model(
    dataset    = env$ds$dataset,
    design     = env$mvdes,
    repnet_des = rn_des,
    distfun    = cordist("pearson"),
    simfun     = "pearson"
  )

  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  perf <- res$performance[[1]]

  expect_false(res$error)
  expect_equal(unname(perf["conn_raw"]), -1.0, tolerance = 1e-10)
})


test_that("repnet: n_items and n_pairs are correct", {
  set.seed(47)
  K <- 7
  items <- paste0("it", seq_len(K))

  E <- matrix(rnorm(K * NVOX), nrow = K, ncol = NVOX)
  rownames(E) <- items

  D_seed <- as.matrix(dist(seq_len(K)))
  rownames(D_seed) <- colnames(D_seed) <- items

  env <- build_rena_roi(E, items, reps_per_item = 2L)

  rn_des <- repnet_design(
    design    = env$mvdes,
    key_var   = ~ .rownum,
    seed_rdm  = D_seed
  )
  mspec <- repnet_model(
    dataset    = env$ds$dataset,
    design     = env$mvdes,
    repnet_des = rn_des,
    distfun    = cordist("pearson"),
    simfun     = "pearson"
  )

  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  perf <- res$performance[[1]]

  expect_false(res$error)
  expect_equal(unname(perf["n_items"]), K)
  expect_equal(unname(perf["n_pairs"]), K * (K - 1) / 2)
})


test_that("repnet: no-confound beta_seed matches OLS coefficient exactly", {
  set.seed(48)
  K <- 10
  items <- paste0("it", seq_len(K))

  E <- matrix(rnorm(K * NVOX), nrow = K, ncol = NVOX)
  rownames(E) <- items

  D_roi <- 1 - cor(t(E))
  d_roi_vec <- D_roi[lower.tri(D_roi)]

  seed_features <- matrix(rnorm(K * 6), nrow = K, ncol = 6)
  D_seed <- 1 - cor(t(seed_features))
  rownames(D_seed) <- colnames(D_seed) <- items
  d_seed_vec <- D_seed[lower.tri(D_seed)]

  ## Without confounds, beta_seed should be the simple OLS slope
  df_hand <- data.frame(d_roi = d_roi_vec, seed = d_seed_vec)
  fit_hand <- lm(d_roi ~ seed, data = df_hand)
  expected_beta <- unname(coef(fit_hand)["seed"])

  env <- build_rena_roi(E, items, reps_per_item = 2L)

  rn_des <- repnet_design(
    design    = env$mvdes,
    key_var   = ~ .rownum,
    seed_rdm  = D_seed
  )
  mspec <- repnet_model(
    dataset    = env$ds$dataset,
    design     = env$mvdes,
    repnet_des = rn_des,
    distfun    = cordist("pearson"),
    simfun     = "pearson"
  )

  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  perf <- res$performance[[1]]

  expect_false(res$error)
  expect_equal(unname(perf["beta_seed"]), expected_beta, tolerance = 1e-10)
})


## ===========================================================================
## SECTION 2: repmed_model analytical tests
## ===========================================================================

test_that("repmed: path decomposition c = c' + ab holds exactly", {
  set.seed(50)
  K <- 20
  items <- paste0("it", seq_len(K))

  ## Deterministic latent structure: x -> m -> y with known coefficients
  z_x <- seq(-2, 2, length.out = K)
  z_m <- 0.8 * z_x + rnorm(K, sd = 0.01)
  z_y <- 0.6 * z_m + 0.3 * z_x + rnorm(K, sd = 0.01)

  X_rdm <- as.matrix(dist(z_x))
  Y_rdm <- as.matrix(dist(z_y))
  rownames(X_rdm) <- colnames(X_rdm) <- items
  rownames(Y_rdm) <- colnames(Y_rdm) <- items

  ## Mediator patterns: each voxel is z_m + tiny noise
  E_m <- matrix(z_m, nrow = K, ncol = NVOX) +
    matrix(rnorm(K * NVOX, sd = 0.001), nrow = K, ncol = NVOX)
  rownames(E_m) <- items

  env <- build_rena_roi(E_m, items, reps_per_item = 2L)

  md_des <- repmed_design(items = items, X_rdm = X_rdm, Y_rdm = Y_rdm)
  mspec <- repmed_model(
    dataset    = env$ds$dataset,
    design     = env$mvdes,
    repmed_des = md_des,
    key_var    = ~ .rownum,
    distfun    = "euclidean"
  )

  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  perf <- res$performance[[1]]

  expect_false(res$error)

  a  <- unname(perf["med_a"])
  b  <- unname(perf["med_b"])
  cp <- unname(perf["med_cprime"])
  c_tot <- unname(perf["med_c"])
  ind <- unname(perf["med_indirect"])

  ## Verify indirect = a * b
  expect_equal(ind, a * b, tolerance = 1e-12)

  ## Verify c ≈ c' + ab (total effect decomposition)
  expect_equal(c_tot, cp + ind, tolerance = 0.05)
})


test_that("repmed: coefficients match hand-computed lm() on vectorized RDMs", {
  set.seed(51)
  K <- 15
  items <- paste0("it", seq_len(K))

  z_x <- rnorm(K)
  z_m <- 0.7 * z_x + rnorm(K, sd = 0.01)
  z_y <- 0.5 * z_m + 0.2 * z_x + rnorm(K, sd = 0.01)

  X_rdm <- as.matrix(dist(z_x))
  Y_rdm <- as.matrix(dist(z_y))
  rownames(X_rdm) <- colnames(X_rdm) <- items
  rownames(Y_rdm) <- colnames(Y_rdm) <- items

  E_m <- matrix(z_m, nrow = K, ncol = NVOX) +
    matrix(rnorm(K * NVOX, sd = 0.001), nrow = K, ncol = NVOX)
  rownames(E_m) <- items

  ## Hand-compute the mediator RDM using euclidean distance
  D_M_hand <- as.matrix(dist(E_m))
  rownames(D_M_hand) <- colnames(D_M_hand) <- items

  ## Vectorize lower triangles (sorted item order, as process_roi does)
  si <- sort(items)
  X_sorted <- X_rdm[si, si]
  Y_sorted <- Y_rdm[si, si]
  M_sorted <- D_M_hand[si, si]

  x_vec <- X_sorted[lower.tri(X_sorted)]
  y_vec <- Y_sorted[lower.tri(Y_sorted)]
  m_vec <- M_sorted[lower.tri(M_sorted)]

  ## Hand-compute path a: m ~ x
  fit_a_hand <- lm(m_vec ~ x_vec)
  a_hand <- coef(fit_a_hand)["x_vec"]

  ## Hand-compute paths b/c': y ~ m + x
  fit_b_hand <- lm(y_vec ~ m_vec + x_vec)
  b_hand  <- coef(fit_b_hand)["m_vec"]
  cp_hand <- coef(fit_b_hand)["x_vec"]

  ## Hand-compute total effect: y ~ x
  fit_c_hand <- lm(y_vec ~ x_vec)
  c_hand <- coef(fit_c_hand)["x_vec"]

  ## Hand-compute Sobel
  sa <- summary(fit_a_hand)$coefficients["x_vec", "Std. Error"]
  sb <- summary(fit_b_hand)$coefficients["m_vec", "Std. Error"]
  ind_hand <- unname(a_hand * b_hand)
  se_ind <- sqrt(unname(b_hand)^2 * sa^2 + unname(a_hand)^2 * sb^2)
  sobel_z_hand <- ind_hand / se_ind
  sobel_p_hand <- 2 * pnorm(-abs(sobel_z_hand))

  ## Now run the model
  env <- build_rena_roi(E_m, items, reps_per_item = 2L)

  md_des <- repmed_design(items = items, X_rdm = X_rdm, Y_rdm = Y_rdm)
  mspec <- repmed_model(
    dataset    = env$ds$dataset,
    design     = env$mvdes,
    repmed_des = md_des,
    key_var    = ~ .rownum,
    distfun    = "euclidean"
  )

  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  perf <- res$performance[[1]]

  expect_false(res$error)

  expect_equal(unname(perf["med_a"]),        unname(a_hand),  tolerance = 1e-4)
  expect_equal(unname(perf["med_b"]),        unname(b_hand),  tolerance = 1e-4)
  expect_equal(unname(perf["med_cprime"]),   unname(cp_hand), tolerance = 1e-4)
  expect_equal(unname(perf["med_c"]),        unname(c_hand),  tolerance = 1e-4)
  expect_equal(unname(perf["med_indirect"]), ind_hand,        tolerance = 1e-4)
  expect_equal(unname(perf["med_sobel_z"]),  sobel_z_hand,    tolerance = 1e-4)
  expect_equal(unname(perf["med_sobel_p"]),  sobel_p_hand,    tolerance = 1e-4)
})


test_that("repmed: zero mediation when mediator is independent of X and Y", {
  set.seed(52)
  K <- 20
  items <- paste0("it", seq_len(K))

  z_x <- rnorm(K)
  z_y <- 0.9 * z_x + rnorm(K, sd = 0.01)
  z_m <- rnorm(K)  # independent

  X_rdm <- as.matrix(dist(z_x))
  Y_rdm <- as.matrix(dist(z_y))
  rownames(X_rdm) <- colnames(X_rdm) <- items
  rownames(Y_rdm) <- colnames(Y_rdm) <- items

  E_m <- matrix(z_m, nrow = K, ncol = NVOX) +
    matrix(rnorm(K * NVOX, sd = 0.001), nrow = K, ncol = NVOX)
  rownames(E_m) <- items

  env <- build_rena_roi(E_m, items, reps_per_item = 2L)

  md_des <- repmed_design(items = items, X_rdm = X_rdm, Y_rdm = Y_rdm)
  mspec <- repmed_model(
    dataset    = env$ds$dataset,
    design     = env$mvdes,
    repmed_des = md_des,
    key_var    = ~ .rownum,
    distfun    = "euclidean"
  )

  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  perf <- res$performance[[1]]

  expect_false(res$error)

  ## Path a should be much smaller than in the mediation case (mediator doesn't

  ## systematically track X). With K=20 RDM pairs the variance of the null
  ## distribution is substantial, so we use a generous threshold.
  expect_lt(abs(unname(perf["med_a"])), 0.7)
  ## Indirect effect should be small (no systematic mediation)
  expect_lt(abs(unname(perf["med_indirect"])), 0.3)
  ## Total effect c should still be large (X→Y direct relationship exists)
  expect_gt(abs(unname(perf["med_c"])), 0.5)
})


test_that("repmed: confound RDMs correctly reduce spurious mediation", {
  set.seed(53)
  K <- 24
  items <- paste0("it", seq_len(K))

  ## Confound drives everything
  z_conf <- rep(c(-1, 0, 1), length.out = K)
  z_x <- z_conf + rnorm(K, sd = 0.01)
  z_m <- z_conf + rnorm(K, sd = 0.01)
  z_y <- z_conf + rnorm(K, sd = 0.01)

  X_rdm <- as.matrix(dist(z_x))
  Y_rdm <- as.matrix(dist(z_y))
  D_conf <- as.matrix(dist(z_conf))
  rownames(X_rdm) <- colnames(X_rdm) <- items
  rownames(Y_rdm) <- colnames(Y_rdm) <- items
  rownames(D_conf) <- colnames(D_conf) <- items

  E_m <- matrix(z_m, nrow = K, ncol = NVOX) +
    matrix(rnorm(K * NVOX, sd = 0.001), nrow = K, ncol = NVOX)
  rownames(E_m) <- items

  env <- build_rena_roi(E_m, items, reps_per_item = 2L)

  ## Without confound
  md_des_no <- repmed_design(items = items, X_rdm = X_rdm, Y_rdm = Y_rdm)
  mspec_no <- repmed_model(
    dataset    = env$ds$dataset,
    design     = env$mvdes,
    repmed_des = md_des_no,
    key_var    = ~ .rownum,
    distfun    = "euclidean"
  )
  res_no <- process_roi(mspec_no, env$roi, rnum = 1L, center_global_id = NA)
  perf_no <- res_no$performance[[1]]

  ## With confound
  md_des_c <- repmed_design(items = items, X_rdm = X_rdm, Y_rdm = Y_rdm,
                            confound_rdms = list(conf = D_conf))
  mspec_c <- repmed_model(
    dataset    = env$ds$dataset,
    design     = env$mvdes,
    repmed_des = md_des_c,
    key_var    = ~ .rownum,
    distfun    = "euclidean"
  )
  res_c <- process_roi(mspec_c, env$roi, rnum = 1L, center_global_id = NA)
  perf_c <- res_c$performance[[1]]

  expect_false(res_no$error)
  expect_false(res_c$error)

  ## Without confound: strong spurious mediation
  expect_gt(abs(unname(perf_no["med_indirect"])), 0.1)
  ## With confound: mediation collapses
  expect_lt(abs(unname(perf_c["med_indirect"])),
            abs(unname(perf_no["med_indirect"])))
})


test_that("repmed: Sobel formula consistency: ind = a*b, p = 2*pnorm(-|z|)", {
  set.seed(54)
  K <- 20
  items <- paste0("it", seq_len(K))

  z_x <- rnorm(K)
  z_m <- 0.8 * z_x + rnorm(K, sd = 0.1)
  z_y <- 0.6 * z_m + rnorm(K, sd = 0.1)

  X_rdm <- as.matrix(dist(z_x))
  Y_rdm <- as.matrix(dist(z_y))
  rownames(X_rdm) <- colnames(X_rdm) <- items
  rownames(Y_rdm) <- colnames(Y_rdm) <- items

  E_m <- matrix(z_m, nrow = K, ncol = NVOX) +
    matrix(rnorm(K * NVOX, sd = 0.001), nrow = K, ncol = NVOX)
  rownames(E_m) <- items

  env <- build_rena_roi(E_m, items, reps_per_item = 2L)

  md_des <- repmed_design(items = items, X_rdm = X_rdm, Y_rdm = Y_rdm)
  mspec <- repmed_model(
    dataset    = env$ds$dataset,
    design     = env$mvdes,
    repmed_des = md_des,
    key_var    = ~ .rownum,
    distfun    = "euclidean"
  )

  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  perf <- res$performance[[1]]

  expect_false(res$error)

  a  <- unname(perf["med_a"])
  b  <- unname(perf["med_b"])
  ind <- unname(perf["med_indirect"])
  sz <- unname(perf["med_sobel_z"])
  sp <- unname(perf["med_sobel_p"])

  ## ind = a * b exactly
  expect_equal(ind, a * b, tolerance = 1e-12)
  ## p = 2 * pnorm(-|z|)
  expect_equal(sp, 2 * pnorm(-abs(sz)), tolerance = 1e-12)
  ## z and indirect have the same sign
  expect_equal(sign(sz), sign(ind))
})


## ===========================================================================
## SECTION 3: repmap_model analytical tests
## ===========================================================================

test_that("repmap: noiseless rank-1 data yields R² = 1.0", {
  skip_if_not_installed("rrpack")
  set.seed(60)

  K <- 20
  P <- 5
  items <- paste0("it", seq_len(K))

  X <- matrix(rnorm(K * P), nrow = K, ncol = P)
  rownames(X) <- items

  ## Rank-1 mapping: Y = X %*% B with B = b * w^T
  b <- rnorm(P)
  w <- rnorm(NVOX)
  B <- b %*% t(w)  # P x NVOX, rank 1
  Y <- X %*% B     # K x NVOX, noiseless
  rownames(Y) <- items

  env <- build_rena_roi(Y, items, reps_per_item = 2L)

  mp_des <- repmap_design(items = items, seed_features = X)
  mspec <- repmap_model(
    dataset    = env$ds$dataset,
    design     = env$mvdes,
    repmap_des = mp_des,
    key_var    = ~ .rownum,
    rank       = 1,
    max_rank   = 5
  )

  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  perf <- res$performance[[1]]

  expect_false(res$error)
  expect_equal(unname(perf["map_r2_mean"]), 1.0, tolerance = 0.01)
  expect_equal(as.integer(perf["map_rank"]), 1L)
})


test_that("repmap: Frobenius norm matches hand-computed sqrt(sum(C^2))", {
  skip_if_not_installed("rrpack")
  set.seed(61)

  K <- 25
  P <- 8
  items <- paste0("it", seq_len(K))

  X <- matrix(rnorm(K * P), nrow = K, ncol = P)
  rownames(X) <- items

  B <- matrix(rnorm(P * NVOX), nrow = P, ncol = NVOX)
  Y <- X %*% B + matrix(rnorm(K * NVOX, sd = 0.1), K, NVOX)
  rownames(Y) <- items

  ## Fit using the internal helper directly
  Xc <- base::scale(X, center = TRUE, scale = FALSE)
  Yc <- base::scale(Y, center = TRUE, scale = FALSE)
  fit <- rMVPA:::.repmap_fit_rrr(Xc, Yc, rank = 3, max_rank = 5)

  expected_frob <- sqrt(sum(fit$C^2))

  ## Now run through process_roi
  env <- build_rena_roi(Y, items, reps_per_item = 2L)

  mp_des <- repmap_design(items = items, seed_features = X)
  mspec <- repmap_model(
    dataset    = env$ds$dataset,
    design     = env$mvdes,
    repmap_des = mp_des,
    key_var    = ~ .rownum,
    rank       = 3,
    max_rank   = 5
  )

  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  perf <- res$performance[[1]]

  expect_false(res$error)
  expect_equal(unname(perf["map_frob_norm"]), expected_frob, tolerance = 0.01)
})


test_that("repmap: rank=0 yields zero R² and zero Frobenius norm", {
  set.seed(62)

  K <- 15
  P <- 6
  items <- paste0("it", seq_len(K))

  X <- matrix(rnorm(K * P), nrow = K, ncol = P)
  rownames(X) <- items
  Y <- matrix(rnorm(K * NVOX), nrow = K, ncol = NVOX)
  rownames(Y) <- items

  env <- build_rena_roi(Y, items, reps_per_item = 2L)

  mp_des <- repmap_design(items = items, seed_features = X)
  mspec <- repmap_model(
    dataset    = env$ds$dataset,
    design     = env$mvdes,
    repmap_des = mp_des,
    key_var    = ~ .rownum,
    rank       = 0
  )

  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  perf <- res$performance[[1]]

  expect_false(res$error)
  expect_equal(as.integer(perf["map_rank"]), 0L)
  expect_equal(unname(perf["map_frob_norm"]), 0.0, tolerance = 1e-15)
  expect_equal(unname(perf["map_r2_mean"]), 0.0, tolerance = 1e-10)
})


test_that("repmap: R² matches hand-computed 1 - SS_res/SS_tot", {
  skip_if_not_installed("rrpack")
  set.seed(63)

  K <- 30
  P <- 8
  items <- paste0("it", seq_len(K))

  X <- matrix(rnorm(K * P), nrow = K, ncol = P)
  rownames(X) <- items

  B <- matrix(rnorm(P * NVOX), nrow = P, ncol = NVOX)
  Y <- X %*% B + matrix(rnorm(K * NVOX, sd = 0.5), K, NVOX)
  rownames(Y) <- items

  ## Hand-compute R²
  Xc <- base::scale(X, center = TRUE, scale = FALSE)
  Yc <- base::scale(Y, center = TRUE, scale = FALSE)
  fit <- rMVPA:::.repmap_fit_rrr(Xc, Yc, rank = 4, max_rank = 8)

  Y_hat <- Xc %*% fit$C
  resid <- Yc - Y_hat
  ss_tot <- colSums(Yc^2)
  ss_res <- colSums(resid^2)
  ss_tot[ss_tot < .Machine$double.eps] <- .Machine$double.eps
  r2_vox <- 1 - ss_res / ss_tot
  expected_r2_mean <- mean(r2_vox)
  expected_r2_med  <- median(r2_vox)

  ## Run through process_roi
  env <- build_rena_roi(Y, items, reps_per_item = 2L)

  mp_des <- repmap_design(items = items, seed_features = X)
  mspec <- repmap_model(
    dataset    = env$ds$dataset,
    design     = env$mvdes,
    repmap_des = mp_des,
    key_var    = ~ .rownum,
    rank       = 4,
    max_rank   = 8
  )

  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  perf <- res$performance[[1]]

  expect_false(res$error)
  expect_equal(unname(perf["map_r2_mean"]),   expected_r2_mean, tolerance = 0.01)
  expect_equal(unname(perf["map_r2_median"]), expected_r2_med,  tolerance = 0.01)
})


## ===========================================================================
## SECTION 4: Edge cases and degenerate inputs
## ===========================================================================

test_that("repnet: constant voxels produce error", {
  set.seed(70)
  K <- 6
  items <- paste0("it", seq_len(K))

  E <- matrix(5.0, nrow = K, ncol = NVOX)
  rownames(E) <- items

  D_seed <- as.matrix(dist(seq_len(K)))
  rownames(D_seed) <- colnames(D_seed) <- items

  env <- build_rena_roi(E, items, reps_per_item = 2L)

  rn_des <- repnet_design(
    design    = env$mvdes,
    key_var   = ~ .rownum,
    seed_rdm  = D_seed
  )
  mspec <- repnet_model(
    dataset    = env$ds$dataset,
    design     = env$mvdes,
    repnet_des = rn_des,
    distfun    = cordist("pearson"),
    simfun     = "pearson"
  )

  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  expect_true(res$error)
})


test_that("repmed: constant mediator produces error", {
  set.seed(71)
  K <- 6
  items <- paste0("it", seq_len(K))

  E_m <- matrix(3.0, nrow = K, ncol = NVOX)
  rownames(E_m) <- items

  X_rdm <- as.matrix(dist(seq_len(K)))
  Y_rdm <- as.matrix(dist(rev(seq_len(K))))
  rownames(X_rdm) <- colnames(X_rdm) <- items
  rownames(Y_rdm) <- colnames(Y_rdm) <- items

  env <- build_rena_roi(E_m, items, reps_per_item = 2L)

  md_des <- repmed_design(items = items, X_rdm = X_rdm, Y_rdm = Y_rdm)
  mspec <- repmed_model(
    dataset    = env$ds$dataset,
    design     = env$mvdes,
    repmed_des = md_des,
    key_var    = ~ .rownum,
    distfun    = "euclidean"
  )

  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  expect_true(res$error)
})


test_that("repnet: fewer than 3 common items produces error", {
  set.seed(72)
  K <- 6
  items <- paste0("it", seq_len(K))

  E <- matrix(rnorm(K * NVOX), nrow = K, ncol = NVOX)
  rownames(E) <- items

  seed_items <- c("it1", "it2")
  D_seed <- matrix(c(0, 1, 1, 0), 2, 2)
  rownames(D_seed) <- colnames(D_seed) <- seed_items

  env <- build_rena_roi(E, items, reps_per_item = 2L)

  rn_des <- repnet_design(
    design    = env$mvdes,
    key_var   = ~ .rownum,
    seed_rdm  = D_seed
  )
  mspec <- repnet_model(
    dataset    = env$ds$dataset,
    design     = env$mvdes,
    repnet_des = rn_des,
    distfun    = cordist("pearson"),
    simfun     = "pearson"
  )

  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  expect_true(res$error)
  expect_match(res$error_message, "at least 3 common items")
})


## ===========================================================================
## SECTION 5: group_means correctness (used by all ReNA models)
## ===========================================================================

test_that("group_means correctly averages repeated observations per item", {
  set.seed(80)
  K <- 5
  V <- 4
  reps <- 3

  items <- paste0("it", seq_len(K))
  key <- rep(items, each = reps)

  E_true <- matrix(seq_len(K * V), nrow = K, ncol = V)
  rownames(E_true) <- items

  X <- E_true[key, , drop = FALSE]

  E_recovered <- group_means(X, margin = 1, group = key)
  E_recovered <- E_recovered[items, , drop = FALSE]

  expect_equal(E_recovered, E_true, tolerance = 1e-15)
})


test_that("group_means with noisy replicates converges to true mean", {
  set.seed(81)
  K <- 5
  V <- 4
  reps <- 100

  items <- paste0("it", seq_len(K))
  key <- rep(items, each = reps)

  E_true <- matrix(rnorm(K * V), nrow = K, ncol = V)
  rownames(E_true) <- items

  noise <- matrix(rnorm(K * reps * V, sd = 0.1), nrow = K * reps, ncol = V)
  X <- E_true[key, , drop = FALSE] + noise

  E_recovered <- group_means(X, margin = 1, group = key)
  E_recovered <- E_recovered[items, , drop = FALSE]

  expect_equal(E_recovered, E_true, tolerance = 0.05)
})
