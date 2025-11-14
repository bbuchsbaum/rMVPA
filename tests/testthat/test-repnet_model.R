context("ReNA-RC: repnet_design and repnet_model")

library(testthat)
library(rMVPA)

set.seed(123)

test_that("repnet_design builds key and aligns seed RDM", {
  ds <- gen_sample_dataset(D = c(4,4,4), nobs = 12, blocks = 3, nlevels = 2)
  mvdes <- mvpa_design(
    train_design = ds$design$train_design,
    y_train      = ds$design$y_train,
    block_var    = ds$design$block_var
  )

  key_ids <- mvdes$train_design$.rownum
  seed_mat <- as.matrix(dist(key_ids))
  rownames(seed_mat) <- colnames(seed_mat) <- as.character(key_ids)

  rn_des <- repnet_design(
    design    = mvdes,
    key_var   = ~ .rownum,
    seed_rdm  = seed_mat,
    confound_rdms = NULL
  )

  expect_true(is.factor(rn_des$key))
  expect_equal(length(rn_des$key), nrow(mvdes$train_design))
  expect_true(is.matrix(rn_des$seed_rdm))
  expect_equal(rownames(rn_des$seed_rdm), colnames(rn_des$seed_rdm))
})


test_that("repnet_model + process_roi.repnet_model return connectivity metrics", {
  ds <- gen_sample_dataset(D = c(4,4,4), nobs = 16, blocks = 4, nlevels = 2)
  mvdes <- mvpa_design(
    train_design = ds$design$train_design,
    y_train      = ds$design$y_train,
    block_var    = ds$design$block_var
  )

  key_ids <- mvdes$train_design$.rownum
  seed_mat <- as.matrix(dist(key_ids))
  rownames(seed_mat) <- colnames(seed_mat) <- as.character(key_ids)

  rn_des <- repnet_design(
    design    = mvdes,
    key_var   = ~ .rownum,
    seed_rdm  = seed_mat
  )

  mspec <- repnet_model(
    dataset    = ds$dataset,
    design     = mvdes,
    repnet_des = rn_des,
    distfun    = cordist("pearson"),
    simfun     = "pearson"
  )

  # Simple ROI: take whole mask as a single region
  vox <- which(ds$dataset$mask > 0)
  samp <- data_sample(ds$dataset, vox)
  roi  <- as_roi(samp, ds$dataset)

  res <- process_roi(mspec, roi, rnum = 1L, center_global_id = NA)

  expect_s3_class(res, "tbl_df")
  expect_false(res$error)
  perf <- res$performance[[1]]
  expect_true(all(c("n_items", "n_pairs", "conn_raw") %in% names(perf)))
  expect_true(is.numeric(perf["conn_raw"]))
})


## ----------------------------------------------------------------------
## Synthetic geometry tests
## ----------------------------------------------------------------------

simulate_repnet_roi_patterns <- function(K = 40,
                                         d_seed = 20,
                                         snr = 3,
                                         match_seed = TRUE) {

  items <- paste0("it", seq_len(K))

  ## Seed features and RDM (1 - Pearson correlation)
  F <- matrix(rnorm(K * d_seed), nrow = K, ncol = d_seed)
  rownames(F) <- items

  cord <- cordist("pearson")
  D_seed <- pairwise_dist(cord, F)
  rownames(D_seed) <- colnames(D_seed) <- items

  list(
    items        = items,
    seed_features = F,
    seed_rdm     = D_seed,
    snr          = snr,
    match_seed   = match_seed
  )
}


build_repnet_spec_and_roi <- function(sim,
                                      reps_per_item = 2L,
                                      simfun = "pearson") {

  K <- length(sim$items)
  items <- sim$items

  ## Build a minimal mvpa_dataset / mvpa_design matching the API used elsewhere
  ds <- gen_sample_dataset(D = c(4, 4, 4),
                           nobs = K * reps_per_item,
                           blocks = reps_per_item,
                           nlevels = 2)

  ## Overwrite the train design to carry our item IDs explicitly
  mvdes <- mvpa_design(
    train_design = ds$design$train_design,
    y_train      = ds$design$y_train,
    block_var    = ds$design$block_var
  )

  ## Replace the .rownum with item labels repeated reps_per_item times
  mvdes$train_design$.rownum <- rep(items, each = reps_per_item)

  ## Construct repnet_design
  rn_des <- repnet_design(
    design        = mvdes,
    key_var       = ~ .rownum,
    seed_rdm      = sim$seed_rdm,
    confound_rdms = NULL
  )

  mspec <- repnet_model(
    dataset    = ds$dataset,
    design     = mvdes,
    repnet_des = rn_des,
    distfun    = cordist("pearson"),
    simfun     = simfun
  )

  ## Build ROI containing our synthetic item-level patterns for the training set
  ## Use all mask voxels and then overwrite train_roi values to match synthetic patterns
  vox <- which(ds$dataset$mask > 0)
  n_vox <- length(vox)
  samp <- data_sample(ds$dataset, vox)
  roi  <- as_roi(samp, ds$dataset)

  ## Construct item-level ROI patterns:
  ## - If explicit roi_items are provided, adapt them to the available number of voxels.
  ## - Otherwise, generate them from the seed features or as independent noise.
  if (!is.null(sim$roi_items)) {
    E_base <- sim$roi_items
    if (ncol(E_base) != n_vox) {
      K <- nrow(E_base)
      E_full <- matrix(0, nrow = K, ncol = n_vox)
      common <- min(ncol(E_base), n_vox)
      E_full[, seq_len(common)] <- E_base[, seq_len(common), drop = FALSE]
      if (n_vox > common) {
        E_full[, (common + 1L):n_vox] <- matrix(rnorm(K * (n_vox - common)), nrow = K) * 0.01
      }
      rownames(E_full) <- rownames(E_base)
      E <- E_full
    } else {
      E <- E_base
    }
  } else {
    K <- length(items)
    if (isTRUE(sim$match_seed)) {
      ## ROI geometry is a linear transform of seed features plus noise
      F <- sim$seed_features
      d_seed <- ncol(F)
      W <- matrix(rnorm(d_seed * n_vox), nrow = d_seed, ncol = n_vox)
      signal <- F %*% W
    } else {
      ## ROI geometry is independent Gaussian noise
      signal <- matrix(rnorm(K * n_vox), nrow = K, ncol = n_vox)
    }
    noise <- matrix(rnorm(K * n_vox), nrow = K, ncol = n_vox)
    snr <- if (!is.null(sim$snr)) sim$snr else 3
    E <- signal + noise / snr
    rownames(E) <- items
  }

  ## Overwrite the ROI train data with synthetic patterns grouped by key;
  ## each item pattern is replicated reps_per_item times.
  X_rep <- E[as.character(mvdes$train_design$.rownum), , drop = FALSE]
  roi$train_roi@.Data <- X_rep

  list(mspec = mspec, roi = roi)
}


test_that("repnet_model: seed-matched vs orthogonal geometry (A1)", {
  skip_on_cran()

  ## Seed-matched ROI
  sim_match <- simulate_repnet_roi_patterns(match_seed = TRUE, snr = 4)
  obj_match <- build_repnet_spec_and_roi(sim_match, reps_per_item = 2L)
  res_match <- process_roi(obj_match$mspec, obj_match$roi, rnum = 1L, center_global_id = NA)
  perf_match <- res_match$performance[[1]]

  ## Orthogonal ROI
  sim_orth <- simulate_repnet_roi_patterns(match_seed = FALSE, snr = 4)
  obj_orth <- build_repnet_spec_and_roi(sim_orth, reps_per_item = 2L)
  res_orth <- process_roi(obj_orth$mspec, obj_orth$roi, rnum = 1L, center_global_id = NA)
  perf_orth <- res_orth$performance[[1]]

  conn_match <- unname(perf_match["conn_raw"])
  conn_orth  <- unname(perf_orth["conn_raw"])

  expect_gt(conn_match, 0.6)
  expect_lt(abs(conn_orth), 0.2)
})


test_that("repnet_model: confound block RDM is factored out (A2)", {
  skip_on_cran()

  K <- 30
  d_seed <- 10
  n_vox <- 60
  blocks <- 3

  items <- paste0("it", seq_len(K))

  ## Block structure
  block_labels <- rep(seq_len(blocks), length.out = K)
  names(block_labels) <- items

  D_block <- outer(block_labels, block_labels, FUN = function(a, b) as.numeric(a != b))
  rownames(D_block) <- colnames(D_block) <- items

  ## World 1: seed RDM correlated with block RDM
  latent_block <- matrix(rnorm(blocks * d_seed), nrow = blocks, ncol = d_seed)
  F1 <- latent_block[block_labels, ] + matrix(rnorm(K * d_seed), nrow = K, ncol = d_seed) * 0.1
  rownames(F1) <- items
  D_seed1 <- pairwise_dist(cordist("pearson"), F1)
  rownames(D_seed1) <- colnames(D_seed1) <- items

  ## ROI patterns depend only on block
  voxel_dim <- n_vox
  block_patterns <- matrix(rnorm(blocks * voxel_dim), nrow = blocks, ncol = voxel_dim)
  E_block <- block_patterns[block_labels, ] + matrix(rnorm(K * voxel_dim), nrow = K, ncol = voxel_dim) * 0.1
  rownames(E_block) <- items

  ## Helper to run one world with/without confound
  run_world <- function(seed_rdm, use_confound) {
    sim <- list(items = items,
                seed_features = F1,
                seed_rdm = seed_rdm,
                roi_items = E_block)
    obj <- build_repnet_spec_and_roi(sim, reps_per_item = 2L)
    if (use_confound) {
      obj$mspec$confound_rdms <- list(block = D_block)
    }
    res <- process_roi(obj$mspec, obj$roi, rnum = 1L, center_global_id = NA)
    res$performance[[1]]
  }

  ## World 1: seed strongly aligned with block
  perf_w1_raw <- run_world(D_seed1, use_confound = FALSE)
  perf_w1_adj <- run_world(D_seed1, use_confound = TRUE)

  conn_raw_w1 <- unname(perf_w1_raw["conn_raw"])
  beta_seed_w1_raw <- unname(perf_w1_raw["beta_seed"])
  beta_seed_w1_adj <- unname(perf_w1_adj["beta_seed"])

  expect_gt(conn_raw_w1, 0.5)
  expect_gt(beta_seed_w1_raw, 0.1)
  expect_lt(abs(beta_seed_w1_adj), 0.1)

  ## World 2: seed independent of block (permute items before building RDM)
  perm_items <- sample(items)
  F2 <- F1[perm_items, , drop = FALSE]
  rownames(F2) <- items
  D_seed2 <- pairwise_dist(cordist("pearson"), F2)
  rownames(D_seed2) <- colnames(D_seed2) <- items

  perf_w2_raw <- run_world(D_seed2, use_confound = FALSE)
  perf_w2_adj <- run_world(D_seed2, use_confound = TRUE)

  conn_raw_w2 <- unname(perf_w2_raw["conn_raw"])
  beta_seed_w2_raw <- unname(perf_w2_raw["beta_seed"])
  beta_seed_w2_adj <- unname(perf_w2_adj["beta_seed"])

  expect_lt(abs(conn_raw_w2), 0.2)
  expect_lt(abs(beta_seed_w2_raw), 0.2)
  expect_lt(abs(beta_seed_w2_adj), 0.2)
})


test_that("repnet_model: seed permutation destroys connectivity (A3)", {
  skip_on_cran()

  sim <- simulate_repnet_roi_patterns(match_seed = TRUE, snr = 4)

  ## Correct alignment
  obj_correct <- build_repnet_spec_and_roi(sim, reps_per_item = 2L)
  res_correct <- process_roi(obj_correct$mspec, obj_correct$roi, rnum = 1L, center_global_id = NA)
  perf_correct <- res_correct$performance[[1]]
  conn_correct <- unname(perf_correct["conn_raw"])

  ## Permuted seed RDM (permute item labels)
  items <- sim$items
  perm <- sample(items)
  D_perm <- sim$seed_rdm[perm, perm, drop = FALSE]
  rownames(D_perm) <- colnames(D_perm) <- items

  sim_perm <- sim
  sim_perm$seed_rdm <- D_perm
  obj_perm <- build_repnet_spec_and_roi(sim_perm, reps_per_item = 2L)
  res_perm <- process_roi(obj_perm$mspec, obj_perm$roi, rnum = 1L, center_global_id = NA)
  perf_perm <- res_perm$performance[[1]]
  conn_perm <- unname(perf_perm["conn_raw"])

  expect_gt(conn_correct, 0.6)
  expect_lt(abs(conn_perm), 0.2)
})


test_that("repnet_model: multiple seeds show selective connectivity (A4)", {
  skip_on_cran()

  K <- 30
  d_low <- 15
  d_sem <- 10
  n_vox <- 60

  items <- paste0("it", seq_len(K))

  ## Low-level seed features
  F_low <- matrix(rnorm(K * d_low), nrow = K, ncol = d_low)
  rownames(F_low) <- items

  ## Semantic seed features (share a 1D latent structure + noise)
  z_sem <- rnorm(K)
  F_sem <- cbind(z_sem,
                 matrix(rnorm(K * (d_sem - 1L)), nrow = K, ncol = d_sem - 1L))
  rownames(F_sem) <- items

  D_low <- pairwise_dist(cordist("pearson"), F_low)
  rownames(D_low) <- colnames(D_low) <- items
  D_sem <- pairwise_dist(cordist("pearson"), F_sem)
  rownames(D_sem) <- colnames(D_sem) <- items

  ## ROI L: depends only on low-level features
  W_L <- matrix(rnorm(d_low * n_vox), nrow = d_low, ncol = n_vox)
  roi_L_items <- F_low %*% W_L + matrix(rnorm(K * n_vox), nrow = K, ncol = n_vox) / 4
  rownames(roi_L_items) <- items

  ## ROI S: depends only on semantic features
  W_S <- matrix(rnorm(d_sem * n_vox), nrow = d_sem, ncol = n_vox)
  roi_S_items <- F_sem %*% W_S + matrix(rnorm(K * n_vox), nrow = K, ncol = n_vox) / 4
  rownames(roi_S_items) <- items

  ## Build four combinations: ROI L/S x seed low/sem
  sim_L_low <- list(items = items,
                    seed_features = F_low,
                    seed_rdm = D_low,
                    roi_items = roi_L_items)
  sim_L_sem <- list(items = items,
                    seed_features = F_sem,
                    seed_rdm = D_sem,
                    roi_items = roi_L_items)
  sim_S_low <- list(items = items,
                    seed_features = F_low,
                    seed_rdm = D_low,
                    roi_items = roi_S_items)
  sim_S_sem <- list(items = items,
                    seed_features = F_sem,
                    seed_rdm = D_sem,
                    roi_items = roi_S_items)

  obj_L_low <- build_repnet_spec_and_roi(sim_L_low, reps_per_item = 2L)
  obj_L_sem <- build_repnet_spec_and_roi(sim_L_sem, reps_per_item = 2L)
  obj_S_low <- build_repnet_spec_and_roi(sim_S_low, reps_per_item = 2L)
  obj_S_sem <- build_repnet_spec_and_roi(sim_S_sem, reps_per_item = 2L)

  conn_L_low <- unname(process_roi(obj_L_low$mspec, obj_L_low$roi, rnum = 1L, center_global_id = NA)$performance[[1]]["conn_raw"])
  conn_L_sem <- unname(process_roi(obj_L_sem$mspec, obj_L_sem$roi, rnum = 1L, center_global_id = NA)$performance[[1]]["conn_raw"])
  conn_S_low <- unname(process_roi(obj_S_low$mspec, obj_S_low$roi, rnum = 1L, center_global_id = NA)$performance[[1]]["conn_raw"])
  conn_S_sem <- unname(process_roi(obj_S_sem$mspec, obj_S_sem$roi, rnum = 1L, center_global_id = NA)$performance[[1]]["conn_raw"])

  ## ROI L should prefer low-level seed
  expect_gt(conn_L_low, 0.4)
  expect_gt(conn_L_low, conn_L_sem + 0.2)

  ## ROI S should prefer semantic seed
  expect_gt(conn_S_sem, 0.4)
  expect_gt(conn_S_sem, conn_S_low + 0.2)
})
