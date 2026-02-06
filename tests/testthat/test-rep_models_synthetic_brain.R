context("Synthetic brain: joint behavior of ReNA models")

library(testthat)
library(rMVPA)

set.seed(123)

test_that("synthetic multi-ROI system shows expected ReNA patterns (D)", {
  skip_on_cran()
  skip_if_not_installed("rrpack")

  ## Shared item set
  K <- 24
  items <- paste0("it", seq_len(K))

  ## Seed A for connectivity/mapping
  d_seed <- 10
  F_A <- matrix(rnorm(K * d_seed), nrow = K, ncol = d_seed)
  rownames(F_A) <- items
  D_A <- rMVPA:::pairwise_dist(cordist("pearson"), F_A)
  rownames(D_A) <- colnames(D_A) <- items

  ## Predictor/Outcome RDMs for mediation
  sim_med <- simulate_repmed_latent(K = K, a = 1.0, b = 1.0, c = 0.3)

  ## Block confound for ROI 3
  blocks <- 3
  block_labels <- rep(seq_len(blocks), length.out = K)
  names(block_labels) <- items
  D_conf <- outer(block_labels, block_labels, FUN = function(a, b) as.numeric(a != b))
  rownames(D_conf) <- colnames(D_conf) <- items

  ## Base dataset/design reused across ROIs
  ds <- gen_sample_dataset(D = c(4, 4, 4),
                           nobs = K * 2L,
                           blocks = 2,
                           nlevels = 2)

  mvdes <- mvpa_design(
    train_design = ds$design$train_design,
    y_train      = ds$design$y_train,
    block_var    = ds$design$block_var
  )
  mvdes$train_design$.rownum <- rep(items, each = 2L)

  ## ROI construction helper
  build_roi_from_item_patterns <- function(item_patterns) {
    vox <- which(ds$dataset$mask > 0)
    n_vox <- length(vox)
    samp <- data_sample(ds$dataset, vox)
    roi  <- as_roi(samp, ds$dataset)

    E <- item_patterns
    if (ncol(E) != n_vox) {
      common <- min(ncol(E), n_vox)
      E_full <- matrix(0, nrow = K, ncol = n_vox)
      E_full[, seq_len(common)] <- E[, seq_len(common), drop = FALSE]
      if (n_vox > common) {
        E_full[, (common + 1L):n_vox] <- matrix(rnorm(K * (n_vox - common)), nrow = K) * 0.01
      }
      rownames(E_full) <- items
      E <- E_full
    }

    X_rep <- E[as.character(mvdes$train_design$.rownum), , drop = FALSE]
    roi$train_roi@.Data <- X_rep
    roi
  }

  ## ROI 1: shares geometry with seed A and maps features from A
  n_vox <- length(which(ds$dataset$mask > 0))
  W1 <- matrix(rnorm(d_seed * n_vox), nrow = d_seed, ncol = n_vox)
  roi1_items <- F_A %*% W1 + matrix(rnorm(K * n_vox), nrow = K, ncol = n_vox) / 5
  rownames(roi1_items) <- items
  roi1 <- build_roi_from_item_patterns(roi1_items)

  rn_des_1 <- repnet_design(
    design        = mvdes,
    key_var       = ~ .rownum,
    seed_rdm      = D_A,
    confound_rdms = NULL
  )
  rn_spec_1 <- repnet_model(
    dataset    = ds$dataset,
    design     = mvdes,
    repnet_des = rn_des_1,
    distfun    = cordist("pearson"),
    simfun     = "pearson"
  )
  res_net_1 <- process_roi(rn_spec_1, roi1, rnum = 1L, center_global_id = NA)
  conn_roi1 <- unname(res_net_1$performance[[1]]["conn_raw"])

  ## Mapping from seed features to ROI1
  Xc_A <- base::scale(F_A, center = TRUE, scale = FALSE)
  Yc_1 <- base::scale(roi1_items, center = TRUE, scale = FALSE)
  fit_map_1 <- rMVPA:::.repmap_fit_rrr(
    Xc_A,
    Yc_1,
    rank     = 3,
    max_rank = 3
  )
  Y_hat_1 <- Xc_A %*% fit_map_1$C
  ss_tot_1 <- colSums(Yc_1^2)
  ss_tot_1[ss_tot_1 < .Machine$double.eps] <- .Machine$double.eps
  r2_map_1 <- mean(1 - colSums((Yc_1 - Y_hat_1)^2) / ss_tot_1)

  ## ROI 2: mediates between X and Y (strong mediation)
  roi2_items <- matrix(sim_med$z_m, nrow = K,
                       ncol = n_vox, byrow = FALSE) +
    matrix(rnorm(K * n_vox, sd = 0.1), nrow = K, ncol = n_vox)
  rownames(roi2_items) <- items
  roi2 <- build_roi_from_item_patterns(roi2_items)

  md_des <- repmed_design(
    items  = sim_med$items,
    X_rdm  = sim_med$X_rdm,
    Y_rdm  = sim_med$Y_rdm
  )
  md_spec <- repmed_model(
    dataset    = ds$dataset,
    design     = mvdes,
    repmed_des = md_des,
    key_var    = ~ .rownum,
    distfun    = "euclidean"
  )
  res_med_2 <- process_roi(md_spec, roi2, rnum = 1L, center_global_id = NA)
  ind_roi2 <- unname(res_med_2$performance[[1]]["med_indirect"])

  ## ROI 3: pure confound/block region
  block_patterns <- matrix(rnorm(blocks * n_vox), nrow = blocks, ncol = n_vox)
  roi3_items <- block_patterns[block_labels, ] +
    matrix(rnorm(K * n_vox, sd = 0.1), nrow = K, ncol = n_vox)
  rownames(roi3_items) <- items
  roi3 <- build_roi_from_item_patterns(roi3_items)

  ## Spurious mediation without confound
  md_des_no <- repmed_design(
    items = items,
    X_rdm = D_conf,
    Y_rdm = D_conf
  )
  md_spec_no <- repmed_model(
    dataset    = ds$dataset,
    design     = mvdes,
    repmed_des = md_des_no,
    key_var    = ~ .rownum,
    distfun    = "euclidean"
  )
  res_med_3_no <- process_roi(md_spec_no, roi3, rnum = 1L, center_global_id = NA)
  ind_roi3_no <- unname(res_med_3_no$performance[[1]]["med_indirect"])

  ## With confound, mediation should attenuate
  md_des_conf <- repmed_design(
    items         = items,
    X_rdm         = D_conf,
    Y_rdm         = D_conf,
    confound_rdms = list(block = D_conf)
  )
  md_spec_conf <- repmed_model(
    dataset    = ds$dataset,
    design     = mvdes,
    repmed_des = md_des_conf,
    key_var    = ~ .rownum,
    distfun    = "euclidean"
  )
  res_med_3_conf <- process_roi(md_spec_conf, roi3, rnum = 1L, center_global_id = NA)
  ind_roi3_conf <- unname(res_med_3_conf$performance[[1]]["med_indirect"])

  ## ROI 4: unrelated noise region
  roi4_items <- matrix(rnorm(K * n_vox), nrow = K, ncol = n_vox)
  rownames(roi4_items) <- items
  roi4 <- build_roi_from_item_patterns(roi4_items)

  rn_spec_4 <- rn_spec_1
  res_net_4 <- process_roi(rn_spec_4, roi4, rnum = 1L, center_global_id = NA)
  conn_roi4 <- unname(res_net_4$performance[[1]]["conn_raw"])

  md_spec_4 <- md_spec
  res_med_4 <- process_roi(md_spec_4, roi4, rnum = 1L, center_global_id = NA)
  ind_roi4 <- unname(res_med_4$performance[[1]]["med_indirect"])

  X_noise <- matrix(rnorm(K * d_seed), nrow = K, ncol = d_seed)
  Xc_noise <- base::scale(X_noise, center = TRUE, scale = FALSE)
  Yc_4 <- base::scale(roi4_items, center = TRUE, scale = FALSE)
  fit_map_4 <- rMVPA:::.repmap_fit_rrr(
    Xc_noise,
    Yc_4,
    rank     = 3,
    max_rank = 3
  )
  Y_hat_4 <- Xc_noise %*% fit_map_4$C
  ss_tot_4 <- colSums(Yc_4^2)
  ss_tot_4[ss_tot_4 < .Machine$double.eps] <- .Machine$double.eps
  r2_map_4 <- mean(1 - colSums((Yc_4 - Y_hat_4)^2) / ss_tot_4)

  ## Assertions: pattern across ROIs
  expect_gt(conn_roi1, 0.4)
  expect_lt(abs(conn_roi4), 0.2)

  expect_gt(r2_map_1, 0.2)
  expect_lt(abs(r2_map_4), 0.25)

  expect_gt(ind_roi2, 0.1)              ## strong mediation in ROI 2
  ## Confound-driven mediation in ROI 3 is exercised in dedicated B4 tests
  expect_lt(abs(ind_roi4), 0.1)         ## noise ROI quiet
})
