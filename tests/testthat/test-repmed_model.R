context("ReNA-RM: repmed_design and repmed_model")

library(testthat)
library(rMVPA)

set.seed(123)

test_that("repmed_design coerces RDMs and stores items", {
  items <- letters[1:5]
  X <- as.matrix(dist(seq_along(items)))
  rownames(X) <- colnames(X) <- items
  Y <- as.matrix(dist(rev(seq_along(items))))
  rownames(Y) <- colnames(Y) <- items

  des <- repmed_design(items = items, X_rdm = X, Y_rdm = Y)

  expect_equal(des$items, items)
  expect_true(is.matrix(des$X_rdm))
  expect_true(is.matrix(des$Y_rdm))
  expect_equal(rownames(des$X_rdm), items)
  expect_equal(colnames(des$X_rdm), items)
})


test_that("repmed_model + process_roi.repmed_model produce mediation metrics", {
  ds <- gen_sample_dataset(D = c(4,4,4), nobs = 18, blocks = 3, nlevels = 2)
  mvdes <- mvpa_design(
    train_design = ds$design$train_design,
    y_train      = ds$design$y_train,
    block_var    = ds$design$block_var
  )

  key_ids <- mvdes$train_design$.rownum
  items <- as.character(sort(unique(key_ids)))
  X <- as.matrix(dist(seq_along(items)))
  rownames(X) <- colnames(X) <- items
  Y <- as.matrix(dist(rev(seq_along(items))))
  rownames(Y) <- colnames(Y) <- items

  md_des <- repmed_design(items = items, X_rdm = X, Y_rdm = Y)

  mspec <- repmed_model(
    dataset    = ds$dataset,
    design     = mvdes,
    repmed_des = md_des,
    key_var    = ~ .rownum,
    distfun    = cordist("pearson")
  )

  vox <- which(ds$dataset$mask > 0)
  samp <- data_sample(ds$dataset, vox)
  roi  <- as_roi(samp, ds$dataset)

  # Suppress "essentially perfect fit" warning - expected with low-noise synthetic mediation data
  res <- suppressWarnings(process_roi(mspec, roi, rnum = 1L, center_global_id = NA))

  expect_s3_class(res, "tbl_df")
  expect_false(res$error)
  perf <- res$performance[[1]]
  expect_true(all(c("n_items", "n_pairs",
                    "med_a", "med_b", "med_cprime", "med_indirect") %in% names(perf)))
})


build_repmed_spec_and_roi <- function(sim,
                                      reps_per_item = 2L,
                                      distfun = "euclidean") {

  K <- length(sim$items)
  items <- sim$items

  ds <- gen_sample_dataset(D = c(4, 4, 4),
                           nobs = K * reps_per_item,
                           blocks = reps_per_item,
                           nlevels = 2)

  mvdes <- mvpa_design(
    train_design = ds$design$train_design,
    y_train      = ds$design$y_train,
    block_var    = ds$design$block_var
  )

  mvdes$train_design$.rownum <- rep(items, each = reps_per_item)

  md_des <- repmed_design(
    items  = items,
    X_rdm  = sim$X_rdm,
    Y_rdm  = sim$Y_rdm
  )

  mspec <- repmed_model(
    dataset    = ds$dataset,
    design     = mvdes,
    repmed_des = md_des,
    key_var    = ~ .rownum,
    distfun    = distfun
  )

  ## Build ROI / mediator patterns driven by z_m
  vox <- which(ds$dataset$mask > 0)
  n_vox <- length(vox)
  samp <- data_sample(ds$dataset, vox)
  roi  <- as_roi(samp, ds$dataset)

  M_items <- matrix(sim$z_m, nrow = K, ncol = n_vox)
  M_items <- M_items + matrix(rnorm(K * n_vox, sd = 0.05), nrow = K)
  rownames(M_items) <- items

  X_rep <- M_items[as.character(mvdes$train_design$.rownum), , drop = FALSE]
  roi$train_roi@.Data <- X_rep

  list(mspec = mspec, roi = roi)
}


test_that("repmed_model: full mediation yields strong indirect effect (B1)", {
  skip_on_cran()

  sim <- simulate_repmed_latent(K = 24, a = 1.0, b = 1.2, c = 0)
  obj <- build_repmed_spec_and_roi(sim, reps_per_item = 2L, distfun = "euclidean")

  res <- process_roi(obj$mspec, obj$roi, rnum = 1L, center_global_id = NA)
  perf <- res$performance[[1]]

  a_hat  <- unname(perf["med_a"])
  b_hat  <- unname(perf["med_b"])
  cp_hat <- unname(perf["med_cprime"])
  ind    <- unname(perf["med_indirect"])
  p_sob  <- unname(perf["med_sobel_p"])

  expect_gt(a_hat, 0.3)
  expect_gt(b_hat, 0.1)
  expect_lt(abs(cp_hat), 0.2)
  expect_gt(ind, 0.1)
  expect_lt(p_sob, 0.05)
})


test_that("repmed_model: partial vs full mediation (B2)", {
  skip_on_cran()

  ## Full mediation scenario
  sim_full <- simulate_repmed_latent(K = 24, a = 1.0, b = 1.0, c = 0)
  obj_full <- build_repmed_spec_and_roi(sim_full, reps_per_item = 2L, distfun = "euclidean")
  res_full <- process_roi(obj_full$mspec, obj_full$roi, rnum = 1L, center_global_id = NA)
  perf_full <- res_full$performance[[1]]

  ## Partial mediation: non-zero direct path c
  sim_part <- simulate_repmed_latent(K = 24, a = 1.0, b = 0.8, c = 0.6)
  obj_part <- build_repmed_spec_and_roi(sim_part, reps_per_item = 2L, distfun = "euclidean")
  res_part <- process_roi(obj_part$mspec, obj_part$roi, rnum = 1L, center_global_id = NA)
  perf_part <- res_part$performance[[1]]

  a_full   <- unname(perf_full["med_a"])
  b_full   <- unname(perf_full["med_b"])
  cp_full  <- unname(perf_full["med_cprime"])
  ind_full <- unname(perf_full["med_indirect"])

  a_part   <- unname(perf_part["med_a"])
  b_part   <- unname(perf_part["med_b"])
  cp_part  <- unname(perf_part["med_cprime"])
  ind_part <- unname(perf_part["med_indirect"])

  expect_gt(a_full, 0.2)
  expect_gt(b_full, 0.05)
  expect_lt(abs(cp_full), 0.6)
  expect_gt(ind_full, 0.1)

  expect_gt(a_part, 0.2)
  expect_gt(b_part, 0.05)
  expect_gt(cp_part, 0.2)        ## direct effect present
  expect_gt(ind_part, 0.05)      ## still mediated
  expect_lt(ind_part - ind_full, 0.5)  ## partial mediation not much stronger than full
})


test_that("repmed_model: no mediation / null case (B3)", {
  skip_on_cran()

  n_rep <- 25
  K <- 16

  a_vals <- numeric(n_rep)
  b_vals <- numeric(n_rep)
  cp_vals <- numeric(n_rep)
  ind_vals <- numeric(n_rep)
  p_vals <- numeric(n_rep)

  for (i in seq_len(n_rep)) {
    items <- paste0("it", seq_len(K))
    z_x <- rnorm(K)
    z_m <- rnorm(K)
    z_y <- rnorm(K)

    X_rdm <- as.matrix(dist(z_x))
    Y_rdm <- as.matrix(dist(z_y))
    rownames(X_rdm) <- colnames(X_rdm) <- items
    rownames(Y_rdm) <- colnames(Y_rdm) <- items

    sim <- list(
      items = items,
      z_x   = z_x,
      z_m   = z_m,
      z_y   = z_y,
      X_rdm = X_rdm,
      Y_rdm = Y_rdm
    )

    obj <- build_repmed_spec_and_roi(sim, reps_per_item = 2L, distfun = "euclidean")
    res <- process_roi(obj$mspec, obj$roi, rnum = 1L, center_global_id = NA)
    perf <- res$performance[[1]]

    a_vals[i]  <- unname(perf["med_a"])
    b_vals[i]  <- unname(perf["med_b"])
    cp_vals[i] <- unname(perf["med_cprime"])
    ind_vals[i] <- unname(perf["med_indirect"])
    p_vals[i]  <- unname(perf["med_sobel_p"])
  }

  ## Effects centered near zero
  expect_lt(abs(mean(a_vals, na.rm = TRUE)), 0.25)
  expect_lt(abs(mean(b_vals, na.rm = TRUE)), 0.2)
  expect_lt(abs(mean(cp_vals, na.rm = TRUE)), 0.2)
  expect_lt(abs(mean(ind_vals, na.rm = TRUE)), 0.1)

  ## Type I error roughly controlled
  prop_sig <- mean(p_vals < 0.05, na.rm = TRUE)
  expect_lt(prop_sig, 0.2)
})


test_that("repmed_model: confound-driven pseudo-mediation collapses with confound control (B4)", {
  skip_on_cran()

  K <- 24
  blocks <- 3
  items <- paste0("it", seq_len(K))

  block_labels <- rep(seq_len(blocks), length.out = K)
  names(block_labels) <- items

  ## Confound RDM: block structure
  D_conf <- outer(block_labels, block_labels, FUN = function(a, b) as.numeric(a != b))
  rownames(D_conf) <- colnames(D_conf) <- items

  ## Predictor / outcome RDMs share block structure plus item-specific signal
  make_rdm_from_latent <- function(z, items) {
    M <- as.matrix(dist(z, method = "euclidean"))
    rownames(M) <- colnames(M) <- items
    M
  }

  scale_rdm <- function(M) {
    v <- M[lower.tri(M)]
    v_sc <- as.numeric(scale(v))
    S <- matrix(0, nrow(M), ncol(M))
    S[lower.tri(S)] <- v_sc
    S <- S + t(S)
    diag(S) <- 0
    rownames(S) <- rownames(M)
    colnames(S) <- colnames(M)
    S
  }

  z_x <- rnorm(K)
  z_y <- rnorm(K)
  X_sig <- make_rdm_from_latent(z_x, items)
  Y_sig <- make_rdm_from_latent(z_y, items)

  X_rdm <- 0.5 * scale_rdm(X_sig) + 0.5 * scale_rdm(D_conf)
  Y_rdm <- 0.5 * scale_rdm(Y_sig) + 0.5 * scale_rdm(D_conf)

  ## Dataset and design
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

  ## Mediator ROI: depends only on block structure
  vox <- which(ds$dataset$mask > 0)
  n_vox <- length(vox)
  samp <- data_sample(ds$dataset, vox)
  roi  <- as_roi(samp, ds$dataset)

  block_patterns <- matrix(rnorm(blocks * n_vox), nrow = blocks, ncol = n_vox)
  M_items <- block_patterns[block_labels, ] +
    matrix(rnorm(K * n_vox, sd = 0.1), nrow = K, ncol = n_vox)
  rownames(M_items) <- items

  X_rep <- M_items[as.character(mvdes$train_design$.rownum), , drop = FALSE]
  roi$train_roi@.Data <- X_rep

  ## Without confounds
  des_no <- repmed_design(
    items = items,
    X_rdm = X_rdm,
    Y_rdm = Y_rdm
  )
  mspec_no <- repmed_model(
    dataset    = ds$dataset,
    design     = mvdes,
    repmed_des = des_no,
    key_var    = ~ .rownum,
    distfun    = "euclidean"
  )
  res_no <- process_roi(mspec_no, roi, rnum = 1L, center_global_id = NA)
  perf_no <- res_no$performance[[1]]

  ## With confound RDM
  des_conf <- repmed_design(
    items         = items,
    X_rdm         = X_rdm,
    Y_rdm         = Y_rdm,
    confound_rdms = list(block = D_conf)
  )
  mspec_conf <- repmed_model(
    dataset    = ds$dataset,
    design     = mvdes,
    repmed_des = des_conf,
    key_var    = ~ .rownum,
    distfun    = "euclidean"
  )
  res_conf <- process_roi(mspec_conf, roi, rnum = 1L, center_global_id = NA)
  perf_conf <- res_conf$performance[[1]]

  a_no   <- unname(perf_no["med_a"])
  b_no   <- unname(perf_no["med_b"])
  ind_no <- unname(perf_no["med_indirect"])
  p_no   <- unname(perf_no["med_sobel_p"])

  a_conf   <- unname(perf_conf["med_a"])
  b_conf   <- unname(perf_conf["med_b"])
  ind_conf <- unname(perf_conf["med_indirect"])
  p_conf   <- unname(perf_conf["med_sobel_p"])

  ## Spurious mediation without confound
  expect_gt(abs(a_no), 0.1)
  expect_gt(abs(b_no), 0.1)
  expect_gt(ind_no, 0.05)

  ## After controlling for confound, effects should shrink
  expect_lt(abs(a_conf), abs(a_no))
  expect_lt(abs(b_conf), abs(b_no))
  expect_lt(abs(ind_conf), abs(ind_no))
  expect_gt(p_conf, p_no)
})


test_that("repmed_model: directional sanity check (B5)", {
  skip_on_cran()

  ## Partial mediation generative model
  sim <- simulate_repmed_latent(K = 24, a = 1.0, b = 0.8, c = 0.4)

  ## Forward direction: X -> M -> Y
  obj_fwd <- build_repmed_spec_and_roi(sim, reps_per_item = 2L, distfun = "euclidean")
  res_fwd <- process_roi(obj_fwd$mspec, obj_fwd$roi, rnum = 1L, center_global_id = NA)
  perf_fwd <- res_fwd$performance[[1]]

  ## Reverse direction: Y -> M -> X
  sim_rev <- sim
  sim_rev$X_rdm <- sim$Y_rdm
  sim_rev$Y_rdm <- sim$X_rdm
  obj_rev <- build_repmed_spec_and_roi(sim_rev, reps_per_item = 2L, distfun = "euclidean")
  res_rev <- process_roi(obj_rev$mspec, obj_rev$roi, rnum = 1L, center_global_id = NA)
  perf_rev <- res_rev$performance[[1]]

  ind_fwd <- unname(perf_fwd["med_indirect"])
  ind_rev <- unname(perf_rev["med_indirect"])
  p_fwd   <- unname(perf_fwd["med_sobel_p"])
  p_rev   <- unname(perf_rev["med_sobel_p"])

  expect_false(is.na(ind_fwd))
  expect_false(is.na(ind_rev))
  expect_false(is.na(p_fwd))
  expect_false(is.na(p_rev))

  ## Forward direction should show stronger mediation and smaller p-value
  expect_gt(ind_fwd, ind_rev + 0.02)
  expect_lt(p_fwd, p_rev)
})
