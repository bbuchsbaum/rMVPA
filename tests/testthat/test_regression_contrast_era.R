context("regression: contrast RSA row order and ERA item alignment")

# ---------------------------------------------------------------------------
# Contrast RSA (msreve_design): rows of the contrast matrix and of nuisance
# RDMs must be aligned to the condition levels used by the cross-validated
# means, not to the order the user happened to type them in.
# ---------------------------------------------------------------------------

.cs_reg_metrics <- c("beta_delta", "beta_only", "delta_only",
                     "recon_score", "beta_delta_norm")

.cs_reg_design <- function(n_cond = 4, n_blocks = 3, per_cell = 2, n_vox = 20) {
  n <- n_cond * n_blocks * per_cell
  Y <- factor(rep(paste0("Cond", seq_len(n_cond)), each = n / n_cond))
  block_var <- factor(rep(seq_len(n_blocks), length.out = n))
  des <- structure(
    list(
      Y = Y,
      block_var = block_var,
      conditions = levels(Y),
      ncond = n_cond,
      design_matrix = matrix(0, n, 1),
      samples = seq_len(n)
    ),
    class = c("mvpa_design", "list")
  )
  set.seed(456)
  means <- matrix(as.numeric(Y), nrow = n, ncol = n_vox) * 10
  data <- means + matrix(rnorm(n * n_vox), n, n_vox)
  mask <- array(1, dim = c(3, 3, 3))
  mask[-(1:n_vox)] <- 0
  structure(
    list(train_data = data, mask = mask, has_test_set = FALSE,
         design = des, nfeatures = n_vox),
    class = c("mvpa_dataset", "list")
  )
}

.cs_reg_run <- function(ms, dset, metrics = .cs_reg_metrics) {
  spec <- contrast_rsa_model(dset, ms, output_metric = metrics,
                             check_collinearity = FALSE)
  des <- dset$design
  cv <- structure(
    list(.folds_val = as.integer(des$block_var),
         .n_folds_val = nlevels(des$block_var)),
    class = c("mock_cv_spec", "cross_validation", "list")
  )
  n_sl <- 7L
  sl_data <- dset$train_data[, seq_len(n_sl), drop = FALSE]
  sl_info <- list(center_local_id = 4, center_global_id = 4,
                  radius = 0, n_voxels = n_sl)
  suppressWarnings(with_mocked_bindings(
    get_nfolds = function(obj, ...) obj$.n_folds_val,
    train_indices = function(obj, fold_num, ...) which(obj$.folds_val != fold_num),
    .package = "rMVPA",
    train_model(spec, sl_data, sl_info, cv)
  ))
}

.cs_reg_C <- function() {
  C <- matrix(c(1, 1, -1, -1,
                1, -1, 1, -1), nrow = 4, ncol = 2,
              dimnames = list(paste0("Cond", 1:4), c("CFC1", "CFC2")))
  C
}

test_that("contrast RSA output does not depend on the row order of named contrasts", {
  dset <- .cs_reg_design()
  C_lvl <- .cs_reg_C()
  perm <- c(3, 1, 4, 2)
  C_perm <- C_lvl[perm, , drop = FALSE]
  expect_identical(rownames(C_perm), rownames(C_lvl)[perm])

  ms_lvl <- msreve_design(dset$design, C_lvl)
  ms_perm <- msreve_design(dset$design, C_perm)
  expect_equal(ms_perm$contrast_matrix, ms_lvl$contrast_matrix)
  expect_identical(rownames(ms_perm$contrast_matrix), levels(dset$design$Y))

  res_lvl <- .cs_reg_run(ms_lvl, dset)
  res_perm <- .cs_reg_run(ms_perm, dset)
  for (m in .cs_reg_metrics) {
    expect_equal(res_perm[[m]], res_lvl[[m]], info = m)
  }
  expect_false(anyNA(res_lvl$beta_only))
  expect_false(anyNA(res_lvl$beta_delta))
})

test_that("contrast RSA rejects named contrast rows that do not match the condition levels", {
  dset <- .cs_reg_design()
  C <- .cs_reg_C()

  C_wrong <- C
  rownames(C_wrong)[4] <- "CondX"
  expect_error(msreve_design(dset$design, C_wrong),
               "Row names of `contrast_matrix` must match the condition labels")

  C_dup <- C
  rownames(C_dup) <- c("Cond1", "Cond1", "Cond3", "Cond4")
  expect_error(msreve_design(dset$design, C_dup),
               "Row names of `contrast_matrix` must match the condition labels")

  C_subset <- C[1:3, , drop = FALSE]
  expect_error(msreve_design(dset$design, C_subset),
               "Row names of `contrast_matrix` must match the condition labels")
})

test_that("unnamed contrast rows are taken in level order", {
  dset <- .cs_reg_design()
  C_lvl <- .cs_reg_C()
  C_unnamed <- C_lvl
  rownames(C_unnamed) <- NULL

  ms_lvl <- msreve_design(dset$design, C_lvl)
  ms_unnamed <- msreve_design(dset$design, C_unnamed)
  expect_identical(rownames(ms_unnamed$contrast_matrix), levels(dset$design$Y))
  expect_equal(ms_unnamed$contrast_matrix, ms_lvl$contrast_matrix)

  res_lvl <- .cs_reg_run(ms_lvl, dset)
  res_unnamed <- .cs_reg_run(ms_unnamed, dset)
  for (m in .cs_reg_metrics) {
    expect_equal(res_unnamed[[m]], res_lvl[[m]], info = m)
  }

  expect_error(msreve_design(dset$design, C_unnamed[1:3, , drop = FALSE]),
               "must match number of conditions")
})

test_that("nuisance RDMs are aligned to condition levels by name", {
  dset <- .cs_reg_design()
  C_lvl <- .cs_reg_C()
  lv <- levels(dset$design$Y)
  set.seed(7)
  N <- as.matrix(stats::dist(matrix(rnorm(8), 4)))
  dimnames(N) <- list(lv, lv)
  perm <- c(2, 4, 1, 3)
  N_perm <- N[perm, perm, drop = FALSE]
  expect_identical(rownames(N_perm), lv[perm])

  ms_ref <- msreve_design(dset$design, C_lvl, nuisance_rdms = list(time = N))
  ms_nuis_perm <- msreve_design(dset$design, C_lvl, nuisance_rdms = list(time = N_perm))
  ms_both_perm <- msreve_design(dset$design, C_lvl[perm, , drop = FALSE],
                                nuisance_rdms = list(time = N_perm))
  ms_unnamed <- msreve_design(dset$design, C_lvl,
                              nuisance_rdms = list(time = unname(N)))
  ms_none <- msreve_design(dset$design, C_lvl)

  expect_identical(rownames(ms_nuis_perm$nuisance_rdms$time), lv)
  expect_equal(ms_nuis_perm$nuisance_rdms$time, N)

  res_ref <- .cs_reg_run(ms_ref, dset)
  for (ms in list(ms_nuis_perm, ms_both_perm, ms_unnamed)) {
    res <- .cs_reg_run(ms, dset)
    for (m in .cs_reg_metrics) {
      expect_equal(res[[m]], res_ref[[m]], info = m)
    }
  }

  # The nuisance predictor is actually used in the regression.
  res_none <- .cs_reg_run(ms_none, dset)
  expect_false(isTRUE(all.equal(res_none$beta_only, res_ref$beta_only)))
})

test_that("nuisance RDMs with unrelated labels are rejected", {
  dset <- .cs_reg_design()
  C_lvl <- .cs_reg_C()
  N <- matrix(c(0, 1, 2, 1, 1, 0, 1, 2, 2, 1, 0, 1, 1, 2, 1, 0), 4, 4)
  rownames(N) <- colnames(N) <- c("a", "b", "c", "d")
  expect_error(suppressWarnings(
    msreve_design(dset$design, C_lvl, nuisance_rdms = list(time = N))
  ), "item labels must match the condition labels")
})

# ---------------------------------------------------------------------------
# ERA models: per-item metadata and confounds must be aligned to item keys,
# with the canonical key order (factor levels, else numeric, else C-sort).
# ---------------------------------------------------------------------------

.era_reg_keys <- as.character(1:12)

# Fixture with one encoding and one retrieval trial per item (12 items, keys
# "1".."12"). Trial rows are shuffled so the design row order differs from the
# key order. `key_type` selects character, numeric or factor keys.
.era_reg_fixture <- function(key_type = "character", levels_override = NULL,
                             seed = 2026) {
  K <- 12L
  p <- 10L
  set.seed(seed)
  toy <- gen_sample_dataset(D = c(3, 3, 3), nobs = K, nlevels = 2, blocks = 2,
                            external_test = TRUE, ntest_obs = K)
  enc_keys <- sample(.era_reg_keys)
  ret_keys <- sample(.era_reg_keys)
  make_key <- function(x) {
    switch(key_type,
           character = x,
           numeric = as.numeric(x),
           factor = factor(x, levels = levels_override %||% .era_reg_keys))
  }
  toy$design$train_design$item <- make_key(enc_keys)
  toy$design$test_design$item <- make_key(ret_keys)

  E_key <- matrix(rnorm(K * p), K, p, dimnames = list(.era_reg_keys, NULL))
  R_key <- E_key + matrix(rnorm(K * p, sd = 0.8), K, p,
                          dimnames = list(.era_reg_keys, NULL))
  list(
    toy = toy,
    E_key = E_key,
    R_key = R_key,
    Xenc = E_key[enc_keys, , drop = FALSE],
    Xret = R_key[ret_keys, , drop = FALSE],
    p = p
  )
}

# Truth table: three blocks of four items, assigned by key.
.era_reg_block_truth <- function(seed = 99) {
  set.seed(seed)
  setNames(sample(rep(c("A", "B", "C"), each = 4)), .era_reg_keys)
}

# Independent computation of the block-specific item specificity, using the
# definitions in the ERA-RSA documentation: for retrieval item i, matched
# similarity minus mean similarity to encoding items in the same (other)
# block; averaged over items.
.era_reg_expected_block_metrics <- function(E_key, R_key, block_by_key) {
  b <- block_by_key[.era_reg_keys]
  S <- stats::cor(t(R_key), t(E_key)) # S[i, j] = cor(retrieval i, encoding j)
  K <- nrow(S)
  same <- diff <- numeric(K)
  for (i in seq_len(K)) {
    same_j <- which(b == b[i] & seq_len(K) != i)
    diff_j <- which(b != b[i])
    same[i] <- S[i, i] - mean(S[i, same_j])
    diff[i] <- S[i, i] - mean(S[i, diff_j])
  }
  c(same = mean(same), diff = mean(diff))
}

.era_reg_rsa_model <- function(fx, item_block, ...) {
  suppressWarnings(era_rsa_model(
    dataset = fx$toy$dataset,
    design = fx$toy$design,
    key_var = ~ item,
    phase_var = ~ block_var,
    item_block = item_block,
    era_components = "item",
    ...
  ))
}

.era_reg_rsa_metrics <- function(model, fx) {
  out <- fit_roi(
    model,
    roi_data = list(train_data = fx$Xenc, test_data = fx$Xret,
                    indices = seq_len(fx$p)),
    context = list(id = 1L)
  )
  expect_false(out$error)
  out$metrics
}

test_that("ERA-RSA block metrics match the truth table for character, numeric and factor keys", {
  truth <- .era_reg_block_truth()
  exp_m <- .era_reg_expected_block_metrics(
    E_key = .era_reg_fixture("character")$E_key,
    R_key = .era_reg_fixture("character")$R_key,
    block_by_key = truth
  )

  for (kt in c("character", "numeric", "factor")) {
    fx <- .era_reg_fixture(kt)
    # Named metadata, in shuffled name order.
    named <- truth[sample(.era_reg_keys)]
    m_named <- .era_reg_rsa_metrics(.era_reg_rsa_model(fx, named), fx)
    # Unnamed metadata, in level order (1..12 for these keys).
    m_unnamed <- .era_reg_rsa_metrics(
      .era_reg_rsa_model(fx, unname(truth[.era_reg_keys])), fx
    )
    for (m in list(m_named, m_unnamed)) {
      expect_equal(m[["era_diag_minus_off_same_block"]], exp_m[["same"]],
                   tolerance = 1e-10, info = kt)
      expect_equal(m[["era_diag_minus_off_diff_block"]], exp_m[["diff"]],
                   tolerance = 1e-10, info = kt)
    }
  }
})

test_that("ERA-RSA metrics do not depend on the order of named item metadata", {
  truth <- .era_reg_block_truth()
  fx <- .era_reg_fixture("character")
  ref <- .era_reg_rsa_metrics(.era_reg_rsa_model(fx, truth), fx)
  shuffled <- truth[sample(.era_reg_keys)]
  got <- .era_reg_rsa_metrics(.era_reg_rsa_model(fx, shuffled), fx)
  expect_identical(got, ref)
})

test_that("unnamed ERA-RSA metadata follows the factor level order of the key", {
  truth <- .era_reg_block_truth()
  rev_keys <- rev(.era_reg_keys)
  fx <- .era_reg_fixture("factor", levels_override = rev_keys)
  exp_m <- .era_reg_expected_block_metrics(fx$E_key, fx$R_key, truth)
  # Levels are reversed, so unnamed values are supplied in reversed order.
  m <- .era_reg_rsa_metrics(
    .era_reg_rsa_model(fx, unname(truth[rev_keys])), fx
  )
  expect_equal(m[["era_diag_minus_off_same_block"]], exp_m[["same"]], tolerance = 1e-10)
  expect_equal(m[["era_diag_minus_off_diff_block"]], exp_m[["diff"]], tolerance = 1e-10)
})

test_that("unnamed ERA-RSA metadata of the wrong length is an error", {
  truth <- .era_reg_block_truth()
  fx <- .era_reg_fixture("character")
  expect_error(
    .era_reg_rsa_model(fx, unname(truth)[1:11]),
    "unnamed with length 11"
  )
})

# ERA-partition: same-block nuisance regressors are built from item metadata.
.era_reg_partition_model <- function(fx, item_block_enc, item_block_ret) {
  suppressWarnings(era_partition_model(
    dataset = fx$toy$dataset,
    design = fx$toy$design,
    key_var = ~ item,
    distfun = eucdist(),
    item_block_enc = item_block_enc,
    item_block_ret = item_block_ret,
    include_procrustes = FALSE,
    compute_xdec_performance = FALSE
  ))
}

.era_reg_partition_metrics <- function(model, fx) {
  out <- fit_roi(
    model,
    roi_data = list(train_data = fx$Xenc, test_data = fx$Xret,
                    indices = seq_len(fx$p)),
    context = list(id = 1L)
  )
  expect_false(out$error)
  out$metrics
}

test_that("ERA-partition same-block nuisance regressors match the truth table", {
  truth <- .era_reg_block_truth()
  keys <- .era_reg_keys
  b <- truth[keys]
  exp_cross <- as.numeric(outer(b, b, "=="))
  exp_enc_lower <- { M <- outer(b, b, "=="); as.numeric(M[lower.tri(M)]) }

  for (kt in c("character", "numeric", "factor")) {
    fx <- .era_reg_fixture(kt)
    unnamed_model <- .era_reg_partition_model(
      fx, unname(truth[keys]), unname(truth[keys])
    )
    named_model <- .era_reg_partition_model(
      fx, truth[sample(keys)], truth[sample(keys)]
    )
    for (model in list(unnamed_model, named_model)) {
      first <- rMVPA:::.era_partition_first_nuisance(model, keys)
      second <- rMVPA:::.era_partition_second_nuisance(model, keys)
      expect_equal(first$same_block_cross, exp_cross, info = kt)
      expect_equal(second$same_block_enc, exp_enc_lower, info = kt)
      expect_equal(second$same_block_ret, exp_enc_lower, info = kt)
    }
  }
})

test_that("ERA-partition block-nuisance delta R2 matches an oracle built from the truth table", {
  truth <- .era_reg_block_truth()
  keys <- .era_reg_keys
  b <- truth[keys]
  same_mat <- outer(b, b, "==")
  for (kt in c("character", "numeric", "factor")) {
    fx <- .era_reg_fixture(kt)
    # Cross-state similarity and encoding/retrieval distances in key order.
    S <- stats::cor(t(fx$R_key[keys, , drop = FALSE]), t(fx$E_key[keys, , drop = FALSE]))
    dE <- as.numeric(stats::dist(fx$E_key[keys, , drop = FALSE]))
    dR <- as.numeric(stats::dist(fx$R_key[keys, , drop = FALSE]))
    # Oracle: same-block regressors from the truth table only.
    first_oracle <- rMVPA:::.era_partition_delta_r2(
      y = as.numeric(S), signal = as.numeric(diag(12)),
      nuisance = list(same_block_cross = as.numeric(same_mat))
    )
    lt <- lower.tri(same_mat)
    second_oracle <- rMVPA:::.era_partition_delta_r2(
      y = dR, signal = dE,
      nuisance = list(same_block_enc = as.numeric(same_mat[lt]),
                      same_block_ret = as.numeric(same_mat[lt]))
    )
    model <- .era_reg_partition_model(fx, unname(truth[keys]), unname(truth[keys]))
    m <- .era_reg_partition_metrics(model, fx)
    expect_equal(m[["first_order_delta_r2"]], first_oracle$delta_r2,
                 tolerance = 1e-8, info = kt)
    expect_equal(m[["second_order_delta_r2"]], second_oracle$delta_r2,
                 tolerance = 1e-8, info = kt)
  }
})

test_that("ERA-partition metrics do not depend on named versus unnamed item metadata", {
  truth <- .era_reg_block_truth()
  fx <- .era_reg_fixture("factor")
  named <- .era_reg_partition_model(fx, truth[sample(.era_reg_keys)],
                                    truth[sample(.era_reg_keys)])
  unnamed <- .era_reg_partition_model(fx, unname(truth[.era_reg_keys]),
                                      unname(truth[.era_reg_keys]))
  expect_identical(.era_reg_partition_metrics(unnamed, fx),
                   .era_reg_partition_metrics(named, fx))
})

test_that("ERA-partition rejects unnamed item metadata of the wrong length", {
  truth <- .era_reg_block_truth()
  fx <- .era_reg_fixture("character")
  expect_error(
    .era_reg_partition_model(fx, unname(truth)[1:11], unname(truth)),
    "unnamed with length 11"
  )
})

# ERA-RSA confounds: partial_against groups must match whole name parts only.
test_that("confound group matching does not match substrings", {
  sel <- rMVPA:::.era_select_confound_names
  expect_identical(sel(c("location", "duplicate_flag"), "category"), character(0))
  expect_identical(sel(c("category", "location"), "category"), "category")
  expect_identical(sel(c("cat_enc", "location"), "category"), "cat_enc")
  expect_identical(sel(c("location", "time_enc"), "time"), "time_enc")
  expect_identical(sel(c("run_enc", "rerun_flag"), "run"), "run_enc")
  expect_identical(sel(c("global_enc", "location"), "global"), "global_enc")
  expect_identical(sel(c("location", "duplicate_flag"), "location"), "location")
})

test_that("ERA-RSA partial geometry ignores confounds named like substrings of 'category'", {
  truth <- .era_reg_block_truth()
  fx <- .era_reg_fixture("character")
  keys <- .era_reg_keys
  set.seed(5)
  M <- as.matrix(stats::dist(matrix(rnorm(12 * 3), 12)))
  dimnames(M) <- list(keys, keys)

  build <- function(name) {
    suppressWarnings(era_rsa_model(
      dataset = fx$toy$dataset,
      design = fx$toy$design,
      key_var = ~ item,
      phase_var = ~ block_var,
      item_block = truth,
      confound_rdms = setNames(list(M), name),
      partial_against = "category"
    ))
  }

  out_loc <- .era_reg_rsa_metrics(build("location"), fx)
  expect_true(is.na(out_loc[["geom_cor_partial"]]))

  out_dup <- .era_reg_rsa_metrics(build("duplicate_flag"), fx)
  expect_true(is.na(out_dup[["geom_cor_partial"]]))

  model_cat <- build("category")
  out_cat <- .era_reg_rsa_metrics(model_cat, fx)
  expect_true(is.finite(out_cat[["geom_cor_partial"]]))

  # Oracle: residual correlation after regressing both geometries on the
  # category confound (lower triangle, keys in canonical order).
  E <- fx$E_key[keys, , drop = FALSE]
  R <- fx$R_key[keys, , drop = FALSE]
  dE <- as.numeric(rMVPA:::pairwise_dist(model_cat$distfun, E)[lower.tri(diag(12))])
  dR <- as.numeric(rMVPA:::pairwise_dist(model_cat$distfun, R)[lower.tri(diag(12))])
  m <- M[lower.tri(M)]
  oracle <- stats::cor(stats::resid(stats::lm(dE ~ m)),
                       stats::resid(stats::lm(dR ~ m)))
  expect_equal(out_cat[["geom_cor_partial"]], oracle, tolerance = 1e-8)
})
