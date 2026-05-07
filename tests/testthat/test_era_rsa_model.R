context("era_rsa_model")

test_that("era_rsa_model output schema documents all scalar metrics", {
  base_schema <- output_schema(structure(
    list(confound_rdms = NULL),
    class = "era_rsa_model"
  ))
  expect_identical(
    names(base_schema),
    c(
      "n_items",
      "era_top1_acc",
      "era_diag_mean",
      "era_diag_minus_off",
      "geom_cor",
      "era_diag_minus_off_same_block",
      "era_diag_minus_off_diff_block",
      "era_lag_cor",
      "geom_cor_partial",
      "geom_cor_run_partial",
      "geom_cor_xrun"
    )
  )
  expect_true(all(base_schema == "scalar"))

  conf_schema <- output_schema(structure(
    list(confound_rdms = list(block = diag(2), time_enc = diag(2))),
    class = "era_rsa_model"
  ))
  expect_true(all(c(
    "beta_enc_geom",
    "beta_block",
    "beta_time_enc",
    "sp_enc_geom",
    "sp_block",
    "sp_time_enc"
  ) %in% names(conf_schema)))
  expect_true(all(conf_schema == "scalar"))
})

test_that("era_rsa_model partial geometry selects nuisance groups and exact names", {
  keys <- paste0("item", seq_len(4))
  time_vec <- c(0, 1, 3, 4, 6, 7)
  signal <- c(-2, -1, 0, 1, 2, 3)
  dE <- signal + time_vec
  dR <- 2 * signal + time_vec

  time_mat <- matrix(0, 4, 4)
  time_mat[lower.tri(time_mat)] <- time_vec
  time_mat <- time_mat + t(time_mat)
  rownames(time_mat) <- colnames(time_mat) <- keys

  got_group <- rMVPA:::.era_partial_geometry_cor(
    dE = dE,
    dR = dR,
    confound_rdms = list(time_enc = time_mat),
    keys = keys,
    partial_against = "time",
    method = "pearson"
  )
  got_exact <- rMVPA:::.era_partial_geometry_cor(
    dE = dE,
    dR = dR,
    confound_rdms = list(time_enc = time_mat),
    keys = keys,
    partial_against = "time_enc",
    method = "pearson"
  )
  oracle <- stats::cor(
    stats::resid(stats::lm(dE ~ time_vec)),
    stats::resid(stats::lm(dR ~ time_vec))
  )

  expect_equal(got_group, oracle, tolerance = 1e-12)
  expect_equal(got_exact, oracle, tolerance = 1e-12)
  expect_true(is.na(rMVPA:::.era_partial_geometry_cor(
    dE = dE,
    dR = dR,
    confound_rdms = list(time_enc = time_mat),
    keys = keys,
    partial_against = "run",
    method = "pearson"
  )))
})

test_that("era_rsa_model derives run partial geometry from item_run metadata", {
  set.seed(11)
  K <- 6
  p <- 5
  toy <- gen_sample_dataset(D = c(3,3,3), nobs = K, nlevels = 2, blocks = 2,
                            external_test = TRUE, ntest_obs = K)
  keys <- paste0("item", seq_len(K))
  toy$design$train_design$item <- factor(keys, levels = keys)
  toy$design$test_design$item  <- factor(keys, levels = keys)

  item_run_enc <- setNames(rep(c("enc_1", "enc_2"), length.out = K), keys)
  item_run_ret <- setNames(rep(c("ret_1", "ret_2", "ret_3"), length.out = K), keys)
  model <- suppressWarnings(era_rsa_model(
    dataset = toy$dataset,
    design = toy$design,
    key_var = ~ item,
    phase_var = ~ block_var,
    item_block = setNames(rep(c("b1", "b2"), length.out = K), keys),
    item_run_enc = item_run_enc,
    item_run_ret = item_run_ret
  ))

  E <- matrix(rnorm(K * p), K, p)
  R <- E + matrix(rnorm(K * p, sd = 0.2), K, p)
  out <- fit_roi(
    model,
    roi_data = list(train_data = E, test_data = R, indices = seq_len(p)),
    context = list(id = 1L)
  )

  expect_false(out$error)
  expect_true(is.finite(out$metrics[["geom_cor_partial"]]))
  expect_true(is.finite(out$metrics[["geom_cor_run_partial"]]))
  expect_equal(out$metrics[["geom_cor_partial"]], out$metrics[["geom_cor_run_partial"]], tolerance = 1e-12)
})

test_that("era_rsa_model merges and selects global nuisance RDMs", {
  set.seed(12)
  K <- 5
  keys <- paste0("item", seq_len(K))
  toy <- gen_sample_dataset(D = c(3, 3, 3), nobs = K, nlevels = 2, blocks = 2,
                            external_test = TRUE, ntest_obs = K)
  toy$design$train_design$item <- factor(keys, levels = keys)
  toy$design$test_design$item <- factor(keys, levels = keys)

  D_enc <- as.matrix(stats::dist(matrix(seq_len(K), ncol = 1)))
  D_ret <- as.matrix(stats::dist(matrix(seq_len(K)^2, ncol = 1)))
  rownames(D_enc) <- colnames(D_enc) <- keys
  rownames(D_ret) <- colnames(D_ret) <- keys

  model <- suppressWarnings(era_rsa_model(
    dataset = toy$dataset,
    design = toy$design,
    key_var = ~ item,
    phase_var = ~ block_var,
    item_block = setNames(rep(c("b1", "b2"), length.out = K), keys),
    global_nuisance = list(D_enc = D_enc, D_ret = D_ret)
  ))

  expect_identical(names(model$confound_rdms), c("global_enc", "global_ret"))
  expect_identical(model$partial_against, c("run", "global"))
  expect_identical(
    rMVPA:::.era_select_confound_names(names(model$confound_rdms), "global"),
    c("global_enc", "global_ret")
  )

  model_null <- suppressWarnings(era_rsa_model(
    dataset = toy$dataset,
    design = toy$design,
    key_var = ~ item,
    phase_var = ~ block_var,
    item_block = setNames(rep(c("b1", "b2"), length.out = K), keys),
    global_nuisance = NULL
  ))
  expect_identical(model_null$partial_against, "run")
})

test_that("global nuisance auto-computes item-level whole-mask RDMs", {
  set.seed(13)
  K <- 6
  keys <- paste0("item", seq_len(K))
  toy <- gen_sample_dataset(D = c(3, 3, 3), nobs = K, nlevels = 2, blocks = 2,
                            external_test = TRUE, ntest_obs = K)
  toy$design$train_design$item <- factor(keys, levels = keys)
  toy$design$test_design$item <- factor(keys, levels = keys)

  got <- rMVPA:::.era_resolve_global_nuisance(
    TRUE,
    dataset = toy$dataset,
    design = toy$design,
    key_var = "item",
    distfun = eucdist()
  )

  expect_identical(got$common_keys, keys)
  expect_equal(dim(got$S_cross), c(K, K))
  expect_equal(dim(got$D_enc), c(K, K))
  expect_equal(dim(got$D_ret), c(K, K))
  expect_identical(rownames(got$S_cross), keys)
  expect_identical(colnames(got$S_cross), keys)
  expect_equal(unname(diag(got$D_enc)), rep(0, K), tolerance = 1e-12)
  expect_equal(unname(diag(got$D_ret)), rep(0, K), tolerance = 1e-12)
})

test_that("ERA models retain evaluated character key_var values", {
  K <- 5
  keys <- paste0("item", seq_len(K))
  toy <- gen_sample_dataset(D = c(3, 3, 3), nobs = K, nlevels = 2, blocks = 2,
                            external_test = TRUE, ntest_obs = K)
  toy$design$train_design$item <- factor(keys, levels = keys)
  toy$design$test_design$item <- factor(keys, levels = keys)
  kv <- "item"

  rsa_model <- suppressWarnings(era_rsa_model(
    dataset = toy$dataset,
    design = toy$design,
    key_var = kv,
    phase_var = ~ block_var,
    item_block = setNames(rep(c("b1", "b2"), length.out = K), keys)
  ))
  partition_model <- suppressWarnings(era_partition_model(
    dataset = toy$dataset,
    design = toy$design,
    key_var = kv,
    include_procrustes = FALSE
  ))

  expect_identical(rsa_model$key_var, "item")
  expect_identical(partition_model$key_var, "item")
})

test_that("era_rsa_model runs regionally and returns expected metrics", {
  skip_on_cran()
  set.seed(123)

  # Create toy encoding/retrieval dataset with external test set
  toy <- gen_sample_dataset(D = c(4,4,4), nobs = 48, nlevels = 3, blocks = 3, external_test = TRUE)

  # Add a shared item key column to both train and test designs
  toy$design$train_design$item <- toy$design$train_design$Y
  toy$design$test_design$item  <- toy$design$test_design$Ytest

  # Simple regional mask with 3 regions
  regionMask <- neuroim2::NeuroVol(sample(1:3, size = length(toy$dataset$mask), replace = TRUE),
                                   neuroim2::space(toy$dataset$mask))

  ms <- era_rsa_model(
    dataset = toy$dataset,
    design  = toy$design,
    key_var = ~ item,
    phase_var = ~ block_var,     # phase label not used in external-test path; keep for interface
    distfun = cordist("pearson"),
    rsa_simfun = "spearman"
  )

  res <- run_regional(ms, regionMask)
  expect_s3_class(res, "regional_mvpa_result")
  expect_true(is.data.frame(res$performance_table))
  # Check a couple of expected metrics exist
  expect_true(all(c("geom_cor", "era_top1_acc") %in% names(res$performance_table)))
})

test_that("era_rsa_model and regional results use dedicated print methods", {
  skip_on_cran()
  set.seed(124)

  toy <- gen_sample_dataset(D = c(4,4,4), nobs = 24, nlevels = 3, blocks = 3, external_test = TRUE)
  toy$design$train_design$item <- toy$design$train_design$Y
  toy$design$test_design$item  <- toy$design$test_design$Ytest

  ms <- era_rsa_model(
    dataset = toy$dataset,
    design  = toy$design,
    key_var = ~ item,
    phase_var = ~ block_var
  )

  regionMask <- neuroim2::NeuroVol(sample(1:2, size = length(toy$dataset$mask), replace = TRUE),
                                   neuroim2::space(toy$dataset$mask))
  res <- run_regional(ms, regionMask)

  ms_out <- paste(capture.output(print(ms)), collapse = "\n")
  res_out <- paste(capture.output(print(res)), collapse = "\n")

  expect_match(ms_out, "Model Specification")
  expect_match(ms_out, "MVPA Dataset")
  expect_match(res_out, "Regional Analysis Results")
  expect_false(grepl("\\$model_spec", res_out))
})


test_that("era_rsa_model searchlight returns metric maps", {
  skip_on_cran()
  set.seed(456)

  toy <- gen_sample_dataset(D = c(4,4,4), nobs = 40, nlevels = 2, blocks = 2, external_test = TRUE)
  toy$design$train_design$item <- toy$design$train_design$Y
  toy$design$test_design$item  <- toy$design$test_design$Ytest

  ms <- era_rsa_model(
    dataset = toy$dataset,
    design  = toy$design,
    key_var = ~ item,
    phase_var = ~ block_var,
    distfun = cordist("pearson"),
    rsa_simfun = "pearson"
  )

  sl <- run_searchlight(ms, radius = 3, method = "standard")
  expect_s3_class(sl, "searchlight_result")
  # At least these metrics should be present as maps
  expect_true(all(c("geom_cor", "era_top1_acc") %in% sl$metrics))
})


test_that("era_rsa_model randomized searchlight returns metric maps", {
  skip_on_cran()
  set.seed(457)

  toy <- gen_sample_dataset(D = c(4,4,4), nobs = 40, nlevels = 2, blocks = 2, external_test = TRUE)
  toy$design$train_design$item <- toy$design$train_design$Y
  toy$design$test_design$item  <- toy$design$test_design$Ytest

  ms <- era_rsa_model(
    dataset = toy$dataset,
    design  = toy$design,
    key_var = ~ item,
    phase_var = ~ block_var,
    distfun = cordist("pearson"),
    rsa_simfun = "pearson"
  )

  sl_rand <- run_searchlight(ms, radius = 3, method = "randomized", niter = 2)
  expect_s3_class(sl_rand, "searchlight_result")
  expect_true(all(c("geom_cor", "era_top1_acc") %in% sl_rand$metrics))
})


test_that("era_rsa_model accepts confound RDMs and emits beta_* metrics", {
  skip_on_cran()
  set.seed(789)

  toy <- gen_sample_dataset(D = c(3,3,3), nobs = 36, nlevels = 3, blocks = 3, external_test = TRUE)
  toy$design$train_design$item <- toy$design$train_design$Y
  toy$design$test_design$item  <- toy$design$test_design$Ytest

  # Build a simple block-based confound RDM at the item level using training block_var
  item_levels <- levels(toy$design$train_design$item)
  # modal block per item (encoding phase)
  Mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }
  enc_block_by_item <- sapply(item_levels, function(it) {
    Mode(toy$design$train_design$block_var[toy$design$train_design$item == it])
  })
  names(enc_block_by_item) <- item_levels
  block_rdm <- outer(enc_block_by_item, enc_block_by_item, FUN = function(a,b) as.numeric(a != b))
  rownames(block_rdm) <- colnames(block_rdm) <- item_levels

  ms <- era_rsa_model(
    dataset = toy$dataset,
    design  = toy$design,
    key_var = ~ item,
    phase_var = ~ block_var,
    confound_rdms = list(block = block_rdm)
  )

  regionMask <- neuroim2::NeuroVol(sample(1:2, size = length(toy$dataset$mask), replace = TRUE),
                                   neuroim2::space(toy$dataset$mask))
  res <- run_regional(ms, regionMask)
  expect_s3_class(res, "regional_mvpa_result")
  # Expect at least one beta_ term from geometry regression
  expect_true(any(grepl("^beta_", names(res$performance_table))))
})


test_that("era_rsa_model computes run-partial, cross-run, block-limited, and lag metrics", {
  skip_on_cran()
  set.seed(101)

  toy <- gen_sample_dataset(D = c(3,3,3), nobs = 36, nlevels = 3, blocks = 3, external_test = TRUE)
  # Shared item keys across phases
  toy$design$train_design$item <- toy$design$train_design$Y
  toy$design$test_design$item  <- toy$design$test_design$Ytest

  # Derive per-item encoding runs from train_design (modal block)
  items <- levels(toy$design$train_design$item)
  Mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }
  item_run_enc <- sapply(items, function(it) Mode(toy$design$train_design$block_var[toy$design$train_design$item == it]))
  names(item_run_enc) <- items
  # Create a retrieval run assignment (independent) for coverage
  item_run_ret <- sample(item_run_enc)
  # Item-level block labels: reuse encoding run as block
  item_block <- factor(item_run_enc, levels = sort(unique(item_run_enc)))
  # Lag vector
  item_lag <- setNames(seq_along(items), items)

  # Build run confound RDMs
  run_enc <- outer(item_run_enc, item_run_enc, FUN = function(a,b) as.numeric(a == b))
  run_ret <- outer(item_run_ret, item_run_ret, FUN = function(a,b) as.numeric(a == b))
  rownames(run_enc) <- colnames(run_enc) <- items
  rownames(run_ret) <- colnames(run_ret) <- items

  ms <- era_rsa_model(
    dataset = toy$dataset,
    design  = toy$design,
    key_var = ~ item,
    phase_var = ~ block_var,
    confound_rdms = list(run_enc = run_enc, run_ret = run_ret),
    item_block = item_block,
    item_lag   = item_lag,
    item_run_enc = factor(item_run_enc),
    item_run_ret = factor(item_run_ret)
  )

  regionMask <- neuroim2::NeuroVol(sample(1:2, size = length(toy$dataset$mask), replace = TRUE),
                                   neuroim2::space(toy$dataset$mask))
  res <- run_regional(ms, regionMask)
  expect_s3_class(res, "regional_mvpa_result")

  nm <- names(res$performance_table)
  expect_true(all(c("geom_cor_run_partial", "geom_cor_xrun", "era_diag_minus_off_same_block", "era_lag_cor") %in% nm))
  # Values should be finite or NA when masking removes all pairs; require not all NA across ROIs
  vals <- res$performance_table$geom_cor_run_partial
  expect_true(any(is.finite(vals)))
  vals2 <- res$performance_table$era_diag_minus_off_same_block
  expect_true(any(is.finite(vals2)))
  vals3 <- res$performance_table$era_lag_cor
  expect_true(any(is.finite(vals3)))
})


test_that("era_rsa_model works with euclidean distfun", {
  skip_on_cran()
  set.seed(202)

  toy <- gen_sample_dataset(D = c(3,3,3), nobs = 30, nlevels = 3, blocks = 3, external_test = TRUE)
  toy$design$train_design$item <- toy$design$train_design$Y
  toy$design$test_design$item  <- toy$design$test_design$Ytest

  ms <- era_rsa_model(
    dataset = toy$dataset,
    design  = toy$design,
    key_var = ~ item,
    phase_var = ~ block_var,
    distfun = eucdist()
  )
  regionMask <- neuroim2::NeuroVol(sample(1:2, size = length(toy$dataset$mask), replace = TRUE),
                                   neuroim2::space(toy$dataset$mask))
  res <- run_regional(suppressWarnings(ms), regionMask)
  expect_s3_class(res, "regional_mvpa_result")
  expect_true("geom_cor" %in% names(res$performance_table))
})

test_that("era_rsa_model warns about missing item-level block/run metadata", {
  set.seed(8)
  toy <- gen_sample_dataset(D = c(3,3,3), nobs = 18, nlevels = 3, blocks = 3, external_test = TRUE)
  toy$design$train_design$item <- toy$design$train_design$Y
  toy$design$test_design$item  <- toy$design$test_design$Ytest

  expect_warning(
    era_rsa_model(toy$dataset, toy$design, key_var = ~ item, phase_var = ~ block_var),
    "item_block.*will be NA"
  )
  expect_warning(
    era_rsa_model(toy$dataset, toy$design, key_var = ~ item, phase_var = ~ block_var),
    "geom_cor_xrun"
  )

  expect_error(
    era_rsa_model(
      toy$dataset, toy$design,
      key_var = ~ item, phase_var = ~ block_var,
      require_run_metadata = TRUE
    ),
    "item_block"
  )
})

test_that("era_rsa_model warns when item_run_enc and item_run_ret share atomic labels", {
  set.seed(9)
  toy <- gen_sample_dataset(D = c(3,3,3), nobs = 18, nlevels = 3, blocks = 3, external_test = TRUE)
  toy$design$train_design$item <- toy$design$train_design$Y
  toy$design$test_design$item  <- toy$design$test_design$Ytest

  items <- as.character(sort(intersect(toy$design$train_design$item, toy$design$test_design$item)))
  ire <- setNames(rep_len(c(1, 2), length(items)), items)
  irr <- setNames(rep_len(c(1, 2), length(items)), items)

  expect_warning(
    era_rsa_model(
      toy$dataset, toy$design,
      key_var = ~ item, phase_var = ~ block_var,
      item_block = setNames(rep_len(c(1, 2), length(items)), items),
      item_run_enc = ire, item_run_ret = irr
    ),
    "phase-scoped labels"
  )
})
