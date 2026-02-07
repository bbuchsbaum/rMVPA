context("era_rsa_model")

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
  res <- run_regional(ms, regionMask)
  expect_s3_class(res, "regional_mvpa_result")
  expect_true("geom_cor" %in% names(res$performance_table))
})
