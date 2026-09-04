# Permutation inference for rsa_model: item-label permutation and design
# diagnostics.  Companion to test_permutation_searchlight.R.

library(testthat)

quiet_logger <- function() {
  futile.logger::flog.threshold(futile.logger::ERROR)
  withr::defer(futile.logger::flog.threshold(futile.logger::INFO), envir = parent.frame())
}

# 24 items in 3 runs of 8; two independent model RDMs.
make_rsa_design <- function(n = 24, blocks = 3, seed = 1, keep_intra_run = FALSE) {
  set.seed(seed)
  a <- dist(matrix(rnorm(n * 4), n, 4))
  b <- dist(matrix(rnorm(n * 4), n, 4))
  block <- rep(seq_len(blocks), each = n / blocks)
  rsa_design(~ a + b, list(a = a, b = b, block = block),
             block_var = ~ block, keep_intra_run = keep_intra_run)
}

# A dataset whose patterns in one corner of the volume are driven by a
# 3-d latent per item, so dist(latent) is a true predictor RDM there.
make_rsa_dataset <- function(n = 24, dims = c(6, 6, 6), seed = 2, snr = 3) {
  set.seed(seed)
  latent <- matrix(rnorm(n * 3), n, 3)
  arr <- array(rnorm(prod(dims) * n), c(dims, n))
  region <- expand.grid(x = 1:3, y = 1:3, z = 1:3)
  W <- matrix(rnorm(3 * nrow(region)), 3, nrow(region))
  signal <- latent %*% W * snr
  for (k in seq_len(nrow(region))) {
    arr[region$x[k], region$y[k], region$z[k], ] <-
      arr[region$x[k], region$y[k], region$z[k], ] + signal[, k]
  }
  vec  <- neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(c(dims, n)))
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, dims), neuroim2::NeuroSpace(dims))
  list(dataset = mvpa_dataset(train_data = vec, mask = mask), latent = latent,
       region_lin = region$x + (region$y - 1) * dims[1] + (region$z - 1) * dims[1] * dims[2])
}

rows_of <- function(ds) {
  d <- dim(ds$train_data)
  d[length(d)]
}

# ---------------------------------------------------------------------------
# permute_labels.rsa_design
# ---------------------------------------------------------------------------

test_that("permute_labels.rsa_design stores a within-block item permutation and nothing else changes", {
  rdes <- make_rsa_design()
  pd   <- permute_labels(rdes, seed = 11L)

  expect_s3_class(pd, "rsa_design")
  expect_equal(sort(pd$item_perm), seq_len(24))
  expect_false(identical(pd$item_perm, seq_len(24)))
  # every item is relabelled to an item of the same run
  expect_equal(rdes$block_var[pd$item_perm], rdes$block_var)

  # design matrix, include mask, and block structure untouched
  expect_identical(pd$model_mat, rdes$model_mat)
  expect_identical(pd$include, rdes$include)
  expect_identical(pd$block_var, rdes$block_var)
  expect_identical(pd$data, rdes$data)

  # reproducible under a seed, different across seeds
  expect_identical(permute_labels(rdes, seed = 11L)$item_perm, pd$item_perm)
  expect_false(identical(permute_labels(rdes, seed = 12L)$item_perm, pd$item_perm))
})

test_that("permute_labels.rsa_design circular_shift shifts within each run", {
  rdes <- make_rsa_design()
  pd   <- permute_labels(rdes, method = "circular_shift", seed = 3L)
  expect_equal(sort(pd$item_perm), seq_len(24))
  expect_equal(rdes$block_var[pd$item_perm], rdes$block_var)
  for (blk in 1:3) {
    pos <- which(rdes$block_var == blk)
    # a cyclic shift of pos: consecutive differences are all 1 except one wrap
    shifted <- pd$item_perm[pos]
    expect_true(all(diff(shifted) %in% c(1L, 1L - length(pos))))
    expect_false(identical(shifted, pos))
  }
})

test_that("permute_labels.rsa_design refuses global shuffling when within-run pairs are excluded", {
  rdes <- make_rsa_design()
  expect_error(permute_labels(rdes, method = "global", seed = 1L),
               "global.*not valid")

  # allowed when the design keeps within-run pairs (no include mask)
  rdes_keep <- make_rsa_design(keep_intra_run = TRUE)
  expect_null(rdes_keep$include)
  pd <- permute_labels(rdes_keep, method = "global", seed = 1L)
  expect_equal(sort(pd$item_perm), seq_len(24))
})

test_that("permute_labels.rsa_design without blocks permutes freely", {
  set.seed(5)
  a <- dist(matrix(rnorm(10 * 3), 10, 3))
  rdes <- rsa_design(~ a, list(a = a))
  pd <- permute_labels(rdes, seed = 2L)
  expect_equal(sort(pd$item_perm), seq_len(10))
})

test_that("permute_labels errors when no within-block shuffle can move anything", {
  set.seed(6)
  a <- dist(matrix(rnorm(6 * 3), 6, 3))
  # six items in six runs: nothing can move within a run, but no pair is
  # excluded either, so items are exchangeable across runs
  rdes <- rsa_design(~ a, list(a = a, block = 1:6), block_var = ~ block)
  expect_true(all(rdes$include))
  expect_error(permute_labels(rdes, seed = 1L), "single item.*use method = 'global'")
  expect_error(permute_labels(rdes, method = "circular_shift", seed = 1L), "single item")
  pd <- permute_labels(rdes, method = "global", seed = 1L)
  expect_equal(sort(pd$item_perm), 1:6)

  # a chance identity draw is a legitimate null sample and must not warn
  b <- dist(matrix(rnorm(4 * 3), 4, 3))
  rdes2 <- rsa_design(~ b, list(b = b, block = c(1, 1, 2, 2)), block_var = ~ block)
  draws <- replicate(60, permute_labels(rdes2)$item_perm, simplify = FALSE)
  expect_true(any(vapply(draws, identical, logical(1), 1:4)))
  expect_no_warning(replicate(60, permute_labels(rdes2)$item_perm))
})

test_that("permute_labels has no method for foreign design classes", {
  expect_error(permute_labels(structure(list(), class = c("vector_rsa_design", "list"))),
               "no method.*vector_rsa_design")
})

# ---------------------------------------------------------------------------
# train_model.rsa_model applies item_perm
# ---------------------------------------------------------------------------

test_that("item_perm on the design equals fitting on row-permuted patterns, every regtype", {
  rdes <- make_rsa_design()
  fx   <- make_rsa_dataset()
  X    <- matrix(rnorm(24 * 30), 24, 30)
  pd   <- permute_labels(rdes, seed = 21L)

  specs <- list(
    list(regtype = "pearson"),
    list(regtype = "spearman"),
    list(regtype = "lm"),
    list(regtype = "lm", semipartial = TRUE),
    list(regtype = "lm", distmethod = "pearson")
  )
  for (sp in specs) {
    args  <- c(list(dataset = fx$dataset, design = rdes, check_collinearity = FALSE), sp)
    model <- do.call(rsa_model, args)
    perm_model        <- model
    perm_model$design <- pd

    via_hook   <- train_model(perm_model, X, y = NULL, indices = NULL)
    via_rows   <- train_model(model, X[pd$item_perm, , drop = FALSE], y = NULL, indices = NULL)
    unpermuted <- train_model(model, X, y = NULL, indices = NULL)

    expect_equal(via_hook, via_rows, info = paste(names(sp), sp, collapse = " "))
    expect_false(isTRUE(all.equal(via_hook, unpermuted)),
                 info = paste("permutation must change the statistic:", sp$regtype))
  }
})

test_that("brain-side item permutation equals inverse permutation of the model RDMs", {
  # This is the equivalence the docs claim: relabelling the neural patterns by
  # pi is the same statistic as permuting rows/cols of every model RDM by
  # pi^{-1}, with the same block_var and include mask.
  rdes <- make_rsa_design()
  fx   <- make_rsa_dataset()
  X    <- matrix(rnorm(24 * 30), 24, 30)
  pd   <- permute_labels(rdes, seed = 31L)
  inv  <- order(pd$item_perm)

  A <- as.matrix(rdes$data$a)[inv, inv]
  B <- as.matrix(rdes$data$b)[inv, inv]
  rdes_inv <- rsa_design(~ a + b, list(a = as.dist(A), b = as.dist(B), block = rdes$data$block),
                         block_var = ~ block)
  expect_identical(rdes_inv$include, rdes$include)

  for (regtype in c("pearson", "lm")) {
    m_hook        <- rsa_model(fx$dataset, rdes, regtype = regtype, check_collinearity = FALSE)
    m_hook$design <- pd
    m_inv         <- rsa_model(fx$dataset, rdes_inv, regtype = regtype, check_collinearity = FALSE)

    expect_equal(train_model(m_hook, X, NULL, NULL),
                 train_model(m_inv, X, NULL, NULL), tolerance = 1e-10, info = regtype)
  }
})

test_that("train_model.rsa_model rejects an item_perm that is not a permutation", {
  rdes <- make_rsa_design()
  fx   <- make_rsa_dataset()
  X    <- matrix(rnorm(24 * 30), 24, 30)
  m    <- rsa_model(fx$dataset, rdes, regtype = "pearson")
  m$design$item_perm <- c(1L, 1L, 3:24)
  expect_error(train_model(m, X, NULL, NULL), "permutation of 1..24")
  m$design$item_perm <- 1:23
  expect_error(train_model(m, X, NULL, NULL), "length 23")
})

# ---------------------------------------------------------------------------
# pair_rsa_design
# ---------------------------------------------------------------------------

test_that("pair_rsa_design within mode with a row subset applies item_perm to that subset", {
  set.seed(41)
  items   <- paste0("it", 1:8)
  row_idx <- c(2L, 4L, 6L, 8L, 10L, 12L, 14L, 16L)
  blocks  <- rep(1:2, each = 4)
  R1 <- as.matrix(dist(matrix(rnorm(8 * 3), 8, 3)))
  R2 <- as.matrix(dist(matrix(rnorm(8 * 3), 8, 3)))
  des <- pair_rsa_design(items, model = list(m1 = R1, m2 = R2),
                         block_var_a = blocks, row_idx_a = row_idx)
  expect_identical(des$block_var_b, blocks)

  pd <- permute_labels(des, seed = 7L)
  expect_s3_class(pd, "pair_rsa_design")
  expect_equal(sort(pd$item_perm), 1:8)
  expect_equal(blocks[pd$item_perm], blocks)
  expect_null(pd$item_perm_b)

  arr  <- array(rnorm(3 * 3 * 2 * 16), c(3, 3, 2, 16))
  vec  <- neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(c(3, 3, 2, 16)))
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, c(3, 3, 2)), neuroim2::NeuroSpace(c(3, 3, 2)))
  ds   <- mvpa_dataset(train_data = vec, mask = mask)
  X    <- matrix(rnorm(16 * 12), 16, 12)

  m        <- rsa_model(ds, des, regtype = "pearson")
  m_perm   <- m
  m_perm$design <- pd
  X_rows   <- X
  X_rows[row_idx, ] <- X[row_idx[pd$item_perm], ]
  expect_equal(train_model(m_perm, X, NULL, NULL), train_model(m, X_rows, NULL, NULL))
})

test_that("pair_rsa_design between mode permutes both item sets within their own blocks", {
  set.seed(42)
  items_a <- paste0("a", 1:8); items_b <- paste0("b", 1:6)
  rows_a  <- 1:8;              rows_b  <- 9:14
  blk_a   <- rep(1:2, each = 4); blk_b <- rep(1:2, each = 3)
  M <- matrix(rnorm(8 * 6), 8, 6)
  des <- pair_rsa_design(items_a, items_b, model = list(m = M), pairs = "between",
                         block_var_a = blk_a, block_var_b = blk_b,
                         row_idx_a = rows_a, row_idx_b = rows_b)
  expect_identical(des$block_var_b, blk_b)
  expect_equal(sum(des$include), 8 * 6 - 4 * 3 - 4 * 3)

  pd <- permute_labels(des, seed = 9L)
  expect_equal(sort(pd$item_perm), 1:8)
  expect_equal(sort(pd$item_perm_b), 1:6)
  expect_equal(blk_a[pd$item_perm], blk_a)
  expect_equal(blk_b[pd$item_perm_b], blk_b)

  arr  <- array(rnorm(3 * 3 * 2 * 14), c(3, 3, 2, 14))
  vec  <- neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(c(3, 3, 2, 14)))
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, c(3, 3, 2)), neuroim2::NeuroSpace(c(3, 3, 2)))
  ds   <- mvpa_dataset(train_data = vec, mask = mask)
  X    <- matrix(rnorm(14 * 12), 14, 12)

  m      <- rsa_model(ds, des, regtype = "pearson")
  m_perm <- m
  m_perm$design <- pd
  X_rows <- X
  X_rows[rows_a, ] <- X[rows_a[pd$item_perm], ]
  X_rows[rows_b, ] <- X[rows_b[pd$item_perm_b], ]
  expect_equal(train_model(m_perm, X, NULL, NULL), train_model(m, X_rows, NULL, NULL))

  # global shuffling is refused here too
  expect_error(permute_labels(des, method = "global"), "global.*not valid")

  # only one side needs to be permutable
  des_one <- pair_rsa_design(items_a, items_b, model = list(m = M), pairs = "between",
                             block_var_a = 1:8, block_var_b = blk_b,
                             row_idx_a = rows_a, row_idx_b = rows_b)
  pd_one <- permute_labels(des_one, seed = 2L)
  expect_equal(pd_one$item_perm, 1:8)
  expect_equal(blk_b[pd_one$item_perm_b], blk_b)
  # any single draw may be the identity by chance (1 in 36 here); across seeds it moves
  perms_b <- lapply(1:12, function(sd) permute_labels(des_one, seed = sd)$item_perm_b)
  expect_true(any(!vapply(perms_b, identical, logical(1), 1:6)))
  des_none <- pair_rsa_design(items_a, items_b, model = list(m = M), pairs = "between",
                              block_var_a = 1:8, block_var_b = 1:6,
                              row_idx_a = rows_a, row_idx_b = rows_b)
  expect_error(permute_labels(des_none, seed = 2L), "No within-block item permutation")
})

# ---------------------------------------------------------------------------
# performance extraction with one-row matrices
# ---------------------------------------------------------------------------

test_that(".extract_perf_values reads the requested column of a one-row performance matrix", {
  perf <- function(a, b) matrix(c(a, b), nrow = 1, dimnames = list(NULL, c("a", "b")))
  tbl  <- tibble::tibble(
    id          = c(10L, 20L, 30L),
    performance = list(perf(0.1, 0.9), perf(0.2, 0.8), NULL),
    error       = c(FALSE, FALSE, TRUE),
    error_message = c("~", "~", "boom")
  )
  expect_equal(unname(rMVPA:::.extract_perf_values(tbl, metric = "b")), c(0.9, 0.8))
  expect_equal(unname(rMVPA:::.extract_perf_values(tbl, metric = "a")), c(0.1, 0.2))
  expect_equal(names(rMVPA:::.extract_perf_values(tbl, metric = "b")), c("10", "20"))
  # inferred metric falls back to the first column; an explicit missing one does not
  expect_equal(unname(rMVPA:::.extract_perf_values(tbl, metric = NULL)), c(0.1, 0.2))
  expect_equal(unname(rMVPA:::.extract_perf_values(tbl, metric = "zzz")), c(0.1, 0.2))
  expect_true(all(is.na(rMVPA:::.extract_perf_values(tbl, metric = "zzz", strict = TRUE))))
  # the positional fallback returns one value even for a list-valued entry
  expect_equal(rMVPA:::.perf_entry_value(list(x = c(1, 2), y = 3), NULL), 1)

  mat <- rMVPA:::.extract_perf_matrix(tbl, metrics = c("b", "a"))
  expect_equal(dim(mat), c(2L, 2L))
  expect_equal(rownames(mat), c("10", "20"))
  expect_equal(colnames(mat), c("b", "a"))
  expect_equal(unname(mat[, "b"]), c(0.9, 0.8))
  expect_equal(unname(mat[, "a"]), c(0.1, 0.2))
})

test_that(".extract_perf_matrix keeps named-vector performance working and handles empties", {
  tbl <- tibble::tibble(
    id = 1:2,
    performance = list(c(Accuracy = 0.8, AUC = 0.9), c(Accuracy = 0.6, AUC = 0.7)),
    error = c(FALSE, FALSE), error_message = c("~", "~")
  )
  mat <- rMVPA:::.extract_perf_matrix(tbl, metrics = c("AUC", "Accuracy"))
  expect_equal(unname(mat[, "AUC"]), c(0.9, 0.7))
  expect_equal(unname(mat[, "Accuracy"]), c(0.8, 0.6))
  empty <- rMVPA:::.extract_perf_matrix(tbl[0, ], metrics = "AUC")
  expect_equal(dim(empty), c(0L, 1L))
  expect_equal(dim(rMVPA:::.extract_perf_matrix(NULL, metrics = c("x", "y"))), c(0L, 2L))
})

# ---------------------------------------------------------------------------
# end to end
# ---------------------------------------------------------------------------

test_that("run_permutation_searchlight on an rsa_model scores every model predictor against one null", {
  skip_on_cran()
  quiet_logger()

  fx <- make_rsa_dataset()
  n  <- rows_of(fx$dataset)
  set.seed(3)
  a <- dist(fx$latent)                               # true structure in the signal region
  b <- dist(matrix(rnorm(n * 3), n, 3))              # unrelated
  block <- rep(1:3, each = n / 3)
  rdes  <- rsa_design(~ a + b, list(a = a, b = b, block = block), block_var = ~ block)
  mspec <- rsa_model(fx$dataset, rdes, regtype = "pearson")

  pc  <- permutation_control(n_perm = 12, subsample = 0.5, seed = 100L,
                             null_method = "global", diagnose = FALSE)
  res <- run_permutation_searchlight(mspec, radius = 2, perm_ctrl = pc)

  expect_s3_class(res, "permutation_result_set")
  expect_named(res, c("a", "b"))
  for (m in c("a", "b")) {
    r <- res[[m]]
    expect_s3_class(r, "permutation_result")
    expect_equal(r$metric, m)
    expect_equal(length(r$p_values), length(r$all_ids))
    expect_true(all(r$p_values >= 0 & r$p_values <= 1, na.rm = TRUE))
    expect_equal(r$n_perm_used, 12L)
    expect_true(r$n_null_vals > 0)
    expect_false(is.null(r$p_map))
  }
  # both metrics were scored against the same permutations
  expect_equal(res$a$n_null_vals, res$b$n_null_vals)

  # the true predictor is detected in the signal corner; the unrelated one is not
  in_region <- res$a$all_ids %in% fx$region_lin
  expect_true(median(res$a$p_values[in_region]) < 0.05)
  expect_true(mean(res$a$p_values[in_region] < 0.05) > mean(res$b$p_values[in_region] < 0.05))
  # under the null, raw p-values for the unrelated predictor are not concentrated
  expect_true(mean(res$b$p_values < 0.05, na.rm = TRUE) < 0.25)

  expect_output(print(res), "Permutation Searchlight Results \\(2 metrics\\)")
  expect_output(print(res), "adj p<0.05")
})

test_that("run_permutation_searchlight on an rsa_model honours a single named metric and rejects unknown ones", {
  skip_on_cran()
  quiet_logger()

  fx <- make_rsa_dataset(dims = c(5, 5, 5))
  n  <- rows_of(fx$dataset)
  set.seed(4)
  a <- dist(fx$latent); b <- dist(matrix(rnorm(n * 3), n, 3))
  rdes  <- rsa_design(~ a + b, list(a = a, b = b))
  mspec <- rsa_model(fx$dataset, rdes, regtype = "lm", check_collinearity = FALSE)

  pc  <- permutation_control(n_perm = 4, subsample = 0.4, seed = 5L,
                             null_method = "global", diagnose = FALSE)
  res <- run_permutation_searchlight(mspec, radius = 2, perm_ctrl = pc, metric = "b")
  expect_s3_class(res, "permutation_result")
  expect_equal(res$metric, "b")

  expect_error(
    run_permutation_searchlight(mspec, radius = 2, perm_ctrl = pc, metric = "nope"),
    "nope"
  )

  # an observed result that lacks one of the spec's predictors must not be
  # scored by silently reusing the first map
  rdes_a <- rsa_design(~ a, list(a = a))
  obs_a  <- run_searchlight(rsa_model(fx$dataset, rdes_a, regtype = "lm"), radius = 2)
  expect_named(obs_a$results, "a")
  expect_error(
    run_permutation_searchlight(mspec, observed = obs_a, radius = 2, perm_ctrl = pc),
    "metric 'b' not found"
  )
  expect_error(
    run_permutation_searchlight(mspec, observed = setNames(rnorm(10), 1:10),
                                radius = 2, perm_ctrl = pc),
    "single metric"
  )
})

test_that("run_permutation_searchlight scores several metrics of an mvpa_model against one null", {
  skip_on_cran()
  quiet_logger()

  dset  <- gen_sample_dataset(c(5, 5, 5), 60, blocks = 3, nlevels = 2)
  mdl   <- load_model("corclass")
  cval  <- blocked_cross_validation(dset$design$block_var)
  mspec <- mvpa_model(mdl, dset$dataset, dset$design, "classification", crossval = cval)
  pc    <- permutation_control(n_perm = 2, subsample = 0.3, seed = 42L,
                               null_method = "global", diagnose = FALSE)

  res <- suppressWarnings(
    run_permutation_searchlight(mspec, radius = 2, perm_ctrl = pc, metric = c("Accuracy", "AUC"))
  )
  expect_s3_class(res, "permutation_result_set")
  expect_named(res, c("Accuracy", "AUC"))
  expect_equal(res$Accuracy$n_null_vals, res$AUC$n_null_vals)
  expect_false(isTRUE(all.equal(res$Accuracy$p_values, res$AUC$p_values)))
  expect_output(summary(res), "\\[AUC\\]")

  # default for an mvpa_model is still the single first metric
  res1 <- suppressWarnings(run_permutation_searchlight(mspec, radius = 2, perm_ctrl = pc))
  expect_s3_class(res1, "permutation_result")
  expect_equal(res1$metric, "Accuracy")
})

test_that("run_permutation_searchlight on an rsa_model works with the full-brain strategy", {
  skip_on_cran()
  quiet_logger()

  fx <- make_rsa_dataset(dims = c(5, 5, 5))
  n  <- rows_of(fx$dataset)
  set.seed(8)
  a <- dist(fx$latent); b <- dist(matrix(rnorm(n * 3), n, 3))
  block <- rep(1:3, each = n / 3)
  rdes  <- rsa_design(~ a + b, list(a = a, b = b, block = block), block_var = ~ block)
  mspec <- rsa_model(fx$dataset, rdes, regtype = "spearman")

  pc  <- permutation_control(n_perm = 3, perm_strategy = "searchlight", seed = 9L,
                             null_method = "global", diagnose = FALSE)
  res <- run_permutation_searchlight(mspec, radius = 2, perm_ctrl = pc)
  expect_s3_class(res, "permutation_result_set")
  expect_named(res, c("a", "b"))
  expect_equal(res$a$perm_strategy, "searchlight")
  expect_equal(res$a$n_null_vals, res$b$n_null_vals)
  expect_true(res$a$n_null_vals >= length(res$a$all_ids))

  # the full-brain pass must actually have applied the item permutation
  perm_spec <- mspec
  perm_spec$design <- permute_labels(rdes, seed = 10L)
  obs  <- res$a$observed
  prm  <- run_searchlight(perm_spec, radius = 2)
  ids  <- res$a$all_ids
  expect_false(isTRUE(all.equal(as.numeric(obs$results$a)[ids], as.numeric(prm$results$a)[ids])))
})

# ---------------------------------------------------------------------------
# design diagnostics
# ---------------------------------------------------------------------------

test_that("rsa_design_diagnostics reports items, pairs, VIF, and effective items", {
  set.seed(12)
  n <- 12
  shared <- matrix(rnorm(n * 4), n, 4)
  a <- dist(shared)
  b <- dist(shared + matrix(rnorm(n * 4, sd = 0.3), n, 4))   # highly collinear with a
  c <- dist(matrix(rnorm(n * 4), n, 4))                       # independent
  des <- rsa_design(~ a + b + c, list(a = a, b = b, c = c))
  d <- rsa_design_diagnostics(des)

  expect_s3_class(d, "rsa_design_diagnostics")
  expect_equal(d$n_items, 12L)
  expect_equal(d$n_pairs, choose(12, 2))
  expect_equal(d$n_predictors, 3L)
  expect_equal(names(d$vif), c("a", "b", "c"))
  expect_equal(sort(d$max_abs_cor_pair), c("a", "b"))
  expect_true(d$max_abs_cor > 0.9)
  expect_true(all(d$vif >= 1))
  expect_true(d$vif[["a"]] > 3 && d$vif[["b"]] > 3)
  expect_true(d$vif[["c"]] < 1.5)
  expect_equal(d$items_per_predictor, d$n_items / d$vif)
  expect_true(d$items_per_predictor[["c"]] > d$threshold)
  expect_true(d$items_per_predictor[["a"]] < d$threshold)
  expect_output(print(d), "below threshold")

  # VIF matches the textbook definition 1 / (1 - R^2_j)
  X <- do.call(cbind, des$model_mat)
  r2_a <- summary(lm(X[, "a"] ~ X[, "b"] + X[, "c"]))$r.squared
  expect_equal(d$vif[["a"]], 1 / (1 - r2_a), tolerance = 1e-8)

  # single predictor: VIF is 1 and no correlation is reported
  d1 <- rsa_design_diagnostics(rsa_design(~ a, list(a = a)))
  expect_equal(unname(d1$vif), 1)
  expect_true(is.na(d1$max_abs_cor))

  # a constant predictor gets NA on its own; the others keep finite VIFs
  cv <- as.dist(matrix(0, n, n))
  d_const <- rsa_design_diagnostics(rsa_design(~ a + b + cv, list(a = a, b = b, cv = cv)))
  expect_true(is.na(d_const$vif[["cv"]]))
  expect_true(all(is.finite(d_const$vif[c("a", "b")])))
  expect_equal(d_const$vif[c("a", "b")], rsa_design_diagnostics(rsa_design(~ a + b, list(a = a, b = b)))$vif,
               tolerance = 1e-8)

  # one NA cell in one RDM does not blank the other predictors
  a_na <- a; a_na[3] <- NA
  d_na <- rsa_design_diagnostics(rsa_design(~ a_na + c, list(a_na = a_na, c = c)))
  expect_true(all(is.finite(d_na$vif)))
})

test_that("rsa_design_diagnostics counts items for pair designs and excluded pairs", {
  set.seed(13)
  items <- paste0("i", 1:10)
  R1 <- as.matrix(dist(matrix(rnorm(10 * 3), 10, 3)))
  des_w <- pair_rsa_design(items, model = list(m = R1), block_var_a = rep(1:2, each = 5))
  d_w <- rsa_design_diagnostics(des_w)
  expect_equal(d_w$n_items, 10L)
  expect_equal(d_w$n_pairs, sum(des_w$include))

  M <- matrix(rnorm(10 * 6), 10, 6)
  des_b <- pair_rsa_design(items, paste0("j", 1:6), model = list(m = M), pairs = "between",
                           row_idx_a = 1:10, row_idx_b = 11:16)
  d_b <- rsa_design_diagnostics(des_b)
  expect_equal(d_b$n_items, 16L)
  expect_equal(d_b$n_pairs, 60L)
})

test_that("rsa_model warns on weak effective support only for regression regtypes", {
  set.seed(14)
  n <- 12
  shared <- matrix(rnorm(n * 4), n, 4)
  a <- dist(shared)
  b <- dist(shared + matrix(rnorm(n * 4, sd = 0.3), n, 4))
  des <- rsa_design(~ a + b, list(a = a, b = b))
  expect_true(rsa_design_diagnostics(des)$max_abs_cor < 0.99)   # passes the hard check

  arr  <- array(rnorm(3 * 3 * 2 * n), c(3, 3, 2, n))
  vec  <- neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(c(3, 3, 2, n)))
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, c(3, 3, 2)), neuroim2::NeuroSpace(c(3, 3, 2)))
  ds   <- mvpa_dataset(train_data = vec, mask = mask)

  expect_warning(suppressMessages(rsa_model(ds, des, regtype = "lm")), "effective items")
  expect_warning(rsa_model(ds, des, regtype = "rfit"), "effective items")
  expect_no_warning(suppressMessages(rsa_model(ds, des, regtype = "lm", check_collinearity = FALSE)))
  expect_no_warning(rsa_model(ds, des, regtype = "pearson"))
  expect_no_warning(rsa_model(ds, des, regtype = "spearman"))

  m <- rsa_model(ds, des, regtype = "pearson")
  expect_s3_class(m$design_diagnostics, "rsa_design_diagnostics")
  expect_equal(m$design_diagnostics$n_items, n)
  expect_output(print(m), "Effective items per predictor")
})

test_that("rsa_model's warning names the actual limit: items or collinearity", {
  set.seed(16)
  n <- 8
  a <- dist(matrix(rnorm(n * 4), n, 4))
  des1 <- rsa_design(~ a, list(a = a))
  arr  <- array(rnorm(3 * 3 * 2 * n), c(3, 3, 2, n))
  vec  <- neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(c(3, 3, 2, n)))
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, c(3, 3, 2)), neuroim2::NeuroSpace(c(3, 3, 2)))
  ds   <- mvpa_dataset(train_data = vec, mask = mask)
  w1 <- capture_warnings(suppressMessages(rsa_model(ds, des1, regtype = "lm")))
  expect_match(w1, "only 8 items")
  expect_match(w1, "Add items or pool runs")
  expect_no_match(w1, "collinear")

  n <- 12
  shared <- matrix(rnorm(n * 4), n, 4)
  a2 <- dist(shared); b2 <- dist(shared + matrix(rnorm(n * 4, sd = 0.3), n, 4))
  des2 <- rsa_design(~ a2 + b2, list(a2 = a2, b2 = b2))
  arr2 <- array(rnorm(3 * 3 * 2 * n), c(3, 3, 2, n))
  ds2  <- mvpa_dataset(train_data = neuroim2::NeuroVec(arr2, neuroim2::NeuroSpace(c(3, 3, 2, n))), mask = mask)
  w2 <- capture_warnings(suppressMessages(rsa_model(ds2, des2, regtype = "lm")))
  expect_match(w2, "less collinear predictors")
  expect_match(w2, "max \\|r\\| between predictors")
})

test_that("rsa_model stays silent when the design supports its predictors", {
  set.seed(15)
  n <- 30
  a <- dist(matrix(rnorm(n * 4), n, 4))
  b <- dist(matrix(rnorm(n * 4), n, 4))
  des <- rsa_design(~ a + b, list(a = a, b = b))
  arr  <- array(rnorm(3 * 3 * 2 * n), c(3, 3, 2, n))
  vec  <- neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(c(3, 3, 2, n)))
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, c(3, 3, 2)), neuroim2::NeuroSpace(c(3, 3, 2)))
  ds   <- mvpa_dataset(train_data = vec, mask = mask)
  expect_no_warning(suppressMessages(rsa_model(ds, des, regtype = "lm")))
  d <- rsa_model(ds, des, regtype = "lm")$design_diagnostics
  expect_true(all(d$items_per_predictor > d$threshold))
})
