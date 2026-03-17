# Tests for R/permutation_searchlight.R
# Functions tested:
#   permutation_control(), permute_labels(), subsample_centers(),
#   compute_local_redundancy(), run_permutation_searchlight()
# Internal helpers (via rMVPA:::):
#   build_adjusted_null(), score_observed(), diagnose_null(),
#   .extract_perf_values(), .align_to_ids()

# ---------------------------------------------------------------------------
# 1. permutation_control validates inputs
# ---------------------------------------------------------------------------

test_that("permutation_control validates inputs", {
  # default works
  pc <- permutation_control()
  expect_s3_class(pc, "permutation_control")
  expect_equal(pc$n_perm, 5L)
  expect_equal(pc$shuffle, "within_block")
  expect_equal(pc$perm_strategy, "iterate")

  pc2 <- permutation_control(perm_strategy = "searchlight")
  expect_equal(pc2$perm_strategy, "searchlight")

  # bad n_perm
  expect_error(permutation_control(n_perm = 0))
  expect_error(permutation_control(n_perm = -1))

  # bad n_bins
  expect_error(permutation_control(n_bins = 1))

  # bad subsample
  expect_error(permutation_control(subsample = 0))
  expect_error(permutation_control(subsample = -1))
  expect_error(permutation_control(perm_strategy = "not-a-strategy"))

  # print method works (call method directly; S3 dispatch requires NAMESPACE registration)
  expect_output(rMVPA:::print.permutation_control(pc), "Permutation Control")
  expect_output(rMVPA:::print.permutation_control(pc2),
                "ignored - full brain per perm")
})


# ---------------------------------------------------------------------------
# 2. permute_labels within_block
# ---------------------------------------------------------------------------

test_that("permute_labels within_block preserves structure", {
  dset   <- gen_sample_dataset(c(5, 5, 5), 60, blocks = 3, nlevels = 3)
  design <- dset$design

  perm <- suppressWarnings(
    permute_labels(design, method = "within_block", seed = 42)
  )

  # Class preserved
  expect_s3_class(perm, "mvpa_design")

  # block_var unchanged
  expect_identical(perm$block_var, design$block_var)

  # labels changed (with high probability given seed)
  expect_false(identical(perm$cv_labels, design$cv_labels))

  # Overall label frequencies preserved
  expect_equal(sort(table(perm$cv_labels)), sort(table(design$cv_labels)))

  # Within each block, label frequencies preserved
  for (blk in unique(design$block_var)) {
    idx <- which(design$block_var == blk)
    expect_equal(sort(table(perm$cv_labels[idx])),
                 sort(table(design$cv_labels[idx])))
  }

  # Reproducibility with same seed
  perm2 <- suppressWarnings(
    permute_labels(design, method = "within_block", seed = 42)
  )
  expect_identical(perm2$cv_labels, perm$cv_labels)
})


# ---------------------------------------------------------------------------
# 3. permute_labels circular_shift
# ---------------------------------------------------------------------------

test_that("permute_labels circular_shift shifts within blocks", {
  dset   <- gen_sample_dataset(c(5, 5, 5), 60, blocks = 3, nlevels = 3)
  design <- dset$design

  perm <- suppressWarnings(
    permute_labels(design, method = "circular_shift", seed = 123)
  )

  expect_s3_class(perm, "mvpa_design")
  expect_identical(perm$block_var, design$block_var)

  # Within each block, the labels are a circular rotation
  for (blk in unique(design$block_var)) {
    idx     <- which(design$block_var == blk)
    orig    <- design$cv_labels[idx]
    shifted <- perm$cv_labels[idx]
    n       <- length(orig)
    if (n > 1) {
      doubled <- c(as.character(orig), as.character(orig))
      found <- FALSE
      for (s in 1:(n - 1)) {
        candidate <- doubled[(s + 1):(s + n)]
        if (identical(candidate, as.character(shifted))) {
          found <- TRUE
          break
        }
      }
      expect_true(found,
                  info = paste("Block", blk, "should be a circular shift"))
    }
  }
})


# ---------------------------------------------------------------------------
# 4. permute_labels global
# ---------------------------------------------------------------------------

test_that("permute_labels global shuffles all labels", {
  dset   <- gen_sample_dataset(c(5, 5, 5), 60, blocks = 3, nlevels = 3)
  design <- dset$design

  perm <- suppressWarnings(
    permute_labels(design, method = "global", seed = 99)
  )

  expect_s3_class(perm, "mvpa_design")
  expect_identical(perm$block_var, design$block_var)
  expect_false(identical(perm$cv_labels, design$cv_labels))

  # Overall frequencies preserved
  expect_equal(sort(table(perm$cv_labels)), sort(table(design$cv_labels)))
})


# ---------------------------------------------------------------------------
# 5. subsample_centers returns correct count and structure
# ---------------------------------------------------------------------------

test_that("subsample_centers returns correct structure", {
  dset <- gen_sample_dataset(c(5, 5, 5), 30, blocks = 3)
  sl   <- get_searchlight(dset$dataset, "standard", 3)

  sub <- suppressWarnings(
    subsample_centers(dset$dataset, sl, fraction = 0.3, seed = 1)
  )

  expect_type(sub, "list")
  expect_true(all(c("center_ids", "vox_list", "covariates") %in% names(sub)))
  expect_equal(length(sub$center_ids), length(sub$vox_list))
  expect_equal(nrow(sub$covariates), length(sub$center_ids))
  expect_true("nfeatures" %in% names(sub$covariates))

  # Check that n_centers is approximately right
  all_ids    <- get_center_ids(dset$dataset)
  expected_n <- max(1, round(length(all_ids) * 0.3))
  # Allow some tolerance due to stratification
  expect_true(
    abs(length(sub$center_ids) - expected_n) <= max(5, expected_n * 0.3)
  )

  # All selected IDs should be valid center IDs
  expect_true(all(sub$center_ids %in% all_ids))

  # When fraction >= 1, return all
  sub_all <- suppressWarnings(
    subsample_centers(dset$dataset, sl, fraction = 1.0, seed = 1)
  )
  expect_equal(length(sub_all$center_ids), length(all_ids))
})


# ---------------------------------------------------------------------------
# 6. build_adjusted_null bins correctly
# ---------------------------------------------------------------------------

test_that("build_adjusted_null creates valid bins", {
  set.seed(42)
  n        <- 500
  nf       <- sample(10:100, n, replace = TRUE)
  null_vals <- rnorm(n, mean = 0.5, sd = 0.1)
  cov      <- data.frame(nfeatures = nf)

  adj <- rMVPA:::build_adjusted_null(null_vals, cov, n_bins = 5,
                                     method = "adjusted")

  expect_s3_class(adj, "adjusted_null")
  expect_equal(adj$method, "adjusted")
  expect_true(adj$n_bins >= 2)

  # Each bin should have some values
  for (b in seq_len(adj$n_bins)) {
    expect_true(adj$bin_stats[[b]]$n > 0)
  }

  # Total values should match
  total_in_bins <- sum(vapply(adj$bin_stats, function(x) x$n, integer(1)))
  expect_equal(total_in_bins, n)

  # Global method produces 1 bin
  adj_global <- rMVPA:::build_adjusted_null(null_vals, cov, method = "global")
  expect_equal(adj_global$n_bins, 1L)
  expect_equal(adj_global$bin_stats[[1]]$n, n)

  # print works (call method directly; S3 dispatch requires NAMESPACE registration)
  expect_output(rMVPA:::print.adjusted_null(adj), "Adjusted Null")
})


# ---------------------------------------------------------------------------
# 7. score_observed produces valid p-values
# ---------------------------------------------------------------------------

test_that("score_observed produces valid p-values", {
  set.seed(42)
  n_null   <- 200
  nf_null  <- sample(10:50, n_null, replace = TRUE)
  null_vals <- rnorm(n_null, mean = 0.5, sd = 0.1)
  null_cov  <- data.frame(nfeatures = nf_null)

  adj <- rMVPA:::build_adjusted_null(null_vals, null_cov, n_bins = 3,
                                     method = "adjusted")

  # Test p-values for some observed values
  n_obs   <- 50
  obs_vals <- rnorm(n_obs, mean = 0.5, sd = 0.15)
  obs_cov  <- data.frame(nfeatures = sample(10:50, n_obs, replace = TRUE))

  p <- rMVPA:::score_observed(obs_vals, adj, obs_cov)

  expect_length(p, n_obs)
  expect_true(all(p >= 0 & p <= 1))

  # Extreme observed value should give small p
  extreme_obs <- c(100)
  extreme_cov <- data.frame(nfeatures = 30L)
  p_extreme   <- rMVPA:::score_observed(extreme_obs, adj, extreme_cov)
  expect_true(p_extreme[1] < 0.1)

  # Value well below null mean should give p near 1
  low_obs <- c(-100)
  p_low   <- rMVPA:::score_observed(low_obs, adj, extreme_cov)
  expect_true(p_low[1] > 0.9)

  # NA handling
  na_obs <- c(NA_real_)
  p_na   <- rMVPA:::score_observed(na_obs, adj, extreme_cov)
  expect_true(is.na(p_na[1]))
})


# ---------------------------------------------------------------------------
# 8. diagnose_null returns correct structure and detects bias
# ---------------------------------------------------------------------------

test_that("diagnose_null returns correct structure and detects bias", {
  set.seed(42)
  n  <- 300
  nf <- seq(10, 100, length.out = n)

  # Case 1: null correlated with nfeatures (bias)
  null_biased <- 0.5 + 0.005 * nf + rnorm(n, 0, 0.01)
  cov         <- data.frame(nfeatures = nf)

  diag <- rMVPA:::diagnose_null(null_biased, cov, n_perm = 5)

  expect_s3_class(diag, "null_diagnostics")
  expect_true("nfeatures_cor" %in% names(diag$checks))
  expect_true(diag$checks$nfeatures_cor$flagged)

  # Case 2: null NOT correlated
  null_clean <- rnorm(n, 0.5, 0.1)
  diag2      <- rMVPA:::diagnose_null(null_clean, cov, n_perm = 5)
  expect_false(diag2$checks$nfeatures_cor$flagged)

  # print works (call method directly; S3 dispatch requires NAMESPACE registration)
  expect_output(rMVPA:::print.null_diagnostics(diag), "Null Distribution Diagnostics")
})


# ---------------------------------------------------------------------------
# 9. Integration: run_permutation_searchlight end-to-end
# ---------------------------------------------------------------------------

test_that("run_permutation_searchlight runs end-to-end", {
  skip_on_cran()

  futile.logger::flog.threshold(futile.logger::ERROR)
  on.exit(futile.logger::flog.threshold(futile.logger::INFO), add = TRUE)

  dset  <- gen_sample_dataset(c(5, 5, 5), 100, blocks = 3, nlevels = 2)
  mdl   <- load_model("corclass")
  cval  <- blocked_cross_validation(dset$design$block_var)
  mspec <- mvpa_model(mdl, dset$dataset, dset$design, "classification",
                      crossval = cval)

  pc <- permutation_control(n_perm = 2, subsample = 0.3, seed = 42L,
                             null_method = "global", correction = "fdr",
                             diagnose = TRUE)

  res <- suppressWarnings(
    run_permutation_searchlight(mspec, radius = 3, perm_ctrl = pc)
  )

  expect_s3_class(res, "permutation_result")
  expect_true(length(res$p_values) > 0)
  expect_true(all(res$p_values >= 0 & res$p_values <= 1, na.rm = TRUE))
  expect_true(length(res$p_adjusted) == length(res$p_values))
  expect_equal(res$n_perm_used, 2L)
  expect_true(res$n_null_vals > 0)
  expect_true(!is.null(res$observed))
  expect_true(!is.null(res$perm_ctrl))

  # print works
  expect_output(print(res), "Permutation Searchlight Result")
  # summary works
  expect_output(summary(res), "Permutation Searchlight Summary")
})


# ---------------------------------------------------------------------------
# 10. Integration with pre-computed observed result
# ---------------------------------------------------------------------------

test_that("run_permutation_searchlight accepts pre-computed observed", {
  skip_on_cran()

  futile.logger::flog.threshold(futile.logger::ERROR)
  on.exit(futile.logger::flog.threshold(futile.logger::INFO), add = TRUE)

  dset  <- gen_sample_dataset(c(5, 5, 5), 100, blocks = 3, nlevels = 2)
  mdl   <- load_model("corclass")
  cval  <- blocked_cross_validation(dset$design$block_var)
  mspec <- mvpa_model(mdl, dset$dataset, dset$design, "classification",
                      crossval = cval)

  # First run the observed searchlight
  obs <- suppressWarnings(
    run_searchlight(mspec, radius = 3, method = "standard")
  )

  # Then pass it in
  pc <- permutation_control(n_perm = 2, subsample = 0.3, seed = 1L,
                             null_method = "adjusted", n_bins = 3,
                             correction = "fdr", diagnose = FALSE)

  res <- suppressWarnings(
    run_permutation_searchlight(mspec, observed = obs, radius = 3,
                                perm_ctrl = pc)
  )

  expect_s3_class(res, "permutation_result")
  expect_true(length(res$p_values) > 0)
  expect_true(all(res$p_values >= 0 & res$p_values <= 1, na.rm = TRUE))
})


# ---------------------------------------------------------------------------
# 11. compute_local_redundancy returns valid values
# ---------------------------------------------------------------------------

test_that("compute_local_redundancy returns valid values", {
  dset <- gen_sample_dataset(c(5, 5, 5), 30, blocks = 3)
  sl   <- get_searchlight(dset$dataset, "standard", 3)

  red <- suppressWarnings(
    compute_local_redundancy(dset$dataset, sl)
  )

  # Should be a numeric vector
  expect_type(red, "double")

  # Length should match number of searchlight centers
  expect_equal(length(red), length(sl))

  # Values should be in [0, 1] (mean absolute correlation)
  expect_true(all(red >= 0 & red <= 1, na.rm = TRUE))

  # Should have names (one per center id)
  expect_false(is.null(names(red)))
  expect_equal(length(names(red)), length(sl))
})


# ---------------------------------------------------------------------------
# 12. run_permutation_searchlight with redundancy adjustment
# ---------------------------------------------------------------------------

test_that("run_permutation_searchlight works with adjust_by = redundancy", {
  skip_on_cran()

  futile.logger::flog.threshold(futile.logger::ERROR)
  on.exit(futile.logger::flog.threshold(futile.logger::INFO), add = TRUE)

  dset  <- gen_sample_dataset(c(5, 5, 5), 100, blocks = 3, nlevels = 2)
  mdl   <- load_model("corclass")
  cval  <- blocked_cross_validation(dset$design$block_var)
  mspec <- mvpa_model(mdl, dset$dataset, dset$design, "classification",
                      crossval = cval)

  # adjust_by = "redundancy" causes compute_local_redundancy to be called
  # inside run_permutation_searchlight before building the null
  pc <- permutation_control(n_perm = 2, subsample = 0.3, seed = 7L,
                             null_method = "global",
                             adjust_by   = "redundancy",
                             correction  = "fdr",
                             diagnose    = FALSE)

  res <- suppressWarnings(
    run_permutation_searchlight(mspec, radius = 3, perm_ctrl = pc)
  )

  expect_s3_class(res, "permutation_result")
  expect_true(length(res$p_values) > 0)
  expect_true(all(res$p_values >= 0 & res$p_values <= 1, na.rm = TRUE))
  expect_equal(res$n_perm_used, 2L)
})


# ---------------------------------------------------------------------------
# 13. .extract_perf_values and .align_to_ids internal helpers
# ---------------------------------------------------------------------------

test_that(".extract_perf_values extracts named metric values correctly", {
  mock_tbl <- tibble::tibble(
    id          = c(1L, 2L, 3L),
    performance = list(
      c(Accuracy = 0.8, AUC = 0.9),
      c(Accuracy = 0.6, AUC = 0.7),
      NULL
    ),
    error         = c(FALSE, FALSE, TRUE),
    error_message = c("~", "~", "some error")
  )

  # Named metric extraction
  vals <- rMVPA:::.extract_perf_values(mock_tbl, metric = "Accuracy")
  expect_equal(length(vals), 2L)
  expect_equal(unname(vals), c(0.8, 0.6))
  expect_equal(names(vals), c("1", "2"))

  # First-metric default (metric = NULL)
  vals2 <- rMVPA:::.extract_perf_values(mock_tbl, metric = NULL)
  expect_equal(unname(vals2), c(0.8, 0.6))

  # Second named metric
  vals3 <- rMVPA:::.extract_perf_values(mock_tbl, metric = "AUC")
  expect_equal(unname(vals3), c(0.9, 0.7))
})

test_that(".extract_perf_values returns empty vector for empty/NULL input", {
  expect_equal(length(rMVPA:::.extract_perf_values(NULL)), 0L)

  empty_tbl <- tibble::tibble(
    id          = integer(0),
    performance = list(),
    error       = logical(0),
    error_message = character(0)
  )
  expect_equal(length(rMVPA:::.extract_perf_values(empty_tbl)), 0L)
})

test_that(".align_to_ids fills NA for missing centers and maps values by name", {
  all_ids <- c(10L, 20L, 30L, 40L)
  vals    <- c("10" = 0.8, "30" = 0.6)   # 20 and 40 missing

  aligned <- rMVPA:::.align_to_ids(vals, all_ids)

  expect_equal(length(aligned), 4L)
  expect_equal(names(aligned), as.character(all_ids))
  expect_equal(unname(aligned[["10"]]), 0.8)
  expect_true(is.na(aligned[["20"]]))
  expect_equal(unname(aligned[["30"]]), 0.6)
  expect_true(is.na(aligned[["40"]]))
})

test_that(".align_to_ids returns all-NA vector when no names match", {
  all_ids <- c(1L, 2L, 3L)
  vals    <- c("99" = 0.5, "100" = 0.9)

  aligned <- rMVPA:::.align_to_ids(vals, all_ids)
  expect_equal(length(aligned), 3L)
  expect_true(all(is.na(aligned)))
})


# ---------------------------------------------------------------------------
# 14. engine passthrough via ... does not break run_permutation_searchlight
# ---------------------------------------------------------------------------

test_that("run_permutation_searchlight forwards engine argument without error", {
  skip_on_cran()

  futile.logger::flog.threshold(futile.logger::ERROR)
  on.exit(futile.logger::flog.threshold(futile.logger::INFO), add = TRUE)

  dset  <- gen_sample_dataset(c(5, 5, 5), 100, blocks = 3, nlevels = 2)
  mdl   <- load_model("corclass")
  cval  <- blocked_cross_validation(dset$design$block_var)
  mspec <- mvpa_model(mdl, dset$dataset, dset$design, "classification",
                      crossval = cval)

  pc <- permutation_control(n_perm = 2, subsample = 0.3, seed = 1L,
                             null_method = "global", diagnose = FALSE)

  # engine = "legacy" is the safest choice that always works
  res <- suppressWarnings(
    run_permutation_searchlight(mspec, radius = 3, perm_ctrl = pc,
                                engine = "legacy")
  )

  expect_s3_class(res, "permutation_result")
  expect_true(length(res$p_values) > 0)
})


# ---------------------------------------------------------------------------
# 15. subsample_centers works with explicit integer n_centers
# ---------------------------------------------------------------------------

test_that("subsample_centers selects exactly n_centers when given an integer", {
  dset <- gen_sample_dataset(c(5, 5, 5), 30, blocks = 3)
  sl   <- get_searchlight(dset$dataset, "standard", 3)

  all_ids <- get_center_ids(dset$dataset)
  # Request fewer than total so stratified path is exercised
  target  <- min(10L, length(all_ids) - 1L)

  sub <- suppressWarnings(
    subsample_centers(dset$dataset, sl, n_centers = target, seed = 42L)
  )

  expect_type(sub, "list")
  expect_true(all(c("center_ids", "vox_list", "covariates") %in% names(sub)))

  # Exact count (after stratified rounding, may be off by at most a few; but
  # the function trims or pads to exactly n_centers)
  expect_equal(length(sub$center_ids), target)
  expect_equal(length(sub$vox_list), target)
  expect_equal(nrow(sub$covariates), target)

  # All selected IDs are valid
  expect_true(all(sub$center_ids %in% all_ids))
})

test_that("subsample_centers returns all centers when n_centers >= total", {
  dset <- gen_sample_dataset(c(5, 5, 5), 30, blocks = 3)
  sl   <- get_searchlight(dset$dataset, "standard", 3)

  all_ids <- get_center_ids(dset$dataset)
  n_total <- length(all_ids)

  sub <- suppressWarnings(
    subsample_centers(dset$dataset, sl, n_centers = n_total + 100L, seed = 1L)
  )

  expect_equal(length(sub$center_ids), n_total)
})


# ---------------------------------------------------------------------------
# 16. searchlight strategy and custom extraction paths
# ---------------------------------------------------------------------------

test_that("run_permutation_searchlight works with perm_strategy = searchlight", {
  skip_on_cran()

  futile.logger::flog.threshold(futile.logger::ERROR)
  on.exit(futile.logger::flog.threshold(futile.logger::INFO), add = TRUE)

  dset  <- gen_sample_dataset(c(5, 5, 5), 100, blocks = 3, nlevels = 2)
  mdl   <- load_model("corclass")
  cval  <- blocked_cross_validation(dset$design$block_var)
  mspec <- mvpa_model(mdl, dset$dataset, dset$design, "classification",
                      crossval = cval)

  pc <- permutation_control(
    n_perm        = 1L,
    perm_strategy = "searchlight",
    null_method   = "global",
    correction    = "none",
    diagnose      = FALSE,
    seed          = 123L
  )

  res <- suppressWarnings(
    run_permutation_searchlight(mspec, radius = 3, perm_ctrl = pc,
                                engine = "legacy")
  )

  expect_s3_class(res, "permutation_result")
  expect_equal(res$perm_strategy, "searchlight")
  expect_equal(length(res$p_values), length(res$all_ids))
  expect_true(res$n_null_vals > 0)
})

test_that("run_permutation_searchlight accepts unnamed observed vectors", {
  skip_on_cran()

  futile.logger::flog.threshold(futile.logger::ERROR)
  on.exit(futile.logger::flog.threshold(futile.logger::INFO), add = TRUE)

  dset  <- gen_sample_dataset(c(5, 5, 5), 100, blocks = 3, nlevels = 2)
  mdl   <- load_model("corclass")
  cval  <- blocked_cross_validation(dset$design$block_var)
  mspec <- mvpa_model(mdl, dset$dataset, dset$design, "classification",
                      crossval = cval)

  all_ids <- get_center_ids(dset$dataset)
  obs_vec <- rep(0.5, length(all_ids))

  pc <- permutation_control(
    n_perm      = 1L,
    null_method = "global",
    correction  = "none",
    diagnose    = FALSE,
    seed        = 99L
  )

  res <- suppressWarnings(
    run_permutation_searchlight(mspec, observed = obs_vec, radius = 3,
                                perm_ctrl = pc, engine = "legacy")
  )

  expect_s3_class(res, "permutation_result")
  expect_equal(length(res$p_values), length(all_ids))
  expect_false(all(is.na(res$p_values)))
})

test_that(".extract_values_from_searchlight_result handles custom outputs", {
  ids <- c(10L, 20L, 30L)
  map <- setNames(c(0.7, 0.8, 0.9), as.character(ids))

  sl_like <- list(results = list(Accuracy = map))
  vals <- rMVPA:::.extract_values_from_searchlight_result(
    sl_like, ids, metric = "Accuracy"
  )
  expect_equal(unname(vals), c(0.7, 0.8, 0.9))
  expect_equal(names(vals), as.character(ids))

  vals2 <- rMVPA:::.extract_values_from_searchlight_result(map, ids)
  expect_equal(unname(vals2), c(0.7, 0.8, 0.9))
  expect_equal(names(vals2), as.character(ids))
})

test_that(".align_to_ids aligns unnamed full-length vectors by position", {
  all_ids <- c(101L, 102L, 103L)
  vals    <- c(0.1, 0.2, 0.3)

  aligned <- rMVPA:::.align_to_ids(vals, all_ids)
  expect_equal(unname(aligned), vals)
  expect_equal(names(aligned), as.character(all_ids))
})
