# Two regions with different content. The left half of the volume carries the
# two phonological features, the right half the acoustic feature, which is
# itself correlated with one phonological feature. A slow drift is added
# everywhere.
.fdm_fixture <- function(n = 120L, dims = c(6L, 4L, 2L), seed = 8801L) {
  set.seed(seed)
  ar <- function() as.numeric(stats::filter(rnorm(n), 0.5, "recursive"))
  m1 <- ar(); p1 <- ar()
  ac <- 0.8 * m1 + 0.6 * ar()
  X <- cbind(m1 = m1, p1 = p1, ac = ac, lex = ar())
  drift <- cbind(drift = seq(-1, 1, length.out = n))
  grid <- expand.grid(x = seq_len(dims[1]), y = seq_len(dims[2]), z = seq_len(dims[3]))
  left <- grid$x <= dims[1] / 2
  v <- nrow(grid)
  B <- matrix(0, 4L, v, dimnames = list(colnames(X), NULL))
  B[c("m1", "p1"), left] <- rnorm(2L * sum(left))
  B["ac", !left] <- rnorm(sum(!left))
  Y <- X %*% B + 4 * drift %*% t(rnorm(v)) + matrix(rnorm(n * v, sd = 0.7), n, v)
  arr <- array(t(Y), dim = c(dims, n))
  space4 <- neuroim2::NeuroSpace(c(dims, n))
  mask <- neuroim2::NeuroVol(array(1, dims), neuroim2::NeuroSpace(dims))
  dataset <- mvpa_dataset(neuroim2::NeuroVec(arr, space4), mask = mask)
  sets <- c("manner", "place", "acoustic", "lexical")
  fs <- feature_sets(
    X, by_set(sets),
    families = list(phonology = c("manner", "place")),
    contrasts = list(manner_vs_place = c(m1 = 1, p1 = -1))
  )
  block <- rep(1:6, each = n / 6L)
  design <- feature_sets_design(fs, block_var_train = block, time_series = TRUE)
  region <- neuroim2::NeuroVol(array(ifelse(left, 1L, 2L), dims),
                               neuroim2::NeuroSpace(dims))
  list(dataset = dataset, design = design, X = X, Y = Y, drift = drift,
       left = left, region = region, block = block, sets = sets)
}

.fdm_fit <- function(model, Y, columns) {
  fit_roi(model, list(train_data = Y[, columns, drop = FALSE], test_data = NULL,
                      indices = which(columns)),
          list(design = model$design, cv_spec = NULL, id = 1L,
               center_global_id = NA))
}

test_that("feature_decode_model recovers which features a region carries", {
  f <- .fdm_fixture()
  model <- feature_decode_model(f$dataset, f$design, outer_crossval = 6L,
                                tune_crossval = 3L, purge = 2L,
                                nuisance = f$drift, partial = TRUE)
  schema <- output_schema(model)
  labels <- c("all", "manner", "place", "acoustic", "lexical", "phonology",
              "manner_vs_place")
  partial_labels <- c("manner", "place", "acoustic", "lexical", "phonology")
  expect_identical(
    names(schema),
    c(as.vector(rbind(paste0(labels, "_r2"), paste0(labels, "_cor"))),
      as.vector(rbind(paste0(partial_labels, "_partial_r2"),
                      paste0(partial_labels, "_partial_cor"))))
  )
  expect_true(all(unlist(schema) == "scalar"))

  left <- .fdm_fit(model, f$Y, f$left)
  right <- .fdm_fit(model, f$Y, !f$left)
  expect_false(left$error); expect_false(right$error)
  expect_identical(names(left$metrics), names(schema))
  L <- left$metrics; R <- right$metrics

  # left carries phonology; right carries acoustics
  expect_gt(L[["phonology_cor"]], 0.8)
  expect_gt(L[["manner_vs_place_cor"]], 0.8)
  expect_gt(R[["acoustic_cor"]], 0.8)
  expect_lt(abs(L[["lexical_cor"]]), 0.25)
  expect_lt(abs(R[["place_cor"]]), 0.25)
  expect_equal(L[["phonology_cor"]], (L[["manner_cor"]] + L[["place_cor"]]) / 2,
               tolerance = 1e-12)

  # the acoustic feature is correlated with manner, so each region "decodes"
  # the other's feature until the shared part is partialled out
  expect_gt(L[["acoustic_cor"]], 0.4)
  expect_lt(abs(L[["acoustic_partial_cor"]]), 0.25)
  expect_gt(R[["manner_cor"]], 0.4)
  expect_lt(abs(R[["manner_partial_cor"]]), 0.25)
  expect_gt(L[["manner_partial_cor"]], 0.5)
  expect_gt(R[["acoustic_partial_cor"]], 0.5)
})

test_that("feature_decode_model matches a direct nested ridge on swapped roles", {
  f <- .fdm_fixture()
  model <- feature_decode_model(f$dataset, f$design, outer_crossval = 6L,
                                tune_crossval = 3L, purge = 1L, nuisance = f$drift)
  got <- .fdm_fit(model, f$Y, f$left)$metrics
  pattern <- f$Y[, f$left]
  colnames(pattern) <- paste0("v", seq_len(ncol(pattern)))
  reference <- rMVPA:::.banded_ridge_nested_cv(
    pattern, f$X, model$outer_folds, model$candidates,
    feature_groups = rep("pattern", ncol(pattern)), inner_v = 3L,
    blocks = f$block, purge = 1L, nuisance = f$drift, seed = 1L
  )
  expect_equal(got[["all_r2"]], mean(reference$performance$r2), tolerance = 1e-6)
  expect_equal(got[["manner_cor"]], reference$performance$correlation[[1L]],
               tolerance = 1e-6)
})

test_that("feature_decode_model time-shift p-values separate signal from none", {
  f <- .fdm_fixture()
  model <- feature_decode_model(f$dataset, f$design, outer_crossval = 6L,
                                tune_crossval = 3L, purge = 2L,
                                nuisance = f$drift, partial = TRUE,
                                n_perm = 99L, min_shift = 10L)
  expect_length(output_schema(model), 48L)
  m <- .fdm_fit(model, f$Y, f$left)$metrics
  expect_identical(names(m), model$metric_names)
  expect_equal(m[["phonology_cor_p"]], 0.01)
  expect_equal(m[["manner_partial_cor_p"]], 0.01)
  expect_gt(m[["lexical_cor_p"]], 0.05)
  expect_identical(m, .fdm_fit(model, f$Y, f$left)$metrics)
})

test_that("feature_decode_model runs through run_regional and run_searchlight", {
  f <- .fdm_fixture()
  model <- feature_decode_model(f$dataset, f$design, outer_crossval = 6L,
                                tune_crossval = 2L, alphas = c(1, 100),
                                nuisance = f$drift)
  regional <- run_regional(model, f$region)
  expect_s3_class(regional, "regional_mvpa_result")
  table <- as.data.frame(regional$performance_table)
  expect_true(all(model$metric_names %in% names(table)))
  expect_equal(nrow(table), 2L)
  expect_gt(table$phonology_cor[[1L]], table$phonology_cor[[2L]])
  expect_gt(table$acoustic_cor[[2L]], 0.8)

  searchlight <- run_searchlight(model, radius = 3, method = "standard")
  expect_s3_class(searchlight, "searchlight_result")
  expect_identical(names(searchlight$results), model$metric_names)
})

test_that("feature_decode_model validates inputs and reports unusable regions", {
  f <- .fdm_fixture()
  build <- function(...) feature_decode_model(f$dataset, f$design,
                                              outer_crossval = 6L, ...)
  expect_error(feature_decode_model(f$dataset, list()), "feature_sets_design")
  expect_error(build(purge = -1), "purge")
  expect_error(build(n_perm = 1.5), "n_perm")
  expect_error(build(partial = NA), "partial")
  expect_error(build(nuisance = f$drift[-1L, , drop = FALSE]), "one row per observation")
  no_blocks <- feature_sets_design(f$design$X_train, time_series = TRUE)
  expect_error(feature_decode_model(f$dataset, no_blocks, outer_crossval = 6L),
               "random row-wise CV is unsafe")

  model <- build(tune_crossval = 2L, alphas = c(1, 100))
  flat <- f$Y
  flat[, f$left] <- 1
  failed <- .fdm_fit(model, flat, f$left)
  expect_true(failed$error)
  expect_match(failed$error_message, "too few usable voxels")
})
