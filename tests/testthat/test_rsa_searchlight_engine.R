# The rsa_fast engine calls train_model.rsa_model() on the same columns as the
# general path, so its maps must be identical for every RSA variant.
testthat::skip_if_not_installed("neuroim2")

rsa_engine_fixture <- function(seed = 1401) {
  set.seed(seed)
  ds <- gen_sample_dataset(c(5, 5, 5), 32, blocks = 4)
  D1 <- dist(matrix(rnorm(32 * 4), 32)); D2 <- dist(matrix(rnorm(32 * 4), 32))
  rdes <- rsa_design(~ D1 + D2, list(D1 = D1, D2 = D2, block = ds$design$block_var),
                     block_var = "block")
  list(ds = ds, rdes = rdes)
}

expect_rsa_identical <- function(ms, radius = 2) {
  fast <- suppressWarnings(run_searchlight(ms, radius = radius, backend = "default"))
  ref <- suppressWarnings(run_searchlight(ms, radius = radius, engine = "legacy", backend = "default"))
  expect_identical(attr(fast, "searchlight_engine"), "rsa_fast")
  expect_identical(sort(names(fast$results)), sort(names(ref$results)))
  for (m in names(ref$results)) {
    expect_identical(as.numeric(neuroim2::values(fast$results[[m]])),
                     as.numeric(neuroim2::values(ref$results[[m]])), info = m)
  }
}

test_that("rsa_fast is identical to the general path across distance and regression types", {
  fx <- rsa_engine_fixture()
  for (dm in c("pearson", "spearman")) for (rt in c("pearson", "spearman", "lm")) {
    expect_rsa_identical(rsa_model(fx$ds$dataset, fx$rdes, distmethod = dm, regtype = rt,
                                   check_collinearity = FALSE))
  }
  expect_rsa_identical(rsa_model(fx$ds$dataset, fx$rdes, regtype = "lm", semipartial = TRUE,
                                 check_collinearity = FALSE))
})

test_that("rsa_fast declines missing values and fingerprints", {
  fx <- rsa_engine_fixture(1402)
  arr <- as.array(fx$ds$dataset$train_data)
  arr[2, 2, 2, 3] <- NA
  fx$ds$dataset$train_data <- neuroim2::NeuroVec(arr, neuroim2::space(fx$ds$dataset$train_data))
  ms <- rsa_model(fx$ds$dataset, fx$rdes, check_collinearity = FALSE)
  res <- suppressWarnings(run_searchlight(ms, radius = 2, backend = "default"))
  expect_identical(attr(res, "searchlight_engine"), "legacy")

  fx2 <- rsa_engine_fixture(1403)
  fp <- rsa_model(fx2$ds$dataset, fx2$rdes, check_collinearity = FALSE, return_fingerprint = TRUE)
  expect_false(rMVPA:::.is_rsa_fast_path(fp, "standard"))
})

test_that("vectorised row ranks equal apply(rank) for multi-column matrices", {
  set.seed(1404)
  X <- matrix(round(rnorm(30 * 12), 1), 30)  # many ties
  expect_identical(rMVPA:::.rdm_rank_rows(X), t(apply(X, 1L, rank, ties.method = "average")))
})
