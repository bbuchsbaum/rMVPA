# Regression tests for two fixes:
#  1. Per-ROI results must not depend on the future plan. The sequential path
#     used the global RNG stream while furrr (seed = TRUE) used per-element
#     L'Ecuyer streams, so stochastic processors gave different numbers under
#     sequential and parallel plans.
#  2. pca_lda fitted LDA on unit-norm scores but predicted on U D scores, and
#     mgsda returned NULL probabilities, so every fold failed.

testthat::skip_if_not_installed("neuroim2")

# ---- helpers ---------------------------------------------------------------

make_rng_spec <- function(shard = FALSE) {
  ds <- gen_sample_dataset(c(4, 4, 4), 18, blocks = 3, nlevels = 2)
  cv <- blocked_cross_validation(ds$design$block_var)
  mspec <- mvpa_model(
    model = load_model("sda_notune"),
    dataset = ds$dataset,
    design = ds$design,
    model_type = "classification",
    crossval = cv
  )
  if (shard) mspec <- use_shard(mspec)
  list(mspec = mspec, mask_idx = which(ds$dataset$mask > 0))
}

# Processor whose output depends on the global RNG stream.
make_draw_processor <- function() {
  function(obj, roi, rnum, center_global_id = NA) {
    draws <- c(runif(2), rnorm(1))
    tibble::tibble(
      result = list(NULL),
      indices = list(NULL),
      performance = list(draws),
      id = rnum,
      error = FALSE,
      error_message = "~",
      warning = FALSE,
      warning_message = "~"
    )
  }
}

# Runs mvpa_iterate from a fixed global seed under the current future plan.
run_draws <- function(shard = FALSE, min_chunk = 16L) {
  spec <- make_rng_spec(shard = shard)
  mask_idx <- spec$mask_idx
  vox_list <- lapply(1:12, function(i) as.integer(mask_idx[i:(i + 2)]))
  ids <- 101:112

  old_opt <- options(rMVPA.searchlight_min_chunk = min_chunk,
                     rMVPA.searchlight_backend_default = "default")
  on.exit(options(old_opt), add = TRUE)

  set.seed(42)
  res <- mvpa_iterate(
    mod_spec = spec$mspec,
    vox_list = vox_list,
    ids = ids,
    verbose = FALSE,
    analysis_type = "searchlight",
    processor = make_draw_processor(),
    fail_fast = TRUE
  )
  list(
    id = res$id,
    draws = res$performance,
    state = get(".Random.seed", envir = globalenv(), inherits = FALSE)
  )
}

expect_same_draws <- function(actual, reference) {
  expect_equal(actual$id, reference$id)
  expect_equal(actual$draws, reference$draws)
  expect_identical(actual$state, reference$state)
}

# Synthetic two-class data with signal in the first 8 voxels.
make_signal_dataset <- function(nobs = 48, D = c(4, 4, 4), blocks = 4, shift = 2) {
  y <- factor(rep(c("a", "b"), length.out = nobs))
  nvox <- prod(D)
  M <- matrix(stats::rnorm(nvox * nobs), nvox, nobs)
  M[1:8, y == "b"] <- M[1:8, y == "b"] + shift
  bvec <- neuroim2::NeuroVec(array(M, c(D, nobs)),
                             neuroim2::NeuroSpace(c(D, nobs), c(1, 1, 1)))
  mask <- as.logical(neuroim2::NeuroVol(array(1, D),
                                        neuroim2::NeuroSpace(D, c(1, 1, 1))))
  dset <- mvpa_dataset(train_data = bvec, mask = mask)
  block_var <- as.integer(as.character(
    cut(seq_len(nobs), blocks, labels = seq_len(blocks))))
  des <- mvpa_design(data.frame(Y = y, block_var = block_var),
                     block_var = "block_var", y_train = ~ Y)
  list(dataset = dset, design = des, y = y)
}

# ---- 1. RNG parity across future plans ---------------------------------------

test_that("per-ROI draws are identical under sequential and multisession plans", {
  skip_on_cran()
  skip_if_not_installed("furrr")
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  future::plan(future::sequential)
  ref <- run_draws()
  expect_gt(length(unique(unlist(ref$draws))), 1L)

  future::plan(future::multisession, workers = 2)
  for (min_chunk in c(1L, 4L, 16L)) {
    expect_same_draws(run_draws(min_chunk = min_chunk), ref)
  }
})

test_that("per-ROI draws are identical under sequential and multicore plans", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("furrr")
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  future::plan(future::sequential)
  ref <- run_draws()

  future::plan(future::multicore, workers = 2)
  for (min_chunk in c(1L, 4L)) {
    expect_same_draws(run_draws(min_chunk = min_chunk), ref)
  }
})

test_that("shard backend draws are identical under sequential and multisession plans", {
  skip_on_cran()
  skip_if_not_installed("shard")
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  future::plan(future::sequential)
  ref <- run_draws(shard = TRUE)

  future::plan(future::multisession, workers = 2)
  expect_same_draws(run_draws(shard = TRUE, min_chunk = 1L), ref)
})

# ---- 2. pca_lda ----------------------------------------------------------------

test_that("pca_lda posteriors match an independent oracle", {
  set.seed(1)
  n <- 40
  p <- 12
  ncomp <- 3
  y <- factor(rep(c("a", "b"), each = n / 2))
  x <- matrix(stats::rnorm(n * p), n, p)
  x[y == "b", 1:4] <- x[y == "b", 1:4] + 1.5
  xnew <- matrix(stats::rnorm(20 * p), 20, p)
  xnew[1:10, 1:4] <- xnew[1:10, 1:4] + 1.5

  mod <- MVPAModels$pca_lda
  fit <- mod$fit(x, y, wts = NULL, param = list(ncomp = ncomp), lev = levels(y),
                 last = NULL, weights = NULL, classProbs = NULL)

  # Oracle: centre/scale the training data, SVD, LDA on scores U D, then
  # project new data with the training centre, scale and loadings V.
  sc <- scale(x)
  ctr <- attr(sc, "scaled:center")
  scl <- attr(sc, "scaled:scale")
  s <- svd(sc)
  U <- s$u[, seq_len(ncomp), drop = FALSE]
  D <- s$d[seq_len(ncomp)]
  V <- s$v[, seq_len(ncomp), drop = FALSE]
  oracle <- MASS::lda(U %*% diag(D, nrow = ncomp), y)
  oracle_post <- function(z) {
    predict(oracle, scale(z, ctr, scl) %*% V)$posterior
  }

  expect_equal(mod$prob(fit, xnew), oracle_post(xnew), tolerance = 1e-8)
  expect_equal(mod$prob(fit, x), oracle_post(x), tolerance = 1e-8)

  oracle_class <- predict(oracle, scale(xnew, ctr, scl) %*% V)$class
  expect_equal(as.character(mod$predict(fit, xnew)), as.character(oracle_class))

  # Training rows get soft posteriors, rows sum to one.
  post_train <- mod$prob(fit, x)
  expect_true(all(abs(rowSums(post_train) - 1) < 1e-8))
  expect_true(all(post_train > 0 & post_train < 1))
})

test_that("pca_lda through mvpa_model and crossval runs and beats chance", {
  set.seed(3)
  ds <- make_signal_dataset()
  cv <- blocked_cross_validation(ds$design$block_var)
  mspec <- mvpa_model(
    model = load_model("pca_lda"),
    dataset = ds$dataset,
    design = ds$design,
    model_type = "classification",
    crossval = cv
  )
  res <- suppressWarnings(run_searchlight(mspec, radius = 2, method = "standard"))
  expect_true(inherits(res, "searchlight_result"))
  acc <- as.numeric(as.array(res$results$Accuracy))
  acc <- acc[is.finite(acc)]
  expect_gt(length(acc), 0L)
  expect_gt(mean(acc), 0.55)
})

# ---- 3. mgsda ------------------------------------------------------------------

test_that("mgsda one-hot probabilities have the expected shape and argmax", {
  levs <- c("a", "b", "c")
  preds <- c(2L, 1L, 3L, 2L)
  probs <- .mgsda_onehot_probs(preds, levs)

  expect_equal(dim(probs), c(4L, 3L))
  expect_equal(colnames(probs), levs)
  expect_true(all(rowSums(probs) == 1))
  expect_equal(levs[max.col(probs, ties.method = "first")], levs[preds])
})

test_that("mgsda prob uses classifier output and matches predict()", {
  skip_if_not_installed("MGSDA")
  set.seed(2)
  n <- 30
  y <- factor(rep(c("a", "b", "c"), each = n / 3))
  x <- matrix(stats::rnorm(n * 10), n, 10)
  x[y == "b", 1:3] <- x[y == "b", 1:3] + 2
  x[y == "c", 4:6] <- x[y == "c", 4:6] + 2
  mod <- MVPAModels$mgsda
  fit <- mod$fit(x, y, wts = NULL, param = list(lambda = 0.1), lev = levels(y),
                 last = NULL, weights = NULL, classProbs = NULL)
  probs <- mod$prob(fit, x)
  expect_equal(dim(probs), c(n, 3L))
  expect_equal(colnames(probs), levels(y))
  expect_true(all(abs(rowSums(probs) - 1) < 1e-12))
  expect_equal(levels(y)[max.col(probs, ties.method = "first")],
               as.character(mod$predict(fit, x)))
})
