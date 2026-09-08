test_that("fixed-rank tuning compares penalties at the rank that will be fitted", {
  testthat::local_mocked_bindings(
    .pattern_fit = function(X, targets, rank, control, graph, penalty, cap_rank = FALSE) {
      path <- list(list(rank = 1, sparse = penalty$sparse), list(rank = 2, sparse = penalty$sparse))
      if (identical(rank, "path")) path else path[[rank]]
    },
    .pattern_loss = function(fit, X, targets) {
      if (fit$sparse == 0.1) c(1, 0.1)[fit$rank] else c(0.5, 0.4)[fit$rank]
    }
  )
  X <- matrix(seq_len(80), 40, 2); Y <- matrix(seq_len(80), 40, 2)
  fixed <- .pattern_select_config(X, Y, rep(1:4, each = 10), pattern_control(),
                                  penalty = list(sparse = c(.1, .8)), rank = 1)
  auto <- .pattern_select_config(X, Y, rep(1:4, each = 10), pattern_control(),
                                 penalty = list(sparse = c(.1, .8)))
  expect_equal(fixed$rank, 1)
  expect_equal(fixed$penalty$sparse, .8)
  expect_equal(auto$rank, 2)
  expect_equal(auto$penalty$sparse, .1)
})

test_that("clustered ROI graph restriction follows retained column positions", {
  set.seed(771)
  ds <- gen_clustered_sample_dataset(D = c(6, 6, 4), nobs = 40, K = 6, nlevels = 2)
  model <- pattern_model(ds$dataset, ds$design, rank = 1,
                         penalty = list(signed_smooth = 1))
  roi <- as_roi(clustered_data_sample(2, c(2L, 4L, 6L)), ds$dataset)
  expect_equal(roi$feature_positions, c(2L, 4L, 6L))
  filtered <- filter_roi(roi)
  expect_equal(filtered$feature_positions, c(2L, 4L, 6L))
  roi_data <- list(train_data = as.matrix(neuroim2::values(filtered$train_roi)),
                   indices = neuroim2::indices(filtered$train_roi),
                   feature_positions = filtered$feature_positions)
  expect_equal(.pattern_roi_graph(model, roi_data)$feature_ids, model$graph$feature_ids[c(2, 4, 6)])
  out <- fit_roi(model, roi_data, list(id = 2))
  expect_false(out$error)
  expect_true(is.finite(out$metrics["Accuracy"]))
  # A constant first cluster is removed without confusing spatial centroids.
  values <- as.matrix(neuroim2::values(roi$train_roi)); values[, 1] <- 0
  roi$train_roi <- neuroim2::ROIVec(neuroim2::space(roi$train_roi),
    neuroim2::coords(roi$train_roi), data = values)
  filtered <- filter_roi(roi)
  expect_equal(filtered$feature_positions, c(4L, 6L))
})

test_that("optional Haufe diagnostics do not discard an evaluation with missing assessment data", {
  sim <- sim_pattern_data(n = 48, dims = c(3, 3, 2), K = 2, external_test = TRUE, n_test = 24, seed = 772)
  # The evaluation already handles missing predictions. Its optional diagnostic
  # must report why it is unavailable instead of throwing after the fit.
  testdata <- sim$dataset$test_data
  arr <- array(t(sim$X_test), c(3, 3, 2, 24)); arr[1, 1, 1, 1] <- NA_real_
  sim$dataset$test_data <- neuroim2::NeuroVec(arr, neuroim2::space(testdata))
  model <- pattern_model(sim$dataset, sim$design, rank = 1)
  utils::capture.output(out <- suppressMessages(run_global(model, return_fits = TRUE)))
  expect_equal(out$haufe_diagnostics[[1]]$status, "non-finite held-out retained features")
})


test_that("fixed-rank tuning respects each fold's eligible rank", {
  testthat::local_mocked_bindings(
    .pattern_fit = function(X, targets, rank, control, graph, penalty, cap_rank = FALSE) {
      cap <- if (1 %in% X[, 1]) 2L else 1L
      if (!identical(rank, "path") && rank > cap && !cap_rank) stop("exceeds eligible rank")
      path <- lapply(seq_len(cap), function(k) list(rank = k, sparse = penalty$sparse))
      if (identical(rank, "path")) path else path[[min(rank, cap)]]
    },
    .pattern_loss = function(fit, X, targets) {
      if (fit$sparse == .1) c(.2, 20)[fit$rank] else c(2, .5)[fit$rank]
    }
  )
  X <- cbind(1:40, 41:80)
  selected <- .pattern_select_config(X, X, rep(1:4, each = 10), pattern_control(),
    penalty = list(sparse = c(.1, .8)), rank = 2)
  # One inner fold fits eligible rank one; three fit requested rank two.
  # Taking the minimum common path rank would incorrectly choose sparse=.1.
  expect_equal(selected$penalty$sparse, .8)
  expect_equal(selected$loss, (2 + 3*.5)/4)
  expect_equal(selected$rank, 2)
})


test_that("coincident ROI coordinates do not alias retained test columns", {
  set.seed(774)
  sp <- neuroim2::NeuroSpace(c(3, 3, 3), c(1, 1, 1))
  coords <- rbind(c(1, 1, 1), c(2, 2, 2), c(1, 1, 1))
  train <- cbind(0, rnorm(20), rnorm(20))
  test <- matrix(rnorm(60), 20, 3)
  roi <- list(train_roi = neuroim2::ROIVec(sp, coords, train),
    test_roi = neuroim2::ROIVec(sp, coords, test), feature_positions = c(8L, 9L, 10L))
  filtered <- filter_roi(roi)
  expect_equal(as.matrix(neuroim2::values(filtered$test_roi)), test[, 2:3], ignore_attr = TRUE)
  expect_equal(filtered$feature_positions, c(9L, 10L))
})


test_that("real fixed-rank tuning caps feature-limited inner fits", {
  set.seed(775)
  Y <- cbind(signal = rnorm(80), other = rnorm(80))
  X <- matrix(2 * Y[, 1] + rnorm(80, sd = .4), ncol = 1)
  blocks <- rep(1:4, each = 20)
  control <- pattern_control(max_rank = 2, noise = list(type = "diag"))
  candidates <- list(sparse = c(.1, .8))
  grid <- .pattern_penalty_grid(candidates)
  folds <- .pattern_inner_folds(nrow(X), blocks)
  # Direct capped fits provide the oracle. Every inner fit has eligible rank
  # one although the caller requested two; dropping them would choose the
  # fallback (first/strongest) penalty instead of the lowest held-out loss.
  losses <- vapply(grid, function(penalty) mean(vapply(folds, function(f) {
    fit <- .pattern_fit(X[f$train, , drop = FALSE], Y[f$train, , drop = FALSE],
      rank = 2, control = control, penalty = penalty, cap_rank = TRUE)
    .pattern_loss(fit, X[f$test, , drop = FALSE], Y[f$test, , drop = FALSE])
  }, numeric(1))), numeric(1))
  selected <- .pattern_select_config(X, Y, blocks, control, penalty = candidates, rank = 2)
  expect_equal(which.min(losses), 2L)
  expect_equal(selected$penalty, grid[[which.min(losses)]])
  expect_equal(selected$loss, min(losses), tolerance = 1e-10)
})
