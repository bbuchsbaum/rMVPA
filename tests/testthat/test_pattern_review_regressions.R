test_that("fixed-rank tuning compares penalties at the rank that will be fitted", {
  testthat::local_mocked_bindings(
    .pattern_fit = function(X, targets, rank, control, graph, penalty) {
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
    .pattern_fit = function(X, targets, rank, control, graph, penalty) {
      cap <- if (1 %in% X[, 1]) 2L else 1L
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
