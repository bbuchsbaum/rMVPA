context("banded_ridge_model")

.brm_test_dataset <- function(Y, active_ids, dims = c(3L, 3L, 2L)) {
  stopifnot(ncol(Y) == prod(dims))
  n <- nrow(Y)
  data_array <- array(as.vector(t(Y)), c(dims, n))
  data_space <- neuroim2::NeuroSpace(c(dims, n), c(1, 1, 1))
  mask_array <- array(0, dims)
  mask_array[active_ids] <- 1
  mask_space <- neuroim2::NeuroSpace(dims, c(1, 1, 1))
  mvpa_dataset(
    neuroim2::NeuroVec(data_array, data_space),
    mask = neuroim2::NeuroVol(mask_array, mask_space)
  )
}

.brm_test_fixture <- function(seed = 8001L,
                              n = 18L,
                              active_ids = c(1L, 3L, 6L, 11L, 18L),
                              constant_response = FALSE) {
  set.seed(seed)
  p <- 6L
  dims <- c(3L, 3L, 2L)
  X <- matrix(rnorm(n * p), n, p,
              dimnames = list(NULL, paste0("x", seq_len(p))))
  beta <- matrix(c(
    1.1, -0.2, 0, 0.5, 0, 0,
    0, 0.8, -0.6, 0, 0.3, 0,
    -0.4, 0, 0.9, 0, 0, -0.2,
    0.2, 0.2, 0.2, -0.3, -0.3, -0.3,
    0, 0, 0, 0, 0, 0
  ), nrow = p)
  Y <- matrix(rnorm(n * prod(dims), sd = 0.4), n, prod(dims))
  Y[, active_ids] <- X %*% beta + matrix(rnorm(n * length(active_ids), sd = 0.2),
                                         n, length(active_ids))
  if (constant_response) Y[, active_ids[[1L]]] <- 3
  dataset <- .brm_test_dataset(Y, active_ids, dims)
  features <- feature_sets(X, blocks(low = 3L, semantic = 3L))
  block_labels <- cut(seq_len(n), breaks = 3L, labels = FALSE)
  blocks_train <- as.integer(block_labels)
  design <- feature_sets_design(
    features, block_var_train = blocks_train, time_series = TRUE
  )
  candidates <- rMVPA:::.banded_ridge_candidates(
    c("low", "semantic"), alphas = c(0.2, 2), method = "fixed",
    theta = rbind(c(low = 0.8, semantic = 0.2),
                  c(low = 0.25, semantic = 0.75))
  )
  list(X = X, Y = Y, active_ids = active_ids, dataset = dataset,
       design = design, blocks = blocks_train, candidates = candidates)
}

.brm_test_model <- function(fixture, ...) {
  args <- list(
    dataset = fixture$dataset,
    design = fixture$design,
    outer_crossval = 3L,
    tune_crossval = 2L,
    candidates = fixture$candidates,
    solver = "direct",
    target_batch_size = 2L,
    return_predictions = TRUE,
    seed = 8002L
  )
  dots <- list(...)
  args[names(dots)] <- dots
  do.call(banded_ridge_model, args)
}

.brm_active_response <- function(fixture) {
  out <- fixture$Y[, fixture$active_ids, drop = FALSE]
  colnames(out) <- paste0("voxel_", fixture$active_ids)
  out
}

test_that("public image workflow matches the certified nested matrix oracle", {
  fixture <- .brm_test_fixture()
  model <- .brm_test_model(fixture, target_batch_size = Inf)
  result <- run_banded_ridge(model)
  observed <- .brm_active_response(fixture)
  oracle <- rMVPA:::.banded_ridge_nested_cv(
    fixture$X, observed, model$outer_folds, model$candidates,
    feature_groups = as.character(fixture$design$X_train$set),
    inner_v = model$inner_v, blocks = fixture$blocks,
    response_batch_size = ncol(observed), seed = model$seed,
    solver = "direct"
  )

  expect_s3_class(model, "banded_ridge_model")
  expect_s3_class(result, "banded_ridge_result")
  expect_equal(result$predictions, oracle$predictions, tolerance = 1e-12)
  expect_equal(result$metrics[, c("response", "n_obs", "mse", "correlation", "r2")],
               oracle$performance, tolerance = 1e-12)
  oracle_selection <- oracle$selections[order(
    match(oracle$selections$response, colnames(observed)),
    oracle$selections$outer_fold
  ), ]
  expect_identical(result$hyperparameters$candidate_id,
                   oracle_selection$candidate_id)
  expect_equal(result$hyperparameters$inner_score,
               oracle_selection$inner_score, tolerance = 1e-12)
  expect_identical(result$active_response_ids, fixture$active_ids)
  expect_identical(result$provenance$outer_fold_ids,
                   vapply(model$outer_folds, `[[`, character(1), "id"))
})

test_that("global voxel identity and sparse-map holes survive every public output", {
  fixture <- .brm_test_fixture(active_ids = c(2L, 5L, 7L, 13L, 17L))
  result <- run_banded_ridge(.brm_test_model(fixture))
  score <- rMVPA:::.banded_ridge_score(.brm_active_response(fixture),
                                       result$predictions)
  inactive <- setdiff(seq_len(prod(c(3L, 3L, 2L))), fixture$active_ids)

  expect_identical(result$metrics$voxel_index, fixture$active_ids)
  expect_identical(result$metrics$response,
                   paste0("voxel_", fixture$active_ids))
  expect_identical(sort(unique(result$hyperparameters$voxel_index)),
                   fixture$active_ids)
  expect_equal(result$metrics$mse, score$mse, tolerance = 0)
  expect_equal(result$metrics$correlation, score$correlation, tolerance = 0)
  expect_equal(result$metrics$r2, score$r2, tolerance = 0)
  for (map in result$maps) {
    values <- as.numeric(map)
    expect_true(all(is.na(values[inactive])))
    expect_true(all(is.finite(values[fixture$active_ids])))
  }
})

test_that("response chunking, reordered partitions, and regional batching are invariant", {
  fixture <- .brm_test_fixture(seed = 8003L)
  model <- .brm_test_model(fixture, target_batch_size = 2L)
  one <- run_banded_ridge(model, target_batch_size = 1L)
  prime <- run_banded_ridge(model, target_batch_size = 3L)
  all <- run_banded_ridge(model, target_batch_size = Inf)
  ids <- fixture$active_ids
  reordered <- run_banded_ridge(
    model,
    response_partitions = list(ids[c(2L, 5L)], ids[c(1L, 4L)], ids[[3L]])
  )
  region_values <- integer(prod(c(3L, 3L, 2L)))
  region_values[ids] <- c(2L, 1L, 2L, 3L, 1L)
  region_mask <- neuroim2::NeuroVol(
    array(region_values, c(3L, 3L, 2L)), neuroim2::space(fixture$dataset$mask)
  )
  regional <- run_regional(model, region_mask)

  for (candidate in list(prime, all, reordered, regional)) {
    expect_equal(candidate$predictions, one$predictions, tolerance = 1e-12)
    expect_equal(candidate$metrics, one$metrics, tolerance = 1e-12)
    expect_identical(candidate$hyperparameters$candidate_id,
                     one$hyperparameters$candidate_id)
    expect_equal(candidate$hyperparameters$inner_score,
                 one$hyperparameters$inner_score, tolerance = 1e-12)
    for (map_name in names(one$maps)) {
      expect_equal(as.numeric(candidate$maps[[map_name]]),
                   as.numeric(one$maps[[map_name]]), tolerance = 1e-12)
    }
  }
  expect_null(one$diagnostics)
  expect_equal(one$provenance$peak_intermediate_dimensions[["responses"]], 1L)
  expect_equal(all$provenance$peak_intermediate_dimensions[["responses"]],
               length(ids))
  expect_gt(all$provenance$max_chunk_result_bytes,
            one$provenance$max_chunk_result_bytes)
  expect_identical(sort(unlist(lapply(
    reordered$provenance$chunk_manifest, `[[`, "voxel_ids"
  ))), ids)
})

test_that("constant responses and pooled OOF metrics have explicit semantics", {
  fixture <- .brm_test_fixture(seed = 8004L, n = 20L,
                               constant_response = TRUE)
  result <- run_banded_ridge(.brm_test_model(fixture, target_batch_size = Inf))
  observed <- .brm_active_response(fixture)
  oracle <- rMVPA:::.banded_ridge_score(observed, result$predictions)

  expect_equal(unname(result$predictions[, 1L]), rep(3, nrow(observed)),
               tolerance = 1e-12)
  expect_equal(result$metrics$mse, oracle$mse, tolerance = 0)
  expect_true(is.na(result$metrics$correlation[[1L]]))
  expect_true(is.na(result$metrics$r2[[1L]]))
  negative <- rMVPA:::.banded_ridge_score(
    matrix(c(-1, 0, 1), ncol = 1L, dimnames = list(NULL, "y")),
    matrix(c(1, 0, -1), ncol = 1L, dimnames = list(NULL, "y"))
  )
  expect_lt(negative$r2, 0)
})

test_that("outer-test response mutations cannot affect that fold's tuning or fit", {
  fixture <- .brm_test_fixture(seed = 8005L)
  model <- .brm_test_model(fixture, target_batch_size = Inf,
                           retain_diagnostics = TRUE)
  baseline <- run_banded_ridge(model)
  target_fold <- model$outer_folds[[2L]]
  mutated_y <- fixture$Y
  mutated_y[target_fold$test, fixture$active_ids] <-
    mutated_y[target_fold$test, fixture$active_ids] + 1e6
  mutated_fixture <- fixture
  mutated_fixture$Y <- mutated_y
  mutated_fixture$dataset <- .brm_test_dataset(
    mutated_y, fixture$active_ids
  )
  mutated_model <- .brm_test_model(
    mutated_fixture, outer_crossval = model$outer_folds,
    target_batch_size = Inf, retain_diagnostics = TRUE
  )
  mutated <- run_banded_ridge(mutated_model)
  before_diag <- baseline$diagnostics[[1L]][[2L]]
  after_diag <- mutated$diagnostics[[1L]][[2L]]

  expect_equal(before_diag$candidate_scores, after_diag$candidate_scores,
               tolerance = 0)
  expect_identical(before_diag$selected$candidate_id,
                   after_diag$selected$candidate_id)
  expect_equal(before_diag$outer_fits, after_diag$outer_fits, tolerance = 0)
  expect_equal(baseline$predictions[target_fold$test, ],
               mutated$predictions[target_fold$test, ], tolerance = 0)
  expect_false(isTRUE(all.equal(
    baseline$metrics$mse, mutated$metrics$mse, tolerance = 0
  )))
})

test_that("certified public solvers agree and emit cache receipts", {
  fixture <- .brm_test_fixture(seed = 8006L)
  direct <- run_banded_ridge(.brm_test_model(
    fixture, solver = "direct", target_batch_size = 2L
  ))
  svd <- run_banded_ridge(.brm_test_model(
    fixture, solver = "svd_primal", target_batch_size = 2L
  ))
  dual <- run_banded_ridge(.brm_test_model(
    fixture, solver = "dual_kernel", target_batch_size = 2L
  ))

  expect_equal(svd$predictions, direct$predictions, tolerance = 1e-10)
  expect_equal(dual$predictions, direct$predictions, tolerance = 1e-10)
  expect_identical(svd$hyperparameters$candidate_id,
                   direct$hyperparameters$candidate_id)
  expect_identical(dual$hyperparameters$candidate_id,
                   direct$hyperparameters$candidate_id)
  expect_gt(svd$provenance$solver_decomposition_count, 0L)
  expect_gt(svd$provenance$solver_cache_hits, 0L)
  expect_gt(dual$provenance$solver_decomposition_count, 0L)
  expect_gt(dual$provenance$solver_cache_hits, 0L)
})

test_that("weight policies retain only their advertised representations", {
  fixture <- .brm_test_fixture(seed = 8007L)
  none <- run_banded_ridge(.brm_test_model(
    fixture, return_predictions = FALSE, weight_retention = "none"
  ))
  primal <- run_banded_ridge(.brm_test_model(
    fixture, return_predictions = FALSE, weight_retention = "primal"
  ))
  dual <- run_banded_ridge(.brm_test_model(
    fixture, return_predictions = FALSE, weight_retention = "dual"
  ))

  expect_null(none$predictions)
  expect_null(none$weights)
  expect_identical(names(primal$weights), paste0("voxel_", fixture$active_ids))
  expect_true(all(vapply(primal$weights, length, integer(1)) == 3L))
  expect_true(all(vapply(primal$weights, function(response) {
    all(vapply(response, function(fold) length(fold$coefficients) == ncol(fixture$X),
               logical(1)))
  }, logical(1))))
  expect_identical(names(dual$weights), paste0("voxel_", fixture$active_ids))
  expect_true(all(vapply(dual$weights, function(response) {
    all(vapply(response, function(fold) {
      length(fold$dual_weights) == length(fold$train)
    }, logical(1)))
  }, logical(1))))
  expect_gt(dual$provenance$weight_decomposition_count, 0L)
})

test_that("constructor and spatial execution fail early on unsafe contracts", {
  fixture <- .brm_test_fixture(seed = 8008L)
  no_blocks <- feature_sets_design(fixture$design$X_train, time_series = TRUE)
  expect_error(banded_ridge_model(
    fixture$dataset, no_blocks, outer_crossval = 3L, tune_crossval = 2L
  ), "time_series.*requires")
  ordinary <- feature_sets_design(fixture$design$X_train)
  expect_error(banded_ridge_model(
    fixture$dataset, ordinary, tune_crossval = 2L
  ), "supply outer_crossval")
  expect_error(.brm_test_model(fixture, purge = -1), "non-negative")
  expect_error(.brm_test_model(fixture, seed = NA_real_), "finite integer")
  expect_error(.brm_test_model(fixture, target_batch_size = 0), "positive integer")

  wrong_rows <- feature_sets_design(feature_sets(
    fixture$X[-1L, , drop = FALSE], blocks(low = 3L, semantic = 3L)
  ))
  expect_error(banded_ridge_model(
    fixture$dataset, wrong_rows, outer_crossval = 3L, tune_crossval = 2L
  ), "rows must match")
  wrong_candidates <- fixture$candidates
  colnames(wrong_candidates$theta) <- c("low", "wrong")
  expect_error(.brm_test_model(fixture, candidates = wrong_candidates),
               "band names exactly")

  model <- .brm_test_model(fixture)
  ids <- fixture$active_ids
  expect_error(run_banded_ridge(model, response_partitions = list(ids[-1L])),
               "every active voxel exactly once")
  expect_error(run_banded_ridge(model, response_partitions = list(ids, ids[[1L]])),
               "exactly once")
  expect_error(run_searchlight(model, radius = 2), "does not support overlapping")
  expect_error(.brm_test_model(fixture, tune_crossval = 3L),
               "fewer independent blocks")

  empty <- fixture$dataset
  empty$mask[] <- 0
  expect_error(banded_ridge_model(
    empty, fixture$design, outer_crossval = 3L, tune_crossval = 2L
  ), "no active responses")
})

test_that("retained-output limits refuse or explicitly downgrade allocations", {
  fixture <- .brm_test_fixture(seed = 8009L)
  expect_error(.brm_test_model(
    fixture, return_predictions = TRUE, max_retained_mb = 1e-6
  ), "above max_retained_mb")
  fallback <- .brm_test_model(
    fixture, return_predictions = FALSE, weight_retention = "primal",
    max_retained_mb = 5e-4, weight_overflow = "none"
  )
  expect_identical(fallback$weight_retention, "none")
  expect_match(fallback$retention_notice, "weight retention was disabled")
  expect_equal(fallback$retention_estimated_bytes[["total"]], 0)

  guarded <- .brm_test_model(
    fixture, target_batch_size = 1L, memory_limit_mb = 4e-4
  )
  expect_error(run_banded_ridge(guarded, target_batch_size = Inf),
               "exceeds the configured estimated storage contract")
})

test_that("feature-set design owns and validates training fold metadata", {
  fixture <- .brm_test_fixture(seed = 8010L)
  expect_identical(fixture$design$block_var_train, fixture$blocks)
  expect_true(fixture$design$time_series)
  expect_identical(fixture$design$train_design$.block, fixture$blocks)
  expect_error(feature_sets_design(
    fixture$design$X_train, block_var_train = 1:2
  ), "one non-missing value per training row")
  expect_error(feature_sets_design(
    fixture$design$X_train, block_var_train = rep(NA, nrow(fixture$X))
  ), "one non-missing value per training row")
  expect_error(feature_sets_design(
    fixture$design$X_train, time_series = NA
  ), "TRUE or FALSE")
})

test_that("print method reports the estimand-facing summary", {
  fixture <- .brm_test_fixture(seed = 8011L)
  result <- run_banded_ridge(.brm_test_model(
    fixture, return_predictions = FALSE, target_batch_size = Inf
  ))
  printed <- capture.output(print(result))
  expect_match(paste(printed, collapse = "\n"), "banded_ridge_result")
  expect_match(paste(printed, collapse = "\n"), "Responses: 5")
  expect_match(paste(printed, collapse = "\n"), "Solver:    direct")
  expect_match(paste(printed, collapse = "\n"), "OOF MSE")
})

test_that("outer_crossval accepts cross_validation objects and applies purge", {
  fixture <- .brm_test_fixture(seed = 8021L)
  n <- nrow(fixture$X)
  purge_rows <- function(train, test, gap) {
    train[vapply(train, function(i) all(abs(i - test) > gap), logical(1))]
  }
  manual <- lapply(seq_len(3L), function(k) {
    test <- which(fixture$blocks == k)
    train <- setdiff(seq_len(n), test)
    list(id = paste0("fold_", k), train = purge_rows(train, test, 1L),
         test = test)
  })

  blocked <- .brm_test_model(
    fixture, outer_crossval = blocked_cross_validation(fixture$blocks),
    purge = 1L
  )
  expect_equal(blocked$outer_folds, manual)
  explicit <- .brm_test_model(fixture, outer_crossval = manual, purge = 1L)
  expect_equal(explicit$outer_folds, manual)
  expect_equal(run_banded_ridge(blocked)$metrics,
               run_banded_ridge(explicit)$metrics, tolerance = 1e-12)

  unpurged <- lapply(manual, function(fold) {
    list(id = fold$id, train = setdiff(seq_len(n), fold$test), test = fold$test)
  })
  expect_error(.brm_test_model(fixture, outer_crossval = unpurged, purge = 1L),
               "validated as supplied")

  custom <- custom_cross_validation(lapply(seq_len(3L), function(k) {
    list(train = which(fixture$blocks != k), test = which(fixture$blocks == k))
  }))
  from_custom <- .brm_test_model(fixture, outer_crossval = custom)
  expect_equal(from_custom$outer_folds, unpurged)

  ordinary <- feature_sets_design(fixture$design$X_train)
  set.seed(8024L)
  kfold <- kfold_cross_validation(n, nfolds = 3L)
  from_kfold <- banded_ridge_model(
    fixture$dataset, ordinary, outer_crossval = kfold, tune_crossval = 2L,
    candidates = fixture$candidates, solver = "direct"
  )
  for (k in seq_len(3L)) {
    expect_identical(from_kfold$outer_folds[[k]]$test,
                     which(kfold$block_var == k))
  }

  expect_error(.brm_test_model(
    fixture,
    outer_crossval = bootstrap_blocked_cross_validation(fixture$blocks, nreps = 2L)
  ), "cannot be converted")
  expect_error(.brm_test_model(
    fixture,
    outer_crossval = twofold_blocked_cross_validation(fixture$blocks, nreps = 2L)
  ), "cannot be converted")
  expect_error(.brm_test_model(
    fixture, outer_crossval = blocked_cross_validation(c(fixture$blocks, 1L))
  ), "block_var of length")
  expect_error(.brm_test_model(
    fixture, tune_crossval = blocked_cross_validation(fixture$blocks)
  ), "accepted only for outer_crossval")
  expect_error(.brm_test_model(fixture, outer_crossval = "three"),
               "cross_validation object")
  forward <- custom_cross_validation(list(
    list(train = 1:6, test = 7:12), list(train = 1:12, test = 13:18)
  ))
  expect_error(.brm_test_model(fixture, outer_crossval = forward),
               "cannot be converted.*exactly once")
  expect_error(.brm_test_model(
    fixture, outer_crossval = data.frame(train = 1L, test = 2L)
  ), "data.frame")
})

test_that("tune_crossval defaults from available training units and validates early", {
  fixture <- .brm_test_fixture(seed = 8022L)
  model <- banded_ridge_model(
    fixture$dataset, fixture$design, outer_crossval = 3L,
    candidates = fixture$candidates, solver = "direct"
  )
  expect_identical(model$inner_v, 2L)
  expect_identical(model$inner_tuning, "nested")
  expect_null(model$inner_folds)
  expect_s3_class(run_banded_ridge(model), "banded_ridge_result")

  ordinary <- feature_sets_design(fixture$design$X_train)
  plain <- banded_ridge_model(
    fixture$dataset, ordinary, outer_crossval = 3L,
    candidates = fixture$candidates, solver = "direct"
  )
  expect_identical(plain$inner_v, 5L)

  two_blocks <- feature_sets_design(
    fixture$design$X_train, block_var_train = rep(1:2, each = 9L),
    time_series = TRUE
  )
  expect_error(banded_ridge_model(
    fixture$dataset, two_blocks, candidates = fixture$candidates
  ), "no usable default")

  bad_inner <- rep(list(list(
    list(id = "a", train = 1:2, test = 3:4),
    list(id = "b", train = 3:4, test = 1:2)
  )), 3L)
  expect_error(.brm_test_model(fixture, tune_crossval = bad_inner),
               "cannot be applied inside outer fold")
  expect_error(.brm_test_model(fixture, tune_crossval = bad_inner[1:2]),
               "one inner-fold list per outer fold")

  # Contiguous outer folds survive the purge; the default random row-wise
  # inner splits cannot, and the error says why.
  no_blocks <- feature_sets_design(fixture$design$X_train, time_series = TRUE)
  expect_error(banded_ridge_model(
    fixture$dataset, no_blocks,
    outer_crossval = blocked_cross_validation(fixture$blocks), purge = 5L,
    candidates = fixture$candidates
  ), "random row-wise splits")
})

test_that("inner tuning is skipped when only one candidate survives the scopes", {
  fixture <- .brm_test_fixture(seed = 8023L)
  single <- banded_ridge_model(
    fixture$dataset, fixture$design, outer_crossval = 3L,
    alphas = 1, theta_method = "fixed", theta = c(low = 0.5, semantic = 0.5),
    solver = "direct", target_batch_size = 2L, return_predictions = TRUE,
    seed = 8002L
  )
  expect_true(is.na(single$inner_v))
  expect_identical(single$inner_tuning, "none")
  expect_null(single$inner_folds)
  result <- run_banded_ridge(single)
  expect_true(all(is.na(result$hyperparameters$inner_score)))
  expect_true(all(result$hyperparameters$alpha == 1))

  observed <- .brm_active_response(fixture)
  oracle <- rMVPA:::.banded_ridge_nested_cv(
    fixture$X, observed, single$outer_folds, single$candidates,
    feature_groups = as.character(fixture$design$X_train$set),
    inner_v = 2L, blocks = fixture$blocks,
    response_batch_size = ncol(observed), seed = single$seed,
    solver = "direct"
  )
  expect_equal(result$predictions, oracle$predictions, tolerance = 1e-12)

  fixed <- .brm_test_model(
    fixture, tune_crossval = 2L, alpha_scope = "fixed", fixed_alpha = 0.2,
    theta_scope = "fixed", fixed_theta = c(low = 0.8, semantic = 0.2),
    delta_sets = "all"
  )
  expect_true(is.na(fixed$inner_v))
  expect_identical(fixed$inner_tuning, "none")
  fixed_result <- run_banded_ridge(fixed)
  expect_true(all(fixed_result$hyperparameters$alpha == 0.2))
  expect_length(unique(fixed_result$hyperparameters$candidate_id), 1L)
  expect_true(all(is.na(fixed_result$hyperparameters$inner_score)))
  expect_gt(nrow(fixed_result$predictive_leave_one_band_out$effects), 0L)

  expect_error(.brm_test_model(fixture, alpha_scope = "fixed", fixed_alpha = 99),
               "No candidate matches")
  # A malformed tune_crossval is rejected even though it would be ignored.
  expect_error(banded_ridge_model(
    fixture$dataset, fixture$design, outer_crossval = 3L,
    tune_crossval = "garbage", alphas = 1, theta_method = "fixed",
    theta = c(low = 0.5, semantic = 0.5)
  ), "integer of at least two")
  expect_error(.brm_test_model(fixture, alpha_scope = "elsewhere"),
               "should be one of")
})
