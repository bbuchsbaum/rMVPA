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
  # the corrupted fold makes this fit lose to the training mean by design
  mutated <- suppressWarnings(run_banded_ridge(mutated_model))
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

test_that("dual solver reuses band Grams and a capped cache reproduces results", {
  fixture <- .brm_test_fixture(seed = 8030L)
  dual <- run_banded_ridge(.brm_test_model(
    fixture, solver = "dual_kernel", target_batch_size = 2L
  ))
  expect_gt(dual$provenance$solver_band_kernel_builds, 0L)
  expect_gt(dual$provenance$solver_band_kernel_hits, 0L)
  expect_identical(dual$provenance$solver_cache_evictions, 0L)
  expect_identical(dual$provenance$solver_cache_oversize, 0L)
  expect_identical(dual$provenance$solver_cache_limit_mb, 1024)
  expect_gt(dual$provenance$solver_cache_peak_mb, 0)

  # a cap above the solver plan's estimate but below the retained working set
  plan_mib <- dual$provenance$solver_plan$estimated_mib[["dual_kernel"]]
  cap <- max(plan_mib * 1.1, dual$provenance$solver_cache_peak_mb / 4)
  expect_lt(cap, dual$provenance$solver_cache_peak_mb)
  capped <- run_banded_ridge(.brm_test_model(
    fixture, solver = "dual_kernel", target_batch_size = 2L,
    memory_limit_mb = cap
  ))
  expect_gt(capped$provenance$solver_cache_evictions, 0L)
  expect_identical(capped$provenance$solver_cache_limit_mb, cap)
  expect_lte(capped$provenance$solver_cache_peak_mb, cap)
  expect_identical(capped$predictions, dual$predictions)
  expect_identical(capped$metrics, dual$metrics)
  expect_identical(capped$hyperparameters, dual$hyperparameters)

  # the plan estimate already exceeds one n x n entry, so any cap the plan
  # accepts fits single entries: caching degrades to eviction, never to oversize
  tiny <- run_banded_ridge(.brm_test_model(
    fixture, solver = "dual_kernel", target_batch_size = 2L,
    memory_limit_mb = plan_mib * 1.01
  ))
  expect_identical(tiny$provenance$solver_cache_oversize, 0L)
  expect_gt(tiny$provenance$solver_cache_evictions, 0L)
  expect_identical(tiny$predictions, dual$predictions)
  expect_identical(sum(tiny$provenance$work_manifest$cache_oversize_added), 0L)
})

# Issue #87: a fixed alpha grid cannot be right for every design, and a grid
# that stops below the optimum used to pass silently.

.brm_saturation_fixture <- function(seed = 8701L, n = 300L, D = 40L,
                                    dims = c(4L, 4L, 2L), signal = 0.02) {
  set.seed(seed)
  p <- D * 3L
  nv <- prod(dims)
  X <- matrix(rnorm(n * p), n, p,
              dimnames = list(NULL, paste0("x", seq_len(p))))
  Y <- X %*% matrix(rnorm(p * nv, sd = signal), p, nv) +
    matrix(rnorm(n * nv), n, nv)
  data_array <- array(as.vector(t(Y)), c(dims, n))
  dataset <- mvpa_dataset(
    neuroim2::NeuroVec(data_array,
                       neuroim2::NeuroSpace(c(dims, n), c(1, 1, 1))),
    mask = neuroim2::NeuroVol(array(1, dims),
                              neuroim2::NeuroSpace(dims, c(1, 1, 1)))
  )
  features <- feature_sets(X, blocks(a = D, b = D, c = D))
  design <- feature_sets_design(
    features, block_var_train = as.integer(cut(seq_len(n), 4L, labels = FALSE)),
    time_series = TRUE
  )
  list(X = X, dataset = dataset, design = design,
       groups = as.character(features$set))
}

.brm_saturation_model <- function(fixture, alpha_scope = "shared", ...) {
  banded_ridge_model(
    fixture$dataset, fixture$design, outer_crossval = 4L, tune_crossval = 2L,
    theta_method = "fixed", theta = c(a = 1 / 3, b = 1 / 3, c = 1 / 3),
    alpha_scope = alpha_scope, theta_scope = "shared", ...
  )
}

test_that("the default alpha grid is anchored on the design, not on a constant", {
  fixture <- .brm_saturation_fixture()
  n <- nrow(fixture$X)
  p <- ncol(fixture$X)
  expected_anchor <- ((n - 1) * (p / 3L)) / min(n - 1, p)
  expect_equal(rMVPA:::.brm_alpha_anchor(fixture$X, fixture$groups),
               expected_anchor)

  model <- .brm_saturation_model(fixture)
  expect_identical(model$alpha_grid_origin, "auto")
  expect_equal(model$alpha_anchor, expected_anchor)
  expect_equal(model$alpha_grid,
               expected_anchor * 10^seq(-3, 3, length.out = 9L))
  expect_identical(model$selectable_alphas, sort(model$alpha_grid))

  # with more rows than columns the rank bound is the column count, so the
  # anchor is (n - 1) / bands: doubling the rows doubles it, and no constant
  # grid can track that
  expect_equal(expected_anchor, (n - 1) / 3)
  wide_n <- .brm_saturation_fixture(seed = 8702L, n = 600L)
  expect_equal(rMVPA:::.brm_alpha_anchor(wide_n$X, wide_n$groups),
               (600 - 1) / 3)

  # with more columns than rows the rank bound is the row count and the anchor
  # is the mean band width, where degenerate columns are not counted
  short <- .brm_saturation_fixture(seed = 8704L, n = 60L)
  expect_equal(rMVPA:::.brm_alpha_anchor(short$X, short$groups),
               ncol(short$X) / 3)
  flat <- short$X
  flat[, 1:20] <- 1
  expect_equal(rMVPA:::.brm_alpha_anchor(flat, short$groups),
               (ncol(short$X) - 20) / 3)

  explicit <- .brm_saturation_model(fixture, alphas = c(1, 10))
  expect_identical(explicit$alpha_grid_origin, "user")
  expect_identical(explicit$alpha_grid, c(1, 10))
  expect_true(is.na(explicit$alpha_anchor))
  expect_error(.brm_saturation_model(fixture, alphas = "wide"),
               "must be \"auto\" or a numeric vector", fixed = TRUE)
  # a factor is not a numeric grid either, and used to fall through to the
  # solver's less helpful message
  expect_error(.brm_saturation_model(fixture, alphas = factor("wide")),
               "must be \"auto\" or a numeric vector", fixed = TRUE)
  expect_identical(
    .brm_saturation_model(fixture, alphas = factor("auto"))$alpha_grid,
    model$alpha_grid
  )

  # a supplied manifest ignores alphas entirely, and says so
  manifest <- rMVPA:::.banded_ridge_candidates(
    c("a", "b", "c"), alphas = c(2, 20, 200), method = "fixed",
    theta = c(a = 1 / 3, b = 1 / 3, c = 1 / 3)
  )
  from_manifest <- .brm_saturation_model(fixture, alphas = c(1, 10),
                                         candidates = manifest)
  expect_identical(from_manifest$alpha_grid_origin, "candidates")
  expect_identical(from_manifest$alpha_grid, c(2, 20, 200))
  expect_true(is.na(from_manifest$alpha_anchor))
})

test_that("the saturation warning names the end, the direction, and the remedy", {
  fixture <- .brm_saturation_fixture()
  ceiling_grid <- .brm_saturation_model(
    fixture, alphas = 10^seq(-2, 2, length.out = 9L), delta_sets = "all"
  )
  message <- tryCatch(run_banded_ridge(ceiling_grid),
                      warning = conditionMessage)
  expect_match(message, "The full model took the largest available alpha (100)",
               fixed = TRUE)
  expect_match(message, "3 leave-one-band-out models (a, b, c), the worst of which",
               fixed = TRUE)
  expect_match(message, "refits are under-penalized", fixed = TRUE)
  expect_match(message, "Widen `alphas`", fixed = TRUE)

  # the same grid supplied as a manifest must not point at an argument that
  # is ignored
  manifest <- rMVPA:::.banded_ridge_candidates(
    c("a", "b", "c"), alphas = 10^seq(-2, 2, length.out = 9L),
    method = "fixed", theta = c(a = 1 / 3, b = 1 / 3, c = 1 / 3)
  )
  manifest_message <- tryCatch(
    run_banded_ridge(.brm_saturation_model(fixture, candidates = manifest)),
    warning = conditionMessage
  )
  expect_match(manifest_message, "Widen the alpha values in `candidates`",
               fixed = TRUE)

  # a floor pile-up names the other end and the other direction
  floor_message <- tryCatch(
    run_banded_ridge(.brm_saturation_model(
      .brm_saturation_fixture(seed = 8706L, signal = 0.6),
      alphas = 10^seq(2, 6, length.out = 5L)
    )),
    warning = conditionMessage
  )
  expect_match(floor_message, "took the smallest available alpha", fixed = TRUE)
  expect_match(floor_message, "refits are over-penalized", fixed = TRUE)
})

test_that("a truncated alpha grid is reported instead of silently saturating", {
  fixture <- .brm_saturation_fixture()
  narrow <- .brm_saturation_model(fixture,
                                  alphas = 10^seq(-2, 2, length.out = 9L))
  expect_warning(narrow_result <- run_banded_ridge(narrow),
                 "alpha grid is saturated")

  narrow_alpha <- narrow_result$selection_diagnostics$alpha
  expect_identical(narrow_alpha$model, "full")
  expect_identical(narrow_alpha$share_at_grid_max, 1)
  expect_identical(narrow_alpha$share_interior, 0)
  expect_identical(narrow_alpha$modal_alpha, 100)
  expect_true(narrow_alpha$saturated)

  # the same data on the design-scaled default finds an interior optimum and
  # predicts better than the training mean it previously lost to
  auto_result <- .brm_saturation_model(fixture)
  auto_result <- run_banded_ridge(auto_result)
  auto_alpha <- auto_result$selection_diagnostics$alpha
  expect_false(auto_alpha$saturated)
  expect_identical(auto_alpha$share_interior, 1)
  # the modal alpha is one of the grid's own values, not a reformatted copy
  expect_true(any(vapply(auto_result$hyperparameters$alpha, identical,
                         logical(1), auto_alpha$modal_alpha)))
  expect_gt(auto_alpha$modal_alpha, narrow_alpha$alpha_grid_max)
  expect_gt(auto_result$selection_diagnostics$fit$median_r2,
            narrow_result$selection_diagnostics$fit$median_r2)
  expect_gt(auto_result$selection_diagnostics$fit$share_r2_positive,
            narrow_result$selection_diagnostics$fit$share_r2_positive)

  expect_identical(auto_result$selection_diagnostics$alpha_grid_origin, "auto")
  expect_identical(auto_result$selection_diagnostics$saturation_share, 0.95)
  printed <- utils::capture.output(print(auto_result))
  expect_true(any(grepl("^Alpha: ", printed)))
  expect_true(any(grepl("^OOF R2: ", printed)))
})

test_that("saturation is judged over the alphas a response could have taken", {
  fixture <- .brm_saturation_fixture()

  # a fixed alpha scope selects one value by construction; that is not evidence
  # that the grid is too narrow
  fixed <- .brm_saturation_model(fixture, alphas = c(1, 10, 100),
                                 alpha_scope = "fixed", fixed_alpha = 100)
  expect_identical(fixed$selectable_alphas, 100)
  fixed_result <- suppressWarnings(run_banded_ridge(fixed))
  fixed_alpha <- fixed_result$selection_diagnostics$alpha
  expect_identical(fixed_alpha$n_alpha_grid, 1L)
  expect_false(fixed_alpha$saturated)
  expect_true(is.na(fixed_alpha$share_interior))
  # a grid of one offered no choice, so it has no boundary to be at
  expect_true(is.na(fixed_alpha$share_at_grid_min))
  expect_true(is.na(fixed_alpha$share_at_grid_max))
  expect_match(utils::capture.output(print(fixed_result)),
               "^Alpha:     fixed at 100$", all = FALSE)

  # two candidate alphas have no interior, so share_interior is undefined, but
  # a one-sided pile-up is still evidence that the optimum is past that end
  two <- .brm_saturation_model(fixture, alphas = c(1, 10))
  two_alpha <- expect_warning(
    run_banded_ridge(two), "alpha grid is saturated"
  )$selection_diagnostics$alpha
  expect_identical(two_alpha$n_alpha_grid, 2L)
  expect_true(is.na(two_alpha$share_interior))
  expect_identical(two_alpha$share_at_grid_max, 1)
  expect_true(two_alpha$saturated)
})

test_that("boundary mass split across both ends is still a saturated grid", {
  # under the default per-response alpha scope a mask of signal and noise
  # voxels sends its boundary mass to opposite ends, so neither end alone
  # reaches the threshold while almost nothing lands inside
  set.seed(8705L)
  n <- 60L
  D <- 8L
  dims <- c(4L, 4L, 2L)
  nv <- prod(dims)
  X <- matrix(rnorm(n * 2L * D), n, 2L * D,
              dimnames = list(NULL, paste0("x", seq_len(2L * D))))
  beta <- matrix(0, 2L * D, nv)
  signal <- seq_len(nv / 2L)
  beta[, signal] <- matrix(rnorm(2L * D * length(signal), sd = 2),
                           2L * D, length(signal))
  Y <- X %*% beta + matrix(rnorm(n * nv), n, nv)
  dataset <- mvpa_dataset(
    neuroim2::NeuroVec(array(as.vector(t(Y)), c(dims, n)),
                       neuroim2::NeuroSpace(c(dims, n), c(1, 1, 1))),
    mask = neuroim2::NeuroVol(array(1, dims),
                              neuroim2::NeuroSpace(dims, c(1, 1, 1)))
  )
  features <- feature_sets(X, blocks(a = D, b = D))
  design <- feature_sets_design(
    features, block_var_train = as.integer(cut(seq_len(n), 3L, labels = FALSE)),
    time_series = TRUE
  )
  model <- banded_ridge_model(
    dataset, design, outer_crossval = 3L, tune_crossval = 2L,
    theta_method = "fixed", theta = c(a = 0.5, b = 0.5),
    alphas = c(0.5, 1, 2)
  )
  result <- expect_warning(run_banded_ridge(model),
                           "inside the grid")
  alpha <- result$selection_diagnostics$alpha
  expect_lt(max(alpha$share_at_grid_min, alpha$share_at_grid_max),
            rMVPA:::.brm_alpha_saturation_share)
  expect_lt(alpha$share_interior, 1 - rMVPA:::.brm_alpha_saturation_share)
  expect_true(alpha$saturated)
})

test_that("leave-one-band-out models are diagnosed alongside the full model", {
  fixture <- .brm_saturation_fixture()
  model <- .brm_saturation_model(fixture, delta_sets = "all",
                                 alphas = 10^seq(-4, -2, length.out = 5L))
  result <- expect_warning(run_banded_ridge(model), "leave-one-band-out model")

  alpha <- result$selection_diagnostics$alpha
  expect_identical(alpha$model,
                   c("full", "without_a", "without_b", "without_c"))
  expect_true(all(alpha$saturated))
  expect_identical(nrow(result$selection_diagnostics$fit), 4L)
  expect_identical(result$selection_diagnostics$fit$model, alpha$model)
})

test_that("a fit that loses to the training mean says so", {
  fixture <- .brm_saturation_fixture()
  starved <- .brm_saturation_model(fixture,
                                   alphas = 10^seq(-4, -2, length.out = 5L))
  result <- expect_warning(run_banded_ridge(starved),
                           "median outer out-of-fold R2")
  fit <- result$selection_diagnostics$fit
  expect_lt(fit$median_r2, rMVPA:::.brm_r2_floor)
  expect_identical(fit$n_responses, nrow(result$metrics))
  expect_identical(fit$n_scored, sum(is.finite(result$metrics$r2)))
  expect_equal(fit$share_r2_positive, mean(result$metrics$r2 > 0))
  expect_equal(fit$mean_r2, mean(result$metrics$r2))

  # a fit that beats the mean is not warned about
  expect_silent(run_banded_ridge(.brm_saturation_model(
    .brm_saturation_fixture(seed = 8703L, signal = 0.5)
  )))
})
