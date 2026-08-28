context("banded_ridge predictive leave-one-band-out attribution")

.bra_test_dataset <- function(Y_active, active_ids = c(1L, 4L),
                              dims = c(2L, 2L, 1L)) {
  n <- nrow(Y_active)
  Y <- matrix(0, n, prod(dims))
  Y[, active_ids] <- Y_active
  array_data <- array(as.vector(t(Y)), c(dims, n))
  data_space <- neuroim2::NeuroSpace(c(dims, n), c(1, 1, 1))
  mask_data <- array(0, dims)
  mask_data[active_ids] <- 1
  mask_space <- neuroim2::NeuroSpace(dims, c(1, 1, 1))
  mvpa_dataset(
    neuroim2::NeuroVec(array_data, data_space),
    mask = neuroim2::NeuroVol(mask_data, mask_space)
  )
}

.bra_test_problem <- function(X,
                              Y_active,
                              groups,
                              group_order = unique(groups),
                              blocks_train = NULL,
                              candidates = NULL,
                              delta_sets = "all",
                              target_batch_size = Inf,
                              return_predictions = TRUE,
                              retain_diagnostics = FALSE,
                              solver = "direct",
                              seed = 8301L) {
  n <- nrow(X)
  if (is.null(blocks_train)) {
    blocks_train <- as.integer(cut(seq_len(n), 3L, labels = FALSE))
  }
  fs <- feature_sets(X, by_set(groups, order = group_order))
  design <- feature_sets_design(
    fs, block_var_train = blocks_train, time_series = TRUE
  )
  if (is.null(candidates)) {
    theta <- setNames(rep(1 / length(group_order), length(group_order)),
                      group_order)
    candidates <- rMVPA:::.banded_ridge_candidates(
      group_order, 0, method = "fixed", theta = theta
    )
  }
  dataset <- .bra_test_dataset(Y_active)
  model <- banded_ridge_model(
    dataset, design, outer_crossval = 3L, tune_crossval = 2L,
    candidates = candidates, delta_sets = delta_sets,
    target_batch_size = target_batch_size,
    return_predictions = return_predictions,
    retain_diagnostics = retain_diagnostics,
    solver = solver, seed = seed
  )
  list(
    X = X, Y_active = Y_active, groups = groups, blocks = blocks_train,
    dataset = dataset, design = design, model = model,
    result = run_banded_ridge(model)
  )
}

.bra_factorial <- function(blocks = 3L) {
  z <- rep(c(-1, -1, 1, 1), blocks)
  u <- rep(c(-1, 1, -1, 1), blocks)
  list(z = z, u = u, block = rep(seq_len(blocks), each = 4L))
}

test_that("suppressor adversary proves dropout effects are not an additive partition", {
  f <- .bra_factorial()
  x1 <- f$z + sqrt(19) * f$u
  x2 <- f$z - sqrt(19) * f$u
  y <- 2 * f$z
  problem <- .bra_test_problem(
    cbind(x1, x2), cbind(y, y), c("plus", "minus"),
    blocks_train = f$block
  )
  out <- problem$result$predictive_leave_one_band_out
  effects <- out$effects

  expect_equal(problem$result$metrics$r2, c(1, 1), tolerance = 1e-12)
  expect_equal(out$reduced_performance$plus$r2, c(0.05, 0.05),
               tolerance = 1e-12)
  expect_equal(out$reduced_performance$minus$r2, c(0.05, 0.05),
               tolerance = 1e-12)
  expect_equal(effects$delta_cv_r2_plus, c(0.95, 0.95), tolerance = 1e-12)
  expect_equal(effects$delta_cv_r2_minus, c(0.95, 0.95), tolerance = 1e-12)
  expect_true(all(effects$delta_cv_r2_plus + effects$delta_cv_r2_minus >
                    problem$result$metrics$r2))
  expect_match(out$estimand, "predictive leave-one-band-out")
  expect_match(out$estimand, "not additive")
  expect_false(out$provenance$additive_partition)
  expect_false(out$provenance$clipping)
})

test_that("orthogonal informative bands have the analytic dropout ordering", {
  f <- .bra_factorial()
  y <- 2 * f$z + 0.5 * f$u
  problem <- .bra_test_problem(
    cbind(f$z, f$u), cbind(y, -y), c("strong", "weak"),
    blocks_train = f$block
  )
  effects <- problem$result$predictive_leave_one_band_out$effects

  expect_equal(problem$result$metrics$r2, c(1, 1), tolerance = 1e-12)
  expect_equal(effects$delta_cv_r2_strong, rep(16 / 17, 2), tolerance = 1e-12)
  expect_equal(effects$delta_cv_r2_weak, rep(1 / 17, 2), tolerance = 1e-12)
  expect_true(all(effects$delta_cv_r2_strong > effects$delta_cv_r2_weak))
  expect_true(all(effects[, 3:4] >= 0))
})

test_that("redundant bands have zero individual dropout without a shared-variance label", {
  f <- .bra_factorial()
  problem <- .bra_test_problem(
    cbind(f$z, f$z), cbind(f$z, -f$z), c("copy_a", "copy_b"),
    blocks_train = f$block
  )
  out <- problem$result$predictive_leave_one_band_out

  expect_equal(out$effects$delta_cv_r2_copy_a, c(0, 0), tolerance = 1e-12)
  expect_equal(out$effects$delta_cv_r2_copy_b, c(0, 0), tolerance = 1e-12)
  expect_false(any(grepl("shared", names(out), ignore.case = TRUE)))
  expect_false(grepl("shared", out$estimand, ignore.case = TRUE))
})

test_that("finite-sample negative dropout effects are retained without clipping", {
  set.seed(1)
  n <- 24L
  signal <- rnorm(n)
  nuisance <- matrix(rnorm(n * 10L), n, 10L)
  X <- cbind(signal, nuisance)
  y <- signal + rnorm(n, sd = 0.3)
  problem <- .bra_test_problem(
    X, cbind(y, y), c("signal", rep("nuisance", 10L)),
    blocks_train = rep(1:3, each = 8L), delta_sets = "nuisance"
  )
  effect <- problem$result$predictive_leave_one_band_out$
    effects$delta_cv_r2_nuisance

  expect_true(all(is.finite(effect)))
  expect_true(all(effect < -0.1))
  expect_equal(
    effect,
    problem$result$metrics$r2 -
      problem$result$predictive_leave_one_band_out$
        reduced_performance$nuisance$r2,
    tolerance = 0
  )
})

test_that("wrapper predictions, metrics, selections, and deltas match independent nested runs", {
  set.seed(8302L)
  n <- 18L
  X <- matrix(rnorm(n * 6L), n, 6L)
  groups <- rep(c("low", "semantic"), each = 3L)
  Y <- cbind(
    X[, 1] - 0.5 * X[, 4] + rnorm(n, sd = 0.2),
    -X[, 2] + 0.7 * X[, 5] + rnorm(n, sd = 0.2)
  )
  colnames(Y) <- c("voxel_1", "voxel_4")
  candidates <- rMVPA:::.banded_ridge_candidates(
    c("low", "semantic"), c(0.1, 1), method = "grid", grid_points = 3L
  )
  problem <- .bra_test_problem(
    X, Y, groups, candidates = candidates, target_batch_size = 1L
  )
  model <- problem$model
  public <- problem$result
  full <- rMVPA:::.banded_ridge_nested_cv(
    X, Y, model$outer_folds, model$candidates,
    feature_groups = groups, inner_v = model$inner_v,
    blocks = problem$blocks, seed = model$seed, solver = "direct"
  )

  expect_equal(public$predictions, full$predictions, tolerance = 1e-12)
  expect_equal(public$metrics$r2, full$performance$r2, tolerance = 1e-12)
  for (band in model$delta_sets) {
    keep <- groups != band
    reduced <- rMVPA:::.banded_ridge_nested_cv(
      X[, keep, drop = FALSE], Y, model$outer_folds,
      model$reduced_candidates[[band]], feature_groups = groups[keep],
      inner_v = model$inner_v, blocks = problem$blocks,
      seed = model$seed, solver = "direct"
    )
    got <- public$predictive_leave_one_band_out
    expect_equal(got$reduced_predictions[[band]], reduced$predictions,
                 tolerance = 1e-12)
    expect_equal(got$reduced_performance[[band]]$r2,
                 reduced$performance$r2, tolerance = 1e-12)
    reduced_selection <- reduced$selections[order(
      match(reduced$selections$response, colnames(Y)),
      reduced$selections$outer_fold
    ), ]
    expect_identical(got$reduced_hyperparameters[[band]]$candidate_id,
                     reduced_selection$candidate_id)
    expect_equal(
      got$effects[[paste0("delta_cv_r2_", band)]],
      full$performance$r2 - reduced$performance$r2,
      tolerance = 1e-12
    )
    expect_identical(
      vapply(reduced$outer_results, `[[`, character(1), "fold_id"),
      got$provenance$matched_outer_fold_ids
    )
  }
})

test_that("outer-test response mutation cannot change full or reduced tuning", {
  set.seed(8303L)
  n <- 18L
  X <- matrix(rnorm(n * 4L), n, 4L)
  groups <- rep(c("a", "b"), each = 2L)
  Y <- cbind(X[, 1] + X[, 3] + rnorm(n, sd = 0.2),
             X[, 2] - X[, 4] + rnorm(n, sd = 0.2))
  base <- .bra_test_problem(
    X, Y, groups, retain_diagnostics = TRUE, target_batch_size = Inf
  )
  target <- base$model$outer_folds[[2L]]
  mutated_y <- Y
  mutated_y[target$test, ] <- mutated_y[target$test, ] + 1e6
  mutated <- .bra_test_problem(
    X, mutated_y, groups, blocks_train = base$blocks,
    candidates = base$model$candidates, retain_diagnostics = TRUE,
    target_batch_size = Inf
  )

  full_before <- base$result$diagnostics[[1L]][[2L]]
  full_after <- mutated$result$diagnostics[[1L]][[2L]]
  expect_equal(full_before$candidate_scores, full_after$candidate_scores,
               tolerance = 0)
  expect_identical(full_before$selected$candidate_id,
                   full_after$selected$candidate_id)
  for (band in base$model$delta_sets) {
    before <- base$result$predictive_leave_one_band_out$
      diagnostics[[band]][[1L]][[2L]]
    after <- mutated$result$predictive_leave_one_band_out$
      diagnostics[[band]][[1L]][[2L]]
    expect_equal(before$candidate_scores, after$candidate_scores, tolerance = 0)
    expect_identical(before$selected$candidate_id,
                     after$selected$candidate_id)
    expect_equal(before$preprocessing, after$preprocessing, tolerance = 0)
    expect_equal(
      base$result$predictive_leave_one_band_out$
        reduced_predictions[[band]][target$test, ],
      mutated$result$predictive_leave_one_band_out$
        reduced_predictions[[band]][target$test, ],
      tolerance = 0
    )
  }
})

test_that("band order and response batching change order only, not the estimand", {
  f <- .bra_factorial()
  X <- cbind(f$z, f$u)
  Y <- cbind(2 * f$z + f$u, f$z - 2 * f$u)
  candidates_ab <- rMVPA:::.banded_ridge_candidates(
    c("a", "b"), 0.2, method = "fixed", theta = c(a = 0.5, b = 0.5)
  )
  candidates_ba <- rMVPA:::.banded_ridge_candidates(
    c("b", "a"), 0.2, method = "fixed", theta = c(b = 0.5, a = 0.5)
  )
  ab <- .bra_test_problem(
    X, Y, c("a", "b"), c("a", "b"), f$block, candidates_ab,
    target_batch_size = 1L
  )
  ba <- .bra_test_problem(
    X, Y, c("a", "b"), c("b", "a"), f$block, candidates_ba,
    target_batch_size = Inf
  )

  expect_equal(ab$result$predictions, ba$result$predictions, tolerance = 1e-12)
  expect_equal(
    ab$result$predictive_leave_one_band_out$effects$delta_cv_r2_a,
    ba$result$predictive_leave_one_band_out$effects$delta_cv_r2_a,
    tolerance = 1e-12
  )
  expect_equal(
    ab$result$predictive_leave_one_band_out$effects$delta_cv_r2_b,
    ba$result$predictive_leave_one_band_out$effects$delta_cv_r2_b,
    tolerance = 1e-12
  )
  expect_identical(
    names(ab$result$predictive_leave_one_band_out$reduced_performance),
    c("a", "b")
  )
  expect_identical(
    names(ba$result$predictive_leave_one_band_out$reduced_performance),
    c("b", "a")
  )
})

test_that("zero-weight bands have exactly zero predictive dropout value", {
  f <- .bra_factorial()
  candidates <- rMVPA:::.banded_ridge_candidates(
    c("used", "zero"), c(0.1, 1), method = "fixed",
    theta = c(used = 1, zero = 0)
  )
  problem <- .bra_test_problem(
    cbind(f$z, f$u), cbind(f$z, -f$z), c("used", "zero"),
    blocks_train = f$block, candidates = candidates, delta_sets = "zero"
  )
  out <- problem$result$predictive_leave_one_band_out

  expect_equal(out$reduced_predictions$zero, problem$result$predictions,
               tolerance = 0)
  expect_equal(out$effects$delta_cv_r2_zero, c(0, 0), tolerance = 0)
  expect_equal(out$candidate_manifests$zero$theta[, "used"], c(1, 1),
               tolerance = 0)
})

test_that("candidate projection is finite, normalized, canonical, and fuzz-stable", {
  for (seed in 8310:8349) {
    set.seed(seed)
    g <- sample(2:4, 1L)
    names <- paste0("g", seq_len(g))
    candidates <- rMVPA:::.banded_ridge_candidates(
      names, c(0, 0.3, 2), method = "random", n_theta = 8L, seed = seed
    )
    for (band in names) {
      reduced <- rMVPA:::.bra_reduced_candidates(candidates, band)
      expect_true(all(is.finite(reduced$theta)))
      expect_true(all(reduced$theta >= 0))
      expect_equal(rowSums(reduced$theta), rep(1, nrow(reduced$theta)),
                   tolerance = 1e-14)
      expect_equal(anyDuplicated(reduced$candidate_id), 0L)
      expect_identical(reduced$group_names, setdiff(names, band))
    }
  }
})

test_that("invalid delta requests and degenerate responses have explicit behavior", {
  f <- .bra_factorial()
  X <- cbind(f$z, f$u)
  Y <- cbind(f$z, rep(3, length(f$z)))
  base <- .bra_test_problem(X, Y, c("a", "b"), blocks_train = f$block)

  expect_error(.bra_test_problem(
    X, Y, c("a", "b"), blocks_train = f$block, delta_sets = "missing"
  ), "unknown feature bands")
  expect_error(.bra_test_problem(
    X, Y, c("a", "b"), blocks_train = f$block,
    delta_sets = c("a", "a")
  ), "unique feature-band names")
  expect_error(.bra_test_problem(
    X, Y, c("a", "b"), blocks_train = f$block,
    delta_sets = c("all", "a")
  ), "cannot be combined")
  expect_error(.bra_test_problem(
    matrix(f$z, ncol = 1L), Y, "only", blocks_train = f$block,
    delta_sets = "all"
  ), "cannot drop the only feature band")

  only_removed_mass <- rMVPA:::.banded_ridge_candidates(
    c("a", "b"), 1, method = "fixed", theta = c(a = 1, b = 0)
  )
  expect_error(.bra_test_problem(
    X, Y, c("a", "b"), blocks_train = f$block,
    candidates = only_removed_mass, delta_sets = "a"
  ), "no positive weight")

  out <- base$result$predictive_leave_one_band_out
  expect_true(is.na(base$result$metrics$r2[[2L]]))
  expect_true(is.na(out$effects$delta_cv_r2_a[[2L]]))
  expect_true(is.na(out$effects$delta_cv_r2_b[[2L]]))
})

test_that("disabled attribution performs no reduced work and enabled cost is explicit", {
  f <- .bra_factorial()
  X <- cbind(f$z, f$u)
  Y <- cbind(f$z + f$u, f$z - f$u)
  disabled <- .bra_test_problem(
    X, Y, c("a", "b"), blocks_train = f$block, delta_sets = NULL,
    target_batch_size = 1L, return_predictions = FALSE
  )
  enabled <- .bra_test_problem(
    X, Y, c("a", "b"), blocks_train = f$block, delta_sets = "all",
    target_batch_size = 1L, return_predictions = FALSE
  )

  expect_null(disabled$result$predictive_leave_one_band_out)
  expect_equal(disabled$result$provenance$conceptual_model_count, 1L)
  expect_true(all(disabled$result$provenance$work_manifest$model_id == "full"))
  expect_equal(disabled$result$provenance$nested_cv_invocations, 2L)
  expect_equal(enabled$result$provenance$conceptual_model_count, 3L)
  expect_equal(enabled$result$provenance$nested_cv_invocations, 6L)
  expect_equal(nrow(enabled$result$provenance$work_manifest), 6L)
  expect_identical(
    enabled$result$provenance$work_manifest$model_id,
    rep(c("full", "without_a", "without_b"), 2L)
  )
  expect_equal(enabled$result$provenance$work_manifest$inner_candidate_fit_calls,
               rep(6L, 6L))
  expect_equal(enabled$result$provenance$work_manifest$outer_refit_calls,
               rep(3L, 6L))
  expect_equal(enabled$result$provenance$work_manifest$total_fit_calls,
               rep(9L, 6L))
  expect_gt(enabled$model$retention_estimated_bytes[["attribution"]], 0)
  expect_gt(enabled$model$retention_estimated_bytes[["total"]],
            disabled$model$retention_estimated_bytes[["total"]])
})

test_that("optimized solver preserves dropout effects and reports reduced work", {
  f <- .bra_factorial()
  X <- cbind(f$z, f$u)
  Y <- cbind(2 * f$z + f$u, f$z - 2 * f$u)
  candidates <- rMVPA:::.banded_ridge_candidates(
    c("a", "b"), c(0.1, 1), method = "grid", grid_points = 3L
  )
  direct <- .bra_test_problem(
    X, Y, c("a", "b"), blocks_train = f$block, candidates = candidates,
    target_batch_size = 1L, solver = "direct"
  )
  svd <- .bra_test_problem(
    X, Y, c("a", "b"), blocks_train = f$block, candidates = candidates,
    target_batch_size = 1L, solver = "svd_primal"
  )

  expect_equal(
    svd$result$predictive_leave_one_band_out$effects,
    direct$result$predictive_leave_one_band_out$effects,
    tolerance = 1e-10
  )
  expect_gt(svd$result$provenance$solver_decomposition_count, 0L)
  expect_gt(svd$result$provenance$solver_cache_hits, 0L)
  expect_true(any(
    svd$result$provenance$work_manifest$decompositions_added > 0L &
      !is.na(svd$result$provenance$work_manifest$excluded_band)
  ))
})
