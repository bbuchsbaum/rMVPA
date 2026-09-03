context("banded_ridge_tuning")

.brt_test_data <- function(seed = 7201L, n = 30L, p = 6L, v = 3L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  Y <- matrix(rnorm(n * v), n, v,
              dimnames = list(NULL, paste0("voxel_", seq_len(v))))
  feature_groups <- rep(c("a", "b"), length.out = p)
  outer <- rMVPA:::.banded_ridge_make_folds(seq_len(n), 3L, seed = seed + 1L)
  inner <- lapply(seq_along(outer), function(ii) {
    rMVPA:::.banded_ridge_make_folds(
      outer[[ii]]$train, 3L, seed = seed + 10L + ii
    )
  })
  candidates <- rMVPA:::.banded_ridge_candidates(
    c("a", "b"), alphas = c(0.2, 3), method = "fixed",
    theta = rbind(c(a = 0.8, b = 0.2), c(a = 0.25, b = 0.75))
  )
  list(X = X, Y = Y, feature_groups = feature_groups, outer = outer,
       inner = inner, candidates = candidates)
}

.brt_test_exhaustive <- function(data, metric = "mse") {
  n <- nrow(data$Y)
  v <- ncol(data$Y)
  predictions <- matrix(NA_real_, n, v,
                        dimnames = list(as.character(seq_len(n)), colnames(data$Y)))
  results <- vector("list", length(data$outer))
  for (oo in seq_along(data$outer)) {
    outer <- data$outer[[oo]]
    inners <- data$inner[[oo]]
    c_count <- length(data$candidates$alpha)
    fold_scores <- array(NA_real_, c(length(inners), c_count, v))
    candidate_scores <- matrix(NA_real_, c_count, v)
    for (cc in seq_len(c_count)) {
      pred <- matrix(NA_real_, length(outer$train), v,
                     dimnames = list(NULL, colnames(data$Y)))
      lookup <- setNames(seq_along(outer$train), outer$train)
      theta <- setNames(data$candidates$theta[cc, ], c("a", "b"))
      for (ii in seq_along(inners)) {
        inner <- inners[[ii]]
        fitted <- rMVPA:::.banded_ridge_fit_predict_direct(
          data$X, data$Y, inner$train, inner$test,
          groups = data$feature_groups,
          alpha = data$candidates$alpha[[cc]], theta = theta
        )
        scored <- rMVPA:::.banded_ridge_score(
          data$Y[inner$test, , drop = FALSE], fitted$predictions
        )[[metric]]
        fold_scores[ii, cc, ] <- scored
        pred[unname(lookup[as.character(inner$test)]), ] <- fitted$predictions
      }
      candidate_scores[cc, ] <- rMVPA:::.banded_ridge_score(
        data$Y[outer$train, , drop = FALSE], pred
      )[[metric]]
    }
    objective <- if (metric == "mse") candidate_scores else -candidate_scores
    selected <- apply(objective, 2L, which.min)
    outer_pred <- matrix(NA_real_, length(outer$test), v,
                         dimnames = list(NULL, colnames(data$Y)))
    for (j in seq_len(v)) {
      cc <- selected[[j]]
      theta <- setNames(data$candidates$theta[cc, ], c("a", "b"))
      fitted <- rMVPA:::.banded_ridge_fit_predict_direct(
        data$X, data$Y[, j, drop = FALSE], outer$train, outer$test,
        groups = data$feature_groups,
        alpha = data$candidates$alpha[[cc]], theta = theta
      )
      outer_pred[, j] <- fitted$predictions
    }
    predictions[outer$test, ] <- outer_pred
    results[[oo]] <- list(fold_scores = fold_scores,
                          candidate_scores = candidate_scores,
                          selected = selected)
  }
  list(predictions = predictions, results = results)
}

.brt_strip_batch_fits <- function(result) {
  lapply(result$outer_results, function(outer) {
    lapply(outer$outer_fits, function(entry) {
      list(candidate_index = entry$candidate_index,
           responses = entry$responses,
           coefficients = entry$fit$coefficients,
           intercept = entry$fit$intercept,
           x_center = entry$fit$x_center,
           x_scale = entry$fit$x_scale)
    })
  })
}

test_that("candidate manifests cover fixed, grid, and deterministic random simplex paths", {
  grid <- rMVPA:::.banded_ridge_candidates(
    c("a", "b", "c"), c(0, 2), method = "grid", grid_points = 3L
  )
  expect_s3_class(grid, "banded_ridge_candidates")
  expect_equal(length(grid$alpha), 12L)
  expect_equal(rowSums(grid$theta), rep(1, 12), tolerance = 0)
  expect_true(any(grid$theta == 0))
  expect_identical(grid$candidate_id, sort(grid$candidate_id)[order(
    match(sort(grid$candidate_id), grid$candidate_id)
  )])

  set.seed(99)
  before <- .Random.seed
  random_1 <- rMVPA:::.banded_ridge_candidates(
    c("a", "b"), 1, method = "random", n_theta = 8L, seed = 7202L
  )
  expect_identical(.Random.seed, before)
  random_2 <- rMVPA:::.banded_ridge_candidates(
    c("a", "b"), 1, method = "random", n_theta = 8L, seed = 7202L
  )
  random_3 <- rMVPA:::.banded_ridge_candidates(
    c("a", "b"), 1, method = "random", n_theta = 8L, seed = 7203L
  )
  expect_identical(random_1, random_2)
  expect_false(isTRUE(all.equal(random_1$theta, random_3$theta)))

  fixed <- rMVPA:::.banded_ridge_candidates(
    c("a", "b"), c(10, 0.1), method = "fixed", theta = c(b = 0.7, a = 0.3)
  )
  expect_identical(fixed$alpha, c(0.1, 10))
  expect_identical(colnames(fixed$theta), c("a", "b"))
})

test_that("fold construction preserves coverage, blocks, purge gaps, and RNG state", {
  set.seed(100)
  before <- .Random.seed
  ordinary <- rMVPA:::.banded_ridge_make_folds(1:23, 5L, seed = 7204L)
  expect_identical(.Random.seed, before)
  expect_identical(sort(unlist(lapply(ordinary, `[[`, "test"))), 1:23)
  expect_lte(max(vapply(ordinary, function(x) length(x$test), integer(1))) -
               min(vapply(ordinary, function(x) length(x$test), integer(1))), 1L)
  expect_true(all(vapply(ordinary, function(x) !length(intersect(x$train, x$test)),
                         logical(1))))
  expect_identical(ordinary,
                   rMVPA:::.banded_ridge_make_folds(1:23, 5L, seed = 7204L))

  blocks <- rep(paste0("run", 1:6), c(3, 5, 4, 6, 2, 4))
  blocked <- rMVPA:::.banded_ridge_make_folds(
    seq_along(blocks), 3L, blocks = blocks, seed = 7205L
  )
  expect_true(all(vapply(blocked, function(x) {
    !length(intersect(blocks[x$train], blocks[x$test]))
  }, logical(1))))

  purged <- rMVPA:::.banded_ridge_make_folds(1:24, 4L, purge = 2L, seed = 7206L)
  expect_true(all(vapply(purged, function(x) {
    all(vapply(x$train, function(i) all(abs(i - x$test) > 2L), logical(1)))
  }, logical(1))))
  expect_error(rMVPA:::.banded_ridge_make_folds(1:3, 4L), "fewer rows")
  expect_error(rMVPA:::.banded_ridge_make_folds(
    1:6, 3L, blocks = rep(c("a", "b"), each = 3)
  ), "fewer independent blocks")
})

test_that("nested scores, selections, and predictions match an exhaustive oracle", {
  data <- .brt_test_data()
  got <- rMVPA:::.banded_ridge_nested_cv(
    data$X, data$Y, data$outer, data$candidates,
    feature_groups = data$feature_groups, inner_folds = data$inner,
    response_batch_size = 2L, seed = 7207L
  )
  oracle <- .brt_test_exhaustive(data)

  for (oo in seq_along(data$outer)) {
    expect_equal(
      unname(got$outer_results[[oo]]$inner_fold_scores),
      unname(oracle$results[[oo]]$fold_scores), tolerance = 1e-12
    )
    expect_equal(
      unname(got$outer_results[[oo]]$candidate_scores),
      unname(oracle$results[[oo]]$candidate_scores), tolerance = 1e-12
    )
    expect_identical(got$outer_results[[oo]]$selected$candidate_index,
                     as.integer(oracle$results[[oo]]$selected))
  }
  expect_equal(unname(got$predictions), unname(oracle$predictions),
               tolerance = 1e-12)
  expect_equal(got$performance,
               rMVPA:::.banded_ridge_score(data$Y, got$predictions),
               tolerance = 0)
  expect_identical(sort(unique(got$fold_performance$outer_fold)),
                   sort(vapply(data$outer, `[[`, character(1), "id")))
  expect_identical(got$metric, "mse")
  expect_match(got$tie_breaking, "canonical numeric")
  expect_identical(got$provenance$outer_fold_ids, vapply(
    data$outer, `[[`, character(1), "id"
  ))
  expect_equal(got$peak_intermediate_dimensions[["candidates"]], 4L)
})

test_that("correlation and R2 are explicit opt-in selection estimands", {
  data <- .brt_test_data(seed = 7208L)
  for (metric in c("correlation", "r2")) {
    got <- rMVPA:::.banded_ridge_nested_cv(
      data$X, data$Y, data$outer, data$candidates,
      feature_groups = data$feature_groups, inner_folds = data$inner,
      metric = metric, seed = 7209L
    )
    oracle <- .brt_test_exhaustive(data, metric = metric)
    expect_equal(unname(got$predictions), unname(oracle$predictions),
                 tolerance = 1e-12)
    expect_identical(got$metric, metric)
  }
})

test_that("response, shared, mixed, and fixed scopes solve their exact constraints", {
  candidates <- rMVPA:::.banded_ridge_candidates(
    c("a", "b"), c(1, 10), method = "fixed",
    theta = rbind(c(a = 0.2, b = 0.8), c(a = 0.8, b = 0.2))
  )
  scores <- matrix(c(
    1, 8, 5,
    5, 2, 4,
    3, 6, 1,
    4, 3, 7
  ), 4, 3, byrow = TRUE)

  response <- rMVPA:::.brt_select_candidates(
    scores, candidates, "response", "response"
  )
  expect_identical(response, apply(scores, 2, which.min))

  shared <- rMVPA:::.brt_select_candidates(scores, candidates, "shared", "shared")
  expect_identical(shared, rep(which.min(rowMeans(scores)), 3L))

  alpha_values <- sort(unique(candidates$alpha))
  alpha_choices <- lapply(alpha_values, function(alpha) {
    idx <- which(candidates$alpha == alpha)
    chosen <- vapply(seq_len(ncol(scores)), function(j) idx[which.min(scores[idx, j])],
                     integer(1))
    list(chosen = chosen, objective = mean(scores[cbind(chosen, seq_len(ncol(scores)))]))
  })
  alpha_winner <- which.min(vapply(alpha_choices, `[[`, numeric(1), "objective"))
  mixed_alpha <- rMVPA:::.brt_select_candidates(
    scores, candidates, "shared", "response"
  )
  expect_identical(mixed_alpha, alpha_choices[[alpha_winner]]$chosen)

  theta_keys <- rMVPA:::.brt_theta_keys(candidates$theta)
  theta_values <- sort(unique(theta_keys))
  theta_choices <- lapply(theta_values, function(theta) {
    idx <- which(theta_keys == theta)
    chosen <- vapply(seq_len(ncol(scores)), function(j) idx[which.min(scores[idx, j])],
                     integer(1))
    list(chosen = chosen, objective = mean(scores[cbind(chosen, seq_len(ncol(scores)))]))
  })
  theta_winner <- which.min(vapply(theta_choices, `[[`, numeric(1), "objective"))
  mixed_theta <- rMVPA:::.brt_select_candidates(
    scores, candidates, "response", "shared"
  )
  expect_identical(mixed_theta, theta_choices[[theta_winner]]$chosen)

  fixed <- rMVPA:::.brt_select_candidates(
    scores, candidates, "fixed", "response", fixed_alpha = 10
  )
  idx <- which(candidates$alpha == 10)
  expected <- vapply(seq_len(ncol(scores)), function(j) idx[which.min(scores[idx, j])],
                     integer(1))
  expect_identical(fixed, expected)
  expect_error(rMVPA:::.brt_select_candidates(
    scores, candidates, "fixed", "response"
  ), "requires fixed_alpha")
})

test_that("constructive heterogeneous responses recover different informative bands", {
  set.seed(7210)
  n <- 96L
  Xa <- matrix(rnorm(n * 3), n, 3)
  Xb <- matrix(rnorm(n * 3), n, 3)
  X <- cbind(Xa, Xb)
  Y <- cbind(
    response_a = 8 * Xa[, 1] - 5 * Xa[, 2] + rnorm(n, sd = 0.02),
    response_b = -7 * Xb[, 1] + 6 * Xb[, 3] + rnorm(n, sd = 0.02)
  )
  outer <- rMVPA:::.banded_ridge_make_folds(1:n, 4L, seed = 7211L)
  candidates <- rMVPA:::.banded_ridge_candidates(
    c("a", "b"), 0.1, method = "fixed",
    theta = rbind(c(a = 1, b = 0), c(a = 0, b = 1))
  )
  got <- rMVPA:::.banded_ridge_nested_cv(
    X, Y, outer, candidates, feature_groups = rep(c("a", "b"), each = 3),
    inner_v = 3L, seed = 7212L
  )
  expect_true(all(got$selections$theta_a[got$selections$response == "response_a"] == 1))
  expect_true(all(got$selections$theta_b[got$selections$response == "response_b"] == 1))
  expect_true(all(got$selections$inner_score < 1))
})

test_that("outer-test outcome and predictor mutations cannot leak into target-fold fitting", {
  data <- .brt_test_data(seed = 7213L, n = 36L)
  base <- rMVPA:::.banded_ridge_nested_cv(
    data$X, data$Y, data$outer, data$candidates,
    feature_groups = data$feature_groups, inner_folds = data$inner,
    seed = 7214L
  )
  target_test <- data$outer[[1]]$test

  changed_y <- data$Y
  changed_y[target_test, ] <- changed_y[target_test, ] * -1e5 + 7e6
  outcome_mutation <- rMVPA:::.banded_ridge_nested_cv(
    data$X, changed_y, data$outer, data$candidates,
    feature_groups = data$feature_groups, inner_folds = data$inner,
    seed = 7214L
  )
  expect_identical(base$outer_results[[1]]$candidate_scores,
                   outcome_mutation$outer_results[[1]]$candidate_scores)
  expect_identical(base$outer_results[[1]]$selected,
                   outcome_mutation$outer_results[[1]]$selected)
  expect_identical(base$outer_results[[1]]$preprocessing,
                   outcome_mutation$outer_results[[1]]$preprocessing)
  expect_identical(.brt_strip_batch_fits(base)[[1]],
                   .brt_strip_batch_fits(outcome_mutation)[[1]])
  expect_identical(base$outer_results[[1]]$predictions,
                   outcome_mutation$outer_results[[1]]$predictions)

  changed_x <- data$X
  changed_x[target_test, ] <- sweep(changed_x[target_test, ], 2,
                                    10^(seq_len(ncol(changed_x)) + 3), "*") + 1e8
  predictor_mutation <- rMVPA:::.banded_ridge_nested_cv(
    changed_x, data$Y, data$outer, data$candidates,
    feature_groups = data$feature_groups, inner_folds = data$inner,
    seed = 7214L
  )
  expect_identical(base$outer_results[[1]]$candidate_scores,
                   predictor_mutation$outer_results[[1]]$candidate_scores)
  expect_identical(base$outer_results[[1]]$preprocessing,
                   predictor_mutation$outer_results[[1]]$preprocessing)
  expect_identical(.brt_strip_batch_fits(base)[[1]],
                   .brt_strip_batch_fits(predictor_mutation)[[1]])
  expect_false(isTRUE(all.equal(base$outer_results[[1]]$predictions,
                                predictor_mutation$outer_results[[1]]$predictions)))
})

test_that("candidate order, response order, repeated calls, and batching are invariant", {
  data <- .brt_test_data(seed = 7215L)
  args <- list(
    X = data$X, Y = data$Y, outer_folds = data$outer,
    candidates = data$candidates, feature_groups = data$feature_groups,
    inner_folds = data$inner, seed = 7216L
  )
  base <- do.call(rMVPA:::.banded_ridge_nested_cv, c(args, list(response_batch_size = 1L)))
  repeated <- do.call(rMVPA:::.banded_ridge_nested_cv, c(args, list(response_batch_size = 1L)))
  expect_identical(base, repeated)

  batched <- do.call(rMVPA:::.banded_ridge_nested_cv, c(args, list(response_batch_size = Inf)))
  expect_equal(base$predictions, batched$predictions, tolerance = 1e-14)
  expect_identical(base$selections, batched$selections)
  expect_identical(
    lapply(base$outer_results, `[[`, "candidate_scores"),
    lapply(batched$outer_results, `[[`, "candidate_scores")
  )

  order <- rev(seq_along(data$candidates$alpha))
  permuted_candidates <- data$candidates
  permuted_candidates$alpha <- permuted_candidates$alpha[order]
  permuted_candidates$theta <- permuted_candidates$theta[order, , drop = FALSE]
  permuted_candidates$candidate_id <- permuted_candidates$candidate_id[order]
  reordered <- do.call(rMVPA:::.banded_ridge_nested_cv, c(
    args[names(args) != "candidates"],
    list(candidates = permuted_candidates, response_batch_size = 1L)
  ))
  expect_identical(base$predictions, reordered$predictions)
  expect_identical(base$selections, reordered$selections)

  response_order <- c(3, 1, 2)
  response_args <- args
  response_args$Y <- data$Y[, response_order, drop = FALSE]
  response_fit <- do.call(rMVPA:::.banded_ridge_nested_cv,
                          c(response_args, list(response_batch_size = 2L)))
  expect_equal(response_fit$predictions[, order(response_order), drop = FALSE],
               base$predictions, tolerance = 1e-14)
  for (oo in seq_along(data$outer)) {
    expect_equal(
      response_fit$outer_results[[oo]]$candidate_scores[, order(response_order), drop = FALSE],
      base$outer_results[[oo]]$candidate_scores, tolerance = 1e-14
    )
  }
})

test_that("nested fold validation enforces outer containment, blocks, purge, and coverage", {
  data <- .brt_test_data(seed = 7217L, n = 36L)
  bad_inner <- data$inner
  outside <- data$outer[[1]]$test[[1]]
  bad_inner[[1]][[1]]$train[[1]] <- outside
  expect_error(rMVPA:::.banded_ridge_nested_cv(
    data$X, data$Y, data$outer, data$candidates,
    feature_groups = data$feature_groups, inner_folds = bad_inner
  ), "outside its allowed outer-training set")

  duplicated_test <- data$outer
  old_test <- duplicated_test[[1]]$test[[1]]
  repeated_test <- duplicated_test[[2]]$test[[1]]
  duplicated_test[[1]]$test[[1]] <- repeated_test
  duplicated_test[[1]]$train <- sort(c(
    setdiff(duplicated_test[[1]]$train, repeated_test), old_test
  ))
  expect_error(rMVPA:::.banded_ridge_nested_cv(
    data$X, data$Y, duplicated_test, data$candidates,
    feature_groups = data$feature_groups
  ), "cover every allowed row exactly once")

  blocks <- rep(paste0("run", 1:12), each = 3)
  blocked_outer <- rMVPA:::.banded_ridge_make_folds(
    1:36, 3L, blocks = blocks, seed = 7218L
  )
  blocked <- rMVPA:::.banded_ridge_nested_cv(
    data$X, data$Y, blocked_outer, data$candidates,
    feature_groups = data$feature_groups, blocks = blocks, inner_v = 3L,
    seed = 7219L
  )
  expect_true(all(vapply(blocked$outer_results, function(outer) {
    all(vapply(outer$inner_folds, function(inner) {
      !length(intersect(blocks[inner$train], blocks[inner$test]))
    }, logical(1)))
  }, logical(1))))

  purged_outer <- rMVPA:::.banded_ridge_make_folds(1:36, 3L, purge = 1L, seed = 7220L)
  purged <- rMVPA:::.banded_ridge_nested_cv(
    data$X, data$Y, purged_outer, data$candidates,
    feature_groups = data$feature_groups, purge = 1L, inner_v = 2L,
    seed = 7221L
  )
  expect_true(all(vapply(purged$outer_results, function(outer) {
    all(vapply(outer$inner_folds, function(inner) {
      all(vapply(inner$train, function(i) all(abs(i - inner$test) > 1), logical(1)))
    }, logical(1)))
  }, logical(1))))

  split_blocks <- rMVPA:::.banded_ridge_make_folds(1:36, 3L, seed = 7222L)
  expect_error(rMVPA:::.banded_ridge_nested_cv(
    data$X, data$Y, split_blocks, data$candidates,
    feature_groups = data$feature_groups, blocks = blocks
  ), "splits at least one declared block")
})

test_that("ties, constants, invalid candidates, boundaries, and minimum folds are explicit", {
  X <- cbind(a = rep(1, 12), b = rep(2, 12))
  Y <- cbind(constant = rep(5, 12), also_constant = rep(-2, 12))
  outer <- rMVPA:::.banded_ridge_make_folds(1:12, 3L, seed = 7223L)
  candidates <- rMVPA:::.banded_ridge_candidates(
    c("a", "b"), c(0.1, 10), method = "fixed",
    theta = rbind(c(a = 1, b = 0), c(a = 0.5, b = 0.5), c(a = 0, b = 1))
  )
  tied <- rMVPA:::.banded_ridge_nested_cv(
    X, Y, outer, candidates, feature_groups = c("a", "b"),
    inner_v = 2L, seed = 7224L
  )
  expect_true(all(tied$selections$candidate_index == 1L))
  expect_true(all(is.finite(tied$predictions)))
  expect_true(all(is.na(tied$performance$correlation)))
  expect_error(rMVPA:::.banded_ridge_nested_cv(
    X, Y, outer, candidates, feature_groups = c("a", "b"),
    inner_v = 2L, metric = "correlation", seed = 7224L
  ), "no valid candidate")

  theta <- rbind(c(a = 0.5, b = 0.5), c(a = 1e-308, b = 1 - 1e-308))
  invalid <- rMVPA:::.brt_canonical_candidates(c(1, 1e308), theta, c("a", "b"))
  set.seed(7225)
  X2 <- matrix(rnorm(24), 12, 2)
  Y2 <- matrix(rnorm(24), 12, 2,
               dimnames = list(NULL, c("v1", "v2")))
  Y2[, "v2"] <- 4 + 1e-12 * Y2[, "v2"]
  recovered <- rMVPA:::.banded_ridge_nested_cv(
    X2, Y2, outer, invalid, feature_groups = c("a", "b"),
    inner_v = 2L, seed = 7226L
  )
  expect_gt(nrow(recovered$candidate_errors), 0L)
  expect_true(all(recovered$selections$alpha == 1))
  expect_true(all(is.finite(recovered$predictions)))

  min_outer <- rMVPA:::.banded_ridge_make_folds(1:6, 3L, seed = 7227L)
  minimal <- rMVPA:::.banded_ridge_nested_cv(
    X2[1:6, , drop = FALSE], Y2[1:6, , drop = FALSE], min_outer,
    rMVPA:::.banded_ridge_candidates(
      c("a", "b"), 1, method = "fixed", theta = c(a = 0.5, b = 0.5)
    ), feature_groups = c("a", "b"), inner_v = 2L, seed = 7228L
  )
  expect_identical(dim(minimal$predictions), c(6L, 2L))

  wrong_bands <- rMVPA:::.banded_ridge_candidates(
    c("a", "missing"), 1, method = "fixed", theta = c(a = 0.5, missing = 0.5)
  )
  expect_error(rMVPA:::.banded_ridge_nested_cv(
    X2, Y2, outer, wrong_bands, feature_groups = c("a", "b"), inner_v = 2L
  ), "column names must match feature-band names")
})

test_that("pure-noise outer predictions remain calibrated near chance", {
  associations <- numeric(24L * 2L)
  cursor <- 1L
  for (simulation in seq_len(24L)) {
    set.seed(7300L + simulation)
    X <- matrix(rnorm(54 * 4), 54, 4)
    Y <- matrix(rnorm(54 * 2), 54, 2,
                dimnames = list(NULL, c("v1", "v2")))
    outer <- rMVPA:::.banded_ridge_make_folds(1:54, 3L, seed = 7400L + simulation)
    candidates <- rMVPA:::.banded_ridge_candidates(
      c("a", "b"), c(0.2, 5), method = "fixed", theta = c(a = 0.5, b = 0.5)
    )
    fit <- rMVPA:::.banded_ridge_nested_cv(
      X, Y, outer, candidates, feature_groups = rep(c("a", "b"), each = 2),
      inner_v = 2L, seed = 7500L + simulation
    )
    fold_correlations <- matrix(
      fit$fold_performance$correlation,
      ncol = 2L, byrow = TRUE
    )
    associations[cursor:(cursor + 1L)] <- colMeans(fold_correlations)
    cursor <- cursor + 2L
  }
  calibration_mean <- mean(associations)
  calibration_se <- stats::sd(associations) / sqrt(length(associations))
  expect_lt(abs(calibration_mean), 0.08)
  expect_lt(calibration_se, 0.04)
})

test_that("nested tuning rejects malformed scopes, candidates, batches, and folds", {
  data <- .brt_test_data(seed = 7229L)
  expect_error(rMVPA:::.banded_ridge_nested_cv(
    data$X, data$Y, data$outer, data$candidates,
    feature_groups = data$feature_groups, alpha_scope = "voxelish"
  ), "one of")
  expect_error(rMVPA:::.banded_ridge_nested_cv(
    data$X, data$Y, data$outer, data$candidates,
    feature_groups = data$feature_groups, response_batch_size = 0
  ), "positive integer or Inf")
  expect_error(rMVPA:::.banded_ridge_nested_cv(
    data$X, data$Y, data$outer, data$candidates,
    feature_groups = data$feature_groups, inner_v = 1L
  ), "at least two")
  expect_error(rMVPA:::.banded_ridge_candidates(
    c("a", "b"), -1, method = "grid"
  ), "non-negative")
  expect_error(rMVPA:::.banded_ridge_candidates(
    c("a", "b"), 1, method = "fixed", theta = c(a = 0.2, b = 0.2)
  ), "sum to one")
  expect_error(rMVPA:::.banded_ridge_candidates(
    c("a", "b"), 1, method = "fixed", theta = c(0.5, 0.5)
  ), "must be named")
})

test_that("nested CV skips inner tuning only when exactly one candidate is available", {
  expect_identical(rMVPA:::.brt_purge_train(1:10, 5L, 2L),
                   c(1L, 2L, 8L, 9L, 10L))
  expect_identical(rMVPA:::.brt_purge_train(1:10, 5L, 0L), 1:10)

  data <- .brt_test_data(seed = 7230L)
  expect_error(rMVPA:::.banded_ridge_nested_cv(
    data$X, data$Y, data$outer, data$candidates,
    feature_groups = data$feature_groups, inner_v = NA
  ), "exactly one candidate")
  expect_error(rMVPA:::.banded_ridge_nested_cv(
    data$X, data$Y, data$outer, data$candidates,
    feature_groups = data$feature_groups, inner_v = NaN
  ), "at least two")

  one <- rMVPA:::.banded_ridge_candidates(
    c("a", "b"), alphas = 0.2, method = "fixed", theta = c(a = 0.8, b = 0.2)
  )
  skipped <- rMVPA:::.banded_ridge_nested_cv(
    data$X, data$Y, data$outer, one,
    feature_groups = data$feature_groups, inner_v = NA, seed = 7231L
  )
  tuned <- rMVPA:::.banded_ridge_nested_cv(
    data$X, data$Y, data$outer, one,
    feature_groups = data$feature_groups, inner_v = 3L, seed = 7231L
  )
  expect_identical(skipped$inner_tuning, "none")
  expect_identical(tuned$inner_tuning, "nested")
  expect_equal(skipped$predictions, tuned$predictions, tolerance = 1e-12)
  expect_true(all(is.na(skipped$selections$inner_score)))
  expect_true(all(is.finite(tuned$selections$inner_score)))
  expect_true(all(vapply(skipped$outer_results, function(x) {
    length(x$inner_folds) == 0L && dim(x$inner_fold_scores)[[1L]] == 0L &&
      length(x$preprocessing) == 0L
  }, logical(1))))
  expect_identical(unname(lengths(skipped$provenance$inner_fold_ids)),
                   rep(0L, length(data$outer)))

  scoped <- rMVPA:::.banded_ridge_nested_cv(
    data$X, data$Y, data$outer, data$candidates,
    feature_groups = data$feature_groups, inner_v = NA,
    alpha_scope = "fixed", fixed_alpha = 0.2,
    theta_scope = "fixed", fixed_theta = c(a = 0.8, b = 0.2), seed = 7231L
  )
  expect_true(all(scoped$selections$alpha == 0.2))
  expect_true(all(scoped$selections$theta_a == 0.8))
  expect_equal(scoped$predictions, skipped$predictions, tolerance = 1e-12)
})

test_that("preprocessing receipts equal the direct-fit standardization without a solve", {
  data <- .brt_test_data(seed = 7230L)
  X <- data$X
  X[, 3] <- 4.5
  fold <- data$inner[[1]][[1]]
  receipt <- rMVPA:::.brt_preprocessing_receipt(X, data$Y, fold, data$feature_groups)
  fit <- rMVPA:::.banded_ridge_fit_direct(
    X[fold$train, , drop = FALSE], data$Y[fold$train, , drop = FALSE],
    groups = data$feature_groups, lambdas = c(a = 1, b = 1)
  )
  expect_identical(receipt$x_center, fit$x_center)
  expect_identical(receipt$x_scale, fit$x_scale)
  expect_identical(receipt$y_center, fit$y_center)
  expect_identical(receipt$constant_x, fit$constant_x)
  expect_true(unname(receipt$constant_x)[[3]])
  expect_identical(unname(receipt$x_scale)[[3]], 1)
  expect_identical(receipt$fold_id, fold$id)
  expect_identical(receipt$train, fold$train)
  expect_identical(receipt$test, fold$test)
  expect_error(rMVPA:::.brt_preprocessing_receipt(
    X, data$Y, list(id = "f", train = 1L, test = 2:3), data$feature_groups
  ), "at least two training rows")
})

test_that("theta-free Gram reuse and cache eviction leave nested dual results unchanged", {
  data <- .brt_test_data(seed = 7231L, n = 36L, p = 8L, v = 3L)
  candidates <- rMVPA:::.banded_ridge_candidates(
    c("a", "b"), alphas = c(0.3, 3), method = "grid", grid_points = 4L
  )
  run <- function(cache, ...) {
    rMVPA:::.banded_ridge_nested_cv(
      data$X, data$Y, data$outer, candidates,
      feature_groups = data$feature_groups, inner_folds = data$inner,
      solver = "dual_kernel", solver_cache = cache, cache_prefix = "public", ...
    )
  }
  same_results <- function(got, reference) {
    expect_identical(got$predictions, reference$predictions)
    expect_identical(got$selections, reference$selections)
    expect_identical(got$performance, reference$performance)
    expect_identical(lapply(got$outer_results, `[[`, "candidate_scores"),
                     lapply(reference$outer_results, `[[`, "candidate_scores"))
    expect_identical(lapply(got$outer_results, `[[`, "preprocessing"),
                     lapply(reference$outer_results, `[[`, "preprocessing"))
  }
  reference <- run(NULL)
  unlimited <- rMVPA:::.banded_ridge_solver_cache()
  full <- run(unlimited)
  same_results(full, reference)

  n_inner <- sum(lengths(data$inner))
  n_folds <- n_inner + length(data$outer)
  n_candidates <- length(candidates$alpha)
  n_theta <- length(unique(rMVPA:::.brt_theta_keys(candidates$theta)))
  refit_calls <- vapply(full$outer_results, function(x) length(x$outer_fits), integer(1))
  refit_thetas <- vapply(full$outer_results, function(x) {
    length(unique(rMVPA:::.brt_theta_keys(
      candidates$theta[x$selected$candidate_index, , drop = FALSE]
    )))
  }, integer(1))
  total_fit_calls <- n_inner * n_candidates + sum(refit_calls)
  expect_identical(unlimited$kernel_build_count, 2L * n_folds)
  expect_identical(unlimited$kernel_hit_count, 2L * (total_fit_calls - n_folds))
  expect_identical(unlimited$decomposition_count,
                   n_inner * n_theta + sum(refit_thetas))
  expect_identical(unlimited$hit_count,
                   total_fit_calls - unlimited$decomposition_count)

  # a leave-one-band-out model on the same folds shares the retained band's Grams
  keep <- data$feature_groups != "b"
  reduced_candidates <- rMVPA:::.banded_ridge_candidates(
    "a", alphas = c(0.3, 3), method = "fixed", theta = rbind(c(a = 1))
  )
  run_reduced <- function(cache, ...) {
    rMVPA:::.banded_ridge_nested_cv(
      data$X[, keep, drop = FALSE], data$Y, data$outer, reduced_candidates,
      feature_groups = data$feature_groups[keep], inner_folds = data$inner,
      solver = "dual_kernel", solver_cache = cache, ...
    )
  }
  builds_before <- unlimited$kernel_build_count
  shared <- run_reduced(unlimited, cache_prefix = "public::without::b",
                        kernel_cache_prefix = "public")
  expect_identical(unlimited$kernel_build_count, builds_before)
  same_results(shared, run_reduced(NULL))
  isolated <- run_reduced(unlimited, cache_prefix = "public::without::b::isolated")
  expect_identical(unlimited$kernel_build_count, builds_before + n_folds)
  same_results(isolated, shared)
  expect_error(run(NULL, kernel_cache_prefix = ""), "kernel_cache_prefix")

  capped <- rMVPA:::.banded_ridge_solver_cache(
    max_bytes = 3 * as.numeric(object.size(matrix(0, 24, 24)))
  )
  evicted <- run(capped)
  same_results(evicted, reference)
  expect_gt(capped$eviction_count, 0L)
  expect_lte(capped$bytes, capped$max_bytes)

  nothing <- rMVPA:::.banded_ridge_solver_cache(max_bytes = 500)
  same_results(run(nothing), reference)
  expect_identical(length(nothing$entries), 0L)
  expect_gt(nothing$oversize_count, 0L)
  expect_identical(nothing$decomposition_count, total_fit_calls)
})
