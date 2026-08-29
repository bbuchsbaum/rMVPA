# Fixed-shape characterization of feature RSA's computational hot paths.
#
# Run from the package root, with native thread limits set before R starts:
#
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
#     VECLIB_MAXIMUM_THREADS=1 \
#     Rscript inst/benchmarks/feature_rsa_hotpaths.R /tmp/feature-rsa.csv 3
#
# The driver checks numerical parity before reporting timings. Wall times are
# descriptive measurements, never package-test thresholds.

frperf_internal <- function(name) {
  get(name, envir = asNamespace("rMVPA"), inherits = FALSE)
}

frperf_standardize <- function(x) {
  x <- as.matrix(x)
  center <- colMeans(x)
  scale <- apply(x, 2L, stats::sd)
  list(
    values = base::scale(x, center = center, scale = scale),
    center = center,
    scale = scale
  )
}

frperf_segment_mse <- function(fit, responses) {
  pred <- fit$validation$pred
  segments <- fit$validation$segments
  out <- matrix(
    NA_real_,
    nrow = length(segments),
    ncol = dim(pred)[3L]
  )
  for (s in seq_along(segments)) {
    idx <- segments[[s]]
    for (k in seq_len(ncol(out))) {
      pred_k <- pred[idx, , k, drop = FALSE][, , 1L]
      out[s, k] <- mean(
        (responses[idx, , drop = FALSE] - pred_k)^2
      )
    }
  }
  out
}

frperf_nested_segment_mse <- function(features,
                                      responses,
                                      ncomp,
                                      segments) {
  ncomp_cv <- min(ncomp, nrow(features) - max(lengths(segments)) - 1L)
  out <- matrix(NA_real_, nrow = length(segments), ncol = ncomp_cv)
  all_rows <- seq_len(nrow(features))
  for (s in seq_along(segments)) {
    test_idx <- segments[[s]]
    train_idx <- all_rows[-test_idx]
    sf <- frperf_standardize(features[train_idx, , drop = FALSE])
    sx <- frperf_standardize(responses[train_idx, , drop = FALSE])
    test_features <- base::scale(
      features[test_idx, , drop = FALSE],
      center = sf$center,
      scale = sf$scale
    )
    test_responses <- base::scale(
      responses[test_idx, , drop = FALSE],
      center = sx$center,
      scale = sx$scale
    )
    fit <- pls::plsr(
      sx$values ~ sf$values,
      ncomp = ncomp_cv,
      scale = FALSE,
      validation = "none"
    )
    for (k in seq_len(ncomp_cv)) {
      pred <- drop(stats::predict(fit, newdata = test_features, ncomp = k))
      if (!is.matrix(pred)) {
        pred <- matrix(pred, nrow = length(test_idx))
      }
      out[s, k] <- mean((test_responses - pred)^2)
    }
  }
  out
}

frperf_ridge_block_mse_materialized <- function(features,
                                                 responses,
                                                 lambdas,
                                                 segments) {
  ridge_decomposition <- frperf_internal(
    ".feature_rsa_ridge_decomposition"
  )
  ridge_fit <- frperf_internal(".feature_rsa_ridge_fit_kernel")
  ridge_predict <- frperf_internal(".feature_rsa_ridge_predict_kernel")
  out <- matrix(NA_real_, nrow = length(segments), ncol = length(lambdas))
  all_rows <- seq_len(nrow(features))
  for (s in seq_along(segments)) {
    test <- segments[[s]]
    train <- all_rows[-test]
    sf <- frperf_standardize(features[train, , drop = FALSE])
    sx <- frperf_standardize(responses[train, , drop = FALSE])
    test_features <- base::scale(
      features[test, , drop = FALSE],
      center = sf$center,
      scale = sf$scale
    )
    observed <- base::scale(
      responses[test, , drop = FALSE],
      center = sx$center,
      scale = sx$scale
    )
    decomposition <- ridge_decomposition(sf$values, sx$values)
    for (j in seq_along(lambdas)) {
      fit <- ridge_fit(
        sf$values,
        sx$values,
        lambdas[[j]],
        decomposition = decomposition
      )
      predicted <- ridge_predict(fit, test_features)
      out[s, j] <- mean((observed - predicted) ^ 2)
    }
  }
  out
}

frperf_onesigma <- function(segment_mse) {
  mean_mse <- colMeans(segment_mse)
  min_idx <- which.min(mean_mse)
  threshold <- mean_mse[[min_idx]] +
    stats::sd(segment_mse[, min_idx]) / sqrt(nrow(segment_mse))
  as.integer(which(mean_mse <= threshold)[1L])
}

frperf_pair <- function(reference,
                        candidate,
                        iterations = 3L,
                        operations_per_timing = 1L) {
  iterations <- as.integer(iterations)
  operations_per_timing <- as.integer(operations_per_timing)
  stopifnot(length(iterations) == 1L, !is.na(iterations), iterations >= 1L)
  stopifnot(
    length(operations_per_timing) == 1L,
    !is.na(operations_per_timing),
    operations_per_timing >= 1L
  )

  reference_output <- reference()
  candidate_output <- candidate()
  reference_times <- candidate_times <- numeric(iterations)

  for (i in seq_len(iterations)) {
    order <- if (i %% 2L) {
      c("reference", "candidate")
    } else {
      c("candidate", "reference")
    }
    for (path in order) {
      gc(verbose = FALSE)
      elapsed <- system.time({
        if (identical(path, "reference")) {
          for (j in seq_len(operations_per_timing)) {
            reference_output <- reference()
          }
        } else {
          for (j in seq_len(operations_per_timing)) {
            candidate_output <- candidate()
          }
        }
      })[["elapsed"]] / operations_per_timing
      if (identical(path, "reference")) {
        reference_times[[i]] <- elapsed
      } else {
        candidate_times[[i]] <- elapsed
      }
    }
  }

  list(
    reference_output = reference_output,
    candidate_output = candidate_output,
    reference_times = reference_times,
    candidate_times = candidate_times,
    reference_median = stats::median(reference_times),
    candidate_median = stats::median(candidate_times),
    operations_per_timing = operations_per_timing
  )
}

frperf_empty_row <- function(section,
                             implementation,
                             reference,
                             n_obs,
                             n_features,
                             n_voxels,
                             n_blocks,
                             max_comps,
                             iterations) {
  data.frame(
    section = section,
    implementation = implementation,
    reference = reference,
    n_obs = as.integer(n_obs),
    n_features = as.integer(n_features),
    n_voxels = as.integer(n_voxels),
    n_blocks = as.integer(n_blocks),
    max_comps = as.integer(max_comps),
    iterations = as.integer(iterations),
    operations_per_timing = 1L,
    elapsed_median_seconds = NA_real_,
    elapsed_min_seconds = NA_real_,
    elapsed_max_seconds = NA_real_,
    reference_median_seconds = NA_real_,
    speedup_vs_reference = NA_real_,
    max_abs_prediction_error = NA_real_,
    max_abs_validation_error = NA_real_,
    selected_ncomp = NA_integer_,
    reference_selected_ncomp = NA_integer_,
    selected_lambda = NA_real_,
    holdout_mse = NA_real_,
    holdout_pattern_correlation = NA_real_,
    holdout_rdm_correlation = NA_real_,
    serialized_bytes = NA_real_,
    reference_serialized_bytes = NA_real_,
    payload_reduction_ratio = NA_real_,
    stringsAsFactors = FALSE
  )
}

frperf_fill_timing <- function(row, pair) {
  row$elapsed_median_seconds <- pair$candidate_median
  row$elapsed_min_seconds <- min(pair$candidate_times)
  row$elapsed_max_seconds <- max(pair$candidate_times)
  row$reference_median_seconds <- pair$reference_median
  row$speedup_vs_reference <- pair$reference_median / pair$candidate_median
  row$operations_per_timing <- pair$operations_per_timing
  row
}

frperf_fixture <- function(seed,
                           n_obs,
                           n_features,
                           n_voxels,
                           n_blocks) {
  stopifnot(
    n_obs >= 8L,
    n_features >= 2L,
    n_voxels >= 2L,
    n_blocks >= 2L,
    n_obs %% n_blocks == 0L
  )
  set.seed(seed)
  latent_dim <- min(6L, n_features)
  latent <- matrix(stats::rnorm(n_obs * latent_dim), n_obs, latent_dim)
  features <- latent %*%
    matrix(stats::rnorm(latent_dim * n_features), latent_dim, n_features) +
    matrix(stats::rnorm(n_obs * n_features, sd = 0.25), n_obs, n_features)
  responses <- features[, seq_len(min(n_features, 8L)), drop = FALSE] %*%
    matrix(
      stats::rnorm(min(n_features, 8L) * n_voxels),
      min(n_features, 8L),
      n_voxels
    ) +
    matrix(stats::rnorm(n_obs * n_voxels, sd = 0.75), n_obs, n_voxels)

  list(
    features = features,
    responses = responses,
    blocks = rep(seq_len(n_blocks), each = n_obs %/% n_blocks)
  )
}

frperf_hotpaths <- function(seed = 2026082861L,
                            n_obs = 120L,
                            n_features = 16L,
                            n_voxels = 100L,
                            n_blocks = 6L,
                            max_comps = 10L,
                            iterations = 3L) {
  if (!requireNamespace("rMVPA", quietly = TRUE) ||
      !requireNamespace("pls", quietly = TRUE)) {
    stop("rMVPA and pls are required for the feature RSA benchmark.",
         call. = FALSE)
  }
  fixture <- frperf_fixture(
    seed, n_obs, n_features, n_voxels, n_blocks
  )
  sf <- frperf_standardize(fixture$features)
  sx <- frperf_standardize(fixture$responses)
  predictors <- sf$values
  responses <- sx$values
  max_comps <- min(as.integer(max_comps), n_features, n_obs - 2L)

  fit_kernel <- frperf_internal(".feature_rsa_fit_kernel")
  predict_kernel <- frperf_internal(".feature_rsa_predict_kernel")
  cv_segment_mse <- frperf_internal(".feature_rsa_cv_segment_mse")
  select_segments <- frperf_internal(".feature_rsa_select_from_segment_mse")
  row_cor <- frperf_internal(".feature_rsa_row_cor")
  ridge_block_mse <- frperf_internal(".feature_rsa_ridge_block_mse")
  ridge_lambdas <- frperf_internal(".feature_rsa_ridge_lambdas")(NULL)

  fit_pair <- frperf_pair(
    reference = function() {
      fit <- pls::plsr(
        responses ~ predictors,
        ncomp = max_comps,
        scale = FALSE,
        validation = "none"
      )
      list(
        fit = fit,
        prediction = drop(stats::predict(
          fit, newdata = predictors, ncomp = max_comps
        ))
      )
    },
    candidate = function() {
      fit <- fit_kernel(
        predictors, responses, ncomp = max_comps, method = "pls"
      )
      list(
        fit = fit,
        prediction = predict_kernel(fit, predictors, ncomp = max_comps)
      )
    },
    iterations = iterations,
    operations_per_timing = 25L
  )
  fit_row <- frperf_empty_row(
    "fit_predict", "compact_matrix_kernel", "pls_formula_mvr",
    n_obs, n_features, n_voxels, n_blocks, max_comps, iterations
  )
  fit_row <- frperf_fill_timing(fit_row, fit_pair)
  fit_row$max_abs_prediction_error <- max(abs(
    fit_pair$candidate_output$prediction -
      fit_pair$reference_output$prediction
  ))
  fit_row$serialized_bytes <- length(serialize(
    fit_pair$candidate_output$fit, NULL
  ))
  fit_row$reference_serialized_bytes <- length(serialize(
    fit_pair$reference_output$fit, NULL
  ))
  fit_row$payload_reduction_ratio <-
    fit_row$reference_serialized_bytes / fit_row$serialized_bytes

  loo_segments <- lapply(seq_len(n_obs), as.integer)
  loo_pair <- frperf_pair(
    reference = function() {
      fit <- pls::plsr(
        responses ~ predictors,
        ncomp = max_comps,
        scale = FALSE,
        validation = "LOO"
      )
      mse <- frperf_segment_mse(fit, responses)
      selected <- frperf_onesigma(mse)
      list(
        fit = fit,
        mse = mse,
        selected = selected,
        prediction = drop(stats::predict(
          fit, newdata = predictors, ncomp = selected
        ))
      )
    },
    candidate = function() {
      mse <- cv_segment_mse(
        predictors,
        responses,
        ncomp = max_comps,
        method = "pls",
        segments = loo_segments
      )
      selected <- select_segments(mse, method = "onesigma")
      fit <- fit_kernel(
        predictors, responses, ncomp = max_comps, method = "pls"
      )
      list(
        fit = fit,
        mse = mse,
        selected = selected,
        prediction = predict_kernel(fit, predictors, ncomp = selected)
      )
    },
    iterations = iterations
  )
  loo_row <- frperf_empty_row(
    "component_selection", "streaming_loo", "pls_formula_loo",
    n_obs, n_features, n_voxels, n_blocks, max_comps, iterations
  )
  loo_row <- frperf_fill_timing(loo_row, loo_pair)
  loo_row$max_abs_prediction_error <- max(abs(
    loo_pair$candidate_output$prediction -
      loo_pair$reference_output$prediction
  ))
  loo_row$max_abs_validation_error <- max(abs(
    loo_pair$candidate_output$mse - loo_pair$reference_output$mse
  ))
  loo_row$selected_ncomp <- loo_pair$candidate_output$selected
  loo_row$reference_selected_ncomp <- loo_pair$reference_output$selected
  loo_row$serialized_bytes <- length(serialize(
    loo_pair$candidate_output$fit, NULL
  ))
  loo_row$reference_serialized_bytes <- length(serialize(
    loo_pair$reference_output$fit, NULL
  ))
  loo_row$payload_reduction_ratio <-
    loo_row$reference_serialized_bytes / loo_row$serialized_bytes

  block_segments <- unname(split(seq_len(n_obs), fixture$blocks))
  blocked_pair <- frperf_pair(
    reference = function() {
      mse <- frperf_nested_segment_mse(
        fixture$features,
        fixture$responses,
        ncomp = max_comps,
        segments = block_segments
      )
      selected <- frperf_onesigma(mse)
      fit <- pls::plsr(
        responses ~ predictors,
        ncomp = max_comps,
        scale = FALSE,
        validation = "none"
      )
      list(
        fit = fit,
        mse = mse,
        selected = selected,
        prediction = drop(stats::predict(
          fit, newdata = predictors, ncomp = selected
        ))
      )
    },
    candidate = function() {
      mse <- cv_segment_mse(
        fixture$features,
        fixture$responses,
        ncomp = max_comps,
        method = "pls",
        segments = block_segments,
        fold_standardize = TRUE
      )
      selected <- select_segments(mse, method = "onesigma")
      fit <- fit_kernel(
        predictors, responses, ncomp = max_comps, method = "pls"
      )
      list(
        fit = fit,
        mse = mse,
        selected = selected,
        prediction = predict_kernel(fit, predictors, ncomp = selected)
      )
    },
    iterations = iterations,
    operations_per_timing = 10L
  )
  blocked_row <- frperf_empty_row(
    "component_selection", "streaming_blocked", "pls_formula_blocked",
    n_obs, n_features, n_voxels, n_blocks, max_comps, iterations
  )
  blocked_row <- frperf_fill_timing(blocked_row, blocked_pair)
  blocked_row$max_abs_prediction_error <- max(abs(
    blocked_pair$candidate_output$prediction -
      blocked_pair$reference_output$prediction
  ))
  blocked_row$max_abs_validation_error <- max(abs(
    blocked_pair$candidate_output$mse - blocked_pair$reference_output$mse
  ))
  blocked_row$selected_ncomp <- blocked_pair$candidate_output$selected
  blocked_row$reference_selected_ncomp <- blocked_pair$reference_output$selected

  ridge_block_pair <- frperf_pair(
    reference = function() {
      frperf_ridge_block_mse_materialized(
        fixture$features,
        fixture$responses,
        ridge_lambdas,
        block_segments
      )
    },
    candidate = function() {
      ridge_block_mse(
        fixture$features,
        fixture$responses,
        ridge_lambdas,
        block_segments
      )
    },
    iterations = iterations,
    operations_per_timing = 5L
  )
  ridge_block_row <- frperf_empty_row(
    "ridge_selection",
    "spectral_block_sse",
    "materialized_block_predictions",
    n_obs,
    n_features,
    n_voxels,
    n_blocks,
    max_comps,
    iterations
  )
  ridge_block_row <- frperf_fill_timing(ridge_block_row, ridge_block_pair)
  ridge_block_row$max_abs_validation_error <- max(abs(
    ridge_block_pair$candidate_output - ridge_block_pair$reference_output
  ))

  # Geometry becomes material for long designs. Use a separately seeded,
  # production-shaped matrix so elapsed time is above the timer resolution
  # even when the fitting fixture is intentionally small.
  set.seed(seed + 1L)
  geometry_n_obs <- max(500L, n_obs)
  geometry_n_voxels <- max(120L, n_voxels)
  geometry <- matrix(
    stats::rnorm(geometry_n_obs * geometry_n_voxels),
    geometry_n_obs,
    geometry_n_voxels
  )
  correlation_pair <- frperf_pair(
    reference = function() {
      stats::cor(t(geometry), use = "pairwise.complete.obs")
    },
    candidate = function() row_cor(geometry),
    iterations = iterations,
    operations_per_timing = 20L
  )
  correlation_row <- frperf_empty_row(
    "trial_geometry", "center_tcrossprod", "stats_cor_pairwise",
    geometry_n_obs, n_features, geometry_n_voxels, n_blocks, max_comps,
    iterations
  )
  correlation_row <- frperf_fill_timing(correlation_row, correlation_pair)
  correlation_row$max_abs_prediction_error <- max(abs(
    correlation_pair$candidate_output - correlation_pair$reference_output
  ))

  rows <- rbind(
    fit_row,
    loo_row,
    blocked_row,
    ridge_block_row,
    correlation_row
  )
  if (any(rows$max_abs_prediction_error > 1e-9, na.rm = TRUE) ||
      any(rows$max_abs_validation_error > 1e-9, na.rm = TRUE) ||
      any(rows$selected_ncomp != rows$reference_selected_ncomp,
          na.rm = TRUE)) {
    stop("Feature RSA hot-path parity validation failed.", call. = FALSE)
  }
  rows
}

frperf_prediction_metrics <- function(predicted, observed, row_cor) {
  predicted <- as.matrix(predicted)
  observed <- as.matrix(observed)
  center_rows <- function(x) {
    sweep(x, 1L, rowMeans(x), "-")
  }
  predicted_centered <- center_rows(predicted)
  observed_centered <- center_rows(observed)
  denominator <- sqrt(
    rowSums(predicted_centered ^ 2) * rowSums(observed_centered ^ 2)
  )
  pattern <- rowSums(predicted_centered * observed_centered) / denominator

  predicted_rdm <- 1 - row_cor(predicted)
  observed_rdm <- 1 - row_cor(observed)
  lower <- lower.tri(predicted_rdm)
  list(
    mse = mean((predicted - observed) ^ 2),
    pattern_correlation = mean(pattern[is.finite(pattern)]),
    rdm_correlation = suppressWarnings(stats::cor(
      predicted_rdm[lower],
      observed_rdm[lower],
      method = "spearman",
      use = "complete.obs"
    ))
  )
}

frperf_ridge_estimators <- function(seed = 2026082864L,
                                    n_obs = 120L,
                                    n_features = 16L,
                                    n_voxels = 100L,
                                    n_blocks = 6L,
                                    max_comps = 10L,
                                    iterations = 3L,
                                    include_glmnet = TRUE) {
  fixture <- frperf_fixture(
    seed,
    n_obs,
    n_features,
    n_voxels,
    n_blocks
  )
  holdout_block <- max(fixture$blocks)
  train <- fixture$blocks != holdout_block
  test <- !train
  train_features <- fixture$features[train, , drop = FALSE]
  train_responses <- fixture$responses[train, , drop = FALSE]
  test_features <- fixture$features[test, , drop = FALSE]
  test_responses <- fixture$responses[test, , drop = FALSE]
  train_blocks <- fixture$blocks[train]
  block_segments <- unname(split(seq_len(sum(train)), train_blocks))
  loo_segments <- lapply(seq_len(sum(train)), as.integer)
  max_comps <- min(
    as.integer(max_comps),
    ncol(train_features),
    nrow(train_features) - 2L
  )

  fit_kernel <- frperf_internal(".feature_rsa_fit_kernel")
  predict_kernel <- frperf_internal(".feature_rsa_predict_kernel")
  cv_segment_mse <- frperf_internal(".feature_rsa_cv_segment_mse")
  select_segments <- frperf_internal(".feature_rsa_select_from_segment_mse")
  ridge_lambdas <- frperf_internal(".feature_rsa_ridge_lambdas")(NULL)
  ridge_decomposition <- frperf_internal(
    ".feature_rsa_ridge_decomposition"
  )
  ridge_gcv <- frperf_internal(".feature_rsa_ridge_gcv_scores")
  ridge_loo <- frperf_internal(".feature_rsa_ridge_loo_mse")
  ridge_blocked <- frperf_internal(".feature_rsa_ridge_block_mse")
  ridge_select <- frperf_internal(".feature_rsa_ridge_select_lambda")
  ridge_fit <- frperf_internal(".feature_rsa_ridge_fit_kernel")
  ridge_predict <- frperf_internal(".feature_rsa_ridge_predict_kernel")
  row_cor <- frperf_internal(".feature_rsa_row_cor")

  standardize_outer <- function() {
    sf <- frperf_standardize(train_features)
    sx <- frperf_standardize(train_responses)
    list(
      sf = sf,
      sx = sx,
      test_features = base::scale(
        test_features,
        center = sf$center,
        scale = sf$scale
      ),
      test_responses = base::scale(
        test_responses,
        center = sx$center,
        scale = sx$scale
      )
    )
  }

  methods <- list(
    pls_streaming_loo = function() {
      scaled <- standardize_outer()
      scores <- cv_segment_mse(
        scaled$sf$values,
        scaled$sx$values,
        ncomp = max_comps,
        method = "pls",
        segments = loo_segments
      )
      selected <- select_segments(scores, method = "onesigma")
      fit <- fit_kernel(
        scaled$sf$values,
        scaled$sx$values,
        ncomp = max_comps,
        method = "pls"
      )
      list(
        prediction = predict_kernel(
          fit,
          scaled$test_features,
          ncomp = selected
        ),
        observed = scaled$test_responses,
        ncomp = selected,
        lambda = NA_real_
      )
    },
    pls_streaming_blocked = function() {
      scaled <- standardize_outer()
      scores <- cv_segment_mse(
        train_features,
        train_responses,
        ncomp = max_comps,
        method = "pls",
        segments = block_segments,
        fold_standardize = TRUE
      )
      selected <- select_segments(scores, method = "onesigma")
      fit <- fit_kernel(
        scaled$sf$values,
        scaled$sx$values,
        ncomp = max_comps,
        method = "pls"
      )
      list(
        prediction = predict_kernel(
          fit,
          scaled$test_features,
          ncomp = selected
        ),
        observed = scaled$test_responses,
        ncomp = selected,
        lambda = NA_real_
      )
    },
    ridge_gcv = function() {
      scaled <- standardize_outer()
      decomposition <- ridge_decomposition(
        scaled$sf$values,
        scaled$sx$values
      )
      selected <- ridge_select(
        ridge_gcv(decomposition, ridge_lambdas),
        ridge_lambdas
      )
      fit <- ridge_fit(
        scaled$sf$values,
        scaled$sx$values,
        ridge_lambdas[[selected]],
        decomposition = decomposition
      )
      list(
        prediction = ridge_predict(fit, scaled$test_features),
        observed = scaled$test_responses,
        ncomp = NA_integer_,
        lambda = ridge_lambdas[[selected]]
      )
    },
    ridge_analytic_loo = function() {
      scaled <- standardize_outer()
      decomposition <- ridge_decomposition(
        scaled$sf$values,
        scaled$sx$values
      )
      selected <- ridge_select(
        ridge_loo(decomposition, ridge_lambdas),
        ridge_lambdas,
        one_se = TRUE
      )
      fit <- ridge_fit(
        scaled$sf$values,
        scaled$sx$values,
        ridge_lambdas[[selected]],
        decomposition = decomposition
      )
      list(
        prediction = ridge_predict(fit, scaled$test_features),
        observed = scaled$test_responses,
        ncomp = NA_integer_,
        lambda = ridge_lambdas[[selected]]
      )
    },
    ridge_blocked = function() {
      scaled <- standardize_outer()
      selected <- ridge_select(
        ridge_blocked(
          train_features,
          train_responses,
          ridge_lambdas,
          block_segments
        ),
        ridge_lambdas,
        one_se = TRUE
      )
      decomposition <- ridge_decomposition(
        scaled$sf$values,
        scaled$sx$values
      )
      fit <- ridge_fit(
        scaled$sf$values,
        scaled$sx$values,
        ridge_lambdas[[selected]],
        decomposition = decomposition
      )
      list(
        prediction = ridge_predict(fit, scaled$test_features),
        observed = scaled$test_responses,
        ncomp = NA_integer_,
        lambda = ridge_lambdas[[selected]]
      )
    }
  )

  if (isTRUE(include_glmnet) && requireNamespace("glmnet", quietly = TRUE)) {
    methods$glmnet_alpha0_blocked_cv <- function() {
      scaled <- standardize_outer()
      fit <- glmnet::cv.glmnet(
        x = scaled$sf$values,
        y = scaled$sx$values,
        family = "mgaussian",
        alpha = 0,
        foldid = as.integer(factor(train_blocks)),
        standardize = FALSE,
        intercept = TRUE,
        parallel = FALSE
      )
      prediction <- drop(stats::predict(
        fit,
        newx = scaled$test_features,
        s = "lambda.min"
      ))
      if (!is.matrix(prediction)) {
        prediction <- matrix(prediction, nrow = nrow(scaled$test_features))
      }
      list(
        prediction = prediction,
        observed = scaled$test_responses,
        ncomp = NA_integer_,
        lambda = fit$lambda.min
      )
    }
  }

  outputs <- vector("list", length(methods))
  names(outputs) <- names(methods)
  operations_per_timing <- integer(length(methods))
  for (j in seq_along(methods)) {
    calibration <- system.time(outputs[[j]] <- methods[[j]]())[["elapsed"]]
    operations_per_timing[[j]] <- if (!is.finite(calibration) ||
                                          calibration <= 0) {
      25L
    } else {
      max(1L, min(25L, as.integer(ceiling(0.02 / calibration))))
    }
  }
  timings <- matrix(
    NA_real_,
    nrow = iterations,
    ncol = length(methods),
    dimnames = list(NULL, names(methods))
  )
  for (i in seq_len(iterations)) {
    order <- seq_along(methods)
    if (i %% 2L == 0L) order <- rev(order)
    for (j in order) {
      gc(verbose = FALSE)
      timings[i, j] <- system.time({
        for (operation in seq_len(operations_per_timing[[j]])) {
          outputs[[j]] <- methods[[j]]()
        }
      })[["elapsed"]] / operations_per_timing[[j]]
    }
  }

  reference_time <- stats::median(timings[, "pls_streaming_loo"])
  do.call(rbind, lapply(seq_along(methods), function(j) {
    output <- outputs[[j]]
    metrics <- frperf_prediction_metrics(
      output$prediction,
      output$observed,
      row_cor
    )
    row <- frperf_empty_row(
      "ridge_estimator_comparison",
      names(methods)[[j]],
      "pls_streaming_loo",
      n_obs,
      n_features,
      n_voxels,
      n_blocks,
      max_comps,
      iterations
    )
    row$elapsed_median_seconds <- stats::median(timings[, j])
    row$elapsed_min_seconds <- min(timings[, j])
    row$elapsed_max_seconds <- max(timings[, j])
    row$reference_median_seconds <- reference_time
    row$speedup_vs_reference <- reference_time / row$elapsed_median_seconds
    row$operations_per_timing <- operations_per_timing[[j]]
    row$selected_ncomp <- output$ncomp
    row$selected_lambda <- output$lambda
    row$holdout_mse <- metrics$mse
    row$holdout_pattern_correlation <- metrics$pattern_correlation
    row$holdout_rdm_correlation <- metrics$rdm_correlation
    row
  }))
}

frperf_ridge_accuracy_suite <- function(
    seeds = 2026082864L + seq_len(20L) - 1L,
    n_obs = 120L,
    n_features = 16L,
    n_voxels = 100L,
    n_blocks = 6L,
    max_comps = 10L,
    include_glmnet = TRUE) {
  seeds <- as.integer(seeds)
  if (!length(seeds) || anyNA(seeds)) {
    stop("seeds must contain at least one finite integer.", call. = FALSE)
  }
  do.call(rbind, lapply(seeds, function(seed) {
    rows <- frperf_ridge_estimators(
      seed = seed,
      n_obs = n_obs,
      n_features = n_features,
      n_voxels = n_voxels,
      n_blocks = n_blocks,
      max_comps = max_comps,
      iterations = 1L,
      include_glmnet = include_glmnet
    )
    rows$seed <- seed
    rows
  }))
}

frperf_ridge_accuracy_summary <- function(rows) {
  required <- c(
    "seed", "implementation", "holdout_mse",
    "holdout_pattern_correlation", "holdout_rdm_correlation"
  )
  if (!is.data.frame(rows) || !all(required %in% names(rows))) {
    stop("rows do not contain ridge accuracy-suite output.", call. = FALSE)
  }
  baseline <- rows[
    rows$implementation == "pls_streaming_loo",
    required,
    drop = FALSE
  ]
  names(baseline)[-(1:2)] <- paste0("baseline_", names(baseline)[-(1:2)])
  baseline$implementation <- NULL
  compared <- merge(rows, baseline, by = "seed", all.x = TRUE, sort = FALSE)
  compared$mse_ratio <- compared$holdout_mse /
    compared$baseline_holdout_mse
  compared$pattern_delta <- compared$holdout_pattern_correlation -
    compared$baseline_holdout_pattern_correlation
  compared$rdm_delta <- compared$holdout_rdm_correlation -
    compared$baseline_holdout_rdm_correlation

  implementations <- unique(compared$implementation)
  do.call(rbind, lapply(implementations, function(implementation) {
    values <- compared[compared$implementation == implementation, , drop = FALSE]
    summary <- data.frame(
      implementation = implementation,
      n_seeds = nrow(values),
      median_mse_ratio = stats::median(values$mse_ratio),
      q90_mse_ratio = as.numeric(stats::quantile(
        values$mse_ratio, 0.9, names = FALSE
      )),
      median_pattern_delta = stats::median(values$pattern_delta),
      q10_pattern_delta = as.numeric(stats::quantile(
        values$pattern_delta, 0.1, names = FALSE
      )),
      median_rdm_delta = stats::median(values$rdm_delta),
      q10_rdm_delta = as.numeric(stats::quantile(
        values$rdm_delta, 0.1, names = FALSE
      )),
      stringsAsFactors = FALSE
    )
    # These are prospective equivalence margins for this dense linear fixture,
    # not universal method guarantees: typical MSE no more than 5% worse,
    # 90th-percentile MSE no more than 20% worse, and correlation losses no
    # larger than .01 typically or .03 at the lower decile.
    summary$within_linear_fixture_margin <- with(summary,
      median_mse_ratio <= 1.05 &&
        q90_mse_ratio <= 1.20 &&
        median_pattern_delta >= -0.01 &&
        q10_pattern_delta >= -0.03 &&
        median_rdm_delta >= -0.01 &&
        q10_rdm_delta >= -0.03
    )
    summary
  }))
}

frperf_ridge_accuracy_receipt <- function(rows, source_root = ".") {
  summary <- frperf_ridge_accuracy_summary(rows)
  timing <- do.call(rbind, lapply(
    split(rows, rows$implementation),
    function(values) {
      data.frame(
        implementation = values$implementation[[1L]],
        median_seconds = stats::median(values$elapsed_median_seconds),
        q10_seconds = as.numeric(stats::quantile(
          values$elapsed_median_seconds, 0.1, names = FALSE
        )),
        q90_seconds = as.numeric(stats::quantile(
          values$elapsed_median_seconds, 0.9, names = FALSE
        )),
        median_speedup_vs_pls_loo = stats::median(
          values$speedup_vs_reference
        ),
        stringsAsFactors = FALSE
      )
    }
  ))
  receipt <- merge(summary, timing, by = "implementation", sort = FALSE)
  receipt <- receipt[match(unique(rows$implementation), receipt$implementation), ]
  receipt$seed_first <- min(rows$seed)
  receipt$seed_last <- max(rows$seed)
  receipt$n_obs <- unique(rows$n_obs)
  receipt$n_features <- unique(rows$n_features)
  receipt$n_voxels <- unique(rows$n_voxels)
  receipt$n_blocks <- unique(rows$n_blocks)
  receipt$max_comps <- unique(rows$max_comps)
  receipt$fixture <- "correlated_dense_linear_last_block_holdout"
  receipt$package_version <- as.character(utils::packageVersion("rMVPA"))
  receipt$pls_version <- as.character(utils::packageVersion("pls"))
  receipt$glmnet_version <- if (requireNamespace("glmnet", quietly = TRUE)) {
    as.character(utils::packageVersion("glmnet"))
  } else {
    NA_character_
  }
  receipt$r_version <- R.version.string
  receipt$platform <- R.version$platform
  receipt$omp_threads <- Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
  receipt$openblas_threads <- Sys.getenv(
    "OPENBLAS_NUM_THREADS", unset = NA_character_
  )
  receipt$mkl_threads <- Sys.getenv("MKL_NUM_THREADS", unset = NA_character_)
  receipt$veclib_threads <- Sys.getenv(
    "VECLIB_MAXIMUM_THREADS", unset = NA_character_
  )
  receipt$git_head <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE,
            stderr = FALSE)[[1L]],
    error = function(e) NA_character_
  )
  receipt$worktree_dirty <- tryCatch(
    length(system2("git", c("status", "--porcelain"), stdout = TRUE,
                   stderr = FALSE)) > 0L,
    error = function(e) NA
  )
  receipt$source_fingerprint <- frperf_source_fingerprint(source_root)
  receipt$benchmark_date <- as.character(Sys.Date())
  rownames(receipt) <- NULL
  receipt
}

frperf_worker_payload <- function(seed = 2026082862L,
                                  n_obs = 500L,
                                  n_features = 20L,
                                  n_blocks = 5L) {
  set.seed(seed)
  latent <- matrix(stats::rnorm(n_obs * n_features), n_obs, n_features)
  similarity <- tcrossprod(base::scale(latent))
  design <- rMVPA::feature_rsa_design(
    S = similarity,
    labels = seq_len(n_obs),
    k = n_features,
    max_comps = min(10L, n_features),
    block_var = rep(seq_len(n_blocks), length.out = n_obs)
  )
  generated <- rMVPA::gen_sample_dataset(
    D = c(3L, 3L, 3L), nobs = n_obs, blocks = n_blocks
  )
  model <- rMVPA::feature_rsa_model(
    generated$dataset,
    design,
    method = "pls",
    ncomp_selection = "max",
    crossval = rMVPA::blocked_cross_validation(design$block_var)
  )
  model$.cv_fold_cache <- frperf_internal(".build_cv_fold_cache")(model)

  legacy_worker <- model
  legacy_worker$dataset <- NULL
  compact_worker <- frperf_internal("as_worker_spec")(model)
  reference_bytes <- length(serialize(legacy_worker, NULL))
  compact_bytes <- length(serialize(compact_worker, NULL))

  row <- frperf_empty_row(
    "worker_payload", "compact_feature_design", "dataset_only_removed",
    n_obs, n_features, NA_integer_, n_blocks, design$max_comps, 0L
  )
  row$serialized_bytes <- compact_bytes
  row$reference_serialized_bytes <- reference_bytes
  row$payload_reduction_ratio <- reference_bytes / compact_bytes
  row
}

frperf_end_to_end <- function(seed = 2026082863L,
                              n_obs = 120L,
                              n_features = 16L,
                              n_regions = 10L,
                              n_blocks = 6L,
                              max_comps = 10L,
                              dims = c(10L, 10L, 10L),
                              iterations = 3L) {
  generated <- rMVPA::gen_sample_dataset(
    D = dims,
    nobs = n_obs,
    blocks = n_blocks
  )
  fixture <- frperf_fixture(
    seed, n_obs, n_features, n_voxels = 12L, n_blocks
  )
  design <- rMVPA::feature_rsa_design(
    F = fixture$features,
    labels = seq_len(n_obs),
    max_comps = max_comps,
    block_var = generated$design$block_var
  )
  crossval <- rMVPA::blocked_cross_validation(generated$design$block_var)
  mask_values <- as.vector(generated$dataset$mask)
  region_values <- integer(length(mask_values))
  inside <- which(mask_values > 0)
  region_values[inside] <- rep_len(seq_len(n_regions), length(inside))
  region_mask <- neuroim2::NeuroVol(
    region_values,
    neuroim2::space(generated$dataset$mask)
  )

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::sequential)

  modes <- c("max", "pve", "blocked", "loo")
  mode_results <- lapply(modes, function(mode) {
    model <- rMVPA::feature_rsa_model(
      generated$dataset,
      design,
      method = "pls",
      ncomp_selection = mode,
      crossval = crossval
    )
    run <- function() suppressWarnings(rMVPA::run_regional(
      model,
      region_mask,
      backend = "default",
      preflight = "off"
    ))
    output <- run()
    elapsed <- numeric(iterations)
    for (i in seq_len(iterations)) {
      gc(verbose = FALSE)
      elapsed[[i]] <- system.time(output <- run())[["elapsed"]]
    }
    if (any(output$fits$error)) {
      stop("Feature RSA end-to-end benchmark produced a failed ROI.",
           call. = FALSE)
    }
    list(
      mode = mode,
      elapsed = elapsed,
      ncomp = output$performance_table$ncomp
    )
  })

  loo_median <- stats::median(
    mode_results[[match("loo", modes)]]$elapsed
  )
  do.call(rbind, lapply(mode_results, function(result) {
    row <- frperf_empty_row(
      "regional_end_to_end",
      paste0("ncomp_", result$mode),
      "ncomp_loo",
      n_obs, n_features, length(inside), n_blocks, max_comps, iterations
    )
    row$elapsed_median_seconds <- stats::median(result$elapsed)
    row$elapsed_min_seconds <- min(result$elapsed)
    row$elapsed_max_seconds <- max(result$elapsed)
    row$reference_median_seconds <- loo_median
    row$speedup_vs_reference <- loo_median / row$elapsed_median_seconds
    row$selected_ncomp <- as.integer(stats::median(result$ncomp, na.rm = TRUE))
    row
  }))
}

frperf_source_fingerprint <- function(source_root = ".") {
  paths <- file.path(source_root, c(
    "R/feature_rsa_model.R",
    "R/mvpa_iterate.R",
    "inst/benchmarks/feature_rsa_hotpaths.R"
  ))
  if (!all(file.exists(paths))) return(NA_character_)
  manifest <- tempfile("rmvpa-feature-rsa-source-", fileext = ".txt")
  on.exit(unlink(manifest, force = TRUE), add = TRUE)
  writeLines(
    paste(basename(paths), unname(tools::md5sum(paths)), sep = "="),
    manifest
  )
  unname(as.character(tools::md5sum(manifest)))
}

frperf_run <- function(output_file = NULL,
                       iterations = 3L,
                       include_end_to_end = TRUE,
                       source_root = ".") {
  rows <- rbind(
    frperf_hotpaths(iterations = iterations),
    frperf_ridge_estimators(iterations = iterations),
    frperf_worker_payload()
  )
  if (isTRUE(include_end_to_end)) {
    rows <- rbind(rows, frperf_end_to_end(iterations = iterations))
  }

  rows$package_version <- as.character(utils::packageVersion("rMVPA"))
  rows$pls_version <- as.character(utils::packageVersion("pls"))
  rows$glmnet_version <- if (requireNamespace("glmnet", quietly = TRUE)) {
    as.character(utils::packageVersion("glmnet"))
  } else {
    NA_character_
  }
  rows$r_version <- R.version.string
  rows$platform <- R.version$platform
  rows$omp_threads <- Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
  rows$openblas_threads <- Sys.getenv(
    "OPENBLAS_NUM_THREADS", unset = NA_character_
  )
  rows$mkl_threads <- Sys.getenv("MKL_NUM_THREADS", unset = NA_character_)
  rows$veclib_threads <- Sys.getenv(
    "VECLIB_MAXIMUM_THREADS", unset = NA_character_
  )
  rows$git_head <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE,
            stderr = FALSE)[[1L]],
    error = function(e) NA_character_
  )
  rows$worktree_dirty <- tryCatch(
    length(system2("git", c("status", "--porcelain"), stdout = TRUE,
                   stderr = FALSE)) > 0L,
    error = function(e) NA
  )
  rows$source_fingerprint <- frperf_source_fingerprint(source_root)
  rows$benchmark_date <- as.character(Sys.Date())

  if (!is.null(output_file) && nzchar(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(rows, output_file, row.names = FALSE)
  }
  rows
}

if (sys.nframe() == 0L) {
  if (file.exists("DESCRIPTION")) {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("pkgload is required to benchmark a source checkout.",
           call. = FALSE)
    }
    pkgload::load_all(".", quiet = TRUE)
  }
  args <- commandArgs(trailingOnly = TRUE)
  output_file <- if (length(args) >= 1L && nzchar(args[[1L]])) {
    args[[1L]]
  } else {
    NULL
  }
  iterations <- if (length(args) >= 2L) as.integer(args[[2L]]) else 3L
  result <- frperf_run(output_file, iterations = iterations)
  print(result[, c(
    "section", "implementation", "elapsed_median_seconds",
    "speedup_vs_reference", "max_abs_prediction_error",
    "max_abs_validation_error", "selected_ncomp",
    "payload_reduction_ratio"
  )], row.names = FALSE)
}
