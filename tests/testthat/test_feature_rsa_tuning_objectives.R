library(testthat)
library(rMVPA)

tuning_reference_standardize <- function(x) {
  x <- as.matrix(x)
  center <- colMeans(x)
  scale <- apply(x, 2L, stats::sd)
  list(
    values = base::scale(x, center = center, scale = scale),
    center = center,
    scale = scale
  )
}

tuning_reference_discrimination <- function(predicted, observed) {
  cross_cor <- suppressWarnings(stats::cor(
    t(predicted),
    t(observed),
    use = "pairwise.complete.obs"
  ))
  mean(diag(cross_cor), na.rm = TRUE) -
    mean(cross_cor[row(cross_cor) != col(cross_cor)], na.rm = TRUE)
}

tuning_reference_rank <- function(predicted, observed) {
  cross_cor <- suppressWarnings(stats::cor(
    t(predicted),
    t(observed),
    use = "pairwise.complete.obs"
  ))
  ranks <- vapply(seq_len(nrow(cross_cor)), function(i) {
    values <- cross_cor[i, ]
    denominator <- sum(is.finite(values)) - 1L
    if (denominator > 0L && is.finite(values[[i]])) {
      (sum(values <= values[[i]], na.rm = TRUE) - 1L) / denominator
    } else {
      NA_real_
    }
  }, numeric(1L))
  mean(ranks, na.rm = TRUE)
}

tuning_reference_pattern_score <- function(predicted, observed, objective) {
  switch(
    objective,
    pattern_discrimination =
      tuning_reference_discrimination(predicted, observed),
    pattern_rank_percentile = tuning_reference_rank(predicted, observed),
    stop("unknown reference objective")
  )
}

tuning_reference_component_scores <- function(features,
                                              responses,
                                              ncomp,
                                              method,
                                              segments,
                                              objective,
                                              raw_space = TRUE) {
  ncomp_cv <- min(ncomp, nrow(features) - max(lengths(segments)) - 1L)
  out <- matrix(NA_real_, nrow = length(segments), ncol = ncomp_cv)
  all_rows <- seq_len(nrow(features))
  fit_fun <- if (identical(method, "pls")) pls::plsr else pls::pcr

  for (s in seq_along(segments)) {
    test <- segments[[s]]
    train <- all_rows[-test]
    sf <- tuning_reference_standardize(features[train, , drop = FALSE])
    sx <- tuning_reference_standardize(responses[train, , drop = FALSE])
    test_features <- base::scale(
      features[test, , drop = FALSE],
      center = sf$center,
      scale = sf$scale
    )
    fit <- fit_fun(
      sx$values ~ sf$values,
      ncomp = ncomp_cv,
      scale = FALSE,
      validation = "none"
    )
    for (k in seq_len(ncomp_cv)) {
      predicted_standardized <- drop(stats::predict(
        fit,
        newdata = test_features,
        ncomp = k
      ))
      if (!is.matrix(predicted_standardized)) {
        predicted_standardized <- matrix(
          predicted_standardized,
          nrow = length(test)
        )
      }
      predicted <- if (isTRUE(raw_space)) {
        sweep(
          sweep(predicted_standardized, 2L, sx$scale, "*"),
          2L,
          sx$center,
          "+"
        )
      } else {
        predicted_standardized
      }
      observed <- if (isTRUE(raw_space)) {
        responses[test, , drop = FALSE]
      } else {
        base::scale(
          responses[test, , drop = FALSE],
          center = sx$center,
          scale = sx$scale
        )
      }
      out[s, k] <- tuning_reference_pattern_score(
        predicted,
        observed,
        objective
      )
    }
  }
  out
}

tuning_reference_ridge_scores <- function(features,
                                          responses,
                                          lambdas,
                                          segments,
                                          objective,
                                          raw_space = TRUE) {
  out <- matrix(NA_real_, nrow = length(segments), ncol = length(lambdas))
  all_rows <- seq_len(nrow(features))
  for (s in seq_along(segments)) {
    test <- segments[[s]]
    train <- all_rows[-test]
    sf <- tuning_reference_standardize(features[train, , drop = FALSE])
    sx <- tuning_reference_standardize(responses[train, , drop = FALSE])
    test_features <- base::scale(
      features[test, , drop = FALSE],
      center = sf$center,
      scale = sf$scale
    )
    for (j in seq_along(lambdas)) {
      penalty <- nrow(sf$values) * lambdas[[j]]
      beta <- solve(
        crossprod(sf$values) + diag(penalty, ncol(sf$values)),
        crossprod(sf$values, sx$values)
      )
      predicted_standardized <- test_features %*% beta
      predicted <- if (isTRUE(raw_space)) {
        sweep(
          sweep(predicted_standardized, 2L, sx$scale, "*"),
          2L,
          sx$center,
          "+"
        )
      } else {
        predicted_standardized
      }
      observed <- if (isTRUE(raw_space)) {
        responses[test, , drop = FALSE]
      } else {
        base::scale(
          responses[test, , drop = FALSE],
          center = sx$center,
          scale = sx$scale
        )
      }
      out[s, j] <- tuning_reference_pattern_score(
        predicted,
        observed,
        objective
      )
    }
  }
  out
}

test_that("linear pattern discrimination matches a dense correlation oracle", {
  set.seed(4521)
  predicted <- matrix(rnorm(17 * 31), 17, 31)
  observed <- predicted * 0.35 + matrix(rnorm(17 * 31), 17, 31)

  expected <- tuning_reference_discrimination(predicted, observed)
  got <- rMVPA:::.feature_rsa_pattern_discrimination_linear(
    predicted,
    observed
  )
  expect_equal(got, expected, tolerance = 2e-12)

  row_order <- sample(seq_len(nrow(predicted)))
  column_order <- sample(seq_len(ncol(predicted)))
  expect_equal(
    rMVPA:::.feature_rsa_pattern_discrimination_linear(
      predicted[row_order, column_order],
      observed[row_order, column_order]
    ),
    got,
    tolerance = 2e-12
  )

  repeated <- matrix(rep(seq_len(31L), each = 8L), 8L, 31L)
  expect_equal(
    rMVPA:::.feature_rsa_pattern_discrimination_linear(repeated, repeated),
    0,
    tolerance = 2e-12
  )
})

test_that("linear pattern discrimination fails closed on undefined inputs", {
  expect_error(
    rMVPA:::.feature_rsa_pattern_discrimination_linear(
      matrix(1:12, 3, 4),
      matrix(1:8, 2, 4)
    ),
    "matching matrix dimensions"
  )
  expect_true(is.na(rMVPA:::.feature_rsa_pattern_discrimination_linear(
    matrix(1:4, 1, 4),
    matrix(1:4, 1, 4)
  )))
  nonfinite <- matrix(rnorm(12), 3, 4)
  nonfinite[1, 1] <- NA_real_
  expect_true(is.na(rMVPA:::.feature_rsa_pattern_discrimination_linear(
    nonfinite,
    matrix(rnorm(12), 3, 4)
  )))

  implementation_symbols <- all.names(
    body(rMVPA:::.feature_rsa_pattern_discrimination_linear)
  )
  expect_false(any(implementation_symbols %in% c(
    "tcrossprod",
    ".feature_rsa_row_cor",
    ".feature_rsa_pattern_metrics_blockwise"
  )))
})

test_that("pattern discrimination is centered at chance under independence", {
  scores <- vapply(seq_len(120L), function(seed) {
    set.seed(7000L + seed)
    rMVPA:::.feature_rsa_pattern_discrimination_linear(
      matrix(rnorm(14 * 48), 14, 48),
      matrix(rnorm(14 * 48), 14, 48)
    )
  }, numeric(1L))

  expect_true(all(is.finite(scores)))
  expect_lt(abs(mean(scores)), 0.01)
})

test_that("component score selection handles maximizing and one-SE simplicity", {
  scores <- cbind(
    c(0.75, 0.76, 0.75, 0.76),
    c(0.76, 0.80, 0.72, 0.76),
    c(0.70, 0.71, 0.69, 0.70)
  )
  expect_identical(
    rMVPA:::.feature_rsa_select_from_segment_scores(
      scores,
      maximize = TRUE,
      one_se = FALSE
    ),
    2L
  )
  expect_identical(
    rMVPA:::.feature_rsa_select_from_segment_scores(
      scores,
      maximize = TRUE,
      one_se = TRUE
    ),
    1L
  )
  expect_identical(
    rMVPA:::.feature_rsa_select_from_segment_scores(
      matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2L),
      maximize = TRUE,
      one_se = FALSE
    ),
    1L
  )
})

test_that("blocked component pattern objectives match independent nested oracles", {
  set.seed(4522)
  n <- 32L
  blocks <- rep(seq_len(4L), each = 8L)
  segments <- unname(split(seq_len(n), blocks))
  features <- matrix(rnorm(n * 6), n, 6)
  beta <- matrix(rnorm(6 * 12), 6, 12)
  response_scale <- exp(seq(log(0.15), log(25), length.out = 12L))
  responses <- sweep(
    features %*% beta + matrix(rnorm(n * 12, sd = 0.9), n, 12),
    2L,
    response_scale,
    "*"
  )
  responses <- sweep(responses, 2L, seq(-30, 30, length.out = 12L), "+")
  dset <- gen_sample_dataset(c(3, 3, 3), n, blocks = 4L)
  design <- feature_rsa_design(
    F = features,
    labels = seq_len(n),
    max_comps = 5L,
    block_var = blocks
  )

  for (method in c("pls", "pca")) {
    for (objective in c(
      "pattern_discrimination",
      "pattern_rank_percentile"
    )) {
      expected <- tuning_reference_component_scores(
        features, responses, 5L, method, segments, objective
      )
      got <- rMVPA:::.feature_rsa_cv_segment_scores(
        features,
        responses,
        ncomp = 5L,
        method = method,
        segments = segments,
        fold_standardize = TRUE,
        objective = objective
      )
      expect_equal(unname(got), unname(expected), tolerance = 5e-10)
      if (method == "pca" && objective == "pattern_discrimination") {
        standardized_space <- tuning_reference_component_scores(
          features,
          responses,
          5L,
          method,
          segments,
          objective,
          raw_space = FALSE
        )
        expect_gt(max(abs(expected - standardized_space)), 0.01)
      }

      expected_ncomp <- which.max(colMeans(expected))
      model <- feature_rsa_model(
        dset$dataset,
        design,
        method = method,
        ncomp_selection = "blocked",
        ncomp_objective = objective,
        crossval = blocked_cross_validation(blocks)
      )
      fit <- train_model(
        model,
        responses,
        features,
        indices = seq_len(ncol(responses)),
        observation_indices = seq_len(n)
      )
      expect_null(fit$error)
      expect_identical(fit$ncomp, as.integer(expected_ncomp))
      expect_identical(fit$ncomp_selection, "blocked")
      expect_identical(fit$ncomp_objective, objective)
      expect_false(fit$ncomp_one_se)
    }
  }
})

test_that("blocked ridge discrimination matches an independent raw-space oracle", {
  set.seed(4523)
  n <- 36L
  blocks <- rep(seq_len(4L), each = 9L)
  segments <- unname(split(seq_len(n), blocks))
  features <- matrix(rnorm(n * 7), n, 7)
  beta <- matrix(rnorm(7 * 14), 7, 14)
  response_scale <- exp(seq(log(0.1), log(30), length.out = 14L))
  responses <- sweep(
    features %*% beta + matrix(rnorm(n * 14, sd = 1.1), n, 14),
    2L,
    response_scale,
    "*"
  )
  responses <- sweep(responses, 2L, seq(-50, 50, length.out = 14L), "+")
  lambdas <- c(0.001, 0.02, 0.3, 4)

  expected <- tuning_reference_ridge_scores(
    features,
    responses,
    lambdas,
    segments,
    "pattern_discrimination"
  )
  got <- rMVPA:::.feature_rsa_ridge_block_pattern_discrimination(
    features,
    responses,
    lambdas,
    segments
  )
  expect_equal(unname(got), unname(expected), tolerance = 5e-10)
  standardized_space <- tuning_reference_ridge_scores(
    features,
    responses,
    lambdas,
    segments,
    "pattern_discrimination",
    raw_space = FALSE
  )
  expect_gt(max(abs(expected - standardized_space)), 0.01)
  expect_identical(
    rMVPA:::.feature_rsa_ridge_select_lambda(-got, lambdas),
    which.max(colMeans(expected))
  )

  dset <- gen_sample_dataset(c(3, 3, 3), n, blocks = 4L)
  design <- feature_rsa_design(
    F = features,
    labels = seq_len(n),
    max_comps = 5L,
    block_var = blocks
  )
  model <- feature_rsa_model(
    dset$dataset,
    design,
    method = "ridge",
    lambda = lambdas,
    lambda_selection = "blocked",
    lambda_objective = "pattern_discrimination",
    crossval = blocked_cross_validation(blocks)
  )
  fit <- train_model(
    model,
    responses,
    features,
    indices = seq_len(ncol(responses)),
    observation_indices = seq_len(n)
  )
  expect_null(fit$error)
  expect_equal(fit$lambda, lambdas[[which.max(colMeans(expected))]])
  expect_identical(fit$lambda_objective, "pattern_discrimination")
  expect_false(fit$lambda_one_se)
})

test_that("MSE remains the public tuning default for this release", {
  set.seed(4524)
  n <- 24L
  blocks <- rep(seq_len(3L), each = 8L)
  dset <- gen_sample_dataset(c(3, 3, 3), n, blocks = 3L)
  design <- feature_rsa_design(
    F = matrix(rnorm(n * 5), n, 5),
    labels = seq_len(n),
    max_comps = 4L,
    block_var = blocks
  )

  component_default <- feature_rsa_model(
    dset$dataset,
    design,
    method = "pca",
    crossval = blocked_cross_validation(blocks)
  )
  expect_identical(component_default$ncomp_objective, "mse")
  expect_true(component_default$ncomp_one_se)

  ridge_default <- feature_rsa_model(
    dset$dataset,
    design,
    method = "ridge",
    crossval = blocked_cross_validation(blocks)
  )
  expect_identical(ridge_default$lambda_objective, "mse")
  expect_false(ridge_default$lambda_one_se)

  component_discrimination <- feature_rsa_model(
    dset$dataset,
    design,
    method = "pls",
    ncomp_selection = "blocked",
    ncomp_objective = "pattern_discrimination",
    crossval = blocked_cross_validation(blocks)
  )
  expect_false(component_discrimination$ncomp_one_se)

  expect_error(feature_rsa_model(
    dset$dataset,
    design,
    method = "pls",
    ncomp_selection = "loo",
    ncomp_objective = "pattern_discrimination",
    crossval = blocked_cross_validation(blocks)
  ), "requires ncomp_selection='blocked'")
  expect_error(feature_rsa_model(
    dset$dataset,
    design,
    method = "ridge",
    lambda_selection = "gcv",
    lambda_objective = "pattern_discrimination",
    crossval = blocked_cross_validation(blocks)
  ), "requires lambda_selection='blocked'")

  component_rank <- feature_rsa_model(
    dset$dataset,
    design,
    method = "pca",
    ncomp_selection = "blocked",
    ncomp_objective = "pattern_rank_percentile",
    crossval = blocked_cross_validation(blocks)
  )
  expect_identical(component_rank$ncomp_objective, "pattern_rank_percentile")
  expect_false(component_rank$ncomp_one_se)
})

test_that("pattern-relative tuning rejects singleton candidate sets", {
  set.seed(4525)
  features <- matrix(rnorm(12 * 4), 12, 4)
  responses <- matrix(rnorm(12 * 8), 12, 8)
  singleton_segments <- as.list(seq_len(12L))

  expect_error(
    rMVPA:::.feature_rsa_cv_segment_scores(
      features,
      responses,
      ncomp = 3L,
      method = "pca",
      segments = singleton_segments,
      fold_standardize = TRUE,
      objective = "pattern_discrimination"
    ),
    "at least two held-out observations"
  )
  expect_error(
    rMVPA:::.feature_rsa_ridge_block_pattern_discrimination(
      features,
      responses,
      c(0.1, 1),
      singleton_segments
    ),
    "at least two held-out observations"
  )
})

test_that("pattern-relative tuning fails closed on undefined pattern scores", {
  set.seed(4526)
  features <- matrix(rnorm(16 * 4), 16, 4)
  # Every response variable varies over observations, but every spatial pattern
  # is constant and therefore has undefined row correlation.
  responses <- outer(seq_len(16L), rep(1, 8L))
  segments <- split(seq_len(16L), rep(seq_len(4L), each = 4L))

  expect_error(
    rMVPA:::.feature_rsa_cv_segment_scores(
      features,
      responses,
      ncomp = 3L,
      method = "pca",
      segments = unname(segments),
      fold_standardize = TRUE,
      objective = "pattern_discrimination"
    ),
    "produced an undefined score"
  )
  expect_error(
    rMVPA:::.feature_rsa_ridge_block_pattern_discrimination(
      features,
      responses,
      c(0.1, 1),
      unname(segments)
    ),
    "produced an undefined score"
  )
})
