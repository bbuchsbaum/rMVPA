library(testthat)
library(rMVPA)

ridge_reference_standardize <- function(x) {
  x <- as.matrix(x)
  center <- colMeans(x)
  scale <- apply(x, 2L, stats::sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  list(
    values = base::scale(x, center = center, scale = scale),
    center = center,
    scale = scale
  )
}

ridge_reference_fit <- function(predictors, responses, lambda) {
  predictors <- as.matrix(predictors)
  responses <- as.matrix(responses)
  x_mean <- colMeans(predictors)
  y_mean <- colMeans(responses)
  xc <- sweep(predictors, 2L, x_mean, "-")
  yc <- sweep(responses, 2L, y_mean, "-")
  penalty <- nrow(predictors) * lambda
  beta <- solve(
    crossprod(xc) + diag(penalty, ncol(xc)),
    crossprod(xc, yc)
  )
  list(beta = beta, x_mean = x_mean, y_mean = y_mean)
}

ridge_reference_predict <- function(fit, newdata) {
  sweep(
    sweep(as.matrix(newdata), 2L, fit$x_mean, "-") %*% fit$beta,
    2L,
    fit$y_mean,
    "+"
  )
}

ridge_reference_gcv <- function(predictors, responses, lambdas) {
  predictors <- as.matrix(predictors)
  responses <- as.matrix(responses)
  n <- nrow(predictors)
  xc <- sweep(predictors, 2L, colMeans(predictors), "-")
  out <- numeric(length(lambdas))
  for (j in seq_along(lambdas)) {
    penalty <- n * lambdas[[j]]
    smoother_centered <- xc %*% solve(
      crossprod(xc) + diag(penalty, ncol(xc)),
      t(xc)
    )
    smoother <- smoother_centered + matrix(1 / n, n, n)
    residual <- responses - smoother %*% responses
    effective_df <- sum(diag(smoother))
    out[[j]] <- mean(residual ^ 2) /
      (1 - effective_df / n) ^ 2
  }
  out
}

ridge_reference_loo_mse <- function(predictors, responses, lambdas) {
  predictors <- as.matrix(predictors)
  responses <- as.matrix(responses)
  n <- nrow(predictors)
  out <- matrix(NA_real_, nrow = n, ncol = length(lambdas))
  for (j in seq_along(lambdas)) {
    for (i in seq_len(n)) {
      fit <- ridge_reference_fit(
        predictors[-i, , drop = FALSE],
        responses[-i, , drop = FALSE],
        lambdas[[j]]
      )
      predicted <- ridge_reference_predict(
        fit,
        predictors[i, , drop = FALSE]
      )
      out[i, j] <- mean((responses[i, , drop = FALSE] - predicted) ^ 2)
    }
  }
  out
}

ridge_reference_block_mse <- function(features,
                                       responses,
                                       lambdas,
                                       segments) {
  all_rows <- seq_len(nrow(features))
  out <- matrix(
    NA_real_,
    nrow = length(segments),
    ncol = length(lambdas)
  )
  for (s in seq_along(segments)) {
    test <- segments[[s]]
    train <- all_rows[-test]
    sf <- ridge_reference_standardize(features[train, , drop = FALSE])
    sy <- ridge_reference_standardize(responses[train, , drop = FALSE])
    ftest <- base::scale(
      features[test, , drop = FALSE],
      center = sf$center,
      scale = sf$scale
    )
    ytest <- base::scale(
      responses[test, , drop = FALSE],
      center = sy$center,
      scale = sy$scale
    )
    for (j in seq_along(lambdas)) {
      fit <- ridge_reference_fit(sf$values, sy$values, lambdas[[j]])
      predicted <- ridge_reference_predict(fit, ftest)
      out[s, j] <- mean((ytest - predicted) ^ 2)
    }
  }
  out
}

ridge_reference_one_se <- function(score_matrix, lambdas) {
  means <- colMeans(score_matrix)
  best <- which.min(means)
  threshold <- means[[best]] +
    stats::sd(score_matrix[, best]) / sqrt(nrow(score_matrix))
  eligible <- which(means <= threshold)
  eligible[[which.max(lambdas[eligible])]]
}

test_that("fixed feature RSA ridge matches the normalized direct solve", {
  set.seed(4501)
  predictors <- matrix(rnorm(42 * 7, mean = 2), 42, 7)
  responses <- matrix(rnorm(42 * 11, mean = -1), 42, 11)
  newdata <- matrix(rnorm(8 * 7, mean = 2), 8, 7)
  colnames(predictors) <- colnames(newdata) <- paste0("f", seq_len(7))
  colnames(responses) <- paste0("v", seq_len(11))
  lambda <- 0.075

  reference <- ridge_reference_fit(predictors, responses, lambda)
  expected <- ridge_reference_predict(reference, newdata)
  fit <- rMVPA:::.feature_rsa_ridge_fit_kernel(
    predictors,
    responses,
    lambda
  )
  got <- rMVPA:::.feature_rsa_ridge_predict_kernel(fit, newdata)

  expect_s3_class(fit, "feature_rsa_ridge_fit")
  expect_equal(unname(fit$coefficients), unname(reference$beta),
               tolerance = 2e-11)
  expect_equal(unname(got), unname(expected), tolerance = 2e-11)
  expect_identical(colnames(got), colnames(responses))
  expect_equal(fit$lambda, lambda)
  expect_true(is.finite(fit$effective_df))
  expect_true(fit$effective_df >= 1 && fit$effective_df <= 1 + ncol(predictors))
})

test_that("spectral ridge GCV matches an explicit hat-matrix oracle", {
  set.seed(4502)
  predictors <- matrix(rnorm(30 * 6), 30, 6)
  responses <- matrix(rnorm(30 * 9), 30, 9)
  lambdas <- c(0, 1e-4, 0.01, 0.2, 2, 20)

  decomposition <- rMVPA:::.feature_rsa_ridge_decomposition(
    predictors,
    responses
  )
  got <- rMVPA:::.feature_rsa_ridge_gcv_scores(
    decomposition,
    lambdas
  )
  expected <- ridge_reference_gcv(predictors, responses, lambdas)

  expect_equal(got, expected, tolerance = 2e-11)
  expect_identical(
    rMVPA:::.feature_rsa_ridge_select_lambda(got, lambdas),
    which.min(expected)
  )
})

test_that("analytic ridge LOO matches exhaustive refitting", {
  set.seed(4503)
  predictors <- matrix(rnorm(18 * 5), 18, 5)
  responses <- matrix(rnorm(18 * 7), 18, 7)
  lambdas <- c(0, 0.002, 0.03, 0.4, 4)

  decomposition <- rMVPA:::.feature_rsa_ridge_decomposition(
    predictors,
    responses
  )
  got <- rMVPA:::.feature_rsa_ridge_loo_mse(
    decomposition,
    lambdas
  )
  expected <- ridge_reference_loo_mse(predictors, responses, lambdas)

  expect_equal(unname(got), unname(expected), tolerance = 3e-10)
  expect_identical(
    rMVPA:::.feature_rsa_ridge_select_lambda(
      got,
      lambdas,
      one_se = TRUE
    ),
    ridge_reference_one_se(expected, lambdas)
  )
})

test_that("blocked ridge tuning reproduces nested standardization oracle", {
  set.seed(4504)
  features <- matrix(rnorm(36 * 6, mean = 3, sd = 2), 36, 6)
  responses <- matrix(rnorm(36 * 10, mean = -2, sd = 4), 36, 10)
  blocks <- rep(paste0("run", 1:4), each = 9)
  segments <- unname(split(seq_len(36), blocks))
  lambdas <- c(0, 0.001, 0.02, 0.3, 3)

  got <- rMVPA:::.feature_rsa_ridge_block_mse(
    features,
    responses,
    lambdas,
    segments
  )
  expected <- ridge_reference_block_mse(
    features,
    responses,
    lambdas,
    segments
  )

  expect_equal(unname(got), unname(expected), tolerance = 3e-10)
  expect_identical(
    rMVPA:::.feature_rsa_ridge_select_lambda(
      got,
      lambdas,
      one_se = TRUE
    ),
    ridge_reference_one_se(expected, lambdas)
  )
})

test_that("ridge public contract validates selection and lambda", {
  set.seed(4505)
  dset <- gen_sample_dataset(c(3, 3, 3), 24, blocks = 3)
  design <- feature_rsa_design(
    F = matrix(rnorm(24 * 5), 24, 5),
    labels = seq_len(24),
    max_comps = 4L,
    block_var = dset$design$block_var
  )

  model <- feature_rsa_model(
    dset$dataset,
    design,
    method = "ridge",
    lambda_selection = "gcv"
  )
  expect_identical(model$method, "ridge")
  expect_identical(model$lambda_selection, "gcv")
  expect_null(model$lambda)

  expect_error(feature_rsa_model(
    dset$dataset,
    design,
    method = "ridge",
    lambda_selection = "fixed"
  ), "requires exactly one")
  expect_error(feature_rsa_model(
    dset$dataset,
    design,
    method = "ridge",
    lambda_selection = "fixed",
    lambda = c(0.1, 1)
  ), "requires exactly one")
  expect_error(feature_rsa_model(
    dset$dataset,
    design,
    method = "ridge",
    lambda = -1
  ), "non-negative")

  no_blocks <- feature_rsa_design(
    F = design$F,
    labels = seq_len(24),
    max_comps = 4L
  )
  expect_error(feature_rsa_model(
    dset$dataset,
    no_blocks,
    method = "ridge",
    lambda_selection = "blocked",
    crossval = blocked_cross_validation(dset$design$block_var)
  ), "requires design\\$block_var")
})

test_that("ridge Feature RSA training and prediction preserve direct-solve output", {
  set.seed(4506)
  dset <- gen_sample_dataset(c(3, 3, 3), 30, blocks = 3)
  features <- matrix(rnorm(30 * 6), 30, 6)
  responses <- matrix(rnorm(30 * 12), 30, 12)
  design <- feature_rsa_design(
    F = features,
    labels = seq_len(30),
    max_comps = 5L,
    block_var = dset$design$block_var
  )
  model <- feature_rsa_model(
    dset$dataset,
    design,
    method = "ridge",
    lambda_selection = "fixed",
    lambda = 0.15
  )
  fit <- train_model(
    model,
    responses,
    features,
    indices = seq_len(ncol(responses)),
    observation_indices = seq_len(nrow(features))
  )
  got <- predict_model(model, fit, features)

  sf <- ridge_reference_standardize(features)
  sy <- ridge_reference_standardize(responses)
  reference <- ridge_reference_fit(sf$values, sy$values, 0.15)
  expected_z <- ridge_reference_predict(reference, sf$values)
  expected <- sweep(
    sweep(expected_z, 2L, sy$scale, "*"),
    2L,
    sy$center,
    "+"
  )

  expect_false(isTRUE(fit$error))
  expect_s3_class(fit$trained_model, "feature_rsa_ridge_fit")
  expect_equal(fit$lambda, 0.15)
  expect_equal(fit$effective_df, fit$trained_model$effective_df)
  expect_true(is.na(fit$ncomp))
  expect_equal(unname(got), unname(expected), tolerance = 3e-10)
})

test_that("ridge kernel is stable under collinearity and feature permutation", {
  set.seed(4507)
  base <- matrix(rnorm(28 * 4), 28, 4)
  predictors <- cbind(base, base[, 1] + base[, 2], base[, 3])
  responses <- matrix(rnorm(28 * 8), 28, 8)
  newdata <- matrix(rnorm(5 * ncol(predictors)), 5, ncol(predictors))
  order <- c(6, 2, 5, 1, 4, 3)

  fit <- rMVPA:::.feature_rsa_ridge_fit_kernel(
    predictors,
    responses,
    lambda = 0.02
  )
  permuted <- rMVPA:::.feature_rsa_ridge_fit_kernel(
    predictors[, order, drop = FALSE],
    responses,
    lambda = 0.02
  )
  expected <- rMVPA:::.feature_rsa_ridge_predict_kernel(fit, newdata)
  got <- rMVPA:::.feature_rsa_ridge_predict_kernel(
    permuted,
    newdata[, order, drop = FALSE]
  )

  expect_true(all(is.finite(fit$coefficients)))
  expect_equal(got, expected, tolerance = 3e-10)
})

test_that("public ridge selectors reproduce independent selection oracles", {
  set.seed(4508)
  n <- 32L
  dset <- gen_sample_dataset(c(3, 3, 3), n, blocks = 4)
  features <- matrix(rnorm(n * 6, mean = 1.5), n, 6)
  responses <- matrix(rnorm(n * 9, mean = -0.5), n, 9)
  blocks <- dset$design$block_var
  lambdas <- c(0.001, 0.01, 0.1, 1, 10)
  design <- feature_rsa_design(
    F = features,
    labels = seq_len(n),
    max_comps = 5L,
    block_var = blocks
  )

  sf <- ridge_reference_standardize(features)
  sy <- ridge_reference_standardize(responses)
  expected <- c(
    gcv = lambdas[[which.min(ridge_reference_gcv(
      sf$values,
      sy$values,
      lambdas
    ))]],
    loo = lambdas[[ridge_reference_one_se(
      ridge_reference_loo_mse(sf$values, sy$values, lambdas),
      lambdas
    )]],
    blocked = lambdas[[ridge_reference_one_se(
      ridge_reference_block_mse(
        features,
        responses,
        lambdas,
        unname(split(seq_len(n), blocks))
      ),
      lambdas
    )]]
  )

  for (selector in names(expected)) {
    model <- feature_rsa_model(
      dset$dataset,
      design,
      method = "ridge",
      lambda_selection = selector,
      lambda = lambdas
    )
    fit <- train_model(
      model,
      responses,
      features,
      indices = seq_len(ncol(responses)),
      observation_indices = seq_len(n)
    )

    expect_false(isTRUE(fit$error), info = selector)
    expect_equal(fit$lambda, unname(expected[[selector]]), info = selector)
    expect_identical(fit$lambda_selection, selector)
  }
})

test_that("zero-penalty ridge is stable in the p-greater-than-n limit", {
  set.seed(4509)
  predictors <- matrix(rnorm(12 * 25), 12, 25)
  predictors[, 25] <- predictors[, 1] + predictors[, 2]
  responses <- matrix(rnorm(12 * 7), 12, 7)

  fit <- rMVPA:::.feature_rsa_ridge_fit_kernel(
    predictors,
    responses,
    lambda = 0
  )
  predicted <- rMVPA:::.feature_rsa_ridge_predict_kernel(fit, predictors)

  expect_true(all(is.finite(fit$coefficients)))
  expect_equal(fit$rank, 11L)
  expect_equal(fit$effective_df, 12, tolerance = 1e-10)
  expect_equal(predicted, responses, tolerance = 2e-10)
})

test_that("standardized ridge is equivariant to column-wise affine units", {
  set.seed(4510)
  n <- 28L
  dset <- gen_sample_dataset(c(3, 3, 3), n, blocks = 4)
  features <- matrix(rnorm(n * 6), n, 6)
  responses <- matrix(rnorm(n * 8), n, 8)
  feature_scale <- c(-3, 0.25, 2, -0.5, 4, 1.5)
  feature_shift <- seq(-2, 3, length.out = ncol(features))
  response_scale <- rep(c(-2, 0.4), length.out = ncol(responses))
  response_shift <- seq(1, 4, length.out = ncol(responses))
  transformed_features <- sweep(
    sweep(features, 2L, feature_scale, "*"),
    2L,
    feature_shift,
    "+"
  )
  transformed_responses <- sweep(
    sweep(responses, 2L, response_scale, "*"),
    2L,
    response_shift,
    "+"
  )

  fit_once <- function(feature_data, response_data) {
    design <- feature_rsa_design(
      F = feature_data,
      labels = seq_len(n),
      block_var = dset$design$block_var
    )
    model <- feature_rsa_model(
      dset$dataset,
      design,
      method = "ridge",
      lambda_selection = "fixed",
      lambda = 0.12
    )
    fit <- train_model(
      model,
      response_data,
      feature_data,
      indices = seq_len(ncol(response_data)),
      observation_indices = seq_len(n)
    )
    predict_model(model, fit, feature_data)
  }

  original <- fit_once(features, responses)
  transformed <- fit_once(transformed_features, transformed_responses)
  expected <- sweep(
    sweep(original, 2L, response_scale, "*"),
    2L,
    response_shift,
    "+"
  )

  expect_equal(transformed, expected, tolerance = 3e-10)
})

test_that("fixed ridge agrees with glmnet's alpha-zero estimator", {
  skip_if_not_installed("glmnet")
  set.seed(4511)
  predictors <- matrix(rnorm(40 * 7), 40, 7)
  responses <- matrix(rnorm(40 * 10), 40, 10)
  lambda <- 0.08

  ridge_fit <- rMVPA:::.feature_rsa_ridge_fit_kernel(
    predictors,
    responses,
    lambda
  )
  expected <- rMVPA:::.feature_rsa_ridge_predict_kernel(
    ridge_fit,
    predictors
  )
  glmnet_fit <- glmnet::glmnet(
    x = predictors,
    y = responses,
    family = "mgaussian",
    alpha = 0,
    lambda = lambda,
    standardize = FALSE,
    intercept = TRUE,
    control = list(thresh = 1e-14, maxit = 1e6)
  )
  got <- drop(stats::predict(
    glmnet_fit,
    newx = predictors,
    s = lambda
  ))
  if (!is.matrix(got)) {
    got <- matrix(got, nrow = nrow(predictors))
  }

  expect_equal(unname(got), unname(expected), tolerance = 2e-7)
})

test_that("regional ridge reports ridge diagnostics instead of ncomp", {
  set.seed(4512)
  n <- 24L
  dset <- gen_sample_dataset(c(3, 3, 3), n, blocks = 3)
  features <- matrix(rnorm(n * 5), n, 5)
  design <- feature_rsa_design(
    F = features,
    labels = seq_len(n),
    block_var = dset$design$block_var
  )
  model <- feature_rsa_model(
    dset$dataset,
    design,
    method = "ridge",
    lambda_selection = "fixed",
    lambda = 0.2
  )
  region_mask <- neuroim2::NeuroVol(
    rep_len(1:2, length(dset$dataset$mask)),
    neuroim2::space(dset$dataset$mask)
  )

  result <- run_regional(
    model,
    region_mask,
    backend = "default",
    preflight = "off"
  )
  performance <- result$performance_table
  schema <- output_schema(model)

  expect_s3_class(result, "regional_mvpa_result")
  expect_false(any(result$fits$error))
  expect_true(all(c("median_lambda", "mean_effective_df") %in%
                    names(performance)))
  expect_false("ncomp" %in% names(performance))
  expect_equal(performance$median_lambda, rep(0.2, nrow(performance)))
  expect_true(all(is.finite(performance$mean_effective_df)))
  expect_setequal(
    tail(names(schema), 2L),
    c("median_lambda", "mean_effective_df")
  )
})
