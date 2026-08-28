library(testthat)
library(rMVPA)

reference_feature_rsa_standardize <- function(x) {
  x <- as.matrix(x)
  center <- colMeans(x)
  scale <- apply(x, 2L, stats::sd)
  list(
    values = base::scale(x, center = center, scale = scale),
    center = center,
    scale = scale
  )
}

reference_feature_rsa_fit <- function(features, responses, ncomp, method,
                                      validation = "none", segments = NULL) {
  fit_fun <- if (identical(method, "pls")) pls::plsr else pls::pcr
  args <- list(
    responses ~ features,
    ncomp = ncomp,
    scale = FALSE,
    validation = validation
  )
  if (!is.null(segments)) {
    args$segments <- segments
  }
  do.call(fit_fun, args)
}

reference_feature_rsa_segment_mse <- function(model, responses) {
  pred <- model$validation$pred
  segments <- model$validation$segments
  ncomp <- dim(pred)[3L]
  out <- matrix(NA_real_, nrow = length(segments), ncol = ncomp)

  for (s in seq_along(segments)) {
    idx <- segments[[s]]
    for (k in seq_len(ncomp)) {
      pred_k <- pred[idx, , k, drop = FALSE][, , 1L]
      out[s, k] <- mean((responses[idx, , drop = FALSE] - pred_k)^2)
    }
  }
  out
}

reference_feature_rsa_onesigma <- function(segment_mse) {
  mean_mse <- colMeans(segment_mse)
  min_idx <- which.min(mean_mse)
  threshold <- mean_mse[[min_idx]] +
    stats::sd(segment_mse[, min_idx]) / sqrt(nrow(segment_mse))
  as.integer(which(mean_mse <= threshold)[1L])
}

reference_feature_rsa_nested_segment_mse <- function(features,
                                                      responses,
                                                      ncomp,
                                                      method,
                                                      segments) {
  ncomp_cv <- min(ncomp, nrow(features) - max(lengths(segments)) - 1L)
  out <- matrix(NA_real_, nrow = length(segments), ncol = ncomp_cv)
  all_rows <- seq_len(nrow(features))

  for (s in seq_along(segments)) {
    test_idx <- segments[[s]]
    train_idx <- all_rows[-test_idx]
    sf <- reference_feature_rsa_standardize(
      features[train_idx, , drop = FALSE]
    )
    sx <- reference_feature_rsa_standardize(
      responses[train_idx, , drop = FALSE]
    )
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
    fit <- reference_feature_rsa_fit(
      sf$values,
      sx$values,
      ncomp = ncomp_cv,
      method = method
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

test_that("feature RSA matrix kernels reproduce pls and pcr predictions", {
  set.seed(4401)
  features <- matrix(rnorm(48 * 9, mean = 2, sd = 1.5), 48, 9)
  responses <- matrix(rnorm(48 * 13, mean = -1, sd = 2), 48, 13)
  new_features <- matrix(rnorm(7 * 9, mean = 2, sd = 1.5), 7, 9)
  colnames(features) <- colnames(new_features) <- paste0("feature_", 1:9)
  colnames(responses) <- paste0("voxel_", 1:13)

  sf <- reference_feature_rsa_standardize(features)
  sx <- reference_feature_rsa_standardize(responses)
  new_scaled <- base::scale(
    new_features,
    center = sf$center,
    scale = sf$scale
  )

  for (method in c("pls", "pca")) {
    reference <- reference_feature_rsa_fit(
      sf$values, sx$values, ncomp = 6L, method = method
    )
    expected <- drop(stats::predict(reference, newdata = new_scaled, ncomp = 4L))
    if (!is.matrix(expected)) {
      expected <- matrix(expected, nrow = nrow(new_scaled))
    }

    fit <- rMVPA:::.feature_rsa_fit_kernel(
      sf$values, sx$values, ncomp = 6L, method = method
    )
    got <- rMVPA:::.feature_rsa_predict_kernel(fit, new_scaled, ncomp = 4L)

    expect_s3_class(fit, "feature_rsa_kernel_fit")
    expect_identical(colnames(got), colnames(expected))
    expect_equal(unname(got), unname(expected), tolerance = 2e-12)
  }
})

test_that("feature RSA matrix kernels respect supported pls algorithms", {
  set.seed(4413)
  features <- matrix(rnorm(40 * 7), 40, 7)
  responses <- matrix(rnorm(40 * 9), 40, 9)
  new_features <- matrix(rnorm(5 * 7), 5, 7)
  old_options <- pls::pls.options()
  on.exit(do.call(pls::pls.options, old_options), add = TRUE)

  for (algorithm in c(
    "kernelpls", "widekernelpls", "simpls", "oscorespls", "nipalspls"
  )) {
    pls::pls.options(plsralg = algorithm)
    reference <- suppressWarnings(pls::plsr(
      responses ~ features,
      ncomp = 4L,
      validation = "none",
      scale = FALSE
    ))
    fit <- suppressWarnings(rMVPA:::.feature_rsa_fit_kernel(
      features, responses, ncomp = 4L, method = "pls"
    ))

    expect_identical(fit$algorithm, algorithm)
    expect_equal(
      unname(rMVPA:::.feature_rsa_predict_kernel(
        fit, new_features, ncomp = 3L
      )),
      unname(drop(stats::predict(
        reference, newdata = new_features, ncomp = 3L
      ))),
      tolerance = 2e-10,
      info = algorithm
    )
  }
})

test_that("streaming feature RSA LOO errors reproduce pls validation", {
  set.seed(4402)
  features <- scale(matrix(rnorm(28 * 7), 28, 7))
  responses <- scale(matrix(rnorm(28 * 11), 28, 11))

  for (method in c("pls", "pca")) {
    reference <- reference_feature_rsa_fit(
      features,
      responses,
      ncomp = 5L,
      method = method,
      validation = "LOO"
    )
    expected_mse <- reference_feature_rsa_segment_mse(reference, responses)
    got_mse <- rMVPA:::.feature_rsa_cv_segment_mse(
      features,
      responses,
      ncomp = 5L,
      method = method,
      segments = reference$validation$segments
    )

    expect_equal(got_mse, expected_mse, tolerance = 2e-12)
    expect_identical(
      rMVPA:::.feature_rsa_select_from_segment_mse(got_mse),
      reference_feature_rsa_onesigma(expected_mse)
    )
  }
})

test_that("blocked feature RSA selection reproduces explicit pls CV segments", {
  set.seed(4410)
  n <- 48L
  blocks <- rep(seq_len(4L), each = n / 4L)
  features <- matrix(rnorm(n * 8), n, 8)
  responses <- matrix(rnorm(n * 12), n, 12)
  dset <- gen_sample_dataset(c(3, 3, 3), n, blocks = 4)
  design <- feature_rsa_design(
    F = features,
    labels = seq_len(n),
    max_comps = 6L,
    block_var = blocks
  )
  segments <- split(seq_len(n), blocks)
  sf <- reference_feature_rsa_standardize(features)
  sx <- reference_feature_rsa_standardize(responses)

  for (method in c("pls", "pca")) {
    model <- feature_rsa_model(
      dset$dataset,
      design,
      method = method,
      ncomp_selection = "blocked",
      crossval = blocked_cross_validation(blocks)
    )
    fit <- train_model(
      model,
      responses,
      features,
      indices = seq_len(ncol(responses))
    )
    got <- predict_model(model, fit, features)

    expected_mse <- reference_feature_rsa_nested_segment_mse(
      features, responses, ncomp = 6L, method = method, segments = segments
    )
    expected_ncomp <- reference_feature_rsa_onesigma(
      expected_mse
    )
    reference <- reference_feature_rsa_fit(
      sf$values, sx$values, ncomp = 6L, method = method
    )
    expected_scaled <- drop(stats::predict(
      reference,
      newdata = sf$values,
      ncomp = expected_ncomp
    ))
    expected <- sweep(
      sweep(expected_scaled, 2L, sx$scale, "*"),
      2L,
      sx$center,
      "+"
    )

    expect_identical(fit$ncomp, expected_ncomp)
    expect_equal(unname(got), unname(expected), tolerance = 3e-11)
  }
})

test_that("blocked feature RSA segment errors use inner-training scaling", {
  set.seed(4415)
  n <- 40L
  blocks <- rep(seq_len(4L), each = n / 4L)
  features <- matrix(rnorm(n * 7), n, 7)
  responses <- matrix(rnorm(n * 10), n, 10)
  segments <- unname(split(seq_len(n), blocks))

  for (method in c("pls", "pca")) {
    expected <- reference_feature_rsa_nested_segment_mse(
      features, responses, ncomp = 5L, method = method, segments = segments
    )
    got <- rMVPA:::.feature_rsa_cv_segment_mse(
      features,
      responses,
      ncomp = 5L,
      method = method,
      segments = segments,
      fold_standardize = TRUE
    )

    expect_equal(got, expected, tolerance = 3e-11)
  }
})

test_that("blocked feature RSA selection is invariant to block labels", {
  blocks <- rep(c("run-a", "run-b", "run-c"), each = 5L)
  relabelled <- c(a = "z", b = "x", c = "y")[sub("run-", "", blocks)]

  original <- rMVPA:::.feature_rsa_block_segments(blocks, n_rows = length(blocks))
  renamed <- rMVPA:::.feature_rsa_block_segments(relabelled, n_rows = length(blocks))

  expect_equal(
    sort(vapply(original, paste, collapse = ",", character(1))),
    sort(vapply(renamed, paste, collapse = ",", character(1)))
  )
})

test_that("blocked feature RSA selection fails closed without repeated blocks", {
  set.seed(4411)
  n <- 20L
  features <- matrix(rnorm(n * 5), n, 5)
  responses <- matrix(rnorm(n * 7), n, 7)
  dset <- gen_sample_dataset(c(3, 3, 3), n, blocks = 2)
  design <- feature_rsa_design(
    F = features,
    labels = seq_len(n),
    max_comps = 4L,
    block_var = rep(1L, n)
  )
  model <- feature_rsa_model(
    dset$dataset,
    design,
    method = "pls",
    ncomp_selection = "blocked",
    crossval = blocked_cross_validation(dset$design$block_var)
  )

  fit <- train_model(
    model,
    responses,
    features,
    indices = seq_len(ncol(responses))
  )

  expect_match(fit$error, "at least two non-empty blocks")
})

test_that("blocked feature RSA selection requires block metadata", {
  dset <- gen_sample_dataset(c(3, 3, 3), 20L, blocks = 2L)
  design <- feature_rsa_design(
    F = matrix(rnorm(20 * 5), 20, 5),
    labels = seq_len(20),
    max_comps = 4L
  )

  expect_error(
    feature_rsa_model(
      dset$dataset,
      design,
      method = "pls",
      ncomp_selection = "blocked",
      crossval = blocked_cross_validation(dset$design$block_var)
    ),
    "requires design\\$block_var"
  )
})

test_that("component selection settings are ignored by glmnet as documented", {
  skip_if_not_installed("glmnet")
  dset <- gen_sample_dataset(c(3, 3, 3), 20L, blocks = 2L)
  design <- feature_rsa_design(
    F = matrix(rnorm(20 * 5), 20, 5),
    labels = seq_len(20),
    max_comps = 4L
  )

  expect_no_error(feature_rsa_model(
    dset$dataset,
    design,
    method = "glmnet",
    ncomp_selection = "blocked",
    pve_threshold = 2,
    crossval = blocked_cross_validation(dset$design$block_var)
  ))
})

test_that("feature RSA reuses the full-data cv.glmnet path", {
  skip_if_not_installed("glmnet")
  set.seed(4416)
  n <- 36L
  features <- matrix(rnorm(n * 7), n, 7)
  responses <- matrix(rnorm(n * 9), n, 9)
  dset <- gen_sample_dataset(c(3, 3, 3), n, blocks = 3L)
  design <- feature_rsa_design(
    F = features,
    labels = seq_len(n),
    max_comps = 5L
  )
  model <- feature_rsa_model(
    dset$dataset,
    design,
    method = "glmnet",
    cv_glmnet = TRUE,
    crossval = blocked_cross_validation(dset$design$block_var)
  )

  set.seed(9416)
  fit <- train_model(
    model, responses, features, indices = seq_len(ncol(responses))
  )
  got <- predict_model(model, fit, features[1:6, , drop = FALSE])

  scaled_features <- base::scale(
    features[1:6, , drop = FALSE],
    center = fit$glmnet_f_mean,
    scale = fit$glmnet_f_sd
  )
  expected_scaled <- drop(stats::predict(
    fit$cv_results$glmnet.fit,
    newx = scaled_features,
    s = fit$cv_results$lambda.min
  ))
  expected <- sweep(
    sweep(expected_scaled, 2L, fit$glmnet_x_sd, "*"),
    2L,
    fit$glmnet_x_mean,
    "+"
  )

  expect_true(fit$cv_glmnet)
  expect_identical(fit$trained_model, fit$cv_results$glmnet.fit)
  expect_equal(fit$lambda_used, fit$cv_results$lambda.min)
  expect_true(all(is.finite(got)))
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("blocked feature RSA selection receives outer-fold observation indices", {
  set.seed(4412)
  n <- 36L
  dset <- gen_sample_dataset(c(3, 3, 3), n, blocks = 4)
  features <- matrix(rnorm(n * 6), n, 6)
  design <- feature_rsa_design(
    F = features,
    labels = seq_len(n),
    max_comps = 4L,
    block_var = dset$design$block_var
  )
  model <- feature_rsa_model(
    dset$dataset,
    design,
    method = "pls",
    ncomp_selection = "blocked",
    crossval = blocked_cross_validation(dset$design$block_var)
  )
  region_mask <- neuroim2::NeuroVol(
    rep(1:2, length.out = length(dset$dataset$mask)),
    neuroim2::space(dset$dataset$mask)
  )

  result <- suppressWarnings(run_regional(
    model,
    region_mask,
    backend = "default",
    preflight = "off"
  ))

  expect_s3_class(result, "regional_mvpa_result")
  expect_equal(nrow(result$performance_table), 2L)
  expect_true(all(is.finite(result$performance_table$ncomp)))
})

test_that("blocked feature RSA selection has default and shard parity", {
  skip_if_not_installed("shard")
  set.seed(4414)
  n <- 36L
  dset <- gen_sample_dataset(c(3, 3, 3), n, blocks = 4)
  features <- matrix(rnorm(n * 6), n, 6)
  design <- feature_rsa_design(
    F = features,
    labels = seq_len(n),
    max_comps = 4L,
    block_var = dset$design$block_var
  )
  model <- feature_rsa_model(
    dset$dataset,
    design,
    method = "pls",
    ncomp_selection = "blocked",
    crossval = blocked_cross_validation(dset$design$block_var)
  )
  region_mask <- neuroim2::NeuroVol(
    rep(1:2, length.out = length(dset$dataset$mask)),
    neuroim2::space(dset$dataset$mask)
  )
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::sequential)

  default <- suppressWarnings(run_regional(
    model, region_mask, backend = "default", preflight = "off"
  ))
  shared <- suppressWarnings(run_regional(
    model, region_mask, backend = "shard", preflight = "off"
  ))

  expect_identical(
    names(shared$performance_table), names(default$performance_table)
  )
  expect_equal(
    shared$performance_table,
    default$performance_table,
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
})

test_that("feature RSA kernel explained variance reproduces pls::explvar", {
  set.seed(4403)
  features <- scale(matrix(rnorm(42 * 8), 42, 8))
  responses <- scale(matrix(rnorm(42 * 10), 42, 10))

  for (method in c("pls", "pca")) {
    reference <- reference_feature_rsa_fit(
      features, responses, ncomp = 6L, method = method
    )
    fit <- rMVPA:::.feature_rsa_fit_kernel(
      features,
      responses,
      ncomp = 6L,
      method = method,
      keep_explained_variance = TRUE
    )
    got <- 100 * fit$Xvar / fit$Xtotvar

    expect_equal(got, as.numeric(pls::explvar(reference)), tolerance = 2e-12)
  }
})

test_that("feature RSA kernel predictions are invariant to feature ordering", {
  set.seed(4404)
  features <- matrix(rnorm(36 * 8), 36, 8)
  responses <- matrix(rnorm(36 * 12), 36, 12)
  new_features <- matrix(rnorm(5 * 8), 5, 8)
  order <- c(8L, 3L, 1L, 6L, 2L, 7L, 5L, 4L)

  for (method in c("pls", "pca")) {
    base_fit <- rMVPA:::.feature_rsa_fit_kernel(
      features, responses, ncomp = 5L, method = method
    )
    perm_fit <- rMVPA:::.feature_rsa_fit_kernel(
      features[, order], responses, ncomp = 5L, method = method
    )
    base_pred <- rMVPA:::.feature_rsa_predict_kernel(
      base_fit, new_features, ncomp = 4L
    )
    perm_pred <- rMVPA:::.feature_rsa_predict_kernel(
      perm_fit, new_features[, order], ncomp = 4L
    )

    expect_equal(perm_pred, base_pred, tolerance = 2e-11)
  }
})

test_that("feature RSA fast row correlation matches stats::cor", {
  set.seed(4405)
  x <- matrix(rnorm(24 * 17), 24, 17)
  y <- matrix(rnorm(19 * 17), 19, 17)

  expect_equal(
    rMVPA:::.feature_rsa_row_cor(x),
    stats::cor(t(x), use = "pairwise.complete.obs"),
    tolerance = 2e-12
  )
  expect_equal(
    rMVPA:::.feature_rsa_row_cor(x, y),
    stats::cor(t(x), t(y), use = "pairwise.complete.obs"),
    tolerance = 2e-12
  )

  x_shifted <- x * 7 + 1e4
  y_shifted <- y * 0.25 - 500
  expect_equal(
    rMVPA:::.feature_rsa_row_cor(x_shifted, y_shifted),
    rMVPA:::.feature_rsa_row_cor(x, y),
    tolerance = 2e-11
  )
})

test_that("feature RSA row correlation preserves constant and missing-data semantics", {
  x <- rbind(
    c(1, 2, 3, 4, 5),
    c(7, 7, 7, 7, 7),
    c(2, NA, 4, 8, 16),
    c(-1, 0, 2, 5, 9)
  )

  expect_equal(
    suppressWarnings(rMVPA:::.feature_rsa_row_cor(x)),
    suppressWarnings(stats::cor(t(x), use = "pairwise.complete.obs")),
    tolerance = 2e-12
  )
})

test_that("feature RSA training uses compact kernels with legacy prediction parity", {
  set.seed(4406)
  dset <- gen_sample_dataset(c(3, 3, 3), 36, blocks = 3)
  features <- matrix(rnorm(36 * 7), 36, 7)
  responses <- matrix(rnorm(36 * 14), 36, 14)
  design <- feature_rsa_design(
    F = features,
    labels = seq_len(36),
    max_comps = 5L
  )

  sf <- reference_feature_rsa_standardize(features)
  sx <- reference_feature_rsa_standardize(responses)

  for (method in c("pls", "pca")) {
    model <- feature_rsa_model(
      dset$dataset,
      design,
      method = method,
      ncomp_selection = "loo",
      crossval = blocked_cross_validation(dset$design$block_var)
    )
    fit <- train_model(model, responses, features, indices = seq_len(ncol(responses)))
    got <- predict_model(model, fit, features)

    reference <- reference_feature_rsa_fit(
      sf$values,
      sx$values,
      ncomp = 5L,
      method = method,
      validation = "LOO"
    )
    expected_ncomp <- reference_feature_rsa_onesigma(
      reference_feature_rsa_segment_mse(reference, sx$values)
    )
    expected_scaled <- drop(stats::predict(
      reference,
      newdata = sf$values,
      ncomp = expected_ncomp
    ))
    expected <- sweep(
      sweep(expected_scaled, 2L, sx$scale, "*"),
      2L,
      sx$center,
      "+"
    )

    expect_s3_class(fit$trained_model, "feature_rsa_kernel_fit")
    expect_null(fit$trained_model$validation)
    expect_identical(fit$ncomp, expected_ncomp)
    expect_equal(unname(got), unname(expected), tolerance = 3e-11)
  }
})

test_that("feature RSA compact kernels handle collinear predictors", {
  set.seed(4407)
  base <- matrix(rnorm(32 * 4), 32, 4)
  features <- cbind(base, base[, 1L] + base[, 2L], base[, 3L])
  responses <- matrix(rnorm(32 * 9), 32, 9)

  for (method in c("pls", "pca")) {
    reference <- reference_feature_rsa_fit(
      features, responses, ncomp = 4L, method = method
    )
    fit <- rMVPA:::.feature_rsa_fit_kernel(
      features, responses, ncomp = 4L, method = method
    )
    expected <- drop(stats::predict(reference, newdata = features, ncomp = 3L))
    got <- rMVPA:::.feature_rsa_predict_kernel(fit, features, ncomp = 3L)

    expect_true(all(is.finite(got)))
    expect_equal(unname(got), unname(expected), tolerance = 5e-10)
  }
})

test_that("feature RSA worker specs omit unused quadratic design payload", {
  set.seed(4408)
  n <- 80L
  latent <- matrix(rnorm(n * 8), n, 8)
  similarity <- tcrossprod(scale(latent))
  design <- feature_rsa_design(
    S = similarity,
    labels = seq_len(n),
    k = 8L,
    max_comps = 5L,
    block_var = rep(seq_len(4L), length.out = n)
  )
  dset <- gen_sample_dataset(c(3, 3, 3), n, blocks = 4)
  model <- feature_rsa_model(
    dset$dataset,
    design,
    method = "pls",
    ncomp_selection = "max",
    crossval = blocked_cross_validation(design$block_var)
  )
  uncached_worker <- rMVPA:::as_worker_spec(model)
  model$.cv_fold_cache <- rMVPA:::.build_cv_fold_cache(model)

  worker <- rMVPA:::as_worker_spec(model)

  expect_null(worker$dataset)
  expect_null(worker$design$S)
  expect_null(worker$design$F)
  expect_null(worker$design$labels)
  expect_null(worker$design$cv_labels)
  expect_null(worker$design$targets)
  expect_equal(worker$.cv_fold_cache$y, model$design$targets)
  expect_null(worker$.cv_fold_cache$ytrain)
  expect_null(worker$.cv_fold_cache$ytest)
  expect_equal(uncached_worker$design$targets, model$design$targets)
  expect_equal(model$design$S, similarity)
  expect_lt(length(serialize(worker, NULL)), 0.5 * length(serialize(model, NULL)))
})

test_that("feature RSA fold formatting retains only merge inputs", {
  set.seed(4409)
  dset <- gen_sample_dataset(c(3, 3, 3), 24, blocks = 3)
  features <- matrix(rnorm(24 * 5), 24, 5)
  responses <- matrix(rnorm(24 * 8), 24, 8)
  design <- feature_rsa_design(F = features, labels = seq_len(24), max_comps = 4L)
  model <- feature_rsa_model(
    dset$dataset,
    design,
    method = "pca",
    ncomp_selection = "max",
    crossval = blocked_cross_validation(dset$design$block_var)
  )
  fit <- train_model(
    model,
    responses[1:16, , drop = FALSE],
    features[1:16, , drop = FALSE],
    indices = seq_len(ncol(responses))
  )
  out <- format_result(
    model,
    fit,
    context = list(
      test = responses[17:24, , drop = FALSE],
      ytest = features[17:24, , drop = FALSE],
      test_ind = 17:24
    )
  )

  expect_named(out$result[[1L]], "ncomp")
  expect_identical(out$result[[1L]]$ncomp, fit$ncomp)
  expect_null(out$performance[[1L]])
})
