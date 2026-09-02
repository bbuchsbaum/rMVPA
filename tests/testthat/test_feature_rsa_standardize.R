context("feature_rsa standardization")

.frs_fixture <- function(seed = 4101L, n = 60L, p_raw = 40L, k = 12L, v = 30L) {
  set.seed(seed)
  Z <- matrix(rnorm(n * 4L), n, 4L)
  raw <- Z %*% matrix(rnorm(4L * p_raw), 4L, p_raw) +
    0.5 * matrix(rnorm(n * p_raw), n, p_raw)
  colnames(raw) <- paste0("f", seq_len(p_raw))
  sv <- svd(base::scale(raw))
  US <- sv$u[, seq_len(k)] %*% diag(sv$d[seq_len(k)])
  U <- sv$u[, seq_len(k)]
  colnames(US) <- colnames(U) <- paste0("pc", seq_len(k))
  responses <- Z[, 1:2] %*% matrix(rnorm(2L * v), 2L, v) +
    matrix(rnorm(n * v), n, v)
  colnames(responses) <- paste0("v", seq_len(v))
  dset <- gen_sample_dataset(c(4L, 4L, 4L), nobs = n, blocks = 3L)
  list(raw = raw, US = US, U = U, responses = responses, dset = dset,
       n = n, v = v)
}

.frs_design <- function(fx, F, max_comps = 6L) {
  feature_rsa_design(
    F = F, labels = paste0("t", seq_len(fx$n)), max_comps = max_comps,
    block_var = fx$dset$design$block_var
  )
}

test_that(".standardize supports center-only standardization", {
  set.seed(4102L)
  X <- matrix(rnorm(80), 20, 4) %*% diag(c(1, 2, 3, 4))
  colnames(X) <- letters[1:4]
  sc <- rMVPA:::.standardize(X)
  ce <- rMVPA:::.standardize(X, "center")
  expect_equal(unname(apply(sc$X_sc, 2, sd)), rep(1, 4))
  expect_equal(unname(colMeans(ce$X_sc)), rep(0, 4))
  expect_equal(unname(apply(ce$X_sc, 2, sd)), unname(apply(X, 2, sd)))
  expect_equal(ce$sd, setNames(rep(1, 4), letters[1:4]))
  expect_equal(ce$mean, sc$mean)
  expect_error(rMVPA:::.standardize(X, "whiten"), "should be one of")
})

test_that("spectrum diagnostic flags exactly orthogonal inputs only", {
  fx <- .frs_fixture()
  diagnose <- rMVPA:::.feature_rsa_spectrum_diagnostic
  whitened <- diagnose(fx$U)
  scores <- diagnose(fx$US)
  structured <- diagnose(fx$raw)
  expect_true(whitened$flat)
  expect_true(scores$flat)
  expect_equal(scores$condition_number, 1, tolerance = 1e-8)
  expect_equal(whitened$condition_number, 1, tolerance = 1e-8)
  expect_false(structured$flat)
  expect_gt(structured$condition_number, 10)
  expect_identical(c(structured$n_rows, structured$n_cols), c(60L, 40L))

  # Chunked accumulation and column thinning do not change the verdict.
  chunked <- diagnose(fx$raw, chunk_rows = 7L)
  expect_equal(chunked$condition_number, structured$condition_number,
               tolerance = 1e-10)
  thinned <- diagnose(fx$US, max_cols = 5L)
  expect_identical(thinned$n_cols, 5L)
  expect_true(thinned$flat)

  # Noise is not flat, even at large n with few columns; orthogonal
  # equal-variance designs are (component ordering really is arbitrary).
  set.seed(4103L)
  noise <- diagnose(sweep(matrix(rnorm(5000 * 5), 5000, 5), 2L, 1:5, "*"))
  expect_false(noise$flat)
  expect_gt(noise$condition_number, 1.02)
  g <- factor(rep(1:8, 8))
  contrasts(g) <- contr.helmert(8)
  helmert <- diagnose(model.matrix(~ g)[, -1])
  expect_true(helmert$flat)

  # Documented limitation: scores computed on a superset of rows are missed.
  superset <- svd(base::scale(matrix(rnorm(300 * 40), 300, 40)))
  subset_scores <- (superset$u[, 1:12] %*% diag(superset$d[1:12]))[1:150, ]
  expect_false(diagnose(subset_scores)$flat)

  expect_null(diagnose(matrix(rnorm(10), 10, 1)))
  expect_null(diagnose(matrix(1, 10, 3)))
  expect_null(diagnose(cbind(fx$U, NA_real_)))
})

test_that("feature_rsa_model resolves feature scaling and warns on a flat spectrum", {
  fx <- .frs_fixture()
  expect_warning(
    model <- feature_rsa_model(
      fx$dset$dataset, .frs_design(fx, fx$US), method = "pca",
      ncomp_selection = "max"
    ),
    "spectrum of F is flat"
  )
  expect_identical(model$feature_standardize, "scale")
  expect_true(model$feature_spectrum$flat)
  expect_warning(
    centered <- feature_rsa_model(
      fx$dset$dataset, .frs_design(fx, fx$US), method = "pca",
      ncomp_selection = "max", feature_standardize = "center"
    ),
    NA
  )
  expect_identical(centered$feature_standardize, "center")
  expect_null(centered$feature_spectrum)
  expect_warning(
    raw_model <- feature_rsa_model(
      fx$dset$dataset, .frs_design(fx, fx$raw), method = "pca",
      ncomp_selection = "max"
    ),
    NA
  )
  expect_false(raw_model$feature_spectrum$flat)
  expect_warning(
    pls_model <- feature_rsa_model(
      fx$dset$dataset, .frs_design(fx, fx$US), method = "pls",
      ncomp_selection = "max"
    ),
    NA
  )
  expect_true(pls_model$feature_spectrum$flat)
  ridge_model <- feature_rsa_model(
    fx$dset$dataset, .frs_design(fx, fx$US), method = "ridge"
  )
  expect_identical(ridge_model$feature_standardize, "scale")
  expect_null(ridge_model$feature_spectrum)
  expect_error(feature_rsa_model(
    fx$dset$dataset, .frs_design(fx, fx$raw), method = "pca",
    feature_standardize = "whiten"
  ), "should be one of")
  expect_output(print(centered), "Feature standardization")

  # A design built from S holds PCA-like scores and defaults to centering.
  S <- tcrossprod(base::scale(fx$raw))
  s_design <- feature_rsa_design(
    S = S, labels = paste0("t", seq_len(fx$n)), k = 6L,
    block_var = fx$dset$design$block_var
  )
  expect_warning(
    s_model <- feature_rsa_model(fx$dset$dataset, s_design, method = "pca",
                                 ncomp_selection = "max"),
    NA
  )
  expect_identical(s_model$feature_standardize, "center")
  expect_null(s_model$feature_spectrum)
  expect_warning(
    s_scaled <- feature_rsa_model(fx$dset$dataset, s_design, method = "pca",
                                  ncomp_selection = "max",
                                  feature_standardize = "scale"),
    "spectrum of F is flat"
  )
  expect_identical(s_scaled$feature_standardize, "scale")
})

test_that("center-only standardization preserves the variance profile of F", {
  fx <- .frs_fixture()
  train <- 1:45
  test <- 46:60
  fit_predict <- function(F, standardize, ncomp = 3L) {
    model <- suppressWarnings(feature_rsa_model(
      fx$dset$dataset, .frs_design(fx, F, max_comps = ncomp), method = "pca",
      ncomp_selection = "max", feature_standardize = standardize
    ))
    fit <- train_model(model, fx$responses[train, ], F[train, ],
                       indices = seq_len(fx$v))
    list(fit = fit, pred = predict_model(model, fit, F[test, ]))
  }
  scaled_US <- fit_predict(fx$US, "scale")
  scaled_U <- fit_predict(fx$U, "scale")
  # Column scaling makes U %*% S and U the same input (the issue's demonstration).
  expect_equal(scaled_US$pred, scaled_U$pred, tolerance = 1e-10)

  centered_US <- fit_predict(fx$US, "center")
  centered_U <- fit_predict(fx$U, "center")
  expect_false(isTRUE(all.equal(centered_US$pred, centered_U$pred,
                                tolerance = 1e-6)))
  expect_identical(centered_US$fit$feature_standardize, "center")
  expect_equal(unname(centered_US$fit$f_sd), rep(1, ncol(fx$US)))

  # Oracle: covariance-PCR on the centered training features.
  pcr_oracle <- function(Ftr, Ytr, Fte, ncomp) {
    mu <- colMeans(Ftr)
    Zc <- sweep(Ftr, 2L, mu)
    V <- svd(Zc)$v[, seq_len(ncomp), drop = FALSE]
    scores <- Zc %*% V
    ymu <- colMeans(Ytr)
    B <- solve(crossprod(scores), crossprod(scores, sweep(Ytr, 2L, ymu)))
    sweep((sweep(Fte, 2L, mu) %*% V) %*% B, 2L, ymu, "+")
  }
  expect_equal(
    unname(centered_US$pred),
    unname(pcr_oracle(fx$US[train, ], fx$responses[train, ], fx$US[test, ], 3L)),
    tolerance = 1e-8
  )
})

test_that("center-only standardization threads through every tuning path", {
  fx <- .frs_fixture()
  cv <- blocked_cross_validation(fx$dset$design$block_var)
  ids <- seq_len(fx$v)
  fit_center <- function(method, ...) {
    model <- feature_rsa_model(
      fx$dset$dataset, .frs_design(fx, fx$US), method = method,
      feature_standardize = "center", crossval = cv, ...
    )
    fit <- train_model(model, fx$responses, fx$US, indices = ids)
    expect_null(fit$error)
    expect_identical(fit$feature_standardize, "center")
    list(model = model, fit = fit)
  }

  for (selection in c("max", "loo", "blocked", "pve")) {
    got <- fit_center("pca", ncomp_selection = selection)
    expect_equal(unname(got$fit$f_sd), rep(1, ncol(fx$US)), info = selection)
    expect_gte(got$fit$ncomp, 1L)
    expect_identical(dim(predict_model(got$model, got$fit, fx$US)),
                     c(fx$n, fx$v))
  }
  for (selection in c("gcv", "loo", "blocked")) {
    got <- fit_center("ridge", lambda_selection = selection)
    expect_equal(unname(got$fit$f_sd), rep(1, ncol(fx$US)), info = selection)
    expect_identical(dim(predict_model(got$model, got$fit, fx$US)),
                     c(fx$n, fx$v))
  }

  fixed <- function(standardize) {
    model <- feature_rsa_model(
      fx$dset$dataset, .frs_design(fx, fx$US), method = "ridge",
      lambda_selection = "fixed", lambda = 1,
      feature_standardize = standardize, crossval = cv
    )
    predict_model(model, train_model(model, fx$responses, fx$US, indices = ids),
                  fx$US)
  }
  expect_false(isTRUE(all.equal(fixed("scale"), fixed("center"))))

  segments <- split(seq_len(fx$n), fx$dset$design$block_var)
  s_scale <- rMVPA:::.feature_rsa_cv_segment_mse(
    fx$US, fx$responses, ncomp = 6L, method = "pca", segments = segments,
    fold_standardize = TRUE
  )
  s_center <- rMVPA:::.feature_rsa_cv_segment_mse(
    fx$US, fx$responses, ncomp = 6L, method = "pca", segments = segments,
    fold_standardize = TRUE, feature_standardize = "center"
  )
  expect_false(isTRUE(all.equal(s_scale, s_center)))

  skip_if_not_installed("glmnet")
  got <- fit_center("glmnet", lambda = 0.1)
  expect_equal(unname(got$fit$glmnet_f_sd), rep(1, ncol(fx$US)))
  expect_null(got$model$feature_spectrum)
  expect_identical(dim(predict_model(got$model, got$fit, fx$US)),
                   c(fx$n, fx$v))
})
