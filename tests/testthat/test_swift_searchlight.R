testthat::skip_if_not_installed("neuroim2")

build_swift_mspec <- function(D = c(4, 4, 4), nobs = 54, blocks = 6) {
  ds <- gen_sample_dataset(D = D, nobs = nobs, nlevels = 3, blocks = blocks)

  y <- factor(rep(letters[1:3], length.out = nobs), levels = letters[1:3])
  block <- rep(seq_len(blocks), each = ceiling(nobs / blocks))[seq_len(nobs)]

  des <- mvpa_design(
    data.frame(y = y, block = block),
    y_train = ~ y,
    block_var = ~ block
  )

  cval <- blocked_cross_validation(des$block_var)
  mspec <- mvpa_model(
    model = load_model("sda_notune"),
    dataset = ds$dataset,
    design = des,
    model_type = "classification",
    crossval = cval
  )

  list(mspec = mspec, dataset = ds)
}

build_swift_ineligible_mspec <- function(D = c(4, 4, 4), nobs = 40, blocks = 5) {
  ds <- gen_sample_dataset(D = D, nobs = nobs, nlevels = 2, blocks = blocks)

  y <- factor(rep(letters[1:2], length.out = nobs), levels = letters[1:2])
  block <- rep(seq_len(blocks), each = ceiling(nobs / blocks))[seq_len(nobs)]

  des <- mvpa_design(
    data.frame(y = y, block = block),
    y_train = ~ y,
    block_var = ~ block
  )

  cval <- blocked_cross_validation(des$block_var)
  mvpa_model(
    model = load_model("sda_notune"),
    dataset = ds$dataset,
    design = des,
    model_type = "classification",
    crossval = cval
  )
}

with_swift_options <- function(code) {
  keys <- c(
    "rMVPA.searchlight_mode",
    "rMVPA.searchlight_profile",
    "rMVPA.swift_searchlight",
    "rMVPA.swift_whitening",
    "rMVPA.warn_legacy_options"
  )
  old <- options()[keys]
  on.exit(options(old), add = TRUE)
  options(
    rMVPA.searchlight_mode = NULL,
    rMVPA.searchlight_profile = NULL,
    rMVPA.warn_legacy_options = FALSE
  )
  force(code)
}

test_that("swift gate follows mode defaults", {
  with_swift_options({
    options(rMVPA.searchlight_mode = "fast")
    options(rMVPA.swift_searchlight = NULL)
    expect_true(rMVPA:::.swift_searchlight_enabled())

    options(rMVPA.searchlight_mode = "legacy")
    options(rMVPA.swift_searchlight = NULL)
    expect_false(rMVPA:::.swift_searchlight_enabled())
  })
})

test_that("swift multiclass evidence matches direct definition", {
  set.seed(4001)
  classes <- c("a", "b", "c")
  y_train <- factor(rep(classes, each = 4), levels = classes)
  y_test <- factor(rep(classes, each = 3), levels = classes)

  x_train <- matrix(rnorm(length(y_train) * 7), nrow = length(y_train), ncol = 7)
  x_test <- matrix(rnorm(length(y_test) * 7), nrow = length(y_test), ncol = 7)

  got <- rMVPA:::.swift_multiclass_evidence(
    x_train = x_train,
    x_test = x_test,
    y_train = y_train,
    y_test = y_test,
    classes = classes
  )

  mu_tr <- sapply(classes, function(cl) colMeans(x_train[y_train == cl, , drop = FALSE]))
  mu_te <- sapply(classes, function(cl) colMeans(x_test[y_test == cl, , drop = FALSE]))
  mu_tr <- t(mu_tr)
  mu_te <- t(mu_te)

  n_tr <- as.numeric(table(y_train)[classes])
  n_te <- as.numeric(table(y_test)[classes])

  d_tr <- sweep(mu_tr, 2, colMeans(x_train), "-")
  d_te <- sweep(mu_te, 2, colMeans(x_test), "-")
  d_tr <- sweep(d_tr, 1, sqrt(n_tr), "*")
  d_te <- sweep(d_te, 1, sqrt(n_te), "*")

  expected <- colSums(d_tr * d_te)
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("swift searchlight returns multiclass metrics and SWIFT_Info map when enabled", {
  skip_if_not_installed("neuroim2")

  with_swift_options({
    options(rMVPA.searchlight_mode = "fast")
    options(rMVPA.swift_searchlight = TRUE)
    options(rMVPA.swift_whitening = "zscore")

    set.seed(4002)
    built <- build_swift_mspec(D = c(4, 4, 4), nobs = 54, blocks = 6)
    res <- run_searchlight(
      built$mspec,
      radius = 2,
      method = "standard",
      backend = "default"
    )

    expect_s3_class(res, "searchlight_result")
    expect_true("Accuracy" %in% names(res$results))
    expect_true("AUC" %in% names(res$results))
    expect_true("SWIFT_Info" %in% names(res$results))
    expect_identical(attr(res, "searchlight_engine"), "swift")

    acc <- as.numeric(neuroim2::values(res$results$Accuracy))
    auc <- as.numeric(neuroim2::values(res$results$AUC))
    vals <- as.numeric(neuroim2::values(res$results$SWIFT_Info))
    expect_true(any(is.finite(acc)))
    expect_true(any(is.finite(auc)))
    expect_true(any(is.finite(vals)))
    expect_true(all(acc[is.finite(acc)] >= 0 & acc[is.finite(acc)] <= 1))
    expect_true(all(auc[is.finite(auc)] >= -1 & auc[is.finite(auc)] <= 1))
  })
})

test_that("swift-disabled path falls back to legacy iterator", {
  skip_if_not_installed("neuroim2")

  with_swift_options({
    options(rMVPA.searchlight_mode = "fast")
    options(rMVPA.swift_searchlight = FALSE)

    set.seed(4003)
    built <- build_swift_mspec(D = c(4, 4, 4), nobs = 54, blocks = 6)
    res <- run_searchlight(
      built$mspec,
      radius = 2,
      method = "standard",
      backend = "default"
    )

    expect_s3_class(res, "searchlight_result")
    expect_false("SWIFT_Info" %in% names(res$results))
    expect_identical(attr(res, "searchlight_engine"), "legacy")
  })
})

test_that("engine override can force legacy even when swift is enabled", {
  skip_if_not_installed("neuroim2")

  with_swift_options({
    options(rMVPA.searchlight_mode = "fast")
    options(rMVPA.swift_searchlight = TRUE)

    set.seed(4004)
    built <- build_swift_mspec(D = c(4, 4, 4), nobs = 54, blocks = 6)
    res <- run_searchlight(
      built$mspec,
      radius = 2,
      method = "standard",
      backend = "default",
      engine = "legacy"
    )

    expect_s3_class(res, "searchlight_result")
    expect_false("SWIFT_Info" %in% names(res$results))
    expect_identical(attr(res, "searchlight_engine"), "legacy")
  })
})

test_that("swift engine runs randomized searchlight without fallback", {
  skip_if_not_installed("neuroim2")

  with_swift_options({
    options(rMVPA.searchlight_mode = "fast")
    options(rMVPA.swift_searchlight = TRUE)

    set.seed(4007)
    built <- build_swift_mspec(D = c(4, 4, 4), nobs = 54, blocks = 6)
    res <- run_searchlight(
      built$mspec,
      radius = 2,
      method = "randomized",
      niter = 1,
      backend = "default",
      engine = "swift"
    )

    expect_s3_class(res, "searchlight_result")
    expect_identical(attr(res, "searchlight_engine"), "swift")
    expect_true("SWIFT_Info" %in% names(res$results))
  })
})

test_that("swift engine runs resampled searchlight without fallback", {
  skip_if_not_installed("neuroim2")

  with_swift_options({
    options(rMVPA.searchlight_mode = "fast")
    options(rMVPA.swift_searchlight = TRUE)

    set.seed(4008)
    built <- build_swift_mspec(D = c(4, 4, 4), nobs = 54, blocks = 6)
    res <- run_searchlight(
      built$mspec,
      radius = 2,
      method = "resampled",
      niter = 12,
      backend = "default",
      engine = "swift"
    )

    expect_s3_class(res, "searchlight_result")
    expect_identical(attr(res, "searchlight_engine"), "swift")
    expect_true("SWIFT_Info" %in% names(res$results))
  })
})

test_that("explicit swift request errors when engine is ineligible", {
  skip_if_not_installed("neuroim2")

  with_swift_options({
    options(rMVPA.searchlight_mode = "fast")
    options(rMVPA.swift_searchlight = TRUE)

    set.seed(4009)
    mspec <- build_swift_ineligible_mspec(D = c(4, 4, 4), nobs = 40, blocks = 5)
    expect_error(
      run_searchlight(
        mspec,
        radius = 2,
        method = "randomized",
        niter = 1,
        backend = "default",
        engine = "swift"
      ),
      regexp = "not eligible"
    )
  })
})

test_that("run_searchlight_base rejects non-legacy engine argument", {
  skip_if_not_installed("neuroim2")

  with_swift_options({
    set.seed(4010)
    built <- build_swift_mspec(D = c(4, 4, 4), nobs = 54, blocks = 6)
    expect_error(
      run_searchlight_base(
        built$mspec,
        radius = 2,
        method = "randomized",
        niter = 1,
        engine = "swift"
      ),
      regexp = "legacy iterator only"
    )
  })
})

test_that("invalid engine override errors clearly", {
  skip_if_not_installed("neuroim2")

  with_swift_options({
    options(rMVPA.searchlight_mode = "fast")
    options(rMVPA.swift_searchlight = TRUE)

    set.seed(4005)
    built <- build_swift_mspec(D = c(4, 4, 4), nobs = 54, blocks = 6)
    expect_error(
      run_searchlight(
        built$mspec,
        radius = 2,
        method = "standard",
        engine = "not_an_engine"
      ),
      regexp = "arg"
    )
  })
})

test_that("swift metric kernel matches generic multiclass performance", {
  set.seed(4006)
  classes <- c("a", "b", "c")
  n <- 30
  b <- 4
  k <- length(classes)

  observed <- factor(sample(classes, n, replace = TRUE), levels = classes)
  raw <- array(rexp(n * b * k), dim = c(n, b, k))
  denom <- apply(raw, c(1, 2), sum)
  probs <- raw
  for (kk in seq_len(k)) {
    probs[, , kk] <- raw[, , kk] / denom
  }
  dimnames(probs)[[3]] <- classes

  test_idx <- seq_len(n)
  split_groups <- list(g1 = seq(1, n, by = 2), g2 = seq(2, n, by = 2))

  got <- rMVPA:::.swift_metric_with_splits_matrix(
    observed = observed,
    probs = probs,
    test_idx = test_idx,
    classes = classes,
    split_groups = split_groups,
    class_metrics = TRUE
  )

  expected <- vapply(seq_len(b), function(j) {
    p <- probs[, j, , drop = TRUE]
    pred <- factor(classes[max.col(p)], levels = classes)
    cres <- multiway_classification_result(
      observed = observed,
      predicted = pred,
      probs = p,
      testind = test_idx
    )
    performance(cres, split_list = split_groups, class_metrics = TRUE)
  }, numeric(ncol(got)))
  expected <- t(expected)
  colnames(expected) <- names(performance(
    multiway_classification_result(
      observed = observed,
      predicted = factor(classes[max.col(probs[, 1, , drop = TRUE])], levels = classes),
      probs = probs[, 1, , drop = TRUE],
      testind = test_idx
    ),
    split_list = split_groups,
    class_metrics = TRUE
  ))

  expect_equal(colnames(got), colnames(expected))
  expect_equal(got, expected, tolerance = 1e-10)
})
