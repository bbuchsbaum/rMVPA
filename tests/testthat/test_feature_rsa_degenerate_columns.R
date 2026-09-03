context("feature_rsa degenerate column handling")

.frdeg_fixture <- function(seed = 8501L, n = 96L, p = 10L, blocks = 4L) {
  set.seed(seed)
  dset <- gen_sample_dataset(c(4L, 4L, 3L), nobs = n, blocks = blocks)
  F <- matrix(rnorm(n * p), n, p)
  colnames(F) <- paste0("f", seq_len(p))
  region_mask <- neuroim2::NeuroVol(
    sample(1:2, size = length(dset$dataset$mask), replace = TRUE),
    neuroim2::space(dset$dataset$mask)
  )
  list(dset = dset, F = F, n = n, region_mask = region_mask,
       cv = blocked_cross_validation(dset$design$block_var))
}

.frdeg_run <- function(fx, F, standardize = "center", max_comps = 5L,
                       method = "pca") {
  des <- feature_rsa_design(
    F = F, labels = paste0("t", seq_len(fx$n)), max_comps = max_comps,
    block_var = fx$dset$design$block_var
  )
  model <- suppressWarnings(feature_rsa_model(
    fx$dset$dataset, des, method = method, crossval = fx$cv,
    ncomp_selection = "max", feature_standardize = standardize
  ))
  suppressMessages(run_regional(model, region_mask = fx$region_mask))
}

.frdeg_metrics <- function(res) {
  as.matrix(res$performance_table[, -1, drop = FALSE])
}

.feature_rsa_kernel_fit_uncapped <- function(F, X, ncomp) {
  centered <- scale(F, center = TRUE, scale = FALSE)
  fit <- pls::svdpc.fit(centered, scale(X, TRUE, TRUE), ncomp = ncomp,
                        center = TRUE, stripped = TRUE)
  fit$coefficients
}

test_that("column screen judges each column against its own magnitude", {
  assert_usable <- rMVPA:::.feature_rsa_assert_columns_usable
  set.seed(8502L)
  M <- matrix(rnorm(200), 50, 4)
  colnames(M) <- letters[1:4]

  expect_silent(assert_usable(M, "M"))
  # An absolute variance threshold fails here; a relative one must not.
  expect_silent(assert_usable(M * 1e-8, "M"))
  # Nor may a column be judged against the widest column in the matrix: F is
  # allowed to mix units, which is what "scale" is for.
  mixed <- M
  mixed[, 2] <- mixed[, 2] * 2e5
  expect_silent(assert_usable(mixed, "M"))
  # A large common offset is not degeneracy either: centering a column of
  # timestamps still leaves most of its significant digits.
  offset <- M
  offset[, 3] <- 1.7e12 + offset[, 3] * 1000
  expect_silent(assert_usable(offset, "M"))

  constant <- M
  constant[, 2] <- 7
  expect_error(assert_usable(constant, "Feature matrix F"),
               "Feature matrix F has 1 of 4 columns")
  expect_error(assert_usable(constant, "Feature matrix F"), "column 2 \\('b'\\)")
  expect_error(assert_usable(constant, "F", remedy = "Drop them."), "Drop them\\.")
  # Rescaling the whole matrix does not change the verdict either way.
  expect_error(assert_usable(constant * 1e-8, "F"), "1 of 4 columns")
  expect_silent(assert_usable(constant, "F", check_constant = FALSE))

  zeros <- M
  zeros[, 1] <- 0
  expect_error(assert_usable(zeros, "F"), "1 of 4 columns")

  # Constant up to rounding is still constant.
  rounding <- M
  rounding[, 4] <- 1e6 + rnorm(50, sd = 1e-9)
  expect_error(assert_usable(rounding, "F"), "column 4 \\('d'\\)")

  expect_error(assert_usable(matrix(1, 20, 3), "F"), "3 of 3 columns")

  missing <- M
  missing[3, 4] <- NA_real_
  expect_error(assert_usable(missing, "Brain data X"),
               "non-finite values \\(first: column 4")
  expect_error(assert_usable(missing, "X", check_constant = FALSE),
               "non-finite values")
})

test_that("component rank tracks numerical, not structural, dimension", {
  rank_of <- rMVPA:::.feature_rsa_component_rank
  set.seed(8503L)
  M <- scale(matrix(rnorm(60 * 8), 60, 8), TRUE, FALSE)
  expect_identical(rank_of(M), 8L)
  expect_identical(rank_of(M * 1e-8), 8L)

  constant <- M
  constant[, 3] <- 0
  expect_identical(rank_of(constant), 7L)

  collinear <- M
  collinear[, 5] <- collinear[, 1] + collinear[, 2]
  expect_identical(rank_of(collinear), 7L)

  # Merely correlated columns are still full rank (r about 0.9999).
  correlated <- M
  correlated[, 5] <- correlated[, 1] + rnorm(60, sd = 0.01)
  expect_identical(rank_of(correlated), 8L)
})

test_that("centered feature RSA accepts constant and small-unit F columns", {
  fx <- .frdeg_fixture()
  plain <- .frdeg_run(fx, fx$F)
  expect_true(all(is.finite(.frdeg_metrics(plain))))

  # Issue #85: neither a constant column nor a change of units is a reason to
  # refuse an analysis that only centers F.
  with_constant <- cbind(fx$F, const = 3)
  expect_true(all(is.finite(.frdeg_metrics(.frdeg_run(fx, with_constant)))))

  rescaled <- .frdeg_run(fx, fx$F * 1e-6)
  expect_true(all(is.finite(.frdeg_metrics(rescaled))))
  # Rescaling F is a no-op for centered PCR, so the results must match.
  expect_equal(.frdeg_metrics(rescaled), .frdeg_metrics(plain), tolerance = 1e-8)
})

test_that("scaled feature RSA refuses constant F columns, with detail", {
  fx <- .frdeg_fixture()
  with_constant <- cbind(fx$F, const = 3)

  scaled_units <- .frdeg_run(fx, fx$F * 1e-6, standardize = "scale")
  expect_true(all(is.finite(.frdeg_metrics(scaled_units))))

  design <- feature_rsa_design(
    F = with_constant, labels = paste0("t", seq_len(fx$n)), max_comps = 5L,
    block_var = fx$dset$design$block_var
  )
  methods <- c("pca", "pls", "ridge")
  if (requireNamespace("glmnet", quietly = TRUE)) methods <- c(methods, "glmnet")

  for (method in methods) {
    # A column that is constant over every row fails in every fold, so it is
    # named once at construction rather than 400 times inside the ROI loop.
    err <- tryCatch(
      suppressWarnings(feature_rsa_model(
        fx$dset$dataset, design, method = method, crossval = fx$cv,
        ncomp_selection = "max", feature_standardize = "scale"
      )),
      error = function(e) conditionMessage(e)
    )
    expect_type(err, "character")
    expect_match(err, "1 of 11 columns that are constant")
    expect_match(err, "column 11 \\('const'\\)")
    expect_true(grepl("feature_standardize = \"center\"", err, fixed = TRUE))

    expect_s3_class(
      suppressWarnings(feature_rsa_model(
        fx$dset$dataset, design, method = method, crossval = fx$cv,
        ncomp_selection = "max", feature_standardize = "center"
      )),
      "feature_rsa_model"
    )
  }
})

test_that("an F column constant only within a fold is refused at fit time", {
  fx <- .frdeg_fixture()
  # Overall this column varies, so construction succeeds; within the first
  # half of the rows it is constant, as a rare binary annotation would be in
  # some training folds.
  fold_F <- fx$F
  fold_F[, 3] <- rep(c(0, 1), each = fx$n / 2)
  train <- seq_len(fx$n / 2)
  design <- feature_rsa_design(
    F = fold_F, labels = paste0("t", seq_len(fx$n)), max_comps = 5L,
    block_var = fx$dset$design$block_var
  )
  X <- matrix(rnorm(fx$n * 6), fx$n, 6)
  methods <- c("pca", "pls", "ridge")
  if (requireNamespace("glmnet", quietly = TRUE)) methods <- c(methods, "glmnet")

  for (method in methods) {
    scaled <- suppressWarnings(feature_rsa_model(
      fx$dset$dataset, design, method = method, crossval = fx$cv,
      ncomp_selection = "max", feature_standardize = "scale"
    ))
    fit <- train_model(scaled, X[train, , drop = FALSE],
                       fold_F[train, , drop = FALSE],
                       indices = seq_len(ncol(X)))
    expect_match(fit$error, "1 of 10 columns that are constant")
    expect_match(fit$error, "column 3 \\('f3'\\)")

    centered <- suppressWarnings(feature_rsa_model(
      fx$dset$dataset, design, method = method, crossval = fx$cv,
      ncomp_selection = "max", feature_standardize = "center"
    ))
    expect_null(train_model(centered, X[train, , drop = FALSE],
                            fold_F[train, , drop = FALSE],
                            indices = seq_len(ncol(X)))$error)
  }
})

test_that("constant brain-data columns are still refused", {
  fx <- .frdeg_fixture()
  design <- feature_rsa_design(
    F = fx$F, labels = paste0("t", seq_len(fx$n)), max_comps = 5L,
    block_var = fx$dset$design$block_var
  )
  model <- suppressWarnings(feature_rsa_model(
    fx$dset$dataset, design, method = "pca", crossval = fx$cv,
    ncomp_selection = "max", feature_standardize = "center"
  ))
  X <- matrix(rnorm(fx$n * 6), fx$n, 6)

  constant_voxel <- X
  constant_voxel[, 2] <- 1
  fit <- train_model(model, constant_voxel, fx$F, indices = seq_len(ncol(X)))
  expect_match(fit$error, "Brain data X has 1 of 6 columns that are constant")

  # X is z-scored, so its units are irrelevant to whether it can be fitted.
  expect_null(train_model(model, X * 1e-6, fx$F,
                          indices = seq_len(ncol(X)))$error)

  missing <- X
  missing[4, 5] <- NA_real_
  expect_match(train_model(model, missing, fx$F,
                           indices = seq_len(ncol(X)))$error,
               "non-finite values")
})

test_that("components are capped at the numerical rank of F", {
  fx <- .frdeg_fixture()
  rank_deficient <- cbind(fx$F, const = 3, sum12 = fx$F[, 1] + fx$F[, 2])
  design <- feature_rsa_design(
    F = rank_deficient, labels = paste0("t", seq_len(fx$n)),
    max_comps = ncol(rank_deficient), block_var = fx$dset$design$block_var
  )
  model <- suppressWarnings(feature_rsa_model(
    fx$dset$dataset, design, method = "pca", crossval = fx$cv,
    ncomp_selection = "max", feature_standardize = "center"
  ))
  X <- matrix(rnorm(fx$n * 6), fx$n, 6)
  fit <- train_model(model, X, rank_deficient, indices = seq_len(ncol(X)))

  expect_null(fit$error)
  expect_equal(fit$ncomp, 10)
  expect_true(all(is.finite(fit$trained_model$coefficients)))

  # Without the cap the eleventh and twelfth components are undefined: pls
  # returns either NaN coefficients or coefficients around 1e27, and both
  # reach the caller as an all-NA performance table rather than as an error.
  uncapped <- .feature_rsa_kernel_fit_uncapped(
    rank_deficient, X, ncol(rank_deficient)
  )
  capped_scale <- max(abs(fit$trained_model$coefficients))
  expect_true(
    anyNA(uncapped) || !all(is.finite(uncapped)) ||
      max(abs(uncapped)) > 1e6 * capped_scale
  )

  res <- .frdeg_run(fx, rank_deficient, max_comps = ncol(rank_deficient))
  expect_true(all(is.finite(.frdeg_metrics(res))))
})

test_that("inner CV scores only component counts every segment defines", {
  set.seed(8504L)
  n <- 60L
  blocks <- rep(1:3, each = n / 3L)
  F <- matrix(rnorm(n * 6), n, 6)
  colnames(F) <- paste0("f", 1:6)
  # This column carries information only inside block 3, so the segment that
  # holds block 3 out trains on a constant column and cannot define six
  # components.
  F[, 6] <- ifelse(blocks == 3L, rnorm(sum(blocks == 3L)), 0)
  X <- matrix(rnorm(n * 5), n, 5)
  segments <- unname(split(seq_len(n), blocks))

  scores <- rMVPA:::.feature_rsa_cv_segment_scores(
    F, X, ncomp = 6L, method = "pca", segments = segments,
    fold_standardize = TRUE, feature_standardize = "center"
  )
  degenerate_segment <- which(vapply(segments, function(idx) {
    train <- setdiff(seq_len(n), idx)
    rMVPA:::.feature_rsa_component_rank(
      scale(F[train, , drop = FALSE], TRUE, FALSE)
    ) < 6L
  }, logical(1)))
  expect_length(degenerate_segment, 1L)

  # The undefined component is left unscored rather than scored with the
  # NaN or ~1e30 numbers pls returns for it.
  expect_true(all(is.finite(scores[-degenerate_segment, ])))
  expect_true(is.na(scores[degenerate_segment, 6L]))
  expect_true(all(is.finite(scores[degenerate_segment, 1:5])))

  # And a component one segment could not score is not a candidate.
  selected <- rMVPA:::.feature_rsa_select_from_segment_scores(
    scores, maximize = FALSE, one_se = FALSE
  )
  expect_true(selected >= 1L && selected <= 5L)

  # The fallback still returns a component when no count is complete.
  patchy <- matrix(c(1, NA, 2, NA), nrow = 2L)
  expect_identical(
    rMVPA:::.feature_rsa_select_from_segment_scores(
      patchy, maximize = FALSE, one_se = FALSE
    ),
    1L
  )
})
