.relational_fixture <- function(n = 12L) {
  set.seed(930)
  E <- matrix(rnorm(n * 32L), n)
  R <- E + matrix(rnorm(n * 32L), n)
  data <- data.frame(precision = rnorm(n), vividness = rnorm(n),
                     task = factor(rep(c("colour", "location"), length.out = n)))
  X <- rbind(E, R)
  dataset <- structure(list(train_data = t(X), mask = seq_len(ncol(X))),
                       class = "mvpa_dataset")
  list(E = E, R = R, X = X, data = data, dataset = dataset,
       Y = cor(t(E), t(R)), n = n)
}

.relational_design <- function(fx, ...) {
  pair_rsa_design(
    seq_len(fx$n), seq_len(fx$n), pairs = "between",
    row_idx_a = seq_len(fx$n), row_idx_b = fx$n + seq_len(fx$n),
    features_b = fx$data,
    model = list(item = function(a, b) as.numeric(a == b)), ...
  )
}

.relational_fit <- function(fx, design, ...) {
  suppressWarnings(rsa_model(fx$dataset, design, distmethod = "pearson",
                            regtype = "lm", measure = "similarity",
                            statistic = "beta", ...))
}

test_that("relational coefficients reproduce classical and modulated ERA", {
  fx <- .relational_fixture()
  simple <- .relational_fit(fx, .relational_design(fx))
  result <- train_model(simple, fx$X, NULL, NULL)
  expect_equal(unname(result["item"]), mean(diag(fx$Y)) - mean(fx$Y[row(fx$Y) != col(fx$Y)]),
               tolerance = 1e-12)

  design <- .relational_design(fx, modulation = list(item = ~ b.precision * b.vividness),
                               nuisance = ~ factor(b.row))
  model <- .relational_fit(fx, design,
                           contrasts = list(difference = c(item.b.precision = 1,
                                                           item.b.vividness = -1)))
  result <- train_model(model, fx$X, NULL, NULL)
  specificity <- vapply(seq_len(fx$n), function(j) {
    fx$Y[j, j] - mean(fx$Y[-j, j])
  }, numeric(1))
  reference <- coef(lm(specificity ~ precision * vividness, data = fx$data))
  expect_equal(unname(result[seq_along(reference)]), unname(reference), tolerance = 1e-12)
  expect_equal(unname(result["contrast_difference"]),
               unname(reference["precision"] - reference["vividness"]), tolerance = 1e-12)
  expect_identical(names(result), names(output_schema(model)))
  expect_false(any(grepl("background", names(result))))
  expect_identical(design$modulation$item, ~ b.precision * b.vividness)
  expect_identical(design$formula_terms$item$assign, c(0L, 1L, 2L, 3L))
})

test_that("factor modulation uses ordinary model.matrix coding", {
  fx <- .relational_fixture()
  d <- .relational_design(fx, modulation = list(item = ~ b.precision * b.task),
                          nuisance = ~ factor(b.row))
  m <- .relational_fit(fx, d)
  y <- diag(fx$Y) - (colSums(fx$Y) - diag(fx$Y)) / (fx$n - 1L)
  expected <- coef(lm(y ~ precision * task, fx$data))
  expect_equal(unname(train_model(m, fx$X, NULL, NULL)), unname(expected), tolerance = 1e-12)
  expect_identical(d$formula_terms$item$contrasts, list(b.task = "contr.treatment"))
})

test_that("rectangular repeated correspondence preserves observations and order", {
  fx <- .relational_fixture()
  ia <- c(1:6, 1L, 4L)
  ib <- c(2L, 1L, 4L, 1L, 9L)
  fa <- fx$data[ia, , drop = FALSE]
  fb <- fx$data[ib, , drop = FALSE]
  fb$precision[4] <- fb$precision[4] + 0.5 # repeated attempts have distinct behaviour
  build <- function(pa, pb) pair_rsa_design(
    ia[pa], ib[pb], pairs = "between", row_idx_a = ia[pa], row_idx_b = fx$n + ib[pb],
    features_a = fa[pa, , drop = FALSE], features_b = fb[pb, , drop = FALSE],
    model = list(item = function(a, b) as.numeric(a == b)),
    modulation = list(item = ~ b.precision), nuisance = ~ factor(b.row))
  original <- build(seq_along(ia), seq_along(ib))
  expect_equal(original$n_a, length(ia))
  expect_equal(original$n_b, length(ib))
  expect_equal(original$model_mat$item, as.numeric(outer(ia, ib, `==`)))
  result <- train_model(.relational_fit(fx, original), fx$X, NULL, NULL)

  # Independent dense regression oracle; no averaging of repeated attempts.
  Y <- cor(t(fx$E[ia, ]), t(fx$R[ib, ]))
  grid <- expand.grid(i = seq_along(ia), j = seq_along(ib))
  reference <- lm(as.vector(Y) ~ factor(grid$j) +
                    as.vector(outer(ia, ib, `==`)) +
                    I(as.vector(outer(ia, ib, `==`)) * fb$precision[grid$j]))
  expect_equal(unname(result), unname(tail(coef(reference), 2L)), tolerance = 1e-12)
  reordered <- build(c(8, 2, 5, 1, 7, 3, 4, 6), c(5, 3, 1, 4, 2))
  expect_equal(train_model(.relational_fit(fx, reordered), fx$X, NULL, NULL),
               result, tolerance = 1e-12)
})

test_that("unequal eligible comparisons produce the correct weighted ERA reduction", {
  fx <- .relational_fixture()
  block_a <- rep(1:3, c(3, 4, 5))
  block_b <- block_a %% 3 + 1L
  d <- .relational_design(fx, modulation = list(item = ~ b.precision),
                          nuisance = ~ factor(b.row),
                          block_var_a = block_a, block_var_b = block_b)
  Y <- fx$Y
  Y[!matrix(d$include, fx$n)] <- NA_real_
  r <- colSums(is.finite(Y)) - 1L
  specificity <- diag(Y) - (colSums(Y, na.rm = TRUE) - diag(Y)) / r
  reference <- coef(lm(specificity ~ precision, fx$data, weights = r / (r + 1)))
  result <- train_model(.relational_fit(fx, d), fx$X, NULL, NULL)
  expect_equal(unname(result), unname(reference), tolerance = 1e-12)
  expect_gt(max(abs(result - coef(lm(specificity ~ precision, fx$data)))), 1e-6)
})

test_that("missing measurement cells rebuild queries on the actual eligible design", {
  fx <- .relational_fixture()
  d <- .relational_design(fx, modulation = list(item = ~ b.precision),
                          nuisance = ~ factor(b.row))
  m <- .relational_fit(fx, d)
  y <- as.vector(fx$Y)
  y[c(3, 7, 17, 18, 40)] <- NA_real_
  X <- cbind(1, do.call(cbind, d$model_mat))
  keep <- is.finite(y)
  reference <- lm.fit(X[keep, , drop = FALSE], y[keep])$coefficients[2:3]
  expect_equal(unname(rMVPA:::.rsa_coefficient_metrics(y, m)), unname(reference), tolerance = 1e-12)

  y[seq_len(fx$n)] <- NA_real_ # an entire target loses support
  expect_error(rMVPA:::.rsa_coefficient_metrics(y, m), "aliased or unsupported")
})

test_that("measurement choice and coefficient fallback preserve the declared observable", {
  fx <- .relational_fixture()
  d <- .relational_design(fx, modulation = list(item = ~ b.precision),
                          nuisance = ~ factor(b.row))
  m <- .relational_fit(fx, d)
  beta <- train_model(m, fx$X, NULL, NULL)
  fallback <- m
  fallback$.coefficient_query <- NULL
  expect_equal(train_model(fallback, fx$X, NULL, NULL), beta, tolerance = 1e-12)
  distance <- m
  distance$measure <- "distance"
  expect_equal(train_model(distance, fx$X, NULL, NULL), -beta, tolerance = 1e-12)

  changed <- fx
  changed$data$precision <- rev(changed$data$precision)
  changed_d <- .relational_design(changed, modulation = list(item = ~ b.precision),
                                  nuisance = ~ factor(b.row))
  expect_identical(d$row_idx_a, changed_d$row_idx_a)
  expect_identical(d$row_idx_b, changed_d$row_idx_b)
  # An unchanged template's marginal correlation uses the same neural evidence.
  marginal <- function(design) train_model(rsa_model(
    fx$dataset, design, distmethod = "pearson", measure = "similarity"), fx$X, NULL, NULL)["item"]
  expect_equal(marginal(d), marginal(changed_d), tolerance = 1e-12)
  expect_false(identical(d$model_mat$item.b.precision, changed_d$model_mat$item.b.precision))
})

test_that("symmetric modulation requires an explicit symmetric rule", {
  features <- data.frame(vividness = 1:6)
  model <- list(semantic = function(a, b) abs(a - b))
  expect_error(pair_rsa_design(1:6, features_a = features, model = model,
                               modulation = list(semantic = ~ b.vividness)), "must be symmetric")
  d <- pair_rsa_design(1:6, features_a = features, model = model,
                       modulation = list(semantic = ~ I((a.vividness + b.vividness) / 2)))
  expect_length(d$model_mat, 2L)
  expect_equal(unname(d$model_mat[[2]]),
               abs(d$pair_index$item_a - d$pair_index$item_b) *
                 (features$vividness[d$pair_index$i] + features$vividness[d$pair_index$j]) / 2)
})

test_that("aliases and invalid specifications fail explicitly", {
  fx <- .relational_fixture()
  expect_error(.relational_design(fx, modulation = list(unknown = ~ b.precision)),
                "keyed by unique names")
  expect_error(.relational_design(fx, modulation = list(item = ~ precision)), "unknown pair metadata")
  expect_error(.relational_design(fx, modulation = list(item = precision ~ vividness)),
                "right-hand-side formula")
  expect_error(.relational_design(fx, nuisance = ~ factor(b.row) + b.precision),
                "rank deficient including its intercept")
  expect_error(.relational_design(fx, nuisance = ~ 0 + factor(b.row)), "rank deficient")
  d <- .relational_design(fx)
  expect_error(.relational_fit(fx, d, semipartial = TRUE), "unconstrained")
  expect_error(rsa_model(fx$dataset, d, statistic = "beta"), "unconstrained")
  expect_error(rsa_model(fx$dataset, d, contrasts = list(x = c(item = 1))), "requires statistic")
  expect_error(.relational_fit(fx, d, contrasts = list(x = c(unknown = 1))), "finite named weights")
  expect_error(.relational_fit(fx, d, contrasts = c(item = 1)), "uniquely named list")
  expect_error(.relational_fit(fx, d, contrasts = list(x = c(item = Inf))), "finite named weights")
})

test_that("relational schema and coefficients pass through existing ROI execution", {
  fx <- .relational_fixture()
  d <- .relational_design(fx, modulation = list(item = ~ b.precision),
                          nuisance = ~ factor(b.row))
  m <- .relational_fit(fx, d)
  expected <- train_model(m, fx$X, NULL, NULL)
  roi <- fit_roi(m, list(train_data = fx$X, indices = 1:32), list(id = "region"))
  expect_false(roi$error)
  expect_equal(roi$metrics, expected, tolerance = 1e-12)
  expect_identical(names(output_schema(m)), names(expected))
})

test_that("matrix-free queries match dense measurements under ranking and centering", {
  fx <- .relational_fixture()
  d <- .relational_design(fx, modulation = list(item = ~ b.precision),
                          nuisance = ~ factor(b.row))
  for (method in c("pearson", "spearman")) {
    for (centering in c("none", "stimulus_mean")) {
      m <- suppressWarnings(rsa_model(
        fx$dataset, d, regtype = "lm", statistic = "beta", measure = "similarity",
        distmethod = method, pattern_center = centering,
        contrasts = list(background = c(`(Intercept)` = 1))))
      centered <- rMVPA:::center_patterns(fx$X, centering)
      Y <- cor(t(centered[1:fx$n, ]), t(centered[fx$n + 1:fx$n, ]), method = method)
      X <- cbind(`(Intercept)` = 1, do.call(cbind, d$model_mat))
      reference <- lm.fit(X, as.vector(Y))$coefficients
      expected <- c(reference[d$model_predictors], contrast_background = reference[1])
      expect_equal(unname(train_model(m, fx$X, NULL, NULL)), unname(expected), tolerance = 1e-12)
      # Fingerprints require the explicit pair response and exercise that route.
      m$return_fingerprint <- TRUE
      expect_equal(unname(train_model(m, fx$X, NULL, NULL)), unname(expected), tolerance = 1e-12)
    }
  }

  # The same contraction applies to symmetric lower-triangle observations.
  within <- pair_rsa_design(1:fx$n, row_idx_a = 1:fx$n,
                            model = list(semantic = dist(fx$data$precision)))
  m <- .relational_fit(fx, within)
  Y <- cor(t(fx$X[1:fx$n, ]))
  X <- cbind(1, within$model_mat$semantic)
  reference <- lm.fit(X, Y[lower.tri(Y)])$coefficients[2]
  expect_equal(unname(train_model(m, fx$X, NULL, NULL)), unname(reference), tolerance = 1e-12)
})

test_that("modulated models execute through regional and searchlight interfaces", {
  fx <- .relational_fixture()
  dims <- c(4L, 4L, 2L)
  data <- neuroim2::NeuroVec(array(t(fx$X), c(dims, 2L * fx$n)),
                            neuroim2::NeuroSpace(c(dims, 2L * fx$n)))
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, dims), neuroim2::NeuroSpace(dims))
  fx$dataset <- mvpa_dataset(data, mask = mask)
  d <- .relational_design(fx, modulation = list(item = ~ b.precision),
                          nuisance = ~ factor(b.row))
  m <- .relational_fit(fx, d, contrasts = list(precision = c(item.b.precision = 1)))
  expected <- train_model(m, fx$X, NULL, NULL)
  regions <- neuroim2::NeuroVol(array(1L, dims), neuroim2::NeuroSpace(dims))
  regional <- run_regional(m, regions)
  expect_equal(as.numeric(regional$performance[1, names(expected)]),
               unname(expected), tolerance = 1e-10)
  searchlight <- run_searchlight(m, radius = 1.5, method = "standard")
  expect_true(all(names(expected) %in% names(searchlight$results)))
  for (nm in names(expected)) {
    expect_true(any(is.finite(neuroim2::values(searchlight$results[[nm]]))))
  }
})

test_that("multiple relationships are jointly adjusted and preserve predictor roles", {
  fx <- .relational_fixture()
  d <- pair_rsa_design(
    1:fx$n, 1:fx$n, pairs = "between", row_idx_a = 1:fx$n,
    row_idx_b = fx$n + 1:fx$n, features_a = fx$data, features_b = fx$data,
    model = list(item = function(a, b) as.numeric(a == b),
                 category = function(a, b, fa, fb) as.numeric(fa$task == fb$task)),
    modulation = list(item = ~ b.precision, category = ~ b.vividness),
    nuisance = ~ factor(b.row))
  index <- expand.grid(i = 1:fx$n, j = 1:fx$n)
  item <- as.vector(diag(fx$n))
  category <- as.vector(outer(fx$data$task, fx$data$task, `==`)) * 1
  X <- cbind(model.matrix(~ factor(index$j)), item,
             item * fx$data$precision[index$j], category,
             category * fx$data$vividness[index$j])
  reference <- tail(lm.fit(X, as.vector(fx$Y))$coefficients, 4L)
  expect_equal(unname(train_model(.relational_fit(fx, d), fx$X, NULL, NULL)),
               unname(reference), tolerance = 1e-12)
  expect_true(all(d$predictor_roles[d$model_predictors] == "model"))
  expect_true(all(d$predictor_roles[d$nuisance_predictors] == "nuisance"))
})

test_that("compiled coefficient queries respect item permutations and inference boundaries", {
  fx <- .relational_fixture()
  d <- .relational_design(fx, modulation = list(item = ~ b.precision),
                          nuisance = ~ factor(b.row))
  m <- .relational_fit(fx, d)
  expect_error(run_permutation_searchlight(m, metric = "item.b.precision"),
                "individual regression coefficient")
  permuted <- permute_labels(d, method = "global", seed = 123)
  m$design <- permuted
  Y <- cor(t(fx$E[permuted$item_perm, ]), t(fx$R[permuted$item_perm_b, ]))
  X <- cbind(1, do.call(cbind, d$model_mat))
  reference <- lm.fit(X, as.vector(Y))$coefficients[2:3]
  expect_equal(unname(train_model(m, fx$X, NULL, NULL)), unname(reference), tolerance = 1e-12)
})

test_that("within-domain coefficient fallback preserves missing neural observations", {
  fx <- .relational_fixture()
  d <- pair_rsa_design(
    1:fx$n, row_idx_a = 1:fx$n, features_a = fx$data,
    model = list(semantic = function(a, b) abs(a - b)),
    modulation = list(semantic = ~ I((a.precision + b.precision) / 2)))
  X <- cbind(1, do.call(cbind, d$model_mat))
  for (method in c("pearson", "spearman")) {
    for (missingness in c("cell", "observation")) {
      patterns <- fx$X
      if (missingness == "cell") patterns[2L, 4L] <- NA_real_ else patterns[2L, ] <- NA_real_
      Y <- cor(t(patterns[1:fx$n, ]), method = method)
      y <- Y[lower.tri(Y)]
      eligible <- is.finite(y)
      expect_equal(sum(eligible), choose(fx$n - 1L, 2L))
      reference <- lm.fit(X[eligible, , drop = FALSE], y[eligible])$coefficients[-1L]
      for (fingerprint in c(FALSE, TRUE)) {
        m <- suppressWarnings(rsa_model(
          fx$dataset, d, regtype = "lm", statistic = "beta", measure = "similarity",
          distmethod = method, return_fingerprint = fingerprint))
        actual <- train_model(m, patterns, NULL, NULL)
        expect_equal(as.numeric(actual), unname(reference), tolerance = 1e-10)
      }
    }
  }
  m <- suppressWarnings(rsa_model(fx$dataset, d, regtype = "lm", statistic = "beta",
                                  measure = "similarity", distmethod = "spearman"))
  patterns <- fx$X
  patterns[1:fx$n, ] <- NA_real_
  expect_error(train_model(m, patterns, NULL, NULL), "aliased or unsupported")
})
