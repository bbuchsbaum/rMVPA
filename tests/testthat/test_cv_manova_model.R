.cvm_fixture <- function() {
  set.seed(2411)
  runs <- rep(c("a", "b", "c"), c(24, 27, 30))
  t <- seq_along(runs)
  trials <- data.frame(place = factor(t %% 2), vowel = factor(t %% 3),
                       drift = sin(t / 7))
  X <- model.matrix(~ place + vowel + drift, trials,
                    contrasts.arg = list(place = contr.sum(2), vowel = contr.sum(3)))
  p <- 8L
  noise <- matrix(rnorm(length(t) * p), ncol = p)
  noise[, 2L] <- noise[, 2L] + .7 * noise[, 1L]
  beta <- matrix(seq_len(ncol(X) * p) / 17, ncol = p)
  Y <- X %*% beta + noise
  dims <- c(2L, 2L, 2L)
  space <- neuroim2::NeuroSpace(dims)
  dataset <- mvpa_dataset(
    neuroim2::NeuroVec(array(t(Y), c(dims, length(t))),
                       neuroim2::NeuroSpace(c(dims, length(t)))),
    mask = as.logical(neuroim2::NeuroVol(array(1, dims), space)))
  region <- neuroim2::NeuroVol(array(rep(1:2, each = 4L), dims), space)
  list(Y = Y, X = X, dataset = dataset, region = region, runs = runs,
       design = manova_design(~ place + vowel + drift, trials))
}

.cvm_fit <- function(model, Y) {
  fit_roi(model, list(train_data = Y, indices = seq_len(ncol(Y))), list(id = 1L))
}

test_that("regional and searchlight D match the independent cvManovaCore reference", {
  f <- .cvm_fixture()
  C <- diag(5)
  hypotheses <- list(place = C[2, ], vowel = C[3:4, ], joint = C[2:4, ])
  model <- cv_manova_model(f$dataset, f$design, f$runs, hypotheses)
  # Octave, allefeld/cvmanova@8c2b115, cvManovaCore, lambda=0; unequal runs,
  # correlated voxel noise, one nuisance regressor and a two-dimensional effect.
  expected <- rbind(c(2.918010872534074, 6.061666953583226, 8.679749054982841),
                    c(15.14362701457055, 32.27578637684744, 45.9457954308638))
  expect_equal(validate_analysis(model, verbose = FALSE)$n_warn, 0L)
  regional <- run_regional(model, f$region)
  table <- as.data.frame(regional$performance_table)
  expect_equal(unname(as.matrix(table[order(table$roinum), names(hypotheses)])),
               unname(expected), tolerance = 1e-9)
  searchlight <- run_searchlight(model, radius = 4, method = "standard")
  whole <- c(17.02929443394475, 36.69149150488192, 52.14420514213938)
  expect_named(searchlight$results, names(hypotheses))
  for (i in seq_along(whole)) {
    expect_equal(as.numeric(neuroim2::values(searchlight$results[[i]])), rep(whole[i], 8),
                 tolerance = 1e-9)
  }
  terms <- cv_manova_model(f$dataset, f$design, f$runs, c("place", "vowel"))
  expect_equal(unname(.cvm_fit(terms, f$Y)$metrics), whole[1:2], tolerance = 1e-9)
})

test_that("contrast coordinates, run labels and diagonal precision have defined units", {
  f <- .cvm_fixture()
  build <- function(effects, ...) cv_manova_model(f$dataset, f$design, f$runs, effects, ...)
  C <- diag(5)[3:4, ]
  colnames(C) <- colnames(f$X)
  rotated <- matrix(c(1, 1, 1, -1), 2) %*% C
  a <- build(list(vowel = C))
  b <- build(list(vowel = rotated[, 5:1]))
  expect_equal(.cvm_fit(a, f$Y)$metrics, .cvm_fit(b, f$Y)$metrics)
  diag_model <- build("vowel", precision = "diagonal")
  univariate <- vapply(seq_len(ncol(f$Y)), function(i) {
    unname(.cvm_fit(a, f$Y[, i, drop = FALSE])$metrics)
  }, numeric(1))
  expect_equal(unname(.cvm_fit(diag_model, f$Y)$metrics), sum(univariate))
  adjusted <- build("vowel", residual_df = rev(a$residual_df - 2))
  expect_equal(adjusted$residual_df, a$residual_df - 2)
  expect_lt(.cvm_fit(adjusted, f$Y)$metrics, .cvm_fit(a, f$Y)$metrics)
})

test_that("opposite feature patterns in independent runs yield negative D", {
  f <- .cvm_fixture()
  block <- ifelse(f$runs == "a", 1L, 2L)
  model <- cv_manova_model(f$dataset, f$design, block, "place")
  y <- f$Y[, 1, drop = FALSE] + 12 * f$X[, 2] * ifelse(block == 1, 1, -1)
  expect_lt(unname(.cvm_fit(model, y)$metrics), 0)
})

test_that("invalid designs, contrasts and unavailable precision are explicit", {
  f <- .cvm_fixture()
  build <- function(...) cv_manova_model(f$dataset, f$design, f$runs, ...)
  expect_error(build("missing"), "Unknown effect")
  expect_error(build(list(a = numeric(5))), "rank zero")
  expect_error(build(list(a = c(1, 0))), "one column per coefficient")
  expect_error(build("place", residual_df = c(20, 20)), "residual_df")
  expect_error(cv_manova_model(f$dataset, f$design, rep(1, length(f$runs))),
               "two independent runs")
  d <- f$design$data
  d$copy <- d$drift
  aliased <- manova_design(~ place + drift + copy, d)
  expect_error(cv_manova_model(f$dataset, aliased, f$runs, "drift"), "not estimable")
  ok <- cv_manova_model(f$dataset, aliased, f$runs, "place")
  expect_false(.cvm_fit(ok, f$Y)$error)
  singular <- .cvm_fit(build("place"), f$Y[, c(1, 1)])
  expect_true(singular$error)
  expect_match(singular$error_message, "Singular")
  exact_fit <- matrix(1, nrow(f$Y), 1)
  expect_true(.cvm_fit(build("place"), exact_fit)$error)
  expect_true(.cvm_fit(build("place", precision = "diagonal"), exact_fit)$error)
  large <- matrix(0, nrow(f$Y), 100)
  expect_match(.cvm_fit(build("place"), large)$error_message, "degrees of freedom")
  expect_match(.cvm_fit(build("place", precision = "diagonal"), large)$error_message,
               "Zero training residual variance")
})
