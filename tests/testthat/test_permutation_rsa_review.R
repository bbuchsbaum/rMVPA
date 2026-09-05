# Regression fixtures for the review of the RSA item-permutation null.
review_rsa_fixture <- function(n = 24L) {
  set.seed(7201)
  u <- runif(n, 1, 2)
  A <- outer(u, u, "+"); diag(A) <- 0
  Z <- matrix(0, n, n)
  Z[lower.tri(Z)] <- rnorm(choose(n, 2))
  N <- 1 + 0.03 * (Z + t(Z)); diag(N) <- 0
  e <- rnorm(n)
  Y <- N + 0.001 * outer(e, e, "+"); diag(Y) <- 0
  K <- 1 - Y; diag(K) <- 1
  # Actual centered neural patterns with correlation-distance RDM Y.
  basis <- matrix(rnorm((n + 1L) * n), n + 1L, n)
  basis <- sweep(basis, 2, colMeans(basis), "-")
  patterns <- t(chol(K)) %*% t(qr.Q(qr(basis)))
  dims <- c(n + 1L, 1L, 1L)
  dataset <- mvpa_dataset(
    train_data = neuroim2::NeuroVec(array(t(patterns), c(dims, n)),
                                   neuroim2::NeuroSpace(c(dims, n))),
    mask = neuroim2::LogicalNeuroVol(array(TRUE, dims), neuroim2::NeuroSpace(dims))
  )
  list(dataset = dataset, patterns = patterns, A = A, N = N, Y = Y)
}

test_that("conditional RSA coefficient inference is refused before any fitting", {
  fx <- review_rsa_fixture()
  expect_equal(1 - cor(t(fx$patterns)), fx$Y, tolerance = 1e-14)
  nuisance_des <- rsa_design(~ A, list(A = fx$A), nuisance = list(N = fx$N))
  two_models <- rsa_design(~ A + N, list(A = fx$A, N = fx$N))
  local_mocked_bindings(
    run_searchlight = function(...) stop("observed pass should not run"),
    mvpa_iterate = function(...) stop("permutation pass should not run"),
    .package = "rMVPA"
  )
  for (des in list(nuisance_des, two_models)) {
    for (reg in c("lm", "rfit")) {
      m <- rsa_model(fx$dataset, des, regtype = reg, check_collinearity = FALSE)
      for (strategy in c("iterate", "searchlight")) {
        pc <- permutation_control(n_perm = 1, perm_strategy = strategy)
        # Selecting one output must not hide the other fitted predictors.
        expect_error(run_permutation_searchlight(m, metric = "A", perm_ctrl = pc),
                     "individual.*coefficient.*joint")
        expect_error(run_permutation_searchlight(m, observed = rep(0, 25),
                                                 metric = "A", perm_ctrl = pc),
                     "individual.*coefficient.*joint")
      }
    }
  }
  m <- rsa_model(fx$dataset, nuisance_des, regtype = "lm", check_collinearity = FALSE)
  m$semipartial <- TRUE
  expect_error(run_permutation_searchlight(m), "individual.*coefficient")
  m$nneg <- list(A = TRUE)
  expect_error(run_permutation_searchlight(m), "individual.*coefficient")

  # Old serialized controls also get the safe default.
  pc_old <- permutation_control(); pc_old$rsa_null <- NULL
  expect_error(run_permutation_searchlight(m, perm_ctrl = pc_old), "individual.*coefficient")
})

test_that("joint RSA inference records the whole design and labels its null", {
  fx <- review_rsa_fixture(8L)
  des <- rsa_design(~ A, list(A = fx$A), nuisance = list(N = fx$N))
  m <- rsa_model(fx$dataset, des, regtype = "lm", distmethod = "pearson",
                 check_collinearity = FALSE)
  old_level <- futile.logger::flog.threshold()
  futile.logger::flog.threshold(futile.logger::ERROR)
  withr::defer(futile.logger::flog.threshold(old_level))
  for (strategy in c("iterate", "searchlight")) {
    pc <- permutation_control(n_perm = 2, subsample = 1, seed = 31,
                               null_method = "global", correction = "none",
                               diagnose = FALSE, perm_strategy = strategy,
                               rsa_null = "joint")
    res <- run_permutation_searchlight(m, radius = 20, perm_ctrl = pc)
    expect_s3_class(res, "permutation_result")
    expect_identical(res$rsa_null, "joint")
    expect_identical(res$null_predictors, c("A", "N"))
    expect_match(res$null_hypothesis, "Joint no association.*nuisance")
    expect_true(all(is.finite(res$p_values)))
    expect_output(print(res), "Joint no association")
    expect_output(summary(res), "Joint no association")
  }
  expect_error(permutation_control(rsa_null = "conditional"), "arg")
})

test_that("single-predictor regression and marginal correlations stay available", {
  fx <- review_rsa_fixture(8L)
  single <- rsa_design(~ A, list(A = fx$A))
  multiple <- rsa_design(~ A, list(A = fx$A), nuisance = list(N = fx$N))
  # The guard is deliberately independent of whether an engine can run.
  local_mocked_bindings(run_searchlight = function(...) stop("reached observed pass"),
                       .package = "rMVPA")
  for (reg in c("lm", "rfit")) {
    m <- rsa_model(fx$dataset, single, regtype = reg, check_collinearity = FALSE)
    expect_error(run_permutation_searchlight(m), "reached observed pass")
  }
  for (reg in c("pearson", "spearman")) {
    m <- rsa_model(fx$dataset, multiple, regtype = reg)
    expect_error(run_permutation_searchlight(m), "reached observed pass")
  }
})

test_that("circular shuffling samples every block rotation including zero", {
  set.seed(23)
  # Two blocks of two items have exactly four allowed transformations.
  draws <- replicate(1024, rMVPA:::.permutation_index(4L, c(1, 1, 2, 2),
                                                    "circular_shift"), simplify = FALSE)
  key <- function(p) paste(p, collapse = ",")
  counts <- table(vapply(draws, key, character(1)))
  expect_setequal(names(counts), c("1,2,3,4", "2,1,3,4", "1,2,4,3", "2,1,4,3"))
  expect_true(all(abs(as.numeric(counts) / 1024 - 0.25) < 0.06))
  # The same rule applies with no explicit blocks and with singleton blocks.
  unblocked <- replicate(256, key(rMVPA:::.permutation_index(3L, NULL, "circular_shift")))
  expect_setequal(unique(unblocked), c("1,2,3", "2,3,1", "3,1,2"))
  expect_identical(rMVPA:::.permutation_index(1L, NULL, "circular_shift"), 1L)
  expect_identical(rMVPA:::.permutation_index(3L, 1:3, "circular_shift"), 1:3)
})

test_that("circular RSA p-values agree with the exact four-state null", {
  fx <- review_rsa_fixture(4L)
  des <- rsa_design(~ A, list(A = fx$A, block = c(1, 1, 2, 2)), block_var = ~ block)
  m <- rsa_model(fx$dataset, des, regtype = "pearson", distmethod = "pearson")
  group <- list(1:4, c(2L, 1L, 3L, 4L), c(1L, 2L, 4L, 3L), c(2L, 1L, 4L, 3L))
  statistic <- function(rows) train_model(m, fx$patterns[rows, , drop = FALSE], NULL, NULL)[["A"]]
  exact <- vapply(group, statistic, numeric(1))
  set.seed(392)
  draws <- replicate(512, permute_labels(des, "circular_shift")$item_perm, simplify = FALSE)
  for (g in group) {
    obs <- statistic(g)
    null <- vapply(draws, function(p) statistic(g[p]), numeric(1))
    pool <- rMVPA:::build_adjusted_null(null, data.frame(nfeatures = rep(5L, length(null))),
                                       method = "global")
    p <- rMVPA:::score_observed(obs, pool, data.frame(nfeatures = 5L))
    expect_equal(p, mean(exact >= obs), tolerance = 0.08)
    expect_gte(p, 0.15) # the exact minimum is .25; the old sampler gave 1/513
  }
})

test_that("single predictors without usable variation have no effective support", {
  fx <- review_rsa_fixture()
  for (value in c(0, 7, NA_real_)) {
    constant <- as.dist(matrix(value, 24, 24))
    des <- rsa_design(~ constant, list(constant = constant))
    d <- rsa_design_diagnostics(des)
    expect_true(is.na(d$vif[["constant"]]))
    expect_true(is.na(d$items_per_predictor[["constant"]]))
    expect_identical(d$unsupported_predictors, "constant")
    expect_output(print(d), "no usable variation")
    expect_warning(suppressMessages(rsa_model(fx$dataset, des, regtype = "lm")),
                   "constant.*no usable variation")
    expect_warning(rsa_model(fx$dataset, des, regtype = "rfit"),
                   "constant.*no usable variation")
    expect_no_warning(rsa_model(fx$dataset, des, regtype = "lm", check_collinearity = FALSE))
  }
  valid <- rsa_design_diagnostics(rsa_design(~ A, list(A = fx$A)))
  expect_identical(valid$unsupported_predictors, character())
  expect_equal(valid$vif[["A"]], 1)
  expect_equal(valid$items_per_predictor[["A"]], 24)
})
