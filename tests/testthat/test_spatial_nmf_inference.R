test_that("spatial_nmf_component_test detects one differing component", {
  set.seed(101)
  W_A <- cbind(runif(6, 1.5, 2.0), runif(6, 0.5, 0.8))
  W_B <- cbind(runif(6, 0.1, 0.3), runif(6, 0.5, 0.8))
  W <- rbind(W_A, W_B)
  groups <- rep(c("A", "B"), each = 6)

  res <- spatial_nmf_component_test(
    W = W,
    groups = groups,
    test = "two_group",
    nperm = 200,
    correction = "maxT",
    seed = 42
  )

  expect_equal(res$n_sig, 1)
  expect_lt(res$table$p_fwer[1], 0.05)
  expect_gt(res$table$p_fwer[2], 0.05)
})

test_that("spatial_nmf_component_test handles one-group null_W", {
  set.seed(102)
  W_obs <- cbind(runif(8, 1.0, 1.5), runif(8, 0.1, 0.2))

  null_list <- lapply(seq_len(200), function(i) {
    cbind(runif(8, 0.1, 0.3), runif(8, 0.1, 0.3))
  })

  res <- spatial_nmf_component_test(
    W = W_obs,
    test = "one_group",
    null_W = null_list,
    correction = "maxT",
    alternative = "greater"
  )

  expect_equal(res$n_sig, 1)
  expect_lt(res$table$p_fwer[1], 0.05)
})

test_that("spatial_nmf_stability returns expected shapes and bounds", {
  set.seed(103)
  X <- matrix(rexp(12 * 15), nrow = 12, ncol = 15)

  fit <- rMVPA:::spatial_nmf_fit(X, k = 3, max_iter = 20, check_every = 5)

  stab <- spatial_nmf_stability(
    fit = fit,
    X = X,
    n_boot = 5,
    top_frac = 0.2,
    max_iter = 15,
    check_every = 5
  )

  expect_equal(dim(stab$mean), c(3, 15))
  expect_equal(dim(stab$sd), c(3, 15))
  expect_equal(dim(stab$cv), c(3, 15))
  expect_equal(dim(stab$selection), c(3, 15))
  expect_true(all(stab$selection >= 0 & stab$selection <= 1))
  expect_equal(length(stab$component_similarity), 3)
})

test_that("spatial_nmf_global_test detects group separation", {
  set.seed(110)
  nA <- 10
  nB <- 10
  k <- 1
  p <- 20

  H <- matrix(rexp(k * p), nrow = k, ncol = p)
  W_A <- matrix(rnorm(nA, 2.5, 0.1), ncol = 1)
  W_B <- matrix(rnorm(nB, 0.3, 0.05), ncol = 1)

  X <- rbind(W_A, W_B) %*% H
  X <- pmax(X + matrix(rnorm((nA + nB) * p, 0, 0.01), nA + nB, p), 0)
  groups <- rep(c("A", "B"), each = nA)

  res <- spatial_nmf_global_test(
    X = X,
    groups = groups,
    k = k,
    nfolds = 4,
    nperm = 60,
    seed = 99,
    max_iter = 50,
    check_every = 10,
    project_args = list(max_iter = 50, check_every = 10),
    return_cv = TRUE
  )

  expect_lt(res$p_value, 0.05)
  expect_gt(res$stat, 0.8)
  expect_equal(length(res$cv$pred), length(groups))
  expect_true(all(res$cv$pred >= 0 & res$cv$pred <= 1))
  expect_equal(length(res$cv$fold_id), length(groups))
})

test_that("spatial_nmf_maps runs inference hooks", {
  mask <- neuroim2::NeuroVol(array(1, dim = c(2, 2, 1)), neuroim2::NeuroSpace(c(2, 2, 1)))
  space_obj <- neuroim2::space(mask)

  mapA <- list(
    neuroim2::NeuroVol(array(1, dim = c(2, 2, 1)), space_obj),
    neuroim2::NeuroVol(array(1.1, dim = c(2, 2, 1)), space_obj),
    neuroim2::NeuroVol(array(0.9, dim = c(2, 2, 1)), space_obj)
  )
  mapB <- list(
    neuroim2::NeuroVol(array(0.1, dim = c(2, 2, 1)), space_obj),
    neuroim2::NeuroVol(array(0.2, dim = c(2, 2, 1)), space_obj),
    neuroim2::NeuroVol(array(0.15, dim = c(2, 2, 1)), space_obj)
  )

  res <- spatial_nmf_maps(
    group_A = mapA,
    group_B = mapB,
    mask = mask,
    k = 1,
    lambda = 0,
    max_iter = 10,
    component_test = list(nperm = 20, correction = "none"),
    global_test = list(
      nperm = 20,
      nfolds = 2,
      metric = "accuracy",
      max_iter = 20,
      check_every = 5,
      project_args = list(max_iter = 20, check_every = 5)
    ),
    stability = list(n_boot = 2, max_iter = 10, check_every = 5)
  )

  expect_true(is.list(res$component_test))
  expect_true(is.list(res$global_test))
  expect_true(is.list(res$stability))
})
