test_that("Regression: grouped feature screening preserves duplicate feature_id groups", {
  set.seed(811)
  n_obs <- 90L
  n_groups <- 20L
  group_size <- 5L
  p <- n_groups * group_size

  feature_ids <- rep(seq_len(n_groups), each = group_size)
  y <- rnorm(n_obs)

  X <- matrix(0, nrow = n_obs, ncol = p)
  hot_groups <- 1:4
  scale_vec <- c(1, 0.8, 0.6, 0.4, 0.2)

  for (g in seq_len(n_groups)) {
    cols <- which(feature_ids == g)
    if (g %in% hot_groups) {
      X[, cols] <- outer(y - mean(y), scale_vec)
    } else {
      X[, cols] <- matrix(0, nrow = n_obs, ncol = length(cols))
    }
  }

  support <- rMVPA:::.spacenet_screen_support(
    X = X,
    y = y,
    screening_percentile = 10,
    min_features = 20L,
    feature_ids = feature_ids
  )

  expect_true(length(support) > 0L)
  support_groups <- unique(feature_ids[support])
  expect_equal(sort(support_groups), hot_groups)
  expect_equal(length(support), length(hot_groups) * group_size)
})

test_that("Regression: PDHG accepts duplicated feature_id support maps in binary fit flow", {
  set.seed(812)
  n_obs <- 40L
  n_groups <- 30L
  group_size <- 4L
  p <- n_groups * group_size

  feature_ids <- rep(seq_len(n_groups), each = group_size)
  y <- ifelse(rbinom(n_obs, size = 1, prob = 0.5) == 1, "case", "control")
  y <- factor(y, levels = c("control", "case"))

  X <- matrix(rnorm(n_obs * p), nrow = n_obs, ncol = p)
  support <- rep(1:15, each = 2L)

  edge_info <- rMVPA:::.spacenet_edges_from_feature_ids(feature_ids = support, spatial_mask = NULL)
  fit <- rMVPA:::.spacenet_pdhg_solver(
    X = X[, support, drop = FALSE],
    y = rMVPA:::.spacenet_y01(y),
    alpha = 0.05,
    l1_ratio = 0.7,
    edges = edge_info$edges,
    d_norm2 = edge_info$d_norm2,
    loss = "logistic",
    max_iter = 180L,
    tol = 1e-5
  )

  expect_equal(length(fit$w), length(support))
  expect_true(all(is.finite(fit$w)))
  expect_true(edge_info$d_norm2 >= 0)
  expect_true(all(is.finite(fit$init$w)))
  expect_true(all(is.finite(fit$init$wbar)))
  expect_true(all(is.finite(fit$init$p_dual)))
  expect_true(all(is.finite(fit$rel_history)))
})

test_that("Regression: CV alpha selection returns finite structure for duplicated feature ids", {
  set.seed(813)
  n <- 80L
  p <- 120L
  feature_ids <- rep(rep(seq_len(24), each = 5L))
  y <- rnorm(n)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)

  cv <- rMVPA:::.spacenet_cv_path(
    X = X,
    y = y,
    feature_ids = feature_ids,
    spatial_mask = NULL,
    l1_ratio = 0.6,
    n_alphas = 7L,
    alpha_min_ratio = 1e-3,
    screening_percentile = 20,
    max_iter = 160L,
    tol = 1e-4,
    inner_folds = 3L,
    is_classif = FALSE
  )

  expect_length(cv$alpha_grid, 7L)
  expect_length(cv$cv_loss, 7L)
  expect_true(all(is.finite(cv$alpha_grid)))
  expect_false(anyNA(cv$cv_loss))
  expect_true(all(is.finite(cv$cv_loss) | is.infinite(cv$cv_loss)))
})
