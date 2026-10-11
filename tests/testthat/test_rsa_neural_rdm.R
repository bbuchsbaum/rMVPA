test_that("neural RDMs agree with pairwise correlations including rank ties", {
  X <- rbind(c(1, 1, 4, 3, 2), c(5, 2, 2, 1, 3),
             c(4, 2, 1, 1, 5), c(1, 3, 5, 2, 2))
  rownames(X) <- c("run2_b", "run1_a", "run2_a", "run1_b")
  for (method in c("pearson", "spearman")) {
    for (center in c("none", "stimulus_mean")) {
      prepared <- if (center == "none") X else scale(X, scale = FALSE)
      expected <- matrix(0, nrow(X), nrow(X), dimnames = list(rownames(X), rownames(X)))
      for (i in seq_len(nrow(X))) {
        for (j in seq_len(nrow(X))) {
          expected[i, j] <- 1 - stats::cor(prepared[i, ], prepared[j, ], method = method)
        }
      }
      diag(expected) <- 0
      got <- rsa_neural_rdm(X, method = method, pattern_center = center)
      expect_equal(got, expected, tolerance = 1e-12)
      expect_identical(got, t(got))
      expect_identical(unname(diag(got)), rep(0, nrow(X)))
      expect_true(all(is.finite(got)))
    }
  }
  expect_identical(rsa_neural_rdm(X), rsa_neural_rdm(X, "spearman"))
  expect_false(isTRUE(all.equal(rsa_neural_rdm(X), rsa_neural_rdm(X, "pearson"))))
})

test_that("known identical, reversed and orthogonal patterns have expected distances", {
  X <- rbind(original = c(-3, -1, 1, 3), identical = c(-3, -1, 1, 3),
             reversed = c(3, 1, -1, -3), orthogonal = c(-1, 3, -3, 1))
  for (method in c("pearson", "spearman")) {
    D <- rsa_neural_rdm(X, method)
    expect_equal(unname(D[1, ]), c(0, 0, 2, 1), tolerance = 1e-12)
    # Correlation is invariant to positive rescaling and row-wise offsets.
    shifted <- sweep(sweep(X, 1, c(2, 3, 4, 5), "*"), 1, c(8, -2, 5, 1), "+")
    expect_equal(rsa_neural_rdm(shifted, method), D, tolerance = 1e-12)
    expect_equal(rsa_neural_rdm(X * 1e-100, method), D, tolerance = 1e-12)
  }
})

test_that("observation names, duplicates, order and uncentered subsets are preserved", {
  set.seed(119)
  X <- matrix(rnorm(30), 6, 5, dimnames = list(c("b", "a", "b", "c", "d", "e"), NULL))
  for (method in c("pearson", "spearman")) {
    D <- rsa_neural_rdm(X, method)
    take <- c(5, 2, 1)
    expect_identical(dimnames(D), list(rownames(X), rownames(X)))
    expect_equal(rsa_neural_rdm(X[take, ], method), D[take, take], tolerance = 1e-12)
    expect_null(dimnames(rsa_neural_rdm(unname(X), method)))
  }
  expect_equal(dim(rsa_neural_rdm(matrix(c(1, 2, 2, 1), 2))), c(2L, 2L))
})

test_that("stimulus centering uses the supplied observation set", {
  set.seed(120)
  X <- matrix(rnorm(30), 6, 5)
  before <- X
  for (method in c("pearson", "spearman")) {
    D <- rsa_neural_rdm(X, method, "stimulus_mean")
    shared <- sweep(X, 2, c(10, -4, 3, 2, 8), "+")
    expect_equal(rsa_neural_rdm(shared, method, "stimulus_mean"), D, tolerance = 1e-12)
    expect_false(isTRUE(all.equal(rsa_neural_rdm(X[1:3, ], method, "stimulus_mean"),
                                 D[1:3, 1:3])))
  }
  expect_identical(X, before)
})

test_that("undefined neural RDM inputs fail explicitly", {
  X <- rbind(c(1, 2, 4), c(3, 1, 2), c(2, 4, 1))
  for (bad in list(1:6, as.data.frame(X), matrix(letters[1:9], 3),
                   matrix(TRUE, 3, 3), X + 1i, array(1:8, c(2, 2, 2)))) {
    expect_error(rsa_neural_rdm(bad), "numeric matrix")
  }
  for (bad in list(matrix(numeric(), 0, 3), matrix(numeric(), 3, 0),
                   X[1, , drop = FALSE], X[, 1, drop = FALSE])) {
    expect_error(rsa_neural_rdm(bad), "at least two observations and two features")
  }
  for (value in c(NA_real_, NaN, Inf, -Inf)) {
    bad <- X
    bad[1, 1] <- value
    expect_error(rsa_neural_rdm(bad), "finite values")
  }
  for (method in c("pearson", "spearman")) {
    expect_error(rsa_neural_rdm(rbind(X, c(2, 2, 2)), method), "constant.*rows: 4")
    # Each raw row varies, but centering makes every row constant.
    expect_error(rsa_neural_rdm(rbind(1:3, 2:4), method, "stimulus_mean"), "constant")
  }
  expect_error(rsa_neural_rdm(X, method = "kendall"), "arg")
  expect_error(rsa_neural_rdm(X, pattern_center = "row_mean"), "arg")
})

test_that("public neural RDMs reproduce RSA fits with design-level run exclusions", {
  set.seed(121)
  n <- 6L
  dims <- c(2L, 2L, 1L)
  X <- matrix(rnorm(n * prod(dims)), n)
  mask <- neuroim2::NeuroVol(array(1, dims), neuroim2::NeuroSpace(dims))
  data <- neuroim2::NeuroVec(array(t(X), c(dims, n)), neuroim2::NeuroSpace(c(dims, n)))
  dataset <- mvpa_dataset(data, mask = mask)
  template <- stats::dist(matrix(rnorm(n * 3), n))
  design <- rsa_design(~ template, list(template = template), block_var = rep(1:2, each = 3))

  for (method in c("pearson", "spearman")) {
    for (center in c("none", "stimulus_mean")) {
      model <- rsa_model(dataset, design, distmethod = method, regtype = "pearson",
                         pattern_center = center)
      D <- rsa_neural_rdm(X, method, center)
      neural <- D[lower.tri(D)][design$include]
      # distmethod shapes the neural RDM only. The second-order comparison is
      # controlled by regtype, here "pearson".
      expected <- stats::cor(neural, model$design$model_mat$template, method = "pearson")
      expect_equal(as.numeric(train_model(model, X, NULL, NULL)), as.numeric(expected),
                   tolerance = 1e-12)
      expect_equal(length(D[lower.tri(D)]), choose(n, 2))
      expect_true(length(neural) < choose(n, 2))
    }
  }
})
