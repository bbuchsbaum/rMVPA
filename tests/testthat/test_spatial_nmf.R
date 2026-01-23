library(Matrix)
testthat::skip_if_not_installed("neuroim2")
testthat::skip_if_not_installed("neurosurf")
library(neuroim2)
library(neurosurf)

test_that("spatial_nmf_fit returns nonnegative factors with correct dimensions", {
  set.seed(1)
  X <- matrix(rexp(50 * 30), nrow = 50, ncol = 30)

  fit <- rMVPA:::spatial_nmf_fit(X, k = 5, max_iter = 40, check_every = 10)

  expect_equal(dim(fit$W), c(50, 5))
  expect_equal(dim(fit$H), c(5, 30))
  expect_true(all(fit$W >= 0))
  expect_true(all(fit$H >= 0))
})

test_that("spatial_nmf_fit decreases objective for rank-1 data", {
  set.seed(10)
  a <- runif(20, 0.2, 1)
  b <- runif(15, 0.2, 1)
  X <- a %*% t(b)

  fit <- rMVPA:::spatial_nmf_fit(
    X,
    k = 1,
    init = "nndsvd",
    max_iter = 80,
    check_every = 1
  )

  obj <- fit$objective
  expect_true(length(obj) > 2)
  diffs <- diff(obj)
  tol <- 1e-6 * pmax(1, obj[-length(obj)])
  expect_true(all(diffs <= tol))

  rel_err <- sqrt(sum((X - fit$W %*% fit$H)^2) / sum(X^2))
  expect_lt(rel_err, 1e-3)
})

test_that("spatial_nmf_fit handles normalization of H", {
  set.seed(11)
  X <- matrix(rexp(12 * 9), nrow = 12, ncol = 9)

  fit <- rMVPA:::spatial_nmf_fit(
    X,
    k = 3,
    normalize = "H",
    max_iter = 30,
    check_every = 5
  )

  row_sums <- rowSums(fit$H)
  expect_true(all(abs(row_sums - 1) < 1e-6))
})

test_that("spatial_nmf_fit leaves exact factorization unchanged with fixed init", {
  set.seed(12)
  W_true <- matrix(runif(10 * 2, 0.2, 1), nrow = 10, ncol = 2)
  H_true <- matrix(runif(2 * 7, 0.2, 1), nrow = 2, ncol = 7)
  X <- W_true %*% H_true

  fit <- rMVPA:::spatial_nmf_fit(
    X,
    k = 2,
    W_init = W_true,
    H_init = H_true,
    max_iter = 1,
    min_iter = 0,
    check_every = 1,
    normalize = "none"
  )

  expect_true(isTRUE(all.equal(fit$W, W_true, tolerance = 1e-12)))
  expect_true(isTRUE(all.equal(fit$H, H_true, tolerance = 1e-12)))
})

test_that("spatial_nmf_fit fixed point holds under scaling ambiguity", {
  set.seed(12)
  W_true <- matrix(runif(10 * 2, 0.2, 1), nrow = 10, ncol = 2)
  H_true <- matrix(runif(2 * 7, 0.2, 1), nrow = 2, ncol = 7)
  X <- W_true %*% H_true

  scale_vec <- c(2.5, 0.4)
  S <- diag(scale_vec)
  W_scaled <- W_true %*% S
  H_scaled <- solve(S) %*% H_true

  fit <- rMVPA:::spatial_nmf_fit(
    X,
    k = 2,
    W_init = W_scaled,
    H_init = H_scaled,
    max_iter = 1,
    min_iter = 0,
    check_every = 1,
    normalize = "none"
  )

  expect_true(isTRUE(all.equal(fit$W, W_scaled, tolerance = 1e-12)))
  expect_true(isTRUE(all.equal(fit$H, H_scaled, tolerance = 1e-12)))
})

test_that("spatial_nmf_project is a fixed point when W_init is exact", {
  set.seed(13)
  W_true <- matrix(runif(8 * 3, 0.2, 1), nrow = 8, ncol = 3)
  H_true <- matrix(runif(3 * 5, 0.2, 1), nrow = 3, ncol = 5)
  X <- W_true %*% H_true

  proj <- rMVPA:::spatial_nmf_project(
    X,
    H = H_true,
    W_init = W_true,
    max_iter = 1,
    min_iter = 0,
    check_every = 1
  )

  expect_true(isTRUE(all.equal(proj$W, W_true, tolerance = 1e-12)))
})

test_that("spatial_nmf_fit fixed point holds under component permutation", {
  set.seed(13)
  W_true <- matrix(runif(8 * 3, 0.2, 1), nrow = 8, ncol = 3)
  H_true <- matrix(runif(3 * 5, 0.2, 1), nrow = 3, ncol = 5)
  X <- W_true %*% H_true
  perm <- c(3, 1, 2)

  W_perm <- W_true[, perm, drop = FALSE]
  H_perm <- H_true[perm, , drop = FALSE]

  fit <- rMVPA:::spatial_nmf_fit(
    X,
    k = 3,
    W_init = W_perm,
    H_init = H_perm,
    max_iter = 1,
    min_iter = 0,
    check_every = 1,
    normalize = "none"
  )

  expect_true(isTRUE(all.equal(fit$W, W_perm, tolerance = 1e-12)))
  expect_true(isTRUE(all.equal(fit$H, H_perm, tolerance = 1e-12)))
})

test_that("spatial_nmf_fit is equivariant to column permutations with matched init", {
  set.seed(14)
  X <- matrix(rexp(6 * 8), nrow = 6, ncol = 8)
  k <- 3
  W_init <- matrix(runif(6 * k), nrow = 6, ncol = k)
  H_init <- matrix(runif(k * 8), nrow = k, ncol = 8)
  perm <- sample(seq_len(ncol(X)))

  fit1 <- rMVPA:::spatial_nmf_fit(
    X, k = k,
    W_init = W_init,
    H_init = H_init,
    max_iter = 1,
    min_iter = 0,
    check_every = 1,
    normalize = "none"
  )

  fit2 <- rMVPA:::spatial_nmf_fit(
    X[, perm, drop = FALSE], k = k,
    W_init = W_init,
    H_init = H_init[, perm, drop = FALSE],
    max_iter = 1,
    min_iter = 0,
    check_every = 1,
    normalize = "none"
  )

  expect_true(isTRUE(all.equal(fit1$W, fit2$W, tolerance = 1e-10)))
  expect_true(isTRUE(all.equal(fit1$H[, perm, drop = FALSE], fit2$H, tolerance = 1e-10)))
})

test_that("build_voxel_adjacency creates expected edge count for small grid", {
  mask <- rep(TRUE, 4)
  A <- rMVPA:::build_voxel_adjacency(mask, dims = c(2, 2, 1), neighbors = 6)

  expect_true(inherits(A, "sparseMatrix"))
  expect_equal(Matrix::nnzero(A), 8)
})

test_that("build_voxel_adjacency respects 2D neighborhood degrees", {
  mask <- rep(TRUE, 9)
  A <- rMVPA:::build_voxel_adjacency(mask, dims = c(3, 3, 1), neighbors = 4)
  deg <- Matrix::rowSums(A)

  center_lin <- 2 + (2 - 1) * 3
  corner_lin <- 1 + (1 - 1) * 3
  edge_lin <- 2 + (1 - 1) * 3

  idx <- which(which(mask) == center_lin)
  expect_equal(deg[idx], 4)

  idx <- which(which(mask) == corner_lin)
  expect_equal(deg[idx], 2)

  idx <- which(which(mask) == edge_lin)
  expect_equal(deg[idx], 3)
})

test_that("build_voxel_adjacency returns empty for empty mask", {
  mask <- rep(FALSE, 5)
  A <- rMVPA:::build_voxel_adjacency(mask, dims = c(5, 1, 1), neighbors = 6)
  expect_equal(dim(A), c(0, 0))
})

test_that("build_voxel_adjacency skips masked-out neighbors", {
  mask <- c(TRUE, FALSE, TRUE)
  A <- rMVPA:::build_voxel_adjacency(mask, dims = c(3, 1, 1), neighbors = 6)
  expect_equal(Matrix::nnzero(A), 0)
})

test_that("build_graph_laplacian matches D - A", {
  mask <- rep(TRUE, 3)
  A <- rMVPA:::build_voxel_adjacency(mask, dims = c(3, 1, 1), neighbors = 6)
  graph <- rMVPA:::build_graph_laplacian(A)

  L <- graph$L
  deg <- graph$degree

  expect_true(isTRUE(all.equal(Matrix::diag(L), deg)))
  expect_equal(as.numeric(L[1, 2]), -1)
  expect_equal(as.numeric(L[1, 3]), 0)
  expect_true(isTRUE(all.equal(Matrix::rowSums(L), rep(0, 3))))
})

test_that("graph regularization reduces smoothness penalty for smooth data", {
  set.seed(2)
  base <- runif(6, 0.3, 1)
  X <- base %*% t(rep(1, 10))

  mask <- rep(TRUE, 10)
  A <- rMVPA:::build_voxel_adjacency(mask, dims = c(10, 1, 1), neighbors = 6)
  graph <- rMVPA:::build_graph_laplacian(A)

  set.seed(3)
  W_init <- matrix(runif(6 * 2), nrow = 6, ncol = 2)
  H_init <- matrix(runif(2 * 10), nrow = 2, ncol = 10)
  smooth_init <- sum(H_init * (H_init * graph$degree - H_init %*% graph$A))

  fit <- rMVPA:::spatial_nmf_fit(X, k = 2, graph = graph, lambda = 10,
                                W_init = W_init, H_init = H_init,
                                max_iter = 100, check_every = 20)

  smooth_final <- sum(fit$H * (fit$H * graph$degree - fit$H %*% graph$A))
  expect_lte(smooth_final, smooth_init)
})

test_that("graph smoothness penalty is zero for constant component maps", {
  mask <- rep(TRUE, 6)
  A <- rMVPA:::build_voxel_adjacency(mask, dims = c(6, 1, 1), neighbors = 6)
  graph <- rMVPA:::build_graph_laplacian(A)

  H_const <- matrix(0.7, nrow = 2, ncol = 6)
  smooth <- sum(H_const * (H_const * graph$degree - H_const %*% graph$A))
  expect_lt(abs(smooth), 1e-12)
})

test_that("weighted graph smoothness scales with adjacency weights", {
  A <- Matrix::Matrix(c(0, 1, 0,
                        1, 0, 1,
                        0, 1, 0), nrow = 3, byrow = TRUE, sparse = TRUE)
  graph1 <- rMVPA:::.prepare_graph(list(A = A, weighted = TRUE), lambda = 1, p = 3)
  H <- matrix(runif(2 * 3, 0.2, 1), nrow = 2, ncol = 3)
  smooth1 <- sum(H * (H * graph1$degree - H %*% graph1$A))

  A2 <- 3 * A
  graph2 <- rMVPA:::.prepare_graph(list(A = A2, weighted = TRUE), lambda = 1, p = 3)
  smooth2 <- sum(H * (H * graph2$degree - H %*% graph2$A))
  expect_true(isTRUE(all.equal(smooth2, 3 * smooth1, tolerance = 1e-10)))
})

test_that("graph-regularized objective is non-increasing", {
  set.seed(4)
  X <- matrix(rexp(12 * 8), nrow = 12, ncol = 8)

  mask <- rep(TRUE, 8)
  A <- rMVPA:::build_voxel_adjacency(mask, dims = c(8, 1, 1), neighbors = 6)
  graph <- rMVPA:::build_graph_laplacian(A)

  fit <- rMVPA:::spatial_nmf_fit(
    X,
    k = 3,
    graph = graph,
    lambda = 2,
    init = "random",
    max_iter = 40,
    check_every = 1
  )

  obj <- fit$objective
  diffs <- diff(obj)
  tol <- 1e-6 * pmax(1, obj[-length(obj)])
  expect_true(all(diffs <= tol))
})

test_that("graph-regularized objective matches manual computation", {
  set.seed(7)
  X <- matrix(rexp(6 * 6), nrow = 6, ncol = 6)

  mask <- rep(TRUE, 6)
  A <- rMVPA:::build_voxel_adjacency(mask, dims = c(6, 1, 1), neighbors = 6)
  graph <- rMVPA:::build_graph_laplacian(A)

  fit <- rMVPA:::spatial_nmf_fit(
    X,
    k = 2,
    graph = graph,
    lambda = 1.5,
    max_iter = 30,
    check_every = 1
  )

  W <- fit$W
  H <- fit$H
  resid <- X - W %*% H
  smooth <- sum(H * (H * graph$degree - H %*% graph$A))
  manual_obj <- sum(resid * resid) + 1.5 * smooth

  expect_true(isTRUE(all.equal(tail(fit$objective, 1), manual_obj, tolerance = 1e-6)))
})

test_that("spatial_nmf_project recovers W for exact factorization", {
  set.seed(5)
  W_true <- matrix(runif(10 * 3, 0.2, 1), nrow = 10, ncol = 3)
  H_true <- matrix(runif(3 * 7, 0.2, 1), nrow = 3, ncol = 7)
  X <- W_true %*% H_true

  proj <- rMVPA:::spatial_nmf_project(
    X,
    H = H_true,
    W_init = W_true,
    max_iter = 200,
    check_every = 10
  )

  rel_err <- sqrt(sum((X - proj$W %*% H_true)^2) / sum(X^2))
  expect_lt(rel_err, 1e-8)
})

test_that("spatial_nmf_fit errors on invalid inputs", {
  X <- matrix(1, nrow = 3, ncol = 2)
  expect_error(rMVPA:::spatial_nmf_fit(-X, k = 1))
  expect_error(rMVPA:::spatial_nmf_fit(X, k = 5))
  expect_error(rMVPA:::spatial_nmf_fit(X, k = 1, lambda = -1))
  expect_error(rMVPA:::spatial_nmf_fit(X, k = 1, check_every = 0))
  expect_error(rMVPA:::spatial_nmf_fit(X, k = 1, tol = 0))
  expect_error(rMVPA:::spatial_nmf_fit(X, k = 1, eps = 0))
  expect_error(rMVPA:::spatial_nmf_fit(X, k = 1, min_iter = 5, max_iter = 2))

  H_bad <- matrix(1, nrow = 2, ncol = 3)
  expect_error(rMVPA:::spatial_nmf_project(X, H_bad))

  H_neg <- matrix(-1, nrow = 2, ncol = 2)
  expect_error(rMVPA:::spatial_nmf_project(X, H_neg))
  expect_error(rMVPA:::spatial_nmf_project(X, matrix(1, nrow = 2, ncol = 2), check_every = 0))
})

test_that("spatial_nmf_fit validates graph dimensions", {
  set.seed(6)
  X <- matrix(rexp(5 * 4), nrow = 5, ncol = 4)
  bad_graph <- list(A = Matrix::Matrix(diag(3), sparse = TRUE))

  expect_error(rMVPA:::spatial_nmf_fit(X, k = 2, graph = bad_graph, lambda = 1))
})

test_that("prepare_graph preserves weights when requested", {
  A <- Matrix::Matrix(c(0, 2, 0,
                        2, 0, 3,
                        0, 3, 0), nrow = 3, byrow = TRUE, sparse = TRUE)
  expect_warning(
    info_bin <- rMVPA:::.prepare_graph(list(A = A), lambda = 1, p = 3),
    "binarizing"
  )
  expect_true(all(info_bin$A@x %in% c(0, 1)))

  info_wt <- rMVPA:::.prepare_graph(list(A = A, weighted = TRUE), lambda = 1, p = 3)
  expect_true(any(info_wt$A@x > 1))
  expect_true(isTRUE(info_wt$meta$weighted))
})

test_that("spatial_nmf_maps builds data matrix for volumetric maps", {
  mask <- neuroim2::NeuroVol(array(1, dim = c(2, 2, 1)), neuroim2::NeuroSpace(c(2, 2, 1)))
  space_obj <- neuroim2::space(mask)

  map1 <- neuroim2::NeuroVol(array(1:4, dim = c(2, 2, 1)), space_obj)
  map2 <- neuroim2::NeuroVol(data = c(10, 30, 40), space = space_obj, indices = c(1, 3, 4))

  res <- spatial_nmf_maps(
    group_A = list(map1, map2),
    mask = mask,
    k = 2,
    lambda = 0,
    return_data = TRUE,
    max_iter = 5
  )

  expect_equal(dim(res$data), c(2, 4))
  expect_equal(as.numeric(res$data[1, ]), 1:4)
  expect_equal(as.numeric(res$data[2, ]), c(10, 0, 30, 40))
})

test_that("spatial_nmf_maps enforces consistent indices across maps", {
  geom_file <- rmvpa_test_surface_geom_file()
  if (identical(geom_file, "")) {
    testthat::skip("surface geometry not available (no packaged test surface found)")
  }
  geom <- neurosurf::read_surf_geometry(geom_file)
  nvert <- length(neurosurf::nodes(geom))
  idx1 <- seq_len(min(20, nvert))
  idx2 <- c(idx1[-1], idx1[1])

  map1 <- neurosurf::NeuroSurface(geometry = geom, indices = idx1, data = runif(length(idx1)))
  map2 <- neurosurf::NeuroSurface(geometry = geom, indices = idx2, data = runif(length(idx2)))
  expect_error(
    spatial_nmf_maps(group_A = list(map1, map2), mask = rep(1, nvert), k = 1, lambda = 0, max_iter = 5),
    "inconsistent indices"
  )

  map3 <- neurosurf::NeuroSurface(geometry = geom, indices = seq_len(nvert), data = runif(nvert))
  expect_error(
    spatial_nmf_maps(group_A = list(map1, map3), mask = rep(1, nvert), k = 1, lambda = 0, max_iter = 5)
  )
})

test_that("spatial_nmf_maps handles NA values with na_action", {
  mask <- neuroim2::NeuroVol(array(1, dim = c(2, 2, 1)), neuroim2::NeuroSpace(c(2, 2, 1)))
  space_obj <- neuroim2::space(mask)

  map_na <- neuroim2::NeuroVol(array(c(1, NA, 2, 3), dim = c(2, 2, 1)), space_obj)

  res <- spatial_nmf_maps(
    group_A = list(map_na),
    mask = mask,
    k = 1,
    lambda = 0,
    na_action = "zero",
    return_data = TRUE,
    max_iter = 5
  )

  expect_equal(as.numeric(res$data[1, 2]), 0)
  expect_error(
    spatial_nmf_maps(
      group_A = list(map_na),
      mask = mask,
      k = 1,
      lambda = 0,
      na_action = "error",
      max_iter = 5
    )
  )
})

test_that("spatial_nmf_maps assigns group labels when group_B is provided", {
  mask <- neuroim2::NeuroVol(array(1, dim = c(2, 2, 1)), neuroim2::NeuroSpace(c(2, 2, 1)))
  space_obj <- neuroim2::space(mask)

  map1 <- neuroim2::NeuroVol(array(1:4, dim = c(2, 2, 1)), space_obj)
  map2 <- neuroim2::NeuroVol(array(4:1, dim = c(2, 2, 1)), space_obj)

  res <- spatial_nmf_maps(
    group_A = list(map1),
    group_B = list(map2),
    mask = mask,
    k = 1,
    lambda = 0,
    max_iter = 5
  )

  expect_equal(as.character(res$groups), c("A", "B"))
})

test_that("spatial_nmf_maps supports surface maps and requires graph for lambda > 0", {
  geom_file <- rmvpa_test_surface_geom_file()
  if (identical(geom_file, "")) {
    testthat::skip("surface geometry not available (no packaged test surface found)")
  }
  geom <- neurosurf::read_surf_geometry(geom_file)
  nvert <- length(neurosurf::nodes(geom))

  map1 <- neurosurf::NeuroSurface(geometry = geom, indices = seq_len(nvert), data = runif(nvert))
  map2 <- neurosurf::NeuroSurface(geometry = geom, indices = seq_len(nvert), data = runif(nvert))

  mask <- rep(1, nvert)
  mask[1:5] <- 0

  res <- spatial_nmf_maps(
    group_A = list(map1, map2),
    mask = mask,
    k = 2,
    lambda = 0,
    return_data = TRUE,
    max_iter = 5
  )

  expect_equal(ncol(res$data), sum(mask != 0))
  expect_true(inherits(res$components[[1]], "NeuroSurface"))

  expect_error(
    spatial_nmf_maps(
      group_A = list(map1),
      mask = mask,
      k = 2,
      lambda = 1,
      max_iter = 5
    )
  )
})

test_that("spatial_nmf_maps requires mask for volumetric inputs", {
  map1 <- neuroim2::NeuroVol(array(1:4, dim = c(2, 2, 1)), neuroim2::NeuroSpace(c(2, 2, 1)))
  expect_error(spatial_nmf_maps(group_A = list(map1), k = 1, lambda = 0))
})

test_that("spatial_nmf_maps errors on mismatched volume dimensions", {
  mask <- neuroim2::NeuroVol(array(1, dim = c(2, 2, 1)), neuroim2::NeuroSpace(c(2, 2, 1)))
  map_bad <- neuroim2::NeuroVol(array(1:3, dim = c(3, 1, 1)), neuroim2::NeuroSpace(c(3, 1, 1)))

  expect_error(spatial_nmf_maps(group_A = list(map_bad), mask = mask, k = 1, lambda = 0))
})
