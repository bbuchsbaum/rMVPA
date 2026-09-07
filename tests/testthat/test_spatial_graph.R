library(testthat)

# Reference adjacency: brute-force neighbour test on mask coordinates.
.ref_adjacency <- function(mask_vol, neighbors = 6) {
  dims <- dim(mask_vol)
  idx <- which(as.numeric(neuroim2::values(mask_vol)) > 0)
  coords <- arrayInd(idx, dims)
  n <- length(idx)
  A <- matrix(0, n, n)
  for (a in seq_len(n)) {
    d <- abs(sweep(coords, 2, coords[a, ], "-"))
    within <- apply(d, 1, max) == 1
    nonzero <- rowSums(d)
    ok <- switch(as.character(neighbors),
                 "6"  = within & nonzero == 1,
                 "18" = within & nonzero <= 2,
                 "26" = within)
    A[a, ok] <- 1
  }
  list(A = A, idx = idx)
}

test_that("volume graph vertices align with get_feature_matrix columns", {
  ds <- gen_sample_dataset(c(5, 4, 3), 6)
  g <- spatial_graph(ds$dataset)
  X <- get_feature_matrix(ds$dataset)

  expect_s3_class(g, "spatial_graph")
  expect_equal(g$n_features, ncol(X))
  expect_equal(g$feature_ids, which(ds$dataset$mask > 0))
  expect_equal(g$domain_type, "volume")
  expect_true(grepl("^volume:5x4x3", g$geometry_id))
  expect_true(Matrix::isSymmetric(g$A))
  expect_equal(as.numeric(Matrix::diag(g$A)), rep(0, g$n_features))
  expect_equal(g$degree, as.numeric(Matrix::rowSums(g$A)))
  expect_equal(as.numeric(Matrix::rowSums(g$L)), rep(0, g$n_features), tolerance = 1e-12)
  expect_true(rMVPA:::.assert_graph_aligned(g, ncol(X)))
  expect_error(rMVPA:::.assert_graph_aligned(g, ncol(X) + 1), "columns")
})

test_that("volume graph matches a brute-force reference for 6, 18, and 26 neighbours", {
  set.seed(4)
  ds <- gen_sample_dataset(c(4, 4, 4), 5)
  mask <- ds$dataset$mask
  for (nb in c(6, 18, 26)) {
    g <- spatial_graph(ds$dataset, neighbors = nb)
    ref <- .ref_adjacency(mask, nb)
    expect_equal(as.matrix(g$A), ref$A, ignore_attr = TRUE, info = paste("neighbors", nb))
    expect_equal(g$feature_ids, ref$idx)
  }
})

test_that("permuting features together with the graph preserves adjacency", {
  ds <- gen_sample_dataset(c(4, 4, 4), 5)
  g <- spatial_graph(ds$dataset)
  set.seed(5)
  perm <- sample(g$n_features)
  g_perm <- restrict_graph(g, perm)
  expect_equal(g_perm$feature_ids, g$feature_ids[perm])
  expect_equal(as.matrix(g_perm$A), as.matrix(g$A)[perm, perm], ignore_attr = TRUE)
  expect_equal(g_perm$parent_index, perm)
  # edges relabelled consistently: each permuted edge maps to an original edge
  e <- graph_edges(g_perm)
  orig <- as.matrix(g$A)
  expect_true(all(orig[cbind(perm[e$from], perm[e$to])] == 1))
})

test_that("restrict_graph keeps the induced subgraph and records parent positions", {
  ds <- gen_sample_dataset(c(4, 4, 4), 5)
  g <- spatial_graph(ds$dataset)
  keep <- rep(FALSE, g$n_features); keep[seq(1, g$n_features, by = 2)] <- TRUE
  g2 <- restrict_graph(g, keep)
  expect_equal(g2$n_features, sum(keep))
  expect_equal(g2$feature_ids, g$feature_ids[keep])
  expect_equal(as.matrix(g2$A), as.matrix(g$A)[keep, keep], ignore_attr = TRUE)
  expect_equal(g2$parent_index, which(keep))
  expect_equal(g2$geometry_id, g$geometry_id)
  # nested restriction composes parent positions
  g3 <- restrict_graph(g2, 1:3)
  expect_equal(g3$parent_index, which(keep)[1:3])
  expect_error(restrict_graph(g, c(1, 1)), "unique")
  expect_error(restrict_graph(g, rep(TRUE, 3)), "one entry per feature")
})

test_that("graph_edges returns each undirected edge once", {
  ds <- gen_sample_dataset(c(4, 4, 4), 5)
  g <- spatial_graph(ds$dataset)
  e <- graph_edges(g)
  expect_equal(nrow(e), g$n_edges)
  expect_true(all(e$from < e$to))
  expect_equal(nrow(e), Matrix::nnzero(g$A) / 2)
})

test_that("raw adjacency input is symmetrized, binarized, and self-loop free", {
  A <- matrix(0, 4, 4)
  A[1, 2] <- 3; A[2, 3] <- 1; A[3, 3] <- 5
  g <- spatial_graph(A, feature_ids = c(10L, 20L, 30L, 40L))
  expect_equal(as.matrix(g$A), rbind(c(0, 1, 0, 0), c(1, 0, 1, 0), c(0, 1, 0, 0), c(0, 0, 0, 0)),
               ignore_attr = TRUE)
  expect_equal(g$feature_ids, c(10L, 20L, 30L, 40L))
  expect_equal(g$domain_type, "custom")

  gw <- spatial_graph(list(A = A, weighted = TRUE))
  expect_equal(as.matrix(gw$A)[1, 2], 3)     # one-directional weight kept both ways
  expect_equal(as.matrix(gw$A)[2, 1], 3)
  expect_true(gw$weighted)
  # explicit arguments win over list fields
  gb <- spatial_graph(list(A = A, weighted = TRUE), weighted = FALSE)
  expect_false(gb$weighted)
  expect_equal(as.matrix(gb$A)[1, 2], 1)
  # negative weights are rejected before symmetrization can hide them
  B <- matrix(0, 3, 3); B[1, 2] <- 1; B[2, 1] <- -1
  expect_error(spatial_graph(list(A = B, weighted = TRUE)), "non-negative")
  expect_error(spatial_graph(list(A = A, feature_ids = 1:3)), "feature ids")
  expect_error(spatial_graph("nope"), "no method")
})

test_that("degenerate masks produce well-formed graphs", {
  ds <- gen_sample_dataset(c(4, 4, 4), 5)
  # single active voxel
  m1 <- ds$dataset$mask; m1[] <- 0; m1[2, 2, 2] <- 1
  ds1 <- ds$dataset; ds1$mask <- m1
  g1 <- spatial_graph(ds1)
  expect_equal(g1$n_features, 1L)
  expect_equal(g1$n_edges, 0L)
  expect_equal(g1$feature_ids, which(m1 > 0))
  # empty mask
  m0 <- ds$dataset$mask; m0[] <- 0
  ds0 <- ds$dataset; ds0$mask <- m0
  g0 <- spatial_graph(ds0)
  expect_equal(g0$n_features, 0L)
  expect_equal(nrow(graph_edges(g0)), 0L)
  # 2-D slab: at most 4 neighbours under the 6-neighbour rule
  ds2 <- gen_sample_dataset(c(6, 6, 1), 5)
  g2 <- spatial_graph(ds2$dataset)
  expect_equal(g2$n_features, ncol(get_feature_matrix(ds2$dataset)))
  expect_lte(max(g2$degree), 4)
})

test_that("spatial_graph object is accepted by the NMF graph preparation", {
  ds <- gen_sample_dataset(c(4, 4, 4), 5)
  g <- spatial_graph(ds$dataset)
  info <- rMVPA:::.prepare_graph(g, lambda = 0.5, p = g$n_features)
  expect_true(info$use_graph)
  expect_equal(as.matrix(info$A), as.matrix(g$A), ignore_attr = TRUE)
})

test_that("multibasis graph has one block per basis and aligns with the stacked feature matrix", {
  ds <- gen_sample_dataset(c(4, 4, 3), 6)
  mb <- mvpa_multibasis_dataset(
    train_data = list(ds$dataset$train_data, ds$dataset$train_data, ds$dataset$train_data),
    mask = ds$dataset$mask
  )
  X <- get_feature_matrix(mb)
  g <- spatial_graph(mb)
  n_vox <- sum(ds$dataset$mask > 0)
  expect_equal(g$n_features, ncol(X))
  expect_equal(g$n_features, 3 * n_vox)
  expect_equal(g$feature_ids, rep(which(ds$dataset$mask > 0), 3))
  expect_equal(g$basis, rep(1:3, each = n_vox))
  expect_equal(g$domain_type, "multibasis_volume")

  # column j of the stacked feature matrix is voxel feature_ids[j] of basis[j]
  for (j in c(1L, n_vox, n_vox + 1L, 2L * n_vox + 3L, 3L * n_vox)) {
    expect_equal(
      as.numeric(X[, j]),
      as.numeric(neuroim2::series(mb$train_data[[g$basis[j]]], g$feature_ids[j])),
      info = paste("column", j)
    )
  }

  # no edges cross basis channels by default
  e <- graph_edges(g)
  expect_true(all(g$basis[e$from] == g$basis[e$to]))
  single <- spatial_graph(ds$dataset)
  expect_equal(g$n_edges, 3 * single$n_edges)

  # connect_basis links the same voxel across channels only
  gc <- spatial_graph(mb, connect_basis = TRUE)
  ec <- graph_edges(gc)
  cross <- ec[gc$basis[ec$from] != gc$basis[ec$to], ]
  expect_equal(nrow(cross), n_vox * 3)              # choose(3, 2) pairs per voxel
  expect_true(all(gc$feature_ids[cross$from] == gc$feature_ids[cross$to]))
})

test_that("clustered graph connects parcels whose voxels touch", {
  set.seed(6)
  ds <- gen_clustered_sample_dataset(c(6, 6, 6), nobs = 8, K = 5)
  g <- spatial_graph(ds$dataset)
  K <- ncol(ds$dataset$train_data@ts)
  expect_equal(g$n_features, K)
  expect_equal(g$n_features, ncol(get_feature_matrix(ds$dataset)))
  expect_equal(g$domain_type, "cluster")
  expect_true(Matrix::isSymmetric(g$A))

  # reference: two clusters adjacent iff some voxel pair is adjacent
  cvol <- ds$dataset$train_data@cvol
  vox <- spatial_graph(list(A = rMVPA:::build_voxel_adjacency(cvol@mask, neighbors = 6)))
  cl <- as.integer(cvol@clusters)
  e <- graph_edges(vox)
  ref <- matrix(0, K, K)
  for (r in seq_len(nrow(e))) {
    a <- cl[e$from[r]]; b <- cl[e$to[r]]
    if (a != b) { ref[a, b] <- 1; ref[b, a] <- 1 }
  }
  expect_equal(as.matrix(g$A), ref, ignore_attr = TRUE)
})

test_that("clustered graph uses column positions, not cluster labels, for non-contiguous ids", {
  set.seed(8)
  ds <- gen_sample_dataset(c(6, 6, 6), 5)
  mask_vol <- ds$dataset$mask
  active <- which(mask_vol > 0)
  coords <- arrayInd(active, dim(mask_vol))
  # four spatial quadrants in x/y with deliberately non-contiguous labels
  quadrant <- 1L + (coords[, 1] > 3) + 2L * (coords[, 2] > 3)
  labels <- c(3L, 7L, 9L, 20L)[quadrant]
  lmask <- neuroim2::LogicalNeuroVol(as.array(mask_vol) > 0, neuroim2::space(mask_vol))
  cvol <- neuroim2::ClusteredNeuroVol(lmask, clusters = labels)
  cvec <- neuroim2::ClusteredNeuroVec(ds$dataset$train_data, cvol)
  cds <- structure(list(train_data = cvec, mask = mask_vol),
                   class = c("mvpa_clustered_dataset", "mvpa_dataset", "list"))

  g <- spatial_graph(cds)
  X <- get_feature_matrix(cds)
  expect_equal(g$n_features, ncol(X))
  expect_equal(g$feature_ids, c(3L, 7L, 9L, 20L))
  # column j is the mean series of cluster feature_ids[j]
  series <- neuroim2::series(ds$dataset$train_data, active)
  for (j in seq_len(4)) {
    expect_equal(as.numeric(X[, j]),
                 as.numeric(rowMeans(series[, labels == g$feature_ids[j], drop = FALSE])),
                 tolerance = 1e-10, info = paste("column", j))
  }
  # quadrants sharing a face are adjacent; diagonal quadrants are not
  A <- as.matrix(g$A)
  expect_equal(A[1, 2], 1); expect_equal(A[1, 3], 1); expect_equal(A[2, 4], 1); expect_equal(A[3, 4], 1)
  expect_equal(A[1, 4], 0); expect_equal(A[2, 3], 0)
})

test_that("surface graph restricts mesh adjacency to masked nodes", {
  skip_if_not_installed("neurosurf")
  fname <- rmvpa_test_surface_geom_file()
  skip_if(identical(fname, ""), "no test surface geometry available")
  geom <- neurosurf::read_surf_geometry(fname)
  n_nodes <- length(neurosurf::nodes(geom))
  set.seed(7)
  mat <- matrix(rnorm(n_nodes * 5), n_nodes, 5)
  nsv <- neurosurf::NeuroSurfaceVector(geom, indices = seq_len(n_nodes), mat = mat)
  keep <- sort(sample(n_nodes, 150))          # scattered, non-contiguous node mask
  mask <- numeric(n_nodes); mask[keep] <- 1
  ds <- mvpa_surface_dataset(nsv, mask = mask)

  g <- spatial_graph(ds)
  X <- get_feature_matrix(ds)
  expect_equal(g$n_features, ncol(X))
  expect_equal(g$feature_ids, keep)
  expect_equal(g$domain_type, "surface")
  # column j of the feature matrix is node feature_ids[j]
  expect_equal(unname(X[, 10]), unname(mat[keep[10], ]))
  expect_equal(unname(X[, 150]), unname(mat[keep[150], ]))
  A_full <- neurosurf::adjacency(geom)
  expect_equal(as.matrix(g$A), (as.matrix(A_full[keep, keep]) != 0) * 1, ignore_attr = TRUE)
})

test_that("spatial_graph prints", {
  ds <- gen_sample_dataset(c(4, 4, 4), 5)
  expect_output(print(spatial_graph(ds$dataset)), "spatial_graph \\[volume\\]")
})
