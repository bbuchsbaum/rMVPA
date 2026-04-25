test_that("kmeans-anchor connectivity returns k anchors with shape n_centers x k", {
  set.seed(2026)
  n <- 60L
  rank_d <- 4L
  centers <- matrix(rnorm(3 * rank_d, sd = 2), 3, rank_d)
  cluster_idx <- sample.int(3L, n, replace = TRUE)
  scores <- centers[cluster_idx, , drop = FALSE] +
            matrix(rnorm(n * rank_d, sd = 0.2), n, rank_d)
  rownames(scores) <- as.character(seq_len(n))
  colnames(scores) <- paste0("PC", seq_len(rank_d))

  fake_sl <- structure(list(results = list(), n_voxels = n,
                            active_voxels = n, metrics = character()),
                       class = c("searchlight_result", "list"))
  attr(fake_sl, "fingerprints") <- list(ids = seq_len(n), scores = scores,
                                        dataset = NULL)

  conn <- model_space_connectivity(fake_sl, k = 3L, build_maps = FALSE,
                                   random_seed = 42L)
  expect_s3_class(conn, "model_space_anchor_connectivity")
  expect_equal(conn$method, "kmeans_anchors")
  expect_equal(conn$anchor_method, "kmeans")
  expect_equal(conn$k, 3L)
  expect_equal(length(conn$anchors), 3L)
  expect_equal(dim(conn$similarity), c(n, 3L))
  expect_equal(rownames(conn$similarity), as.character(seq_len(n)))
})

test_that("anchor self-similarity is the maximum in its column for cosine scale", {
  set.seed(1)
  n <- 40L
  scores <- matrix(rnorm(n * 5), n, 5)
  fake_sl <- structure(list(results = list()), class = c("searchlight_result", "list"))
  attr(fake_sl, "fingerprints") <- list(ids = seq_len(n), scores = scores,
                                        dataset = NULL)

  conn <- model_space_connectivity(fake_sl, k = 4L, scale = "norm",
                                   build_maps = FALSE, random_seed = 7L)

  # Each anchor should be the row index in `ids` matching the entry of
  # `anchors`. With cosine similarity (unit-norm rows), self-similarity = 1.
  for (j in seq_along(conn$anchors)) {
    self_pos <- which(as.integer(rownames(conn$similarity)) == conn$anchors[j])
    expect_equal(conn$similarity[self_pos, j], 1, tolerance = 1e-12)
    expect_lte(max(conn$similarity[, j]) - 1, 1e-12)
  }
})

test_that("k larger than n_centers is clamped without crashing", {
  scores <- matrix(rnorm(15), 5, 3)
  fake_sl <- structure(list(results = list()), class = c("searchlight_result", "list"))
  attr(fake_sl, "fingerprints") <- list(ids = seq_len(5), scores = scores,
                                        dataset = NULL)
  conn <- model_space_connectivity(fake_sl, k = 100L, build_maps = FALSE,
                                   random_seed = 1L)
  expect_lte(conn$k, 5L)
  expect_equal(dim(conn$similarity), c(5L, conn$k))
})

test_that("explicit seeds bypass kmeans and produce shape n_centers x length(seeds)", {
  scores <- matrix(rnorm(40), 8, 5)
  ids <- 11:18
  fake_sl <- structure(list(results = list()), class = c("searchlight_result", "list"))
  attr(fake_sl, "fingerprints") <- list(ids = ids, scores = scores, dataset = NULL)

  conn <- model_space_connectivity(fake_sl, seeds = c(12L, 17L),
                                   build_maps = FALSE)
  expect_equal(conn$method, "explicit_anchors")
  expect_equal(conn$anchor_method, "explicit")
  expect_equal(conn$anchors, c(12L, 17L))
  expect_equal(dim(conn$similarity), c(8L, 2L))
})

test_that("explicit anchor similarities match raw and cosine oracle calculations", {
  set.seed(221)
  scores <- matrix(rnorm(35), 7, 5)
  ids <- 101:107
  seeds <- c(102L, 106L)
  anchor_pos <- match(seeds, ids)
  fake_sl <- structure(list(results = list()), class = c("searchlight_result", "list"))
  attr(fake_sl, "fingerprints") <- list(ids = ids, scores = scores, dataset = NULL)

  raw <- model_space_connectivity(fake_sl, seeds = seeds, scale = "raw",
                                  build_maps = FALSE)
  expect_equal(unname(raw$similarity),
               unname(scores %*% t(scores[anchor_pos, , drop = FALSE])),
               tolerance = 1e-12)

  norms <- sqrt(rowSums(scores^2))
  scores_unit <- scores / norms
  norm <- model_space_connectivity(fake_sl, seeds = seeds, scale = "norm",
                                   build_maps = FALSE)
  expect_equal(unname(norm$similarity),
               unname(scores_unit %*% t(scores_unit[anchor_pos, , drop = FALSE])),
               tolerance = 1e-12)
})

test_that("explicit anchor connectivity is equivariant to searchlight row order", {
  set.seed(303)
  scores <- matrix(rnorm(54), 9, 6)
  ids <- 501:509
  seeds <- c(502L, 507L, 509L)
  fake_sl <- structure(list(results = list()), class = c("searchlight_result", "list"))
  attr(fake_sl, "fingerprints") <- list(ids = ids, scores = scores, dataset = NULL)

  perm <- c(4L, 9L, 1L, 7L, 2L, 8L, 3L, 5L, 6L)
  fake_perm <- structure(list(results = list()), class = c("searchlight_result", "list"))
  attr(fake_perm, "fingerprints") <- list(ids = ids[perm],
                                          scores = scores[perm, , drop = FALSE],
                                          dataset = NULL)

  conn <- model_space_connectivity(fake_sl, seeds = seeds, build_maps = FALSE)
  conn_perm <- model_space_connectivity(fake_perm, seeds = seeds, build_maps = FALSE)

  expect_equal(conn_perm$similarity[rownames(conn$similarity),
                                    colnames(conn$similarity)],
               conn$similarity, tolerance = 1e-12)
  expect_equal(conn_perm$anchors, conn$anchors)
})

test_that("kmeans anchor ids are nearest observed fingerprints to centroids", {
  set.seed(404)
  scores <- matrix(rnorm(120), 30, 4)
  fake_sl <- structure(list(results = list()), class = c("searchlight_result", "list"))
  attr(fake_sl, "fingerprints") <- list(ids = seq_len(30), scores = scores,
                                        dataset = NULL)

  conn <- model_space_connectivity(fake_sl, k = 5L, random_seed = 12L,
                                   build_maps = FALSE)

  norms <- sqrt(rowSums(scores^2))
  F_use <- scores / norms
  for (j in seq_len(conn$k)) {
    members <- which(conn$cluster_id == j)
    d <- rowSums((F_use[members, , drop = FALSE] -
                    matrix(conn$centroids[j, ], length(members), ncol(F_use),
                           byrow = TRUE))^2)
    expect_equal(conn$anchor_pos[j], members[which.min(d)])
  }
})

test_that("anchor connectivity handles non-finite and zero fingerprints explicitly", {
  scores <- rbind(
    c(1, 0, 0),
    c(0, 0, 0),
    c(NA, 1, 0),
    c(0, 1, 0),
    c(0, 0, 1)
  )
  fake_sl <- structure(list(results = list()), class = c("searchlight_result", "list"))
  attr(fake_sl, "fingerprints") <- list(ids = 11:15, scores = scores,
                                        dataset = NULL)

  conn <- model_space_connectivity(fake_sl, k = 10L, build_maps = FALSE)
  expect_equal(conn$n_centers, 4L)
  expect_false("13" %in% rownames(conn$similarity))
  expect_true(all(is.finite(conn$similarity)))
  zero_row <- which(rownames(conn$similarity) == "12")
  expect_true(all(conn$similarity[zero_row, ] == 0))
})

test_that("searchlight anchor connectivity validates API edge cases", {
  scores <- matrix(rnorm(40), 8, 5)
  ids <- 11:18
  fake_sl <- structure(list(results = list()), class = c("searchlight_result", "list"))
  attr(fake_sl, "fingerprints") <- list(ids = ids, scores = scores, dataset = NULL)

  expect_error(model_space_connectivity(fake_sl, model_rdms = matrix(1, 2, 2)),
               "`model_rdms` is not used")
  expect_error(model_space_connectivity(fake_sl, k = NA_integer_, build_maps = FALSE),
               "`k` must be a positive integer scalar")
  expect_error(model_space_connectivity(fake_sl, k = 1.5, build_maps = FALSE),
               "`k` must be a positive integer scalar")
  expect_error(model_space_connectivity(fake_sl, nstart = 0L, build_maps = FALSE),
               "`nstart` must be a positive integer scalar")
  expect_error(model_space_connectivity(fake_sl, iter.max = 0L, build_maps = FALSE),
               "`iter.max` must be a positive integer scalar")
  expect_error(model_space_connectivity(fake_sl, seeds = c(12L, 12L),
                                        build_maps = FALSE),
               "must be unique")
  expect_error(model_space_connectivity(fake_sl, seeds = c(12L, NA_integer_),
                                        build_maps = FALSE),
               "non-missing")
})

test_that("random_seed does not perturb caller RNG state", {
  set.seed(20260425)
  scores <- matrix(rnorm(60), 15, 4)
  old_seed <- .Random.seed
  expected_next <- stats::runif(4)
  assign(".Random.seed", old_seed, envir = .GlobalEnv)

  fake_sl <- structure(list(results = list()), class = c("searchlight_result", "list"))
  attr(fake_sl, "fingerprints") <- list(ids = seq_len(15), scores = scores,
                                        dataset = NULL)

  invisible(model_space_connectivity(fake_sl, k = 3L, random_seed = 99L,
                                     build_maps = FALSE))
  expect_equal(stats::runif(4), expected_next)
})

test_that("missing fingerprints produce a clear error", {
  fake_sl <- structure(list(results = list()), class = c("searchlight_result", "list"))
  expect_error(model_space_connectivity(fake_sl), "no stored fingerprints")
})

test_that("end-to-end run_searchlight captures fingerprints when return_fingerprint=TRUE", {
  skip_if_not_installed("neuroim2")
  set.seed(909)
  n <- 14L
  arr  <- array(rnorm(prod(c(4, 4, 4, n))), c(4, 4, 4, n))
  sp   <- neuroim2::NeuroSpace(c(4, 4, 4, n))
  vec  <- neuroim2::NeuroVec(arr, sp)
  mask <- neuroim2::LogicalNeuroVol(array(1, c(4, 4, 4)),
                                    neuroim2::NeuroSpace(c(4, 4, 4)))
  ds   <- mvpa_dataset(train_data = vec, mask = mask)

  R1 <- as.matrix(stats::dist(matrix(rnorm(n * 4), n, 4)))
  R2 <- as.matrix(stats::dist(matrix(rnorm(n * 4), n, 4)))
  des <- pair_rsa_design(items_a = paste0("it", seq_len(n)),
                         model = list(R1 = R1, R2 = R2))
  mdl <- rsa_model(ds, des, distmethod = "pearson", regtype = "pearson",
                   return_fingerprint = TRUE)

  res <- suppressWarnings(run_searchlight(mdl, radius = 6, method = "standard",
                                          verbose = FALSE))
  fp <- attr(res, "fingerprints")
  expect_false(is.null(fp))
  expect_true(nrow(fp$scores) > 0)

  conn <- model_space_connectivity(res, k = 3L, build_maps = FALSE,
                                   random_seed = 11L)
  expect_s3_class(conn, "model_space_anchor_connectivity")
  expect_equal(dim(conn$similarity), c(nrow(fp$scores), conn$k))
})
