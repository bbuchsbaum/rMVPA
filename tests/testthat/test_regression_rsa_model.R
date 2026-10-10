# Regression tests for two RSA correctness fixes.
#
# 1. `regtype` selects how the neural RDM is compared with each model RDM
#    (Pearson or Spearman). `distmethod` only shapes the neural RDM. Before the
#    fix, both regtypes used distmethod, so regtype = "pearson" and
#    regtype = "spearman" gave the same result.
# 2. `rsa_design(block_var = )` accepts character, factor and numeric block
#    labels. Before the fix, a character block made `include` and the model
#    matrix NA, and a factor in `data` was rejected.
#
# The oracles below are computed by hand from the raw patterns with base R
# (cor() and outer()). They do not call the package's RDM, ranking or kernel
# code. Sphere membership for the searchlight oracle is recomputed from voxel
# coordinates.
testthat::skip_if_not_installed("neuroim2")

rsa_regtype_fixture <- function(seed = 2601) {
  set.seed(seed)
  dims <- c(5L, 5L, 5L)
  nobs <- 20L
  block <- rep(1:4, each = 5L)
  # Latent item positions. The model RDM is a monotone but non-linear function
  # of their distances, and the voxel patterns carry the same geometry.
  u <- matrix(stats::rnorm(nobs * 2L), nobs)
  W <- matrix(stats::rnorm(2L * prod(dims)), 2L)
  X <- u %*% W + matrix(stats::rnorm(nobs * prod(dims), sd = 0.5), nobs)
  arr <- array(0, c(dims, nobs))
  for (t in seq_len(nobs)) arr[, , , t] <- array(X[t, ], dims)
  vec <- neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(c(dims, nobs)))
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, dims), neuroim2::NeuroSpace(dims))
  ds <- mvpa_dataset(vec, mask = mask)
  D0 <- as.matrix(stats::dist(u))
  model <- stats::as.dist(exp(D0) - 1)
  list(ds = ds, vec = vec, dims = dims, nobs = nobs, block = block, model = model,
       patterns = as.matrix(neuroim2::series(vec, seq_len(prod(dims)))))
}

# Hand-computed RSA statistic for one set of patterns (rows = observations).
rsa_oracle <- function(X, model, block, distmethod, regtype, keep_intra = FALSE) {
  nd <- 1 - stats::cor(t(X), method = distmethod)
  lt <- lower.tri(nd)
  neural <- nd[lt]
  model_vec <- as.vector(model)
  b <- as.integer(factor(block))
  keep <- if (keep_intra) rep(TRUE, length(neural)) else outer(b, b, "!=")[lt]
  stats::cor(neural[keep], model_vec[keep], method = regtype)
}

# TRUE when two included neural distances are equal to 12 digits. Spearman
# neural distances are discrete (integer sums of squared rank differences), so
# small patterns tie often. Then a Spearman comparison depends on rounding of
# the last bits, and the exact oracle is not defined.
rsa_has_ties <- function(X, distmethod, block) {
  nd <- 1 - stats::cor(t(X), method = distmethod)
  lt <- lower.tri(nd)
  b <- as.integer(factor(block))
  anyDuplicated(round(nd[lt][outer(b, b, "!=")[lt]], 12)) > 0
}

# Sphere of radius r (voxel units) around every grid voxel, within the full mask.
rsa_sphere_oracle <- function(fx, distmethod, regtype, radius = 1, keep_intra = FALSE) {
  coords <- arrayInd(seq_len(prod(fx$dims)), fx$dims)
  vapply(seq_len(prod(fx$dims)), function(cn) {
    d2 <- rowSums((coords - matrix(coords[cn, ], nrow(coords), 3L, byrow = TRUE))^2)
    sel <- which(d2 <= radius^2)
    rsa_oracle(fx$patterns[, sel, drop = FALSE], fx$model, fx$block,
               distmethod, regtype, keep_intra)
  }, numeric(1))
}

rsa_regtype_model <- function(fx, distmethod, regtype, keep_intra_run = FALSE) {
  rdes <- rsa_design(~ model, list(model = fx$model, block = fx$block),
                     block_var = "block", keep_intra_run = keep_intra_run)
  rsa_model(fx$ds, rdes, distmethod = distmethod, regtype = regtype,
            check_collinearity = FALSE)
}

rsa_combos <- expand.grid(distmethod = c("pearson", "spearman"),
                          regtype = c("pearson", "spearman"),
                          stringsAsFactors = FALSE)

test_that("regtype pearson and spearman differ on a non-linear monotone relation", {
  fx <- rsa_regtype_fixture()
  p <- rsa_regtype_model(fx, "spearman", "pearson")
  s <- rsa_regtype_model(fx, "spearman", "spearman")
  rp <- suppressWarnings(train_model(p, fx$patterns, NULL, seq_len(prod(fx$dims))))
  rs <- suppressWarnings(train_model(s, fx$patterns, NULL, seq_len(prod(fx$dims))))
  expect_gt(abs(rp[["model"]] - rs[["model"]]), 1e-3)
})

test_that("train_model equals the hand-computed statistic for every distmethod/regtype", {
  fx <- rsa_regtype_fixture()
  for (i in seq_len(nrow(rsa_combos))) {
    dm <- rsa_combos$distmethod[i]
    rt <- rsa_combos$regtype[i]
    info <- paste(dm, rt)
    ms <- rsa_regtype_model(fx, dm, rt)
    if (rt == "spearman") {
      expect_false(rsa_has_ties(fx$patterns, dm, fx$block), info = info)
    }
    got <- suppressWarnings(train_model(ms, fx$patterns, NULL, seq_len(prod(fx$dims))))
    expected <- rsa_oracle(fx$patterns, fx$model, fx$block, dm, rt)
    expect_equal(unname(got[["model"]]), expected, tolerance = 1e-10, info = info)

    ref <- rMVPA:::.with_reference_paths("rsa_fast_kernel", {
      ms_ref <- rsa_regtype_model(fx, dm, rt)
      suppressWarnings(train_model(ms_ref, fx$patterns, NULL, seq_len(prod(fx$dims))))
    })
    expect_equal(unname(ref[["model"]]), expected, tolerance = 1e-10,
                 info = paste(info, "reference path"))
  }
})

test_that("include masks restrict the oracle to between-block pairs", {
  fx <- rsa_regtype_fixture(2602)
  ms <- rsa_regtype_model(fx, "spearman", "spearman")
  b <- as.integer(factor(fx$block))
  lt <- lower.tri(matrix(0, fx$nobs, fx$nobs))
  expect_identical(ms$design$include, outer(b, b, "!=")[lt])
  expected <- rsa_oracle(fx$patterns, fx$model, fx$block, "spearman", "spearman")
  all_pairs <- rsa_oracle(fx$patterns, fx$model, fx$block, "spearman", "spearman",
                          keep_intra = TRUE)
  expect_false(isTRUE(all.equal(expected, all_pairs)))
  got <- suppressWarnings(train_model(ms, fx$patterns, NULL, seq_len(prod(fx$dims))))
  expect_equal(unname(got[["model"]]), expected, tolerance = 1e-10)
})

test_that("run_regional equals the hand-computed statistic for every combination", {
  fx <- rsa_regtype_fixture(2603)
  roi_idx <- list(seq_len(60L), seq.int(61L, prod(fx$dims)))
  roi_arr <- array(0L, fx$dims)
  roi_arr[roi_idx[[1]]] <- 1L
  roi_arr[roi_idx[[2]]] <- 2L
  rois <- neuroim2::NeuroVol(roi_arr, neuroim2::NeuroSpace(fx$dims))
  for (i in seq_len(nrow(rsa_combos))) {
    dm <- rsa_combos$distmethod[i]
    rt <- rsa_combos$regtype[i]
    info <- paste(dm, rt)
    ms <- rsa_regtype_model(fx, dm, rt)
    res <- suppressWarnings(run_regional(ms, rois, verbose = FALSE))
    tab <- as.data.frame(res$performance_table)
    for (k in 1:2) {
      X_roi <- fx$patterns[, roi_idx[[k]], drop = FALSE]
      if (rt == "spearman") {
        expect_false(rsa_has_ties(X_roi, dm, fx$block), info = paste(info, "roi", k))
      }
      expected <- rsa_oracle(X_roi, fx$model, fx$block, dm, rt)
      row <- tab[tab$roinum == k, ]
      expect_equal(row$model, expected, tolerance = 1e-10, info = paste(info, "roi", k))
    }
    ref <- rMVPA:::.with_reference_paths("rsa_fast_kernel", {
      ms_ref <- rsa_regtype_model(fx, dm, rt)
      suppressWarnings(run_regional(ms_ref, rois, verbose = FALSE))
    })
    expect_equal(as.data.frame(ref$performance_table)$model, tab$model,
                 tolerance = 1e-10, info = paste(info, "reference path"))
  }
})

test_that("searchlight maps equal the hand-computed statistic and agree across engines", {
  fx <- rsa_regtype_fixture(2604)
  for (i in seq_len(nrow(rsa_combos))) {
    dm <- rsa_combos$distmethod[i]
    rt <- rsa_combos$regtype[i]
    info <- paste(dm, rt)
    ms <- rsa_regtype_model(fx, dm, rt)
    # Spearman/Spearman spheres contain tied neural distances (see
    # rsa_has_ties), so the package and the oracle break ties at the last bit
    # differently. That combination gets a 1e-3 tolerance. Every other
    # combination must match to 1e-10.
    tol <- if (dm == "spearman" && rt == "spearman") 1e-3 else 1e-10
    expected <- rsa_sphere_oracle(fx, dm, rt, radius = 2)
    fast <- suppressWarnings(run_searchlight(ms, radius = 2, engine = "rsa_fast",
                                             backend = "default"))
    legacy <- suppressWarnings(run_searchlight(ms, radius = 2, engine = "legacy",
                                               backend = "default"))
    expect_identical(attr(fast, "searchlight_engine"), "rsa_fast")
    map_fast <- as.numeric(neuroim2::values(fast$results$model))
    map_legacy <- as.numeric(neuroim2::values(legacy$results$model))
    expect_equal(map_fast, expected, tolerance = tol, info = paste(info, "rsa_fast"))
    expect_equal(map_legacy, expected, tolerance = tol, info = paste(info, "legacy"))
    expect_identical(map_fast, map_legacy, info = paste(info, "engines identical"))

    ref <- rMVPA:::.with_reference_paths("rsa_fast_kernel", {
      ms_ref <- rsa_regtype_model(fx, dm, rt)
      suppressWarnings(run_searchlight(ms_ref, radius = 2, engine = "rsa_fast",
                                       backend = "default"))
    })
    expect_equal(as.numeric(neuroim2::values(ref$results$model)), map_fast,
                 tolerance = 1e-10, info = paste(info, "reference path"))
  }
})

test_that("block_var character, factor and numeric give the same include mask", {
  lab <- rep(c("run_a", "run_b", "run_c", "run_d"), each = 5L)
  n <- length(lab)
  set.seed(2605)
  model <- stats::dist(stats::rnorm(n))
  b <- as.integer(factor(lab))
  lt <- lower.tri(matrix(0, n, n))
  oracle_inc <- outer(b, b, "!=")[lt]
  num <- match(lab, unique(lab))

  des_chr <- rsa_design(~ model, list(model = model, block = lab), block_var = "block")
  des_fac <- rsa_design(~ model, list(model = model, block = factor(lab)), block_var = "block")
  des_num <- rsa_design(~ model, list(model = model, block = num), block_var = "block")
  des_raw <- rsa_design(~ model, list(model = model), block_var = lab)
  des_fml <- rsa_design(~ model, list(model = model, block = lab), block_var = ~ block)

  expect_identical(des_chr$include, oracle_inc)
  expect_identical(des_fac$include, oracle_inc)
  expect_identical(des_num$include, oracle_inc)
  expect_identical(des_raw$include, oracle_inc)
  expect_identical(des_fml$include, oracle_inc)
  # Numeric block values keep the original rule: include = dist(block) != 0.
  expect_identical(des_num$include, as.vector(stats::dist(num)) != 0)

  for (des in list(des_chr, des_fac, des_num, des_raw, des_fml)) {
    expect_false(anyNA(des$include))
    expect_false(anyNA(unlist(des$model_mat)))
    expect_identical(des$model_mat$model, as.vector(model)[oracle_inc])
  }

  # The block labels still reach the model: fits run without error.
  ms <- rsa_model(mvpa_dataset(neuroim2::NeuroVec(array(stats::rnorm(3 * 3 * 3 * n), c(3, 3, 3, n)),
                                                  neuroim2::NeuroSpace(c(3, 3, 3, n))),
                               mask = neuroim2::LogicalNeuroVol(array(TRUE, c(3, 3, 3)),
                                                                neuroim2::NeuroSpace(c(3, 3, 3)))),
                  des_fac, check_collinearity = FALSE)
  expect_s3_class(ms, "rsa_model")
})

test_that("keep_intra_run = TRUE keeps every pair with a block_var", {
  lab <- rep(c("a", "b"), each = 4L)
  model <- stats::dist(seq_along(lab))
  des <- rsa_design(~ model, list(model = model, block = lab), block_var = "block",
                    keep_intra_run = TRUE)
  expect_null(des$include)
})

test_that("factor predictors and factor block terms in the formula are still rejected", {
  f <- factor(rep(c("x", "y"), 4L))
  expect_error(rsa_design(~ f, list(f = f)), "illegal variable type")
  expect_error(rsa_design(~ model, list(model = stats::dist(seq_len(8)), block = f),
                          nuisance = list(block = f), block_var = "block"),
               "illegal variable type")
})
