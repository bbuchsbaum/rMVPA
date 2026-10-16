## Regression tests for #130 (statistical maps must not read as significant
## outside the mask) and #131 (n_voxels reports the in-mask voxel count).

.oom_setup <- function() {
  set.seed(1)
  ds <- gen_sample_dataset(c(6, 6, 6), 20, blocks = 2, nlevels = 2)
  m <- array(0, c(6, 6, 6))
  m[2:5, 2:5, 2:5] <- 1 # 64 of 216 voxels
  mask <- neuroim2::NeuroVol(m, neuroim2::space(ds$dataset$mask))
  dset <- mvpa_dataset(ds$dataset$train_data, mask = mask)
  list(ds = ds, dset = dset, m = m, out = m == 0, in_ids = which(m > 0))
}

.map_vals <- function(x) as.numeric(as.array(x))

.vrsa_spec <- function(env, nperm = 5) {
  labs <- paste0("L", seq_len(20))
  D <- as.matrix(stats::dist(matrix(stats::rnorm(20 * 20), 20)))
  dimnames(D) <- list(labs, labs)
  rdes <- vector_rsa_design(D = D, labels = labs, block_var = env$ds$design$block_var)
  vector_rsa_model(env$dset, rdes, distfun = cordist(), nperm = nperm)
}

test_that("run_permutation_searchlight p-maps are NA outside the mask (#130)", {
  skip_if_not_installed("sda")
  env <- .oom_setup()
  futile.logger::flog.threshold(futile.logger::ERROR)
  cval <- blocked_cross_validation(env$ds$design$block_var)
  mspec <- mvpa_model(load_model("sda_notune"), env$dset, env$ds$design,
                      "classification", crossval = cval)
  pc <- permutation_control(n_perm = 5, subsample = 0.5, seed = 1L,
                            correction = "fdr", diagnose = FALSE)
  perm <- suppressWarnings(run_permutation_searchlight(mspec, radius = 2, perm_ctrl = pc))

  for (nm in c("p_map", "p_adj_map")) {
    v <- .map_vals(perm[[nm]])
    expect_equal(sum(v[env$out] < 0.05, na.rm = TRUE), 0, info = nm)
    expect_true(all(is.na(v[env$out])), info = nm)
  }
  # In-mask values are exactly the computed p-values.
  ids <- as.integer(perm$all_ids)
  expect_setequal(ids, env$in_ids)
  expect_equal(.map_vals(perm$p_map)[ids], unname(perm$p_values))
  expect_equal(.map_vals(perm$p_adj_map)[ids], unname(perm$p_adjusted))
})

test_that("p_/z_ metric maps from searchlights are NA outside the mask (#130)", {
  env <- .oom_setup()
  futile.logger::flog.threshold(futile.logger::ERROR)
  set.seed(2)
  vs <- .vrsa_spec(env)

  sl <- suppressWarnings(run_searchlight(vs, radius = 2, method = "standard"))
  for (nm in c("p_rsa_score", "z_rsa_score")) {
    v <- .map_vals(sl$results[[nm]])
    expect_true(all(is.na(v[env$out])), info = nm)
    expect_equal(sum(v[env$out] < 0.05, na.rm = TRUE), 0, info = nm)
    expect_true(all(!is.na(v[env$in_ids])), info = nm)
  }
  p_in <- .map_vals(sl$results$p_rsa_score)[env$in_ids]
  expect_true(all(p_in >= 0 & p_in <= 1))
  # Performance maps keep their historical zero fill.
  expect_true(all(.map_vals(sl$results$rsa_score)[env$out] == 0))

  slr <- suppressWarnings(run_searchlight(vs, radius = 2, method = "randomized", niter = 2))
  for (nm in c("p_rsa_score", "z_rsa_score")) {
    v <- .map_vals(slr$results[[nm]])
    expect_true(all(is.na(v[env$out])), info = nm)
    expect_equal(sum(v[env$out] < 0.05, na.rm = TRUE), 0, info = nm)
  }
})

test_that("regional p_/z_ maps are NA outside the regions and keep region values (#130)", {
  env <- .oom_setup()
  futile.logger::flog.threshold(futile.logger::ERROR)
  set.seed(3)
  vs <- .vrsa_spec(env)
  reg <- array(0, c(6, 6, 6))
  reg[2:5, 2:5, 2:3] <- 1
  reg[2:5, 2:5, 4:5] <- 2
  region_mask <- neuroim2::NeuroVol(reg, neuroim2::space(env$dset$mask))
  rr <- suppressWarnings(run_regional(vs, region_mask))

  pt <- rr$performance_table
  for (nm in c("p_rsa_score", "z_rsa_score")) {
    v <- .map_vals(rr$vol_results[[nm]])
    expect_true(all(is.na(v[reg == 0])), info = nm)
    expect_equal(sum(v[reg == 0] < 0.05, na.rm = TRUE), 0, info = nm)
    for (r in 1:2) {
      expect_true(all(v[reg == r] == pt[[nm]][pt$roinum == r]), info = nm)
    }
  }
})

test_that("region maps keep NA region values and apply the fill outside regions (#130)", {
  sp <- neuroim2::NeuroSpace(c(3, 3, 1))
  lab <- neuroim2::NeuroVol(array(c(0, 1, 1, 2, 2, 0, 3, 3, 0), c(3, 3, 1)), sp)
  out <- .map_region_values(lab, c(1, 2, 3), c(0.01, NA, 0.5), fill = NA_real_)
  expect_equal(.map_vals(out), c(NA, 0.01, 0.01, NA, NA, NA, 0.5, 0.5, NA))
  out0 <- .map_region_values(lab, c(1, 2, 3), c(0.7, 0.8, 0.9))
  expect_equal(.map_vals(out0), c(0, 0.7, 0.7, 0.8, 0.8, 0, 0.9, 0.9, 0))

  # Clustered datasets broadcast through the same helper.
  cds <- structure(list(region_mask = lab), class = c("mvpa_clustered_dataset", "mvpa_dataset", "list"))
  pmap <- build_output_map(cds, c(0.01, 0.02, 0.03), 1:3, fill = NA_real_)
  v <- .map_vals(pmap)
  expect_true(all(is.na(v[c(1, 6, 9)])))
  expect_equal(v[c(2, 3, 4, 5, 7, 8)], c(0.01, 0.01, 0.02, 0.02, 0.03, 0.03))
})

test_that("build_volume_map fill applies to sparse and dense inputs (#130)", {
  env <- .oom_setup()
  ids <- env$in_ids
  vals <- seq_along(ids) / 1000
  v0 <- .map_vals(build_volume_map(env$dset, vals, ids))
  expect_true(all(v0[env$out] == 0))
  vna <- .map_vals(build_volume_map(env$dset, vals, ids, fill = NA_real_))
  expect_true(all(is.na(vna[env$out])))
  expect_equal(vna[ids], vals)

  dense <- numeric(216)
  dense[ids] <- vals
  vd <- .map_vals(build_volume_map(env$dset, dense, NULL, fill = NA_real_))
  expect_true(all(is.na(vd[env$out])))
  expect_equal(vd[ids], vals)

  expect_true(is.na(.metric_map_fill("p_rsa_score")))
  expect_true(is.na(.metric_map_fill("z_rsa_score")))
  expect_identical(.metric_map_fill("Accuracy"), 0)
  expect_identical(.metric_map_fill("prob_A"), 0)
})

test_that("spatial_nmf_voxelwise_stats p-maps from sparse maps are NA outside the mask (#130)", {
  sp <- neuroim2::NeuroSpace(c(4, 4, 2))
  idx <- c(2L, 5L, 9L, 17L, 30L)
  mean_v <- c(2, -1, 0.5, 3, 0.1)
  sd_v <- c(1, 1, 1, 1, 1)
  mk <- function(v) neuroim2::SparseNeuroVol(data = v, space = sp, indices = idx)
  stab <- list(maps = list(mean = list(mk(mean_v)), sd = list(mk(sd_v))))

  res <- spatial_nmf_voxelwise_stats(stability = stab, mask_indices = idx)
  p <- .map_vals(res$p[[1]])
  z <- .map_vals(res$z[[1]])
  out <- setdiff(seq_len(32), idx)
  expect_equal(sum(p[out] < 0.05, na.rm = TRUE), 0)
  expect_true(all(is.na(p[out])))
  expect_true(all(is.na(z[out])))
  expect_equal(z[idx], mean_v / sd_v)
  expect_equal(p[idx], 2 * stats::pnorm(-abs(mean_v / sd_v)))

  # Without mask_indices the stored voxels define the mask.
  res2 <- spatial_nmf_voxelwise_stats(stability = stab)
  expect_true(all(is.na(.map_vals(res2$p[[1]])[out])))
})

test_that("searchlight n_voxels counts in-mask voxels, not the grid (#131)", {
  env <- .oom_setup()
  expect_equal(estimate_mask_size(env$dset), 64)
  expect_equal(estimate_mask_size(list(mask = c(1, 0, 1, 1, 0))), 3)

  futile.logger::flog.threshold(futile.logger::ERROR)
  set.seed(4)
  sl <- suppressWarnings(run_searchlight(.vrsa_spec(env, nperm = 0), radius = 2))
  expect_equal(sl$n_voxels, 64)
  txt <- capture.output(print(sl))
  expect_true(any(grepl("in Mask:\\s+64\\b", txt)))
  expect_false(any(grepl("216", txt)))
})
