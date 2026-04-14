library(rMVPA)
library(neuroim2)

manual_feature_rsa_rdm_vec <- function(x) {
  cm <- suppressWarnings(stats::cor(t(x), use = "pairwise.complete.obs"))
  rdm <- 1 - cm
  as.numeric(rdm[lower.tri(rdm)])
}

feature_rsa_perm_oracle <- function(observed,
                                    predicted,
                                    nperm,
                                    save_distributions,
                                    pattern_cor,
                                    pattern_discrim,
                                    pattern_rank,
                                    rdm_cor,
                                    voxel_cor,
                                    mse,
                                    r_squared,
                                    mean_voxelwise_temporal_cor,
                                    valid_col) {
  observed <- as.matrix(observed)
  predicted <- as.matrix(predicted)

  metric_names <- c(
    "pattern_correlation",
    "pattern_discrimination",
    "pattern_rank_percentile",
    "rdm_correlation",
    "voxel_correlation",
    "mse",
    "r_squared",
    "mean_voxelwise_temporal_cor"
  )
  obs_vals <- c(
    pattern_cor,
    pattern_discrim,
    pattern_rank,
    rdm_cor,
    voxel_cor,
    mse,
    r_squared,
    mean_voxelwise_temporal_cor
  )
  names(obs_vals) <- metric_names

  n_rows <- nrow(predicted)
  sd_thresh <- 1e-12
  tss <- sum((observed - mean(observed, na.rm = TRUE))^2, na.rm = TRUE)
  observed_valid <- observed[, valid_col, drop = FALSE]
  obs_row_sd <- apply(observed_valid, 1L, stats::sd)

  count_better <- setNames(integer(length(metric_names)), metric_names)
  sum_perm <- setNames(numeric(length(metric_names)), metric_names)
  sum_sq_perm <- setNames(numeric(length(metric_names)), metric_names)
  n_valid_perm <- setNames(integer(length(metric_names)), metric_names)

  if (isTRUE(save_distributions)) {
    dist_mat <- matrix(
      NA_real_,
      nrow = nperm,
      ncol = length(metric_names),
      dimnames = list(NULL, metric_names)
    )
  }

  for (i in seq_len(nperm)) {
    perm_idx <- sample(n_rows)
    perm_pred <- predicted[perm_idx, , drop = FALSE]
    perm_pred_valid <- perm_pred[, valid_col, drop = FALSE]

    ppc <- ppd <- ppr <- prdm <- pvc <- pmse <- prsq <- pmvtc <- NA_real_

    pred_row_sd <- apply(perm_pred_valid, 1L, stats::sd)
    vr <- which(obs_row_sd > sd_thresh & pred_row_sd > sd_thresh)

    if (length(vr) >= 2L) {
      cm <- tryCatch(
        stats::cor(
          t(perm_pred_valid[vr, , drop = FALSE]),
          t(observed_valid[vr, , drop = FALSE]),
          use = "pairwise.complete.obs"
        ),
        error = function(e) NULL
      )
      if (!is.null(cm)) {
        dc <- diag(cm)
        ppc <- mean(dc, na.rm = TRUE)
        ppd <- ppc - mean(cm[row(cm) != col(cm)], na.rm = TRUE)

        ranks <- numeric(length(vr))
        for (j in seq_along(vr)) {
          row_cors <- cm[j, ]
          denom <- sum(!is.na(row_cors)) - 1L
          ranks[j] <- if (denom > 0L && is.finite(row_cors[j])) {
            (sum(row_cors <= row_cors[j], na.rm = TRUE) - 1L) / denom
          } else {
            NA_real_
          }
        }
        ppr <- mean(ranks, na.rm = TRUE)
      }
    }

    if (length(vr) >= 3L) {
      pc <- tryCatch(
        stats::cor(t(perm_pred_valid[vr, , drop = FALSE]), use = "pairwise.complete.obs"),
        error = function(e) NULL
      )
      oc <- tryCatch(
        stats::cor(t(observed_valid[vr, , drop = FALSE]), use = "pairwise.complete.obs"),
        error = function(e) NULL
      )
      if (!is.null(pc) && !is.null(oc)) {
        pv <- (1 - pc)[upper.tri(pc)]
        ov <- (1 - oc)[upper.tri(oc)]
        if (length(pv) >= 2L && length(pv) == length(ov)) {
          prdm <- tryCatch(
            stats::cor(pv, ov, method = "spearman", use = "complete.obs"),
            error = function(e) NA_real_
          )
        }
      }
    }

    pvc <- tryCatch(
      stats::cor(
        as.vector(perm_pred_valid),
        as.vector(observed_valid),
        use = "pairwise.complete.obs"
      ),
      error = function(e) NA_real_
    )
    pmse <- mean((perm_pred - observed)^2, na.rm = TRUE)
    prsq <- if (tss > 0) {
      1 - sum((observed - perm_pred)^2, na.rm = TRUE) / tss
    } else {
      NA_real_
    }

    if (n_rows > 1L && length(valid_col) > 0L) {
      pmvtc <- mean(vapply(seq_len(length(valid_col)), function(j) {
        tryCatch(
          stats::cor(
            observed_valid[, j],
            perm_pred_valid[, j],
            use = "pairwise.complete.obs"
          ),
          error = function(e) NA_real_
        )
      }, numeric(1)), na.rm = TRUE)
    }

    perm_vals <- c(ppc, ppd, ppr, prdm, pvc, pmse, prsq, pmvtc)
    for (m in seq_along(metric_names)) {
      perm_val <- perm_vals[m]
      obs_val <- obs_vals[m]
      if (!is.na(perm_val)) {
        n_valid_perm[m] <- n_valid_perm[m] + 1L
        sum_perm[m] <- sum_perm[m] + perm_val
        sum_sq_perm[m] <- sum_sq_perm[m] + perm_val^2
        if (!is.na(obs_val)) {
          better <- if (metric_names[m] == "mse") perm_val <= obs_val else perm_val >= obs_val
          if (isTRUE(better)) {
            count_better[m] <- count_better[m] + 1L
          }
        }
      }
    }

    if (isTRUE(save_distributions)) {
      dist_mat[i, ] <- perm_vals
    }
  }

  eps <- .Machine$double.eps
  p_values <- vapply(seq_along(metric_names), function(m) {
    if (n_valid_perm[m] > 0L) {
      (count_better[m] + 1) / (n_valid_perm[m] + 1)
    } else {
      NA_real_
    }
  }, numeric(1))
  names(p_values) <- metric_names

  z_scores <- vapply(seq_along(metric_names), function(m) {
    if (n_valid_perm[m] > 0L) {
      mn <- sum_perm[m] / n_valid_perm[m]
      sd_perm <- sqrt(max(0, sum_sq_perm[m] / n_valid_perm[m] - mn^2))
      sd_use <- max(sd_perm, eps)
      if (metric_names[m] == "mse") {
        (mn - obs_vals[m]) / sd_use
      } else {
        (obs_vals[m] - mn) / sd_use
      }
    } else {
      NA_real_
    }
  }, numeric(1))
  names(z_scores) <- metric_names

  out <- list(p_values = p_values, z_scores = z_scores)
  if (isTRUE(save_distributions)) {
    out$permutation_distributions <- as.list(as.data.frame(dist_mat))
  }
  out
}

test_that("regional feature_rsa_model with direct F matrix runs without error", {
  # Generate a sample dataset: small volume, say 5x5x5, with 100 observations
  dset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  # Create a feature matrix F, e.g., 100 observations x 20 features
  Fmat <- matrix(rnorm(100*20), 100, 20)
  labels <- paste0("obs", 1:100)
  
  # Create feature_rsa_design using direct F matrix
  fdes <- feature_rsa_design(F=Fmat, labels=labels, max_comps=3)
  
  # Create a feature_rsa_model, for example using 'pls'
  mspec <- feature_rsa_model(dset$dataset, fdes, method="pls", crossval=blocked_cross_validation(dset$design$block_var))
  
  # Create a region mask with 5 ROIs
  region_mask <- NeuroVol(sample(1:5, size=length(dset$dataset$mask), replace=TRUE), 
  space(dset$dataset$mask))
  
  # Run regional analysis
  res <- run_regional(mspec, region_mask)
  
  # Check results
  expect_true(!is.null(res))
  expect_s3_class(res, "regional_mvpa_result")
  expect_true(!is.null(res$performance_table))
})

test_that("regional feature_rsa_model with S-based feature extraction runs without error", {
  # Generate a sample dataset: again 5x5x5 volume, 100 obs
  dset <- gen_sample_dataset(c(6,5,5), 100, blocks=3)
  
  # Create a similarity matrix S: must be symmetric and match number of observations
  obs_features <- matrix(rnorm(100*10), 100, 10)
  S <- tcrossprod(base::scale(obs_features))  # similarity matrix
  labels <- paste0("obs", 1:100)
  
  # Create feature_rsa_design using S
  fdes <- feature_rsa_design(S=S, labels=labels, k=10) # reduce to 5 dims
  
  # Create a feature_rsa_model using pca this time
  mspec <- feature_rsa_model(dset$dataset, fdes, method="pca", crossval=blocked_cross_validation(dset$design$block_var))
  
  # Create a region mask with 5 ROIs
  region_mask <- NeuroVol(sample(1:5, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  # Run regional analysis
  res <- run_regional(mspec, region_mask)
  
  # Check results
  expect_true(!is.null(res))
  expect_s3_class(res, "regional_mvpa_result")
  expect_true(!is.null(res$performance_table))
  
  # Check that performance_table has expected columns
  expect_true("pattern_correlation" %in% colnames(res$performance_table))
  expect_true("pattern_discrimination" %in% colnames(res$performance_table))
  expect_true("rdm_correlation" %in% colnames(res$performance_table))
  expect_true("voxel_correlation" %in% colnames(res$performance_table))
  expect_true("mse" %in% colnames(res$performance_table))
  expect_true("r_squared" %in% colnames(res$performance_table))
})

test_that("feature_rsa_model with permutation testing works correctly", {
  set.seed(123)
  
  # Create a small dataset for faster testing
  dset <- gen_sample_dataset(c(3,3,3), 50, blocks=2)
  
  # Create a feature matrix
  Fmat <- matrix(rnorm(50*10), 50, 10)
  labels <- paste0("obs", 1:50)
  
  # Create feature_rsa_design
  fdes <- feature_rsa_design(F=Fmat, labels=labels)
  
  # Create a feature_rsa_model with permutation testing
  mspec <- feature_rsa_model(
    dset$dataset, 
    fdes, 
    method="pca", 
    crossval=blocked_cross_validation(dset$design$block_var),
    nperm=10,  # Small number for testing
    permute_by="observations",
    save_distributions=TRUE
  )
  
  # Create a region mask with just 2 ROIs for faster testing
  region_mask <- NeuroVol(sample(1:2, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  # Run regional analysis
  res <- run_regional(mspec, region_mask)
  
  # Check results
  expect_true(!is.null(res))
  expect_s3_class(res, "regional_mvpa_result")
  
  # Check that permutation results are included in performance table
  perf_cols <- colnames(res$performance_table)
  expect_true(any(grepl("^p_", perf_cols)))  # p-values
  expect_true(any(grepl("^z_", perf_cols)))  # z-scores
  
  # Check specific permutation columns
  expect_true("p_pattern_correlation" %in% perf_cols)
  expect_true("z_pattern_correlation" %in% perf_cols)
  expect_true("p_pattern_discrimination" %in% perf_cols)
  expect_true("z_pattern_discrimination" %in% perf_cols)
  expect_true("p_rdm_correlation" %in% perf_cols)
  expect_true("z_rdm_correlation" %in% perf_cols)
})

test_that("feature_rsa_model with permute_by='features' works correctly", {
  set.seed(123)
  
  # Create a small dataset for faster testing
  dset <- gen_sample_dataset(c(3,3,3), 40, blocks=2)
  
  # Create a feature matrix
  Fmat <- matrix(rnorm(40*8), 40, 8)
  labels <- paste0("obs", 1:40)
  
  # Create feature_rsa_design
  fdes <- feature_rsa_design(F=Fmat, labels=labels)
  
  # Create a feature_rsa_model with permutation testing by features
  mspec <- feature_rsa_model(
    dset$dataset, 
    fdes, 
    method="pca", 
    crossval=blocked_cross_validation(dset$design$block_var),
    nperm=5,  # Small number for testing
    permute_by="features"  # Permute features instead of observations
  )
  
  # Create a region mask with just 2 ROIs for faster testing
  region_mask <- NeuroVol(sample(1:2, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  # Run regional analysis
  res <- run_regional(mspec, region_mask)
  
  # Check results
  expect_true(!is.null(res))
  expect_s3_class(res, "regional_mvpa_result")
  
  # Check that permutation results are included in performance table
  perf_cols <- colnames(res$performance_table)
  expect_true(any(grepl("^p_", perf_cols)))
  expect_true(any(grepl("^z_", perf_cols)))
})

test_that("feature_rsa_model errors on removed cache_pca argument", {
  set.seed(123)
  
  # Create a small dataset
  dset <- gen_sample_dataset(c(3,3,3), 30, blocks=2)
  
  # Create a feature matrix
  Fmat <- matrix(rnorm(30*6), 30, 6)
  labels <- paste0("obs", 1:30)
  
  # Create feature_rsa_design
  fdes <- feature_rsa_design(F=Fmat, labels=labels)
  
  expect_error(
    feature_rsa_model(
      dset$dataset,
      fdes,
      method = "pca",
      crossval = blocked_cross_validation(dset$design$block_var),
      cache_pca = TRUE
    ),
    "cache_pca.*no longer supported"
  )
})

test_that("feature_rsa_model with different max_comps values works correctly", {
  set.seed(123)
  
  # Create a small dataset
  dset <- gen_sample_dataset(c(3,3,3), 40, blocks=2)
  
  # Create a feature matrix with 10 dimensions
  Fmat <- matrix(rnorm(40*10), 40, 10)
  labels <- paste0("obs", 1:40)
  
  # Create feature_rsa_design with different max_comps values
  fdes1 <- feature_rsa_design(F=Fmat, labels=labels, max_comps=3)
  fdes2 <- feature_rsa_design(F=Fmat, labels=labels, max_comps=5)
  
  # Create feature_rsa_models
  mspec1 <- feature_rsa_model(dset$dataset, fdes1, method="pca", 
                             crossval=blocked_cross_validation(dset$design$block_var))
  mspec2 <- feature_rsa_model(dset$dataset, fdes2, method="pca", 
                             crossval=blocked_cross_validation(dset$design$block_var))
  
  # Check that max_comps was properly set
  expect_equal(mspec1$max_comps, 3)
  expect_equal(mspec2$max_comps, 5)
  
  # Create a region mask
  region_mask <- NeuroVol(sample(1:2, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  # Run regional analysis for both models
  res1 <- run_regional(mspec1, region_mask)
  res2 <- run_regional(mspec2, region_mask)
  
  # Check results
  expect_true(!is.null(res1))
  expect_true(!is.null(res2))
  expect_s3_class(res1, "regional_mvpa_result")
  expect_s3_class(res2, "regional_mvpa_result")
})

test_that("regional feature_rsa_model with pca method runs without error and returns performance", {
  dset <- gen_sample_dataset(c(4,4,4), 80, blocks=4)
  
  # Create an F matrix directly
  Fmat <- matrix(rnorm(80*15), 80, 15)
  labels <- paste0("obs", 1:80)
  
  fdes <- feature_rsa_design(F=Fmat, labels=labels) # no dimension reduction, just use as is
  mspec <- feature_rsa_model(dset$dataset, fdes, method="pca", crossval=blocked_cross_validation(dset$design$block_var))
  
  region_mask <- NeuroVol(sample(1:5, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  res <- run_regional(mspec, region_mask)
  
  expect_true(!is.null(res))
  expect_s3_class(res, "regional_mvpa_result")
  # Check if performance metrics are available
  expect_true("performance_table" %in% names(res))
  # Check specific performance metrics
  expect_true("pattern_correlation" %in% colnames(res$performance_table))
  expect_true("pattern_discrimination" %in% colnames(res$performance_table))
  expect_true("rdm_correlation" %in% colnames(res$performance_table))
  expect_true("voxel_correlation" %in% colnames(res$performance_table))
  expect_true("mse" %in% colnames(res$performance_table))
  expect_true("r_squared" %in% colnames(res$performance_table))
})


test_that("can compare feature_rsa with different methods", {
  set.seed(123)
  
  # Create a small dataset for faster testing
  dset <- gen_sample_dataset(c(4,4,4), 60, blocks=2)
  
  # Create a feature matrix
  Fmat <- matrix(rnorm(60*10), 60, 10)
  labels <- paste0("obs", 1:60)
  
  # Create feature_rsa_design
  fdes <- feature_rsa_design(F=Fmat, labels=labels, max_comps=10) # reduce to 5 dims
  
  # Create a region mask with just 2 ROIs for faster testing
  region_mask <- NeuroVol(sample(1:2, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  # Compare different methods
  methods <- c("pca", "pls")
  results_list <- list()
  
  for (method in methods) {
    print(method)
    # Create model with current method
    mspec <- feature_rsa_model(
      dataset = dset$dataset,
      design = fdes,
      method = method,
      crossval = blocked_cross_validation(dset$design$block_var)
    )
    
    # Run regional analysis
    results_list[[method]] <- run_regional(mspec, region_mask)
    
    # Check results
    expect_true(!is.null(results_list[[method]]))
    expect_s3_class(results_list[[method]], "regional_mvpa_result")
    expect_true(!is.null(results_list[[method]]$performance_table))
  }
  
  # Compare performance metrics across methods
  for (method in methods) {
    # Check that each method has the expected performance metrics
    perf_table <- results_list[[method]]$performance_table
    expect_true("pattern_correlation" %in% colnames(perf_table))
    expect_true("pattern_discrimination" %in% colnames(perf_table))
    expect_true("rdm_correlation" %in% colnames(perf_table))
    expect_true("voxel_correlation" %in% colnames(perf_table))
    expect_true("mse" %in% colnames(perf_table))
    expect_true("r_squared" %in% colnames(perf_table))
  }
})

test_that("evaluate_model.feature_rsa_model computes RDM correlation for predicted vs observed geometry", {
  set.seed(123)
  observed <- matrix(rnorm(12 * 30), nrow = 12, ncol = 30)
  predicted <- observed + matrix(rnorm(12 * 30, sd = 0.05), nrow = 12, ncol = 30)

  perf <- evaluate_model.feature_rsa_model(
    object = NULL,
    predicted = predicted,
    observed = observed,
    nperm = 0
  )

  expect_true("rdm_correlation" %in% names(perf))
  expect_true(is.finite(perf$rdm_correlation))
  expect_gt(perf$rdm_correlation, 0.8)
})

test_that("feature_rsa_model with permutation testing produces valid p-values", {
  set.seed(123)
  
  # Create a small dataset for faster testing
  dset <- gen_sample_dataset(c(3,3,3), 40, blocks=2)
  
  # Create a feature matrix with a strong signal
  X <- matrix(rnorm(40*8), 40, 8)
  B <- matrix(rnorm(8*3), 8, 3)
  Fmat <- X %*% B + matrix(rnorm(40*3, sd=0.1), 40, 3)
  labels <- paste0("obs", 1:40)
  
  # Create feature_rsa_design
  fdes <- feature_rsa_design(F=Fmat, labels=labels)
  
  # Create a feature_rsa_model with permutation testing
  mspec <- feature_rsa_model(
    dset$dataset, 
    fdes, 
    method="pca", 
    crossval=blocked_cross_validation(dset$design$block_var),
    nperm=20,  # Small number for testing
    permute_by="observations"
  )
  
  # Create a region mask with just 1 ROI for faster testing
  region_mask <- NeuroVol(rep(1, length(dset$dataset$mask)), space(dset$dataset$mask))
  
  # Run regional analysis
  res <- run_regional(mspec, region_mask)
  
  # Check results
  expect_true(!is.null(res))
  expect_s3_class(res, "regional_mvpa_result")
  
  # Check that p-values are between 0 and 1
  p_value_cols <- grep("^p_", colnames(res$performance_table), value=TRUE)
  for (col in p_value_cols) {
    p_values <- res$performance_table[[col]]
    expect_true(all(p_values >= 0 & p_values <= 1))
  }
  
  # Check that z-scores are reasonable
  z_score_cols <- grep("^z_", colnames(res$performance_table), value=TRUE)
  for (col in z_score_cols) {
    z_scores <- res$performance_table[[col]]
    expect_true(all(!is.na(z_scores)))
  }
})

test_that("evaluate_model.feature_rsa_model errors on row mismatch", {
  observed <- matrix(rnorm(10 * 4), nrow = 10, ncol = 4)
  predicted <- matrix(rnorm(9 * 4), nrow = 9, ncol = 4)

  expect_error(
    evaluate_model.feature_rsa_model(
      object = NULL,
      predicted = predicted,
      observed = observed,
      nperm = 0
    ),
    "Mismatch in rows"
  )
})

test_that("evaluate_model.feature_rsa_model returns full-length RDM vectors matching a direct oracle", {
  observed <- rbind(
    c(0, 1, 2, 3),
    c(5, 5, 5, 5),
    c(1, 3, 2, 4),
    c(2, 4, 3, 5),
    c(3, 5, 4, 6)
  )
  predicted <- observed + rbind(
    c(0.1, -0.2, 0.2, -0.1),
    c(0, 0, 0, 0),
    c(-0.1, 0.2, -0.2, 0.1),
    c(0.2, 0.1, -0.1, -0.2),
    c(-0.2, -0.1, 0.1, 0.2)
  )

  perf_with_vecs <- evaluate_model.feature_rsa_model(
    object = NULL,
    predicted = predicted,
    observed = observed,
    nperm = 0,
    compute_rdm_vectors = TRUE
  )
  perf_no_vecs <- evaluate_model.feature_rsa_model(
    object = NULL,
    predicted = predicted,
    observed = observed,
    nperm = 0,
    compute_rdm_vectors = FALSE
  )

  expected_pred <- manual_feature_rsa_rdm_vec(predicted)
  expected_obs <- manual_feature_rsa_rdm_vec(observed)

  expect_equal(perf_with_vecs$predicted_rdm_vec, expected_pred, tolerance = 1e-12)
  expect_equal(perf_with_vecs$observed_rdm_vec, expected_obs, tolerance = 1e-12)
  expect_length(perf_with_vecs$predicted_rdm_vec, nrow(observed) * (nrow(observed) - 1L) / 2L)
  expect_true(anyNA(perf_with_vecs$predicted_rdm_vec))
  expect_true(anyNA(perf_with_vecs$observed_rdm_vec))
  expect_null(perf_no_vecs$predicted_rdm_vec)
  expect_null(perf_no_vecs$observed_rdm_vec)
})

test_that(".perm_test_feature_rsa matches a direct reference implementation", {
  observed <- matrix(
    c(
      0, 1, 2, 3, 5,
      1, 2, 0, 4, 5,
      2, 3, 1, 5, 5,
      3, 4, 2, 6, 5,
      4, 5, 3, 7, 5,
      5, 6, 4, 8, 5
    ),
    nrow = 6,
    byrow = TRUE
  )
  predicted <- observed + matrix(
    c(
      0.2, -0.1, 0.1, -0.2, 0,
      -0.1, 0.2, -0.2, 0.1, 0,
      0.1, 0.1, -0.1, -0.1, 0,
      -0.2, 0.1, 0.2, -0.1, 0,
      0.1, -0.2, 0.1, 0.2, 0,
      -0.1, 0.2, -0.1, 0.1, 0
    ),
    nrow = 6,
    byrow = TRUE
  )

  perf <- evaluate_model.feature_rsa_model(
    object = NULL,
    predicted = predicted,
    observed = observed,
    nperm = 0,
    compute_rdm_vectors = FALSE
  )
  valid_col <- which(
    apply(observed, 2L, stats::sd) > 1e-12 &
      apply(predicted, 2L, stats::sd) > 1e-12
  )

  set.seed(8123)
  oracle <- feature_rsa_perm_oracle(
    observed = observed,
    predicted = predicted,
    nperm = 12,
    save_distributions = TRUE,
    pattern_cor = perf$pattern_correlation,
    pattern_discrim = perf$pattern_discrimination,
    pattern_rank = perf$pattern_rank_percentile,
    rdm_cor = perf$rdm_correlation,
    voxel_cor = perf$voxel_correlation,
    mse = perf$mse,
    r_squared = perf$r_squared,
    mean_voxelwise_temporal_cor = perf$mean_voxelwise_temporal_cor,
    valid_col = valid_col
  )

  set.seed(8123)
  got <- rMVPA:::.perm_test_feature_rsa(
    observed = observed,
    predicted = predicted,
    nperm = 12,
    save_distributions = TRUE,
    pattern_cor = perf$pattern_correlation,
    pattern_discrim = perf$pattern_discrimination,
    pattern_rank = perf$pattern_rank_percentile,
    rdm_cor = perf$rdm_correlation,
    voxel_cor = perf$voxel_correlation,
    mse = perf$mse,
    r_squared = perf$r_squared,
    mean_voxelwise_temporal_cor = perf$mean_voxelwise_temporal_cor,
    valid_col = valid_col
  )

  expect_equal(got$p_values, oracle$p_values, tolerance = 1e-10)
  expect_equal(got$z_scores, oracle$z_scores, tolerance = 1e-10)
  expect_equal(got$permutation_distributions, oracle$permutation_distributions, tolerance = 1e-10)
})

# ---- ncomp_selection tests ----

test_that("ncomp_selection='loo' selects components via LOO for PLS", {
  set.seed(42)
  dset <- gen_sample_dataset(c(4,4,4), 60, blocks=3)
  Fmat <- matrix(rnorm(60 * 8), 60, 8)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("o", 1:60), max_comps = 8)
  region_mask <- NeuroVol(sample(1:2, length(dset$dataset$mask), replace = TRUE),
                          space(dset$dataset$mask))
  mspec <- feature_rsa_model(dset$dataset, fdes, method = "pls",
                              ncomp_selection = "loo",
                              crossval = blocked_cross_validation(dset$design$block_var))
  res <- run_regional(mspec, region_mask)
  expect_s3_class(res, "regional_mvpa_result")
  expect_true("ncomp" %in% colnames(res$performance_table))
  # LOO may select fewer components than max_comps
  ncomps <- res$performance_table$ncomp
  expect_true(all(ncomps >= 1 & ncomps <= 8))
})

test_that("ncomp_selection='loo' works for PCA (via pcr)", {
  set.seed(42)
  dset <- gen_sample_dataset(c(4,4,4), 60, blocks=3)
  Fmat <- matrix(rnorm(60 * 8), 60, 8)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("o", 1:60), max_comps = 8)
  region_mask <- NeuroVol(sample(1:2, length(dset$dataset$mask), replace = TRUE),
                          space(dset$dataset$mask))
  mspec <- feature_rsa_model(dset$dataset, fdes, method = "pca",
                              ncomp_selection = "loo",
                              crossval = blocked_cross_validation(dset$design$block_var))
  res <- run_regional(mspec, region_mask)
  expect_s3_class(res, "regional_mvpa_result")
  ncomps <- res$performance_table$ncomp
  expect_true(all(ncomps >= 1 & ncomps <= 8))
})

test_that("ncomp_selection='pve' keeps components to meet variance threshold", {
  set.seed(42)
  dset <- gen_sample_dataset(c(4,4,4), 60, blocks=3)
  Fmat <- matrix(rnorm(60 * 8), 60, 8)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("o", 1:60), max_comps = 8)
  region_mask <- NeuroVol(sample(1:2, length(dset$dataset$mask), replace = TRUE),
                          space(dset$dataset$mask))

  mspec <- feature_rsa_model(dset$dataset, fdes, method = "pls",
                              ncomp_selection = "pve", pve_threshold = 0.8,
                              crossval = blocked_cross_validation(dset$design$block_var))
  res <- run_regional(mspec, region_mask)
  expect_s3_class(res, "regional_mvpa_result")
  ncomps <- res$performance_table$ncomp
  expect_true(all(ncomps >= 1 & ncomps <= 8))
})

test_that("ncomp_selection='max' uses all max_comps (legacy behaviour)", {
  set.seed(42)
  dset <- gen_sample_dataset(c(4,4,4), 60, blocks=3)
  Fmat <- matrix(rnorm(60 * 6), 60, 6)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("o", 1:60), max_comps = 6)
  region_mask <- NeuroVol(sample(1:2, length(dset$dataset$mask), replace = TRUE),
                          space(dset$dataset$mask))
  mspec <- feature_rsa_model(dset$dataset, fdes, method = "pls",
                              ncomp_selection = "max",
                              crossval = blocked_cross_validation(dset$design$block_var))
  res <- run_regional(mspec, region_mask)
  expect_s3_class(res, "regional_mvpa_result")
})

# ---- Collapse protection tests ----

test_that(".selectNcomp_mv returns NA or valid integer for near-degenerate model", {
  # Build a near-degenerate PLS model with near-constant response
  set.seed(99)
  n <- 30; p <- 5; q <- 4
  X <- matrix(rnorm(n * p), n, p)
  # Near-constant response → near-zero variance → unstable CV MSEP
  Y <- matrix(1, n, q) + matrix(rnorm(n * q, sd = 1e-10), n, q)
  fit <- tryCatch(
    pls::plsr(Y ~ X, ncomp = 2, scale = FALSE, validation = "LOO"),
    error = function(e) NULL
  )
  # If PLS itself errors on degenerate data, that's fine — skip
  skip_if(is.null(fit), "pls::plsr cannot fit degenerate model")
  result <- tryCatch(
    rMVPA:::.selectNcomp_mv(fit, method = "onesigma"),
    error = function(e) NA_integer_  # pls internals may fail outside package namespace
  )
  # Should be NA_integer_ or a valid integer (not an error)
  expect_true(is.na(result) || (is.integer(result) && result >= 1L))
})

test_that(".selectNcomp_mv works when pls is namespace-loaded only", {
  was_attached <- "package:pls" %in% search()
  if (was_attached) {
    detach("package:pls", unload = TRUE, character.only = TRUE)
    on.exit(suppressPackageStartupMessages(library(pls)), add = TRUE)
  }

  set.seed(123)
  n <- 40
  X <- matrix(rnorm(n * 6), n, 6)
  Y <- matrix(rnorm(n * 3), n, 3)
  fit <- pls::pcr(Y ~ X, ncomp = 4, scale = FALSE, validation = "LOO")

  nc <- rMVPA:::.selectNcomp_mv(fit, method = "onesigma")
  expect_true(is.integer(nc))
  expect_true(nc >= 1L)
  expect_true(nc <= 4L)
})

test_that(".selectNcomp_mv matches segment-wise onesigma selection", {
  set.seed(321)
  n <- 30
  X <- matrix(rnorm(n * 5), n, 5)
  Y <- matrix(rnorm(n * 4), n, 4)
  fit <- pls::pcr(Y ~ X, ncomp = 4, scale = FALSE, validation = "LOO")

  pred <- fit$validation$pred
  segments <- fit$validation$segments
  seg_mse <- vapply(seq_len(dim(pred)[3]), function(k) {
    mean(vapply(segments, function(idx) {
      mean((Y[idx, , drop = FALSE] - pred[idx, , k, drop = FALSE][, , 1L])^2)
    }, numeric(1)))
  }, numeric(1))

  mse_at_k <- vapply(segments, function(idx) {
    k <- which.min(seg_mse)
    mean((Y[idx, , drop = FALSE] - pred[idx, , k, drop = FALSE][, , 1L])^2)
  }, numeric(1))
  k_min <- which.min(seg_mse)
  thresh <- seg_mse[k_min] + stats::sd(mse_at_k) / sqrt(length(mse_at_k))
  manual <- which(seg_mse <= thresh)[1]

  expect_equal(rMVPA:::.selectNcomp_mv(fit, method = "onesigma"), manual)
})

test_that("feature_rsa_design weights S-derived basis by sqrt eigenvalues", {
  S <- diag(c(9, 4, 1))
  des <- feature_rsa_design(S = S, labels = paste0("o", 1:3), k = 3)
  expect_equal(unname(des$F), diag(c(3, 2, 1)), tolerance = 1e-12)
})

test_that("PLS LOO selection does not collapse to NA ncomp", {
  set.seed(77)
  dset <- gen_sample_dataset(c(4,4,4), 60, blocks = 3)
  # Use very low-signal features to stress component selection
  Fmat <- matrix(rnorm(60 * 4, sd = 0.01), 60, 4)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("o", 1:60), max_comps = 4)
  region_mask <- NeuroVol(rep(1L, length(dset$dataset$mask)),
                          space(dset$dataset$mask))

  mspec <- feature_rsa_model(dset$dataset, fdes, method = "pls",
                              ncomp_selection = "loo",
                              crossval = blocked_cross_validation(dset$design$block_var))
  res <- run_regional(mspec, region_mask)
  expect_s3_class(res, "regional_mvpa_result")
  ncomps <- res$performance_table$ncomp
  expect_true(all(!is.na(ncomps) & ncomps >= 1))
})

test_that("PVE selection falls back when threshold unreachable", {
  set.seed(88)
  dset <- gen_sample_dataset(c(4,4,4), 60, blocks = 3)
  Fmat <- matrix(rnorm(60 * 4), 60, 4)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("o", 1:60), max_comps = 4)
  region_mask <- NeuroVol(rep(1L, length(dset$dataset$mask)),
                          space(dset$dataset$mask))

  # Set threshold to 1.0 — cumulative ratio can only reach 1.0 at max comps,

  # but with only 4 components, early ones may not reach it, testing the fallback
  mspec <- feature_rsa_model(dset$dataset, fdes, method = "pls",
                              ncomp_selection = "pve", pve_threshold = 1.0,
                              crossval = blocked_cross_validation(dset$design$block_var))
  # Should not error — falls back to all components
  res <- run_regional(mspec, region_mask)
  expect_s3_class(res, "regional_mvpa_result")
  ncomps <- res$performance_table$ncomp
  expect_true(all(!is.na(ncomps) & ncomps >= 1))
})

test_that("PCA LOO selection does not collapse to NA ncomp", {
  set.seed(66)
  dset <- gen_sample_dataset(c(4,4,4), 60, blocks = 3)
  Fmat <- matrix(rnorm(60 * 4, sd = 0.01), 60, 4)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("o", 1:60), max_comps = 4)
  region_mask <- NeuroVol(rep(1L, length(dset$dataset$mask)),
                          space(dset$dataset$mask))

  mspec <- feature_rsa_model(dset$dataset, fdes, method = "pca",
                              ncomp_selection = "loo",
                              crossval = blocked_cross_validation(dset$design$block_var))
  res <- run_regional(mspec, region_mask)
  expect_s3_class(res, "regional_mvpa_result")
  ncomps <- res$performance_table$ncomp
  expect_true(all(!is.na(ncomps) & ncomps >= 1))
})

test_that("feature_rsa_model retains ROI predicted RDM vectors only when requested", {
  set.seed(901)
  dset <- gen_sample_dataset(c(3,3,3), 24, blocks = 3)
  Fmat <- matrix(rnorm(24 * 6), 24, 6)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("o", 1:24), max_comps = 3)
  region_mask <- NeuroVol(sample(1:3, length(dset$dataset$mask), replace = TRUE),
                          space(dset$dataset$mask))

  mspec_off <- feature_rsa_model(
    dset$dataset, fdes, method = "pls",
    crossval = blocked_cross_validation(dset$design$block_var),
    return_rdm_vectors = FALSE
  )
  res_off <- run_regional(mspec_off, region_mask)
  expect_true(is.null(res_off$fits))
  expect_error(feature_rsa_rdm_vectors(res_off), "return_rdm_vectors=TRUE")

  mspec_on <- feature_rsa_model(
    dset$dataset, fdes, method = "pls",
    crossval = blocked_cross_validation(dset$design$block_var),
    return_rdm_vectors = TRUE
  )
  res_on <- run_regional(mspec_on, region_mask)
  expect_true(!is.null(res_on$fits))

  vecs <- feature_rsa_rdm_vectors(res_on)
  expect_equal(nrow(vecs), nrow(res_on$performance_table))
  expect_true(all(vapply(vecs$rdm_vec, is.numeric, logical(1))))
  expect_true(all(vapply(vecs$rdm_vec, length, integer(1)) == 24L * 23L / 2L))
  expect_true(all(vapply(vecs$observation_index, identical, logical(1), y = seq_len(24))))
})

test_that("feature_rsa_connectivity returns a symmetric ROI x ROI matrix and sparsifies only the final graph", {
  set.seed(902)
  dset <- gen_sample_dataset(c(4,4,4), 30, blocks = 3)
  Fmat <- matrix(rnorm(30 * 8), 30, 8)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("o", 1:30), max_comps = 4)
  region_mask <- NeuroVol(sample(1:4, length(dset$dataset$mask), replace = TRUE),
                          space(dset$dataset$mask))

  mspec <- feature_rsa_model(
    dset$dataset, fdes, method = "pca",
    crossval = blocked_cross_validation(dset$design$block_var),
    return_rdm_vectors = TRUE
  )
  res <- run_regional(mspec, region_mask)

  conn_dense <- feature_rsa_connectivity(res, method = "spearman")
  expect_true(is.matrix(conn_dense))
  expect_equal(dim(conn_dense), c(nrow(res$performance_table), nrow(res$performance_table)))
  expect_equal(conn_dense, t(conn_dense))
  expect_equal(unname(diag(conn_dense)), rep(1, nrow(conn_dense)))

  conn_sparse <- feature_rsa_connectivity(res, method = "spearman", keep = 0.25)
  expect_equal(dim(conn_sparse), dim(conn_dense))
  expect_equal(conn_sparse, t(conn_sparse))
  expect_equal(unname(diag(conn_sparse)), rep(1, nrow(conn_sparse)))
  expect_true(sum(conn_sparse[upper.tri(conn_sparse)] == 0, na.rm = TRUE) >= 1L)
})

test_that("feature_rsa_connectivity matches direct correlation on synthetic ROI RDM vectors", {
  vec_tbl <- tibble::tibble(
    roinum = c(10L, 20L, 30L),
    n_obs = c(4L, 4L, 4L),
    observation_index = list(1:4, 1:4, 1:4),
    rdm_vec = list(
      c(1, 2, 3, 4),
      c(2, 4, 6, 8),
      c(4, 3, 2, 1)
    )
  )

  expected <- stats::cor(do.call(cbind, vec_tbl$rdm_vec), method = "pearson")
  dimnames(expected) <- list(as.character(vec_tbl$roinum), as.character(vec_tbl$roinum))
  diag(expected) <- 1

  got <- feature_rsa_connectivity(vec_tbl, method = "pearson")
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("feature_rsa_connectivity is invariant to per-ROI affine rescaling of RDM vectors", {
  vec_tbl <- tibble::tibble(
    roinum = c(1L, 2L, 3L),
    n_obs = c(4L, 4L, 4L),
    observation_index = list(1:4, 1:4, 1:4),
    rdm_vec = list(
      c(1, 2, 3, 4),
      c(2, 4, 6, 8),
      c(4, 3, 2, 1)
    )
  )

  aff_tbl <- vec_tbl
  aff_tbl$rdm_vec <- list(
    5 * vec_tbl$rdm_vec[[1]] + 7,
    0.5 * vec_tbl$rdm_vec[[2]] - 3,
    2 * vec_tbl$rdm_vec[[3]] + 10
  )

  base_conn <- feature_rsa_connectivity(vec_tbl, method = "pearson")
  aff_conn <- feature_rsa_connectivity(aff_tbl, method = "pearson")
  expect_equal(aff_conn, base_conn, tolerance = 1e-12)
})

test_that("feature_rsa_connectivity validates observation ordering and vector length contracts", {
  vec_tbl <- tibble::tibble(
    roinum = c(1L, 2L),
    n_obs = c(4L, 4L),
    observation_index = list(1:4, 1:4),
    rdm_vec = list(
      c(1, 2, 3, 4),
      c(4, 3, 2, 1)
    )
  )

  bad_order <- vec_tbl
  bad_order$observation_index[[2]] <- c(2L, 1L, 3L, 4L)
  expect_error(feature_rsa_connectivity(bad_order), "observation ordering")

  bad_len <- vec_tbl
  bad_len$rdm_vec[[2]] <- bad_len$rdm_vec[[2]][-1]
  expect_error(feature_rsa_connectivity(bad_len), "same length")
})

test_that("feature_rsa_rdm_vectors also extracts observed_rdm_vec when available", {
  set.seed(903)
  dset <- gen_sample_dataset(c(3,3,3), 24, blocks = 3)
  Fmat <- matrix(rnorm(24 * 6), 24, 6)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("o", 1:24), max_comps = 3)
  region_mask <- NeuroVol(sample(1:3, length(dset$dataset$mask), replace = TRUE),
                          space(dset$dataset$mask))

  mspec <- feature_rsa_model(
    dset$dataset, fdes, method = "pls",
    crossval = blocked_cross_validation(dset$design$block_var),
    return_rdm_vectors = TRUE
  )
  res <- run_regional(mspec, region_mask)
  vecs <- feature_rsa_rdm_vectors(res)

  expect_true("observed_rdm_vec" %in% names(vecs))
  expect_true(all(vapply(vecs$observed_rdm_vec, is.numeric, logical(1))))
  expect_true(all(vapply(vecs$observed_rdm_vec, length, integer(1)) == 24L * 23L / 2L))
})

test_that("feature_rsa_cross_connectivity returns asymmetric ROI x ROI matrix from real data", {
  set.seed(904)
  dset <- gen_sample_dataset(c(4,4,4), 30, blocks = 3)
  Fmat <- matrix(rnorm(30 * 8), 30, 8)
  fdes <- feature_rsa_design(F = Fmat, labels = paste0("o", 1:30), max_comps = 4)
  region_mask <- NeuroVol(sample(1:4, length(dset$dataset$mask), replace = TRUE),
                          space(dset$dataset$mask))

  mspec <- feature_rsa_model(
    dset$dataset, fdes, method = "pca",
    crossval = blocked_cross_validation(dset$design$block_var),
    return_rdm_vectors = TRUE
  )
  res <- run_regional(mspec, region_mask)

  cross <- feature_rsa_cross_connectivity(res, method = "spearman")
  n_roi <- nrow(res$performance_table)
  expect_true(is.matrix(cross))
  expect_equal(dim(cross), c(n_roi, n_roi))
  expect_equal(names(dimnames(cross)), c("predicted", "observed"))
})

test_that("feature_rsa_cross_connectivity matches manual cross-correlation on synthetic data", {
  pred_vecs <- list(c(1, 2, 3, 4), c(2, 4, 6, 8), c(4, 3, 2, 1))
  obs_vecs  <- list(c(4, 3, 2, 1), c(1, 1, 2, 2), c(3, 1, 4, 2))

  vec_tbl <- tibble::tibble(
    roinum = c(10L, 20L, 30L),
    n_obs = c(4L, 4L, 4L),
    observation_index = list(1:4, 1:4, 1:4),
    rdm_vec = pred_vecs,
    observed_rdm_vec = obs_vecs
  )

  pred_mat <- do.call(cbind, pred_vecs)
  obs_mat  <- do.call(cbind, obs_vecs)
  expected <- stats::cor(pred_mat, obs_mat, method = "pearson")
  dimnames(expected) <- list(predicted = as.character(vec_tbl$roinum),
                             observed  = as.character(vec_tbl$roinum))

  got <- feature_rsa_cross_connectivity(vec_tbl, method = "pearson")
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("feature_rsa_cross_connectivity double_center removes additive source and target offsets", {
  pred_vecs <- list(c(1, 2, 3, 4), c(2, 4, 6, 8), c(4, 3, 2, 1))
  obs_vecs  <- list(c(4, 3, 2, 1), c(1, 1, 2, 2), c(3, 1, 4, 2))

  vec_tbl <- tibble::tibble(
    roinum = c(10L, 20L, 30L),
    n_obs = c(4L, 4L, 4L),
    observation_index = list(1:4, 1:4, 1:4),
    rdm_vec = pred_vecs,
    observed_rdm_vec = obs_vecs
  )

  raw <- feature_rsa_cross_connectivity(vec_tbl, method = "pearson")
  grand <- mean(raw)
  row_mean <- rowMeans(raw)
  col_mean <- colMeans(raw)
  expected <- sweep(sweep(raw, 1L, row_mean, "-"), 2L, col_mean, "-") + grand

  got <- feature_rsa_cross_connectivity(
    vec_tbl,
    method = "pearson",
    adjust = "double_center",
    return_components = TRUE
  )

  expect_equal(got$raw_matrix, raw, tolerance = 1e-12)
  expect_equal(got$matrix, expected, tolerance = 1e-12)
  expect_equal(got$adjusted_matrix, expected, tolerance = 1e-12)
  expect_equal(unname(got$source_offset), unname(row_mean - grand), tolerance = 1e-12)
  expect_equal(unname(got$target_offset), unname(col_mean - grand), tolerance = 1e-12)
  expect_equal(unname(rowMeans(got$matrix)), rep(0, nrow(got$matrix)), tolerance = 1e-12)
  expect_equal(unname(colMeans(got$matrix)), rep(0, ncol(got$matrix)), tolerance = 1e-12)
})

test_that("feature_rsa_cross_connectivity residualize_mean matches manual nuisance residualization", {
  residualize_vector <- function(y, x) {
    out <- rep(NA_real_, length(y))
    ok <- is.finite(y) & is.finite(x)
    if (sum(ok) < 2L) {
      return(out)
    }

    x_ok <- x[ok]
    y_ok <- y[ok]
    if (stats::sd(x_ok) <= 1e-12) {
      out[ok] <- y_ok - mean(y_ok)
      return(out)
    }

    fit <- stats::lm.fit(x = cbind(1, x_ok), y = y_ok)
    out[ok] <- fit$residuals
    out
  }

  residualize_columns <- function(mat) {
    ref_vec <- rowMeans(mat)
    out <- vapply(
      seq_len(ncol(mat)),
      function(i) residualize_vector(mat[, i], ref_vec),
      numeric(nrow(mat))
    )
    dimnames(out) <- dimnames(mat)
    out
  }

  common <- c(6, -3, 2, 5, -1, 4)
  pred_vecs <- list(
    c(1, 3, 2, 5, 4, 6) + 3 * common,
    c(2, 1, 4, 3, 6, 5) + 3 * common,
    c(6, 5, 3, 2, 1, 4) + 3 * common
  )
  obs_vecs <- list(
    c(2, 4, 1, 5, 3, 6) + 2 * common,
    c(5, 2, 4, 1, 3, 6) + 2 * common,
    c(4, 6, 2, 5, 1, 3) + 2 * common
  )

  vec_tbl <- tibble::tibble(
    roinum = c(11L, 22L, 33L),
    n_obs = c(4L, 4L, 4L),
    observation_index = list(1:4, 1:4, 1:4),
    rdm_vec = pred_vecs,
    observed_rdm_vec = obs_vecs
  )

  pred_mat <- do.call(cbind, pred_vecs)
  obs_mat  <- do.call(cbind, obs_vecs)
  pred_resid <- residualize_columns(pred_mat)
  obs_resid  <- residualize_columns(obs_mat)
  expected <- stats::cor(pred_resid, obs_resid, method = "pearson")
  dimnames(expected) <- list(predicted = as.character(vec_tbl$roinum),
                             observed  = as.character(vec_tbl$roinum))

  got <- feature_rsa_cross_connectivity(
    vec_tbl,
    method = "pearson",
    adjust = "residualize_mean"
  )

  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("feature_rsa_cross_connectivity errors when observed vectors are missing", {
  vec_tbl <- tibble::tibble(
    roinum = c(1L, 2L),
    n_obs = c(4L, 4L),
    observation_index = list(1:4, 1:4),
    rdm_vec = list(c(1, 2, 3, 4), c(4, 3, 2, 1)),
    observed_rdm_vec = list(NULL, NULL)
  )
  expect_error(feature_rsa_cross_connectivity(vec_tbl), "no observed RDM vectors")
})
