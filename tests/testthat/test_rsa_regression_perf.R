# Parity and oracle tests for the RSA regression performance work:
#  * run_lm_semipartial(): cached lm kernel vs the lm() reference and an lm() oracle
#  * repnet_model / repmed_model: fixed-design QR fits vs lm() oracles
#  * table_to_rdm(): match()-based fill vs the former n^2 list-lookup implementation

## ===========================================================================
## Reference implementation of table_to_rdm (HEAD), kept verbatim for parity.
## ===========================================================================

table_to_rdm_reference <- function(similarity_table,
                        labels,
                        label1_col = "label1",
                        label2_col = "label2",
                        similarity_col = "similarity",
                        default_similarity = 0,
                        symmetric = TRUE,
                        self_similarity = 1,
                        as_dist = TRUE) {
  
  # Input validation
  if (!is.data.frame(similarity_table)) {
    stop("similarity_table must be a data frame")
  }
  
  if (!all(c(label1_col, label2_col, similarity_col) %in% names(similarity_table))) {
    stop(sprintf("similarity_table must contain columns: %s, %s, %s",
                 label1_col, label2_col, similarity_col))
  }
  
  # Extract columns
  label1_vals <- as.character(similarity_table[[label1_col]])
  label2_vals <- as.character(similarity_table[[label2_col]])
  similarity_vals <- as.numeric(similarity_table[[similarity_col]])
  
  # Check similarity values are in valid range
  if (any(!is.na(similarity_vals) & (similarity_vals < 0 | similarity_vals > 1))) {
    warning("Similarity values should be between 0 and 1. Values outside this range will be clipped.")
    similarity_vals <- pmax(0, pmin(1, similarity_vals))
  }
  
  # Check for duplicate entries
  if (symmetric) {
    # For symmetric tables, (A,B) and (B,A) are the same
    pair_keys <- paste(pmin(label1_vals, label2_vals), 
                      pmax(label1_vals, label2_vals), 
                      sep = "___")
  } else {
    # For asymmetric tables, (A,B) and (B,A) are different
    pair_keys <- paste(label1_vals, label2_vals, sep = "___")
  }
  
  if (any(duplicated(pair_keys))) {
    warning("Duplicate label pairs found in similarity_table. Using first occurrence.")
    # Keep only first occurrence
    keep_idx <- !duplicated(pair_keys)
    label1_vals <- label1_vals[keep_idx]
    label2_vals <- label2_vals[keep_idx]
    similarity_vals <- similarity_vals[keep_idx]
  }
  
  # Create lookup structure
  similarity_lookup <- list()
  for (i in seq_along(label1_vals)) {
    key1 <- paste(label1_vals[i], label2_vals[i], sep = "___")
    similarity_lookup[[key1]] <- similarity_vals[i]
    
    if (symmetric && label1_vals[i] != label2_vals[i]) {
      # Add reverse direction for symmetric case
      key2 <- paste(label2_vals[i], label1_vals[i], sep = "___")
      similarity_lookup[[key2]] <- similarity_vals[i]
    }
  }
  
  # Initialize RDM matrix
  n <- length(labels)
  rdm_matrix <- matrix(NA, nrow = n, ncol = n)
  rownames(rdm_matrix) <- labels
  colnames(rdm_matrix) <- labels
  
  # Fill the matrix
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) {
        # Diagonal elements (self-similarity)
        if (is.na(self_similarity)) {
          # Look up self-similarity in table
          key <- paste(labels[i], labels[j], sep = "___")
          sim_val <- similarity_lookup[[key]]
          if (is.null(sim_val)) {
            rdm_matrix[i, j] <- 0  # Default dissimilarity for diagonal
          } else {
            rdm_matrix[i, j] <- 1 - sim_val
          }
        } else {
          # Use specified self-similarity
          rdm_matrix[i, j] <- 1 - self_similarity
        }
      } else {
        # Off-diagonal elements
        key <- paste(labels[i], labels[j], sep = "___")
        sim_val <- similarity_lookup[[key]]
        
        if (is.null(sim_val)) {
          # Use default similarity for missing pairs
          rdm_matrix[i, j] <- 1 - default_similarity
        } else {
          # Convert similarity to dissimilarity
          rdm_matrix[i, j] <- 1 - sim_val
        }
      }
    }
  }
  
  # Ensure matrix is symmetric if specified
  if (symmetric) {
    # Average the upper and lower triangles to ensure perfect symmetry
    rdm_matrix <- (rdm_matrix + t(rdm_matrix)) / 2
  }
  
  # Return in requested format
  if (as_dist) {
    as.dist(rdm_matrix)
  } else {
    rdm_matrix
  }
}

## ===========================================================================
## Semi-partial correlations (rsa_model)
## ===========================================================================

sp_fixture <- function(seed = 2024L, n_obs = 24L, n_vox = 40L,
                       pair_mask = FALSE, nuisance = FALSE, alias = FALSE) {
  set.seed(seed)
  train_mat <- matrix(rnorm(n_obs * n_vox), nrow = n_obs, ncol = n_vox)
  dset <- list(train_data = matrix(rnorm(n_vox * n_obs), nrow = n_vox, ncol = n_obs),
               mask = seq_len(n_vox))
  class(dset) <- "mvpa_dataset"

  D1 <- dist(matrix(rnorm(n_obs * 6L), nrow = n_obs, ncol = 6L))
  D2 <- dist(matrix(rnorm(n_obs * 6L), nrow = n_obs, ncol = 6L))
  data <- list(D1 = D1, D2 = D2)
  formula <- ~ D1 + D2
  if (alias) {
    data$D3 <- D1                      # exact duplicate: aliased by lm()
    formula <- ~ D1 + D2 + D3
  }
  nuis <- if (nuisance) list(nuis = dist(matrix(rnorm(n_obs * 3L), nrow = n_obs, ncol = 3L))) else list()
  pm <- NULL
  if (pair_mask) {
    pm <- sample(c(TRUE, FALSE), choose(n_obs, 2L), replace = TRUE, prob = c(0.75, 0.25))
  }
  rdes <- rsa_design(formula, data, nuisance = nuis, pair_mask = pm)
  build <- function() {
    rsa_model(dataset = dset, design = rdes, distmethod = "pearson", regtype = "lm",
              semipartial = TRUE, check_collinearity = FALSE)
  }
  list(build = build, mspec = build(), train_mat = train_mat, rdes = rdes)
}

# Independent oracle: lm() fitted in the test, sr = sign(t) * sqrt(t^2 * MSE / TSS).
sp_oracle <- function(mspec, train_mat) {
  cm <- 1 - stats::cor(t(train_mat))
  dvec <- cm[lower.tri(cm)]
  if (!is.null(mspec$design$include)) dvec <- dvec[mspec$design$include]
  X <- as.data.frame(mspec$design$model_mat, stringsAsFactors = FALSE)
  fit <- stats::lm(dvec ~ ., data = cbind(data.frame(dvec = dvec), X))
  smry <- summary(fit)
  cf <- stats::coef(smry)
  tv <- setNames(rep(NA_real_, ncol(X)), names(X))
  est <- intersect(names(X), rownames(cf))
  tv[est] <- cf[est, "t value"]
  tss <- sum((dvec - mean(dvec))^2)
  setNames(sign(tv) * sqrt(tv^2 * smry$sigma^2 / tss), names(X))
}

sp_reference_values <- function(fx) {
  rMVPA:::.with_reference_paths("rsa_fast_kernel", {
    ref_spec <- fx$build()
    expect_null(ref_spec$.fast_kernel)
    rMVPA:::train_model.rsa_model(ref_spec, fx$train_mat, y = NULL, indices = NULL)
  })
}

test_that("semipartial kernel matches the lm() reference and the lm() oracle to 1e-12", {
  designs <- list(
    plain = sp_fixture(seed = 11L),
    include_mask = sp_fixture(seed = 12L, pair_mask = TRUE),
    nuisance = sp_fixture(seed = 13L, nuisance = TRUE),
    mask_and_nuisance = sp_fixture(seed = 14L, pair_mask = TRUE, nuisance = TRUE)
  )
  for (nm in names(designs)) {
    fx <- designs[[nm]]
    expect_true(is.list(fx$mspec$.fast_kernel$lm), info = nm)

    fast <- rMVPA:::train_model.rsa_model(fx$mspec, fx$train_mat, y = NULL, indices = NULL)
    ref <- sp_reference_values(fx)
    oracle <- sp_oracle(fx$mspec, fx$train_mat)

    expect_identical(names(fast), names(oracle), info = nm)
    expect_identical(names(ref), names(oracle), info = nm)
    expect_false(any(is.na(fast)), info = nm)
    expect_lt(max(abs(fast - ref)), 1e-12)
    expect_lt(max(abs(fast - oracle)), 1e-12)
    expect_lt(max(abs(ref - oracle)), 1e-12)
  }
})

test_that("semipartial correlations satisfy sr^2 = delta R^2 for each predictor", {
  fx <- sp_fixture(seed = 21L, nuisance = TRUE)
  sr <- rMVPA:::train_model.rsa_model(fx$mspec, fx$train_mat, y = NULL, indices = NULL)
  cm <- 1 - stats::cor(t(fx$train_mat))
  y <- cm[lower.tri(cm)]
  X <- as.data.frame(fx$mspec$design$model_mat)
  full <- summary(stats::lm(y ~ ., data = cbind(data.frame(y = y), X)))$r.squared
  for (term in names(X)) {
    reduced <- summary(stats::lm(y ~ ., data = cbind(data.frame(y = y), X[setdiff(names(X), term)])))$r.squared
    expect_equal(unname(sr[[term]])^2, full - reduced, tolerance = 1e-10, info = term)
  }
})

test_that("rank-deficient semipartial designs keep the lm() NA pattern and skip the kernel", {
  fx <- sp_fixture(seed = 31L, alias = TRUE)
  expect_null(fx$mspec$.fast_kernel$lm)

  fast <- rMVPA:::train_model.rsa_model(fx$mspec, fx$train_mat, y = NULL, indices = NULL)
  ref <- sp_reference_values(fx)
  oracle <- sp_oracle(fx$mspec, fx$train_mat)

  expect_true(is.na(fast[["D3"]]))
  expect_identical(is.na(fast), is.na(ref))
  expect_identical(is.na(fast), is.na(oracle))
  keep <- !is.na(fast)
  expect_lt(max(abs(fast[keep] - oracle[keep])), 1e-12)
})

test_that("semipartial results agree across run_regional and run_searchlight engines", {
  testthat::skip_if_not_installed("neuroim2")
  set.seed(1402)
  ds <- gen_sample_dataset(c(5, 5, 5), 32, blocks = 4)
  D1 <- dist(matrix(rnorm(32 * 4), 32))
  D2 <- dist(matrix(rnorm(32 * 4), 32))
  N1 <- dist(matrix(rnorm(32 * 3), 32))
  rdes <- rsa_design(~ D1 + D2, list(D1 = D1, D2 = D2, block = ds$design$block_var),
                     block_var = "block", nuisance = list(nuis = N1))
  build <- function() {
    rsa_model(ds$dataset, rdes, distmethod = "pearson", regtype = "lm",
              semipartial = TRUE, check_collinearity = FALSE)
  }
  ms <- build()
  expect_true(is.list(ms$.fast_kernel$lm))
  ms_ref <- rMVPA:::.with_reference_paths("rsa_fast_kernel", build())

  region_mask <- neuroim2::NeuroVol(sample(1:5, size = length(ds$dataset$mask), replace = TRUE),
                                    neuroim2::space(ds$dataset$mask))
  reg_fast <- suppressWarnings(run_regional(ms, region_mask))
  reg_ref <- rMVPA:::.with_reference_paths("rsa_fast_kernel",
                                           suppressWarnings(run_regional(ms_ref, region_mask)))
  expect_equal(reg_fast$performance_table, reg_ref$performance_table, tolerance = 1e-12)

  sl_fast <- suppressWarnings(run_searchlight(ms, radius = 2, backend = "default"))
  sl_legacy <- suppressWarnings(run_searchlight(ms, radius = 2, engine = "legacy", backend = "default"))
  sl_ref <- rMVPA:::.with_reference_paths("rsa_fast_kernel",
    suppressWarnings(run_searchlight(ms_ref, radius = 2, engine = "legacy", backend = "default")))
  expect_identical(attr(sl_fast, "searchlight_engine"), "rsa_fast")
  expect_identical(sort(names(sl_fast$results)), sort(names(sl_ref$results)))
  for (m in names(sl_ref$results)) {
    f <- as.numeric(neuroim2::values(sl_fast$results[[m]]))
    l <- as.numeric(neuroim2::values(sl_legacy$results[[m]]))
    r <- as.numeric(neuroim2::values(sl_ref$results[[m]]))
    expect_identical(f, l, info = m)
    expect_lt(max(abs(f - r), na.rm = TRUE), 1e-12)
    expect_identical(is.na(f), is.na(r), info = m)
  }
})

## ===========================================================================
## repnet_model / repmed_model
## ===========================================================================

rp_nvox <- 64L

rp_roi_fixture <- function(E_items, items, reps_per_item = 2L) {
  K <- length(items)
  n_obs <- K * reps_per_item
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = n_obs, blocks = reps_per_item, nlevels = 2)
  mvdes <- mvpa_design(train_design = ds$design$train_design,
                       y_train = ds$design$y_train, block_var = ds$design$block_var)
  mvdes$train_design$.rownum <- rep(items, each = reps_per_item)
  vox <- which(ds$dataset$mask > 0)
  samp <- data_sample(ds$dataset, vox)
  roi <- as_roi(samp, ds$dataset)
  roi$train_roi@.Data <- E_items[as.character(mvdes$train_design$.rownum), , drop = FALSE]
  list(ds = ds, mvdes = mvdes, roi = roi)
}

# Independent oracle: lm() on the given predictors; sp = sign(t) sqrt(t^2 MSE / TSS).
lm_terms_oracle <- function(y, X) {
  d <- cbind(data.frame(.y = y), as.data.frame(X, stringsAsFactors = FALSE))
  fit <- stats::lm(.y ~ ., data = d)
  smry <- summary(fit)
  cf <- stats::coef(smry)
  nm <- names(X)
  yy <- fit$model$.y
  tss <- sum((yy - mean(yy))^2)
  tv <- cf[nm, "t value"]
  list(beta = setNames(cf[nm, "Estimate"], nm),
       sp = setNames(sign(tv) * sqrt(tv^2 * smry$sigma^2 / tss), nm))
}

rp_seed_conf <- function(seed, K, conf_names = "conf1") {
  set.seed(seed)
  items <- paste0("it", seq_len(K))
  E <- matrix(rnorm(K * rp_nvox), nrow = K, ncol = rp_nvox)
  rownames(E) <- items
  F_seed <- matrix(rnorm(K * 6L), nrow = K, ncol = 6L)
  D_seed <- 1 - cor(t(F_seed))
  rownames(D_seed) <- colnames(D_seed) <- items
  confs <- lapply(conf_names, function(nm) {
    F_c <- matrix(rnorm(K * 5L), nrow = K, ncol = 5L)
    D <- 1 - cor(t(F_c))
    rownames(D) <- colnames(D) <- items
    D
  })
  names(confs) <- conf_names
  list(items = items, E = E, D_seed = D_seed, confs = confs)
}

test_that("repnet beta and semipartial terms equal the lm() oracle (two confounds)", {
  sc <- rp_seed_conf(401L, K = 14L, conf_names = c("conf1", "conf2"))
  env <- rp_roi_fixture(sc$E, sc$items)
  rn_des <- repnet_design(design = env$mvdes, key_var = ~ .rownum,
                          seed_rdm = sc$D_seed, confound_rdms = sc$confs)
  mspec <- repnet_model(dataset = env$ds$dataset, design = env$mvdes, repnet_des = rn_des,
                        distfun = cordist("pearson"), simfun = "pearson")
  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  perf <- res$performance[[1]]
  expect_false(res$error)

  D_roi <- 1 - cor(t(sc$E))
  tl <- lower.tri(D_roi)
  X <- data.frame(seed = sc$D_seed[tl], conf1 = sc$confs$conf1[tl], conf2 = sc$confs$conf2[tl])
  oracle <- lm_terms_oracle(D_roi[tl], X)

  for (term in names(X)) {
    expect_equal(unname(perf[[paste0("beta_", term)]]), unname(oracle$beta[[term]]),
                 tolerance = 1e-10, info = paste("beta", term))
    expect_equal(unname(perf[[paste0("sp_", term)]]), unname(oracle$sp[[term]]),
                 tolerance = 1e-10, info = paste("sp", term))
  }
})

test_that("repnet beta_seed without confounds equals the OLS slope and sp_seed equals the seed correlation", {
  sc <- rp_seed_conf(402L, K = 12L, conf_names = character(0))
  env <- rp_roi_fixture(sc$E, sc$items)
  rn_des <- repnet_design(design = env$mvdes, key_var = ~ .rownum, seed_rdm = sc$D_seed)
  mspec <- repnet_model(dataset = env$ds$dataset, design = env$mvdes, repnet_des = rn_des,
                        distfun = cordist("pearson"), simfun = "pearson")
  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  perf <- res$performance[[1]]
  D_roi <- 1 - cor(t(sc$E))
  tl <- lower.tri(D_roi)
  d_seed <- sc$D_seed[tl]
  oracle <- lm_terms_oracle(D_roi[tl], data.frame(seed = d_seed))
  expect_equal(unname(perf[["beta_seed"]]), unname(oracle$beta[["seed"]]), tolerance = 1e-10)
  expect_equal(unname(perf[["sp_seed"]]), unname(oracle$sp[["seed"]]), tolerance = 1e-10)
  # One predictor: sp equals the Pearson correlation with the ROI RDM, signed.
  expect_equal(unname(perf[["sp_seed"]]), unname(stats::cor(d_seed, D_roi[tl])), tolerance = 1e-10)
})

test_that("repnet confound changed after construction is honoured (fixed-design cache)", {
  sc <- rp_seed_conf(403L, K = 13L, conf_names = "conf1")
  env <- rp_roi_fixture(sc$E, sc$items)
  rn_des <- repnet_design(design = env$mvdes, key_var = ~ .rownum, seed_rdm = sc$D_seed)
  mspec <- repnet_model(dataset = env$ds$dataset, design = env$mvdes, repnet_des = rn_des,
                        distfun = cordist("pearson"), simfun = "pearson")
  mspec$confound_rdms <- sc$confs
  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  perf <- res$performance[[1]]
  D_roi <- 1 - cor(t(sc$E))
  tl <- lower.tri(D_roi)
  oracle <- lm_terms_oracle(D_roi[tl], data.frame(seed = sc$D_seed[tl],
                                                   conf1 = sc$confs$conf1[tl]))
  expect_equal(unname(perf[["beta_seed"]]), unname(oracle$beta[["seed"]]), tolerance = 1e-10)
  expect_equal(unname(perf[["beta_conf1"]]), unname(oracle$beta[["conf1"]]), tolerance = 1e-10)
})

test_that("repnet with NA confound entries drops incomplete pairs without erroring", {
  sc <- rp_seed_conf(404L, K = 12L, conf_names = "conf1")
  sc$confs$conf1[2, 5] <- NA
  sc$confs$conf1[5, 2] <- NA
  env <- rp_roi_fixture(sc$E, sc$items)
  rn_des <- repnet_design(design = env$mvdes, key_var = ~ .rownum, seed_rdm = sc$D_seed,
                          confound_rdms = sc$confs)
  mspec <- repnet_model(dataset = env$ds$dataset, design = env$mvdes, repnet_des = rn_des,
                        distfun = cordist("pearson"), simfun = "pearson")
  res <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)
  expect_false(res$error)
  perf <- res$performance[[1]]

  D_roi <- 1 - cor(t(sc$E))
  tl <- lower.tri(D_roi)
  X <- data.frame(seed = sc$D_seed[tl], conf1 = sc$confs$conf1[tl])
  oracle <- lm_terms_oracle(D_roi[tl], X)
  expect_true(all(is.finite(perf[c("beta_seed", "beta_conf1", "sp_seed", "sp_conf1")])))
  expect_equal(unname(perf[["beta_seed"]]), unname(oracle$beta[["seed"]]), tolerance = 1e-10)
  expect_equal(unname(perf[["sp_seed"]]), unname(oracle$sp[["seed"]]), tolerance = 1e-10)
})

test_that("repnet seed_conf_r2 equals the lm() R-squared of seed on confounds", {
  sc <- rp_seed_conf(405L, K = 13L, conf_names = c("conf1", "conf2"))
  # Make the seed partly a confound so the shared share is substantial.
  sc$D_seed <- sc$D_seed + 0.8 * sc$confs$conf1
  env <- rp_roi_fixture(sc$E, sc$items)
  rn_des <- repnet_design(design = env$mvdes, key_var = ~ .rownum,
                          seed_rdm = sc$D_seed, confound_rdms = sc$confs)
  mspec <- repnet_model(dataset = env$ds$dataset, design = env$mvdes, repnet_des = rn_des,
                        distfun = cordist("pearson"), simfun = "pearson")
  tl <- lower.tri(sc$D_seed)
  oracle <- summary(stats::lm(sc$D_seed[tl] ~ sc$confs$conf1[tl] + sc$confs$conf2[tl]))$r.squared
  expect_gt(oracle, 0.2)

  perf <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)$performance[[1]]
  expect_equal(unname(perf[["seed_conf_r2"]]), oracle, tolerance = 1e-10)
  # beta_seed stays the regression coefficient (not the semi-partial)
  D_roi <- 1 - cor(t(sc$E))
  X <- data.frame(seed = sc$D_seed[tl], conf1 = sc$confs$conf1[tl], conf2 = sc$confs$conf2[tl])
  expect_equal(unname(perf[["beta_seed"]]),
               unname(lm_terms_oracle(D_roi[tl], X)$beta[["seed"]]), tolerance = 1e-10)

  # The uncached path (no fixed design) gives the same value.
  mspec$.fixed_design <- NULL
  perf2 <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)$performance[[1]]
  expect_equal(unname(perf2[["seed_conf_r2"]]), oracle, tolerance = 1e-10)
})

test_that("repnet seed_conf_r2 is absent without confounds and 1 for a collinear seed", {
  sc <- rp_seed_conf(406L, K = 12L, conf_names = character(0))
  env <- rp_roi_fixture(sc$E, sc$items)
  rn_des <- repnet_design(design = env$mvdes, key_var = ~ .rownum, seed_rdm = sc$D_seed)
  mspec <- repnet_model(dataset = env$ds$dataset, design = env$mvdes, repnet_des = rn_des)
  perf <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)$performance[[1]]
  expect_false("seed_conf_r2" %in% names(perf))

  sc <- rp_seed_conf(407L, K = 12L, conf_names = "conf1")
  sc$confs$conf1 <- 2 * sc$D_seed + 1
  rownames(sc$confs$conf1) <- colnames(sc$confs$conf1) <- sc$items
  env <- rp_roi_fixture(sc$E, sc$items)
  rn_des <- repnet_design(design = env$mvdes, key_var = ~ .rownum,
                          seed_rdm = sc$D_seed, confound_rdms = sc$confs)
  mspec <- repnet_model(dataset = env$ds$dataset, design = env$mvdes, repnet_des = rn_des)
  perf <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)$performance[[1]]
  expect_equal(unname(perf[["seed_conf_r2"]]), 1, tolerance = 1e-10)
})

test_that("repnet seed_conf_r2 uses only pairs with complete confounds", {
  sc <- rp_seed_conf(408L, K = 12L, conf_names = "conf1")
  sc$confs$conf1[2, 5] <- NA
  sc$confs$conf1[5, 2] <- NA
  env <- rp_roi_fixture(sc$E, sc$items)
  rn_des <- repnet_design(design = env$mvdes, key_var = ~ .rownum,
                          seed_rdm = sc$D_seed, confound_rdms = sc$confs)
  mspec <- repnet_model(dataset = env$ds$dataset, design = env$mvdes, repnet_des = rn_des)
  perf <- process_roi(mspec, env$roi, rnum = 1L, center_global_id = NA)$performance[[1]]
  tl <- lower.tri(sc$D_seed)
  oracle <- summary(stats::lm(sc$D_seed[tl] ~ sc$confs$conf1[tl]))$r.squared
  expect_equal(unname(perf[["seed_conf_r2"]]), oracle, tolerance = 1e-10)
})

# repmed oracle: paths a, b, c' and c by lm(), Sobel test with lm() standard errors
repmed_oracle <- function(m, x, y, conf = NULL) {
  X_a <- data.frame(x = x)
  X_b <- data.frame(m = m, x = x)
  if (!is.null(conf)) {
    X_a <- cbind(X_a, conf)
    X_b <- cbind(X_b, conf)
  }
  X_c <- if (is.null(conf)) data.frame(x = x) else cbind(data.frame(x = x), conf)
  fa <- stats::lm(m ~ ., data = cbind(data.frame(m = m), X_a))
  fb <- stats::lm(y ~ ., data = cbind(data.frame(y = y), X_b))
  fc <- stats::lm(y ~ ., data = cbind(data.frame(y = y), X_c))
  ca <- stats::coef(summary(fa)); cb <- stats::coef(summary(fb)); cc <- stats::coef(summary(fc))
  a <- ca["x", "Estimate"]; sa <- ca["x", "Std. Error"]
  b <- cb["m", "Estimate"]; sb <- cb["m", "Std. Error"]
  indirect <- a * b
  se_ind <- sqrt(b^2 * sa^2 + a^2 * sb^2)
  z <- indirect / se_ind
  list(med_a = a, med_b = b, med_cprime = cb["x", "Estimate"], med_c = cc["x", "Estimate"],
       med_indirect = indirect, med_sobel_z = z, med_sobel_p = 2 * stats::pnorm(-abs(z)))
}

rp_repmed_fixture <- function(seed, K = 16L, with_conf = TRUE) {
  set.seed(seed)
  items <- paste0("it", seq_len(K))
  z_x <- rnorm(K)
  z_m <- 0.7 * z_x + rnorm(K)
  z_y <- 0.5 * z_m + 0.3 * z_x + rnorm(K)
  rdm <- function(v) { D <- as.matrix(dist(v)); rownames(D) <- colnames(D) <- items; D }
  X_rdm <- rdm(z_x); Y_rdm <- rdm(z_y)
  conf_rdms <- if (with_conf) list(conf1 = rdm(rnorm(K))) else NULL
  E_m <- matrix(z_m, nrow = K, ncol = rp_nvox) + matrix(rnorm(K * rp_nvox), nrow = K, ncol = rp_nvox)
  rownames(E_m) <- items
  list(items = items, E = E_m, X_rdm = X_rdm, Y_rdm = Y_rdm, confs = conf_rdms)
}

rp_repmed_model <- function(fx, Y_rdm = fx$Y_rdm) {
  env <- rp_roi_fixture(fx$E, fx$items)
  md_des <- repmed_design(items = fx$items, X_rdm = fx$X_rdm, Y_rdm = Y_rdm,
                          confound_rdms = fx$confs)
  mspec <- repmed_model(dataset = env$ds$dataset, design = env$mvdes, repmed_des = md_des,
                        key_var = ~ .rownum, distfun = "euclidean")
  list(mspec = mspec, roi = env$roi)
}

test_that("repmed paths, total effect and Sobel test equal the lm() oracle", {
  for (with_conf in c(TRUE, FALSE)) {
    fx <- rp_repmed_fixture(501L, with_conf = with_conf)
    rm <- rp_repmed_model(fx)
    res <- process_roi(rm$mspec, rm$roi, rnum = 1L, center_global_id = NA)
    perf <- res$performance[[1]]
    expect_false(res$error, info = paste("conf", with_conf))

    m <- as.vector(stats::dist(fx$E))
    tl <- lower.tri(fx$X_rdm)
    x <- fx$X_rdm[tl]; y <- fx$Y_rdm[tl]
    conf <- if (with_conf) fx$confs$conf1[tl] else NULL
    oracle <- repmed_oracle(m, x, y, conf)
    for (nm in names(oracle)) {
      expect_equal(unname(perf[[nm]]), unname(oracle[[nm]]), tolerance = 1e-10,
                   info = paste(nm, "conf", with_conf))
    }
  }
})

test_that("repmed fits path c once at construction and reuses it per ROI", {
  fx <- rp_repmed_fixture(502L, with_conf = TRUE)
  rm <- rp_repmed_model(fx)
  orig <- rMVPA:::.rep_path_fit
  calls <- 0L
  testthat::local_mocked_bindings(
    .rep_path_fit = function(...) { calls <<- calls + 1L; orig(...) },
    .package = "rMVPA"
  )
  res <- process_roi(rm$mspec, rm$roi, rnum = 1L, center_global_id = NA)
  expect_false(res$error)
  expect_equal(calls, 2L)  # paths a and b only; path c comes from the cache

  # Changing the outcome RDM on the existing spec invalidates the cached path c:
  # it is refitted, and the result follows the new outcome.
  Y_new <- rp_repmed_fixture(503L, with_conf = TRUE)$Y_rdm
  rm$mspec$Y_rdm <- Y_new
  calls <- 0L
  res2 <- process_roi(rm$mspec, rm$roi, rnum = 1L, center_global_id = NA)
  perf2 <- res2$performance[[1]]
  expect_false(res2$error)
  expect_equal(calls, 3L)
  tl <- lower.tri(fx$X_rdm)
  m <- as.vector(stats::dist(fx$E))
  oracle2 <- repmed_oracle(m, fx$X_rdm[tl], Y_new[tl], fx$confs$conf1[tl])
  expect_equal(unname(perf2[["med_c"]]), oracle2$med_c, tolerance = 1e-10)
  expect_equal(unname(perf2[["med_b"]]), oracle2$med_b, tolerance = 1e-10)
})

## ===========================================================================
## table_to_rdm: match()-based fill vs the former list-lookup implementation
## ===========================================================================

capture_table_msgs <- function(expr) {
  msgs <- character(0)
  val <- withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = val, messages = msgs)
}

capture_table_result <- function(expr) {
  tryCatch(capture_table_msgs(expr), error = function(e) list(error = conditionMessage(e)))
}

random_similarity_table <- function(labels, n_rows, seed) {
  set.seed(seed)
  table_labels <- c(labels, paste0("extra", 1:2))
  data.frame(
    label1 = sample(table_labels, n_rows, replace = TRUE),
    label2 = sample(table_labels, n_rows, replace = TRUE),
    similarity = c(stats::runif(n_rows - 3L, -0.2, 1.2), NA, 0.4, 0.9),
    stringsAsFactors = FALSE
  )
}

test_that("table_to_rdm matches the former implementation on random sparse tables", {
  labels <- c("cat", "dog", "bird", "car", "plane", "boat", "fish")
  grid <- expand.grid(symmetric = c(TRUE, FALSE), self = c(1, 0.3, NA),
                      default = c(0, 0.25), as_dist = c(TRUE, FALSE),
                      stringsAsFactors = FALSE)
  for (seed in 1:6) {
    tab <- random_similarity_table(labels, n_rows = 12L + seed * 3L, seed = seed)
    for (g in seq_len(nrow(grid))) {
      args <- c(list(similarity_table = tab, labels = labels), as.list(grid[g, ]))
      old <- capture_table_result(do.call(table_to_rdm_reference, args))
      new <- capture_table_result(do.call(table_to_rdm, args))
      expect_identical(new, old, info = paste("seed", seed, "row", g))
    }
  }
})

test_that("table_to_rdm keeps first-occurrence duplicates and warnings", {
  labels <- c("a", "b", "c", "d")
  tab <- data.frame(label1 = c("a", "b", "a", "c", "a", "d"),
                    label2 = c("b", "a", "b", "d", "c", "d"),
                    similarity = c(0.2, 0.9, 0.7, 0.5, 0.3, 0.1))
  for (sym in c(TRUE, FALSE)) {
    old <- capture_table_result(table_to_rdm_reference(tab, labels, symmetric = sym, as_dist = FALSE))
    new <- capture_table_result(table_to_rdm(tab, labels, symmetric = sym, as_dist = FALSE))
    expect_identical(new, old, info = paste("symmetric", sym))
    expect_true(any(grepl("Duplicate", new$messages)))
  }
})

test_that("table_to_rdm matches the former implementation for repeated labels and key-separator collisions", {
  labels <- c("a", "b", "a", "c")
  tab <- data.frame(label1 = c("a___b", "a", "x"), label2 = c("c", "b___c", "a"),
                    similarity = c(0.6, 0.2, 0.5))
  for (sym in c(TRUE, FALSE)) {
    old <- capture_table_result(table_to_rdm_reference(tab, labels, symmetric = sym, as_dist = FALSE))
    new <- capture_table_result(table_to_rdm(tab, labels, symmetric = sym, as_dist = FALSE))
    expect_identical(new, old, info = paste("symmetric", sym))
  }
})

test_that("table_to_rdm matches the former implementation for NA labels", {
  labels <- c("a", "b", "c")
  tab <- data.frame(label1 = c("a", NA, "b"), label2 = c("b", "c", "c"),
                    similarity = c(0.4, 0.6, 0.8))
  for (sym in c(TRUE, FALSE)) {
    old <- capture_table_result(table_to_rdm_reference(tab, labels, symmetric = sym, as_dist = FALSE))
    new <- capture_table_result(table_to_rdm(tab, labels, symmetric = sym, as_dist = FALSE))
    expect_identical(new, old, info = paste("symmetric", sym))
  }
})

test_that("table_to_rdm matches the former implementation for empty tables and zero labels", {
  empty <- data.frame(label1 = character(0), label2 = character(0), similarity = numeric(0))
  one <- data.frame(label1 = "a", label2 = "a", similarity = 0.5)
  cases <- list(list(empty, c("a", "b")), list(one, "a"),
                list(one, character(0)), list(one, c("a", "b")))
  for (cs in cases) for (sym in c(TRUE, FALSE)) for (ad in c(TRUE, FALSE)) {
    old <- capture_table_result(table_to_rdm_reference(cs[[1]], cs[[2]], symmetric = sym, as_dist = ad))
    new <- capture_table_result(table_to_rdm(cs[[1]], cs[[2]], symmetric = sym, as_dist = ad))
    expect_identical(new, old, info = paste(length(cs[[2]]), sym, ad))
  }
})

test_that("table_to_rdm completes for 300 labels", {
  labels <- paste0("L", seq_len(300))
  tab <- random_similarity_table(labels, n_rows = 2000L, seed = 99L)
  rdm <- suppressWarnings(table_to_rdm(tab, labels, as_dist = FALSE))
  expect_identical(dim(rdm), c(300L, 300L))
  expect_identical(unname(rdm), unname(t(rdm)))
  expect_identical(unname(diag(rdm)), rep(0, 300))
})

