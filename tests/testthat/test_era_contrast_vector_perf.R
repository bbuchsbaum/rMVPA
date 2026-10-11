context("ERA, contrast and vector RSA performance paths")

# ---------------------------------------------------------------------------
# Fixtures and golden outputs.
#
# fixtures/era_contrast_vector_perf_golden.rds holds the outputs of
# perf_golden_outputs() computed by the implementation before the performance
# changes (git HEAD prior to this change). The tests below compare the current
# implementation against that file. Regenerate it only from the reference
# implementation.
# ---------------------------------------------------------------------------

perf_era_fixture <- function(K = 12L, dims = c(4L, 4L, 4L), seed = 731L,
                             na_conf = FALSE, const_conf = FALSE,
                             dup_conf = FALSE) {
  set.seed(seed)
  keys <- sprintf("item_%02d", seq_len(K))
  p <- prod(dims)
  E <- matrix(stats::rnorm(K * p), nrow = K)
  R <- E + matrix(stats::rnorm(K * p, sd = 0.35), nrow = K)
  E[, 1L] <- 0
  R[, 1L] <- seq_len(K)

  enc_vec <- neuroim2::NeuroVec(
    array(as.numeric(t(E)), dim = c(dims, K)),
    neuroim2::NeuroSpace(c(dims, K))
  )
  ret_vec <- neuroim2::NeuroVec(
    array(as.numeric(t(R)), dim = c(dims, K)),
    neuroim2::NeuroSpace(c(dims, K))
  )
  mask <- neuroim2::NeuroVol(array(1, dim = dims), neuroim2::NeuroSpace(dims))

  train_tab <- data.frame(
    item = factor(keys, levels = keys),
    vividness = NA_real_,
    retrieval_run = factor(rep(c("run_1", "run_2"), length.out = K)),
    trial_order = seq_len(K)
  )
  test_tab <- train_tab
  test_tab$vividness <- sample(seq_len(K))

  dataset <- mvpa_dataset(enc_vec, test_data = ret_vec, mask = mask)
  design <- mvpa_design(train_tab, test_design = test_tab,
                        y_train = ~ item, y_test = ~ item)

  sym <- function(v) {
    M <- as.matrix(stats::dist(matrix(v, ncol = 1)))
    dimnames(M) <- list(keys, keys)
    M
  }
  conf <- list(
    time  = sym(seq_len(K) + stats::rnorm(K)),
    block = sym(as.numeric(rep(1:3, length.out = K)))
  )
  if (const_conf) conf$constant <- sym(rep(1, K))
  if (dup_conf) conf$time_copy <- conf$time
  if (na_conf) conf$time[1, 2] <- NA_real_

  list(
    dataset = dataset, design = design, keys = keys, conf = conf,
    run_enc = rep(c("e1", "e2", "e3"), length.out = K),
    run_ret = rep(c("r1", "r1", "r2", "r2", "r3"), length.out = K)
  )
}

perf_era_model <- function(fx, partial_against = "run",
                           components = c("item", "identification", "geometry")) {
  suppressWarnings(era_rsa_model(
    dataset = fx$dataset,
    design = fx$design,
    key_var = ~ item,
    pairing = "one_to_one",
    era_components = components,
    confound_rdms = fx$conf,
    item_run_enc = stats::setNames(fx$run_enc, fx$keys),
    item_run_ret = stats::setNames(fx$run_ret, fx$keys),
    partial_against = partial_against
  ))
}

# Returns list(design = msreve_design, dataset = mvpa_dataset) for a
# 4-condition contrast searchlight.
perf_contrast_fixture <- function(seed = 4242L, keep_intra_run = FALSE,
                                  with_blocks = TRUE, with_nuisance = TRUE) {
  set.seed(seed)
  n_cond <- 4L
  n_blocks <- 3L
  reps <- 2L
  cond_names <- paste0("cond_", seq_len(n_cond))
  cond <- rep(cond_names, each = n_blocks * reps)
  run <- rep(seq_len(n_blocks), times = n_cond * reps)
  n_trials <- length(cond)
  dims <- c(4L, 4L, 4L)

  signal <- matrix(stats::rnorm(n_trials * prod(dims)), nrow = n_trials)
  for (i in seq_len(n_trials)) {
    signal[i, 1:8] <- signal[i, 1:8] + match(cond[i], cond_names)
  }
  vec <- neuroim2::NeuroVec(
    array(as.numeric(t(signal)), dim = c(dims, n_trials)),
    neuroim2::NeuroSpace(c(dims, n_trials))
  )
  mask <- neuroim2::NeuroVol(array(1, dim = dims), neuroim2::NeuroSpace(dims))
  dataset <- mvpa_dataset(vec, mask = mask)

  train_tab <- data.frame(
    cond = factor(cond, levels = cond_names),
    run = factor(run)
  )
  mvpa_des <- if (with_blocks) {
    mvpa_design(train_tab, y_train = ~ cond, block_var = ~ run)
  } else {
    mvpa_design(train_tab, y_train = ~ cond)
  }
  mvpa_des$keep_intra_run <- keep_intra_run

  C <- cbind(c(1, 1, -1, -1) / 2, c(1, -1, 0, 0) / sqrt(2))
  colnames(C) <- c("C1", "C2")
  rownames(C) <- cond_names

  nuisance <- NULL
  if (with_nuisance) {
    M <- matrix(stats::rnorm(n_cond * n_cond), n_cond)
    M <- (M + t(M)) / 2
    diag(M) <- 0
    dimnames(M) <- list(cond_names, cond_names)
    nuisance <- list(nuis_a = M)
  }
  list(
    design = msreve_design(mvpa_des, C, nuisance_rdms = nuisance),
    dataset = dataset
  )
}

perf_vector_fixture <- function(seed = 99L) {
  set.seed(seed)
  n <- 12L
  X <- matrix(stats::rnorm(n * 27L), nrow = n)
  labels <- paste0("L", seq_len(n))
  D <- as.matrix(stats::dist(matrix(stats::rnorm(n * 3L), n)))
  rownames(D) <- colnames(D) <- labels
  block <- rep(1:3, each = 4L)
  vec <- neuroim2::NeuroVec(
    array(as.numeric(t(X)), dim = c(3L, 3L, 3L, n)),
    neuroim2::NeuroSpace(c(3L, 3L, 3L, n))
  )
  mask <- neuroim2::NeuroVol(array(1, dim = c(3L, 3L, 3L)),
                             neuroim2::NeuroSpace(c(3L, 3L, 3L)))
  list(
    X = X, D = D, block = block,
    dataset = mvpa_dataset(vec, mask = mask),
    design = vector_rsa_design(D, factor(labels), block)
  )
}

# Named list of ERA searchlight maps (metric -> numeric vector).
perf_era_maps <- function(fx, partial_against = "run") {
  model <- perf_era_model(fx, partial_against = partial_against)
  res <- suppressWarnings(run_searchlight(
    model, radius = 1, method = "standard",
    engine = "legacy", backend = "default", preflight = "off"
  ))
  lapply(res$results, function(v) as.numeric(neuroim2::values(v)))
}

# Named list of contrast searchlight maps.
perf_contrast_maps <- function(fx) {
  spec <- contrast_rsa_model(
    dataset = fx$dataset,
    design = fx$design,
    output_metric = c("beta_delta", "beta_only", "delta_only"),
    check_collinearity = FALSE
  )
  res <- suppressWarnings(run_searchlight(
    spec, radius = 1, method = "standard",
    engine = "legacy", backend = "default", preflight = "off"
  ))
  lapply(res$results, function(v) as.numeric(neuroim2::values(v)))
}

# Permutation run used both for the golden file and the tests. The RNG is
# seeded after model construction and before trial scoring, so the stream
# consumed by evaluate_model() is fully determined by `seed`.
perf_vector_run <- function(vfx, distfun, rsa_simfun, nperm, seed) {
  vm <- vector_rsa_model(
    vfx$dataset, vfx$design, distfun = distfun, rsa_simfun = rsa_simfun,
    nperm = nperm, save_distributions = TRUE
  )
  set.seed(seed)
  scores <- rMVPA:::compute_trial_scores(vm, vfx$X)
  perm <- evaluate_model.vector_rsa_model(
    vm, predicted = NULL, observed = scores,
    roi_data_for_perm = vfx$X, nperm = nperm, save_distributions = TRUE
  )
  list(scores = scores, perm = perm$permutation_results, rng_after = .Random.seed)
}

# Everything compared against the golden file.
perf_golden_outputs <- function() {
  out <- list(era = list(), contrast = list(), vector = list())

  out$era[["base_run"]] <- perf_era_maps(perf_era_fixture())
  out$era[["base_time_run"]] <- perf_era_maps(perf_era_fixture(), c("time", "run"))
  out$era[["nonfinite_conf"]] <- perf_era_maps(perf_era_fixture(na_conf = TRUE, seed = 17L))
  out$era[["constant_conf"]] <- perf_era_maps(perf_era_fixture(const_conf = TRUE, seed = 99L))
  out$era[["duplicate_conf"]] <- perf_era_maps(perf_era_fixture(dup_conf = TRUE, seed = 5L))

  out$contrast[["blocks_nuisance"]] <- perf_contrast_maps(perf_contrast_fixture())
  out$contrast[["blocks_only"]] <- perf_contrast_maps(
    perf_contrast_fixture(seed = 77L, with_nuisance = FALSE)
  )

  vfx <- perf_vector_fixture()
  # Permutation nulls are checked against a recomputing oracle instead: the
  # null is now drawn within blocks (it was drawn across blocks before).
  out$vector$second_order <- second_order_similarity(
    distfun = cordist(), X = vfx$X, Dref = vfx$D, block_var = vfx$block,
    method = "spearman"
  )
  out
}

perf_golden_path <- function() {
  testthat::test_path("fixtures", "era_contrast_vector_perf_golden.rds")
}

# Maps must agree on NA pattern and values to `tolerance`.
expect_maps_equal <- function(actual, expected, tolerance = 1e-10, info = NULL) {
  expect_identical(names(actual), names(expected), info = info)
  for (nm in names(expected)) {
    expect_identical(is.na(actual[[nm]]), is.na(expected[[nm]]),
                     info = paste(info, nm))
    expect_equal(actual[[nm]], expected[[nm]], tolerance = tolerance,
                 info = paste(info, nm))
  }
}

# ---------------------------------------------------------------------------
# ERA: lm oracles at the helper level
# ---------------------------------------------------------------------------

test_that("ERA confound betas and semi-partials match an lm oracle", {
  fx <- perf_era_fixture()
  model <- perf_era_model(fx)
  cache <- model$era_geometry_cache
  K <- length(fx$keys)
  n <- K * (K - 1L) / 2L
  lower <- function(M) M[lower.tri(M)]

  set.seed(5)
  dE <- stats::rnorm(n)
  dR <- 0.4 * dE + stats::rnorm(n)
  conf <- data.frame(time = lower(fx$conf$time), block = lower(fx$conf$block))
  dat <- cbind(enc_geom = dE, conf)

  fast <- rMVPA:::.era_geometry_fast_confound_fit(cache, dE, dR)
  expect_false(is.null(fast))

  full <- stats::lm(dR ~ ., data = dat)
  cf <- stats::coef(summary(full))
  terms <- c("enc_geom", "time", "block")
  expect_equal(fast$beta, stats::setNames(cf[terms, "Estimate"], paste0("beta_", terms)),
               tolerance = 1e-10)

  # Semi-partial correlation of term j is sign(t_j) * sqrt(R2_full - R2_without_j).
  r2_full <- summary(full)$r.squared
  expected_sp <- vapply(terms, function(term) {
    reduced <- stats::lm(dR ~ ., data = dat[, setdiff(terms, term), drop = FALSE])
    sign(cf[term, "t value"]) * sqrt(r2_full - summary(reduced)$r.squared)
  }, numeric(1L))
  expect_equal(fast$sp, stats::setNames(expected_sp, paste0("sp_", terms)),
               tolerance = 1e-10)
})

test_that("ERA partial geometry correlations match an lm residual oracle", {
  fx <- perf_era_fixture()
  model <- perf_era_model(fx)
  cache <- model$era_geometry_cache
  K <- length(fx$keys)
  n <- K * (K - 1L) / 2L
  lower <- function(M) M[lower.tri(M)]
  all_conf <- rMVPA:::.era_geometry_confound_rdms(
    fx$conf, model$item_run_enc, model$item_run_ret, fx$keys
  )
  Z <- data.frame(run_enc = lower(as.matrix(all_conf$run_enc)),
                  run_ret = lower(as.matrix(all_conf$run_ret)))

  set.seed(6)
  for (rep_i in 1:3) {
    dE <- stats::rnorm(n)
    dR <- 0.5 * dE + stats::rnorm(n)
    res_e <- stats::resid(stats::lm(dE ~ run_enc + run_ret, data = Z))
    res_r <- stats::resid(stats::lm(dR ~ run_enc + run_ret, data = Z))
    for (method in c("pearson", "spearman")) {
      got <- rMVPA:::.era_geometry_partial_cor_cached(
        cache, "run", dE, dR, all_conf, fx$keys, method
      )
      expect_equal(got, stats::cor(res_e, res_r, method = method),
                   tolerance = 1e-10, info = paste(rep_i, method))
    }
  }
})

test_that("ERA confound regressions fall back on rank deficiency and non-finite pairs", {
  K <- 12L
  n <- K * (K - 1L) / 2L
  set.seed(8)
  dE <- stats::rnorm(n)
  dR <- stats::rnorm(n)

  # Constant confound: rank-deficient step-1 design, so the cached path declines.
  fx_const <- perf_era_fixture(const_conf = TRUE, seed = 99L)
  m_const <- perf_era_model(fx_const)
  expect_null(m_const$era_geometry_cache$step1)
  expect_null(rMVPA:::.era_geometry_fast_confound_fit(
    m_const$era_geometry_cache, dE, dR
  ))

  # Duplicated confound: 'time_copy' is identical to 'time', so any design
  # that includes both is rank deficient and must fall back to lm().
  fx_dup <- perf_era_fixture(dup_conf = TRUE, seed = 5L)
  pa_dup <- c("time", "run")
  m_dup <- perf_era_model(fx_dup, partial_against = pa_dup)
  expect_null(m_dup$era_geometry_cache$step1)
  all_dup <- rMVPA:::.era_geometry_confound_rdms(
    fx_dup$conf, m_dup$item_run_enc, m_dup$item_run_ret, fx_dup$keys
  )
  entry_dup <- m_dup$era_geometry_cache$partial[[paste(pa_dup, collapse = "\r")]]
  expect_identical(entry_dup$status, "fallback")
  got_dup <- rMVPA:::.era_geometry_partial_cor_cached(
    m_dup$era_geometry_cache, pa_dup, dE, dR, all_dup, fx_dup$keys, "pearson"
  )
  expected_dup <- rMVPA:::.era_partial_geometry_cor(
    dE, dR, all_dup, fx_dup$keys, pa_dup, "pearson"
  )
  expect_identical(is.na(got_dup), is.na(expected_dup))
  expect_equal(got_dup, expected_dup, tolerance = 1e-10)

  # Non-finite pairs: the cached step-1 path declines and the partial path
  # reproduces the lm() result, including any NA.
  fx <- perf_era_fixture()
  model <- perf_era_model(fx)
  all_conf <- rMVPA:::.era_geometry_confound_rdms(
    fx$conf, model$item_run_enc, model$item_run_ret, fx$keys
  )
  dE_na <- dE
  dE_na[3L] <- NA_real_
  expect_null(rMVPA:::.era_geometry_fast_confound_fit(
    model$era_geometry_cache, dE_na, dR
  ))
  got_na <- rMVPA:::.era_geometry_partial_cor_cached(
    model$era_geometry_cache, "run", dE_na, dR, all_conf, fx$keys, "spearman"
  )
  expected_na <- rMVPA:::.era_partial_geometry_cor(
    dE_na, dR, all_conf, fx$keys, "run", "spearman"
  )
  expect_identical(is.na(got_na), is.na(expected_na))
  expect_equal(got_na, expected_na, tolerance = 1e-10)
})

test_that("ERA searchlight geometry maps equal the pre-change implementation", {
  golden <- readRDS(perf_golden_path())
  cases <- list(
    base_run = list(fx = perf_era_fixture(), pa = "run"),
    base_time_run = list(fx = perf_era_fixture(), pa = c("time", "run")),
    nonfinite_conf = list(fx = perf_era_fixture(na_conf = TRUE, seed = 17L), pa = "run"),
    constant_conf = list(fx = perf_era_fixture(const_conf = TRUE, seed = 99L), pa = "run"),
    duplicate_conf = list(fx = perf_era_fixture(dup_conf = TRUE, seed = 5L), pa = "run")
  )
  for (nm in names(cases)) {
    maps <- perf_era_maps(cases[[nm]]$fx, cases[[nm]]$pa)
    geom <- grep("^(geom_|beta_|sp_)", names(maps), value = TRUE)
    expect_true(length(geom) > 0L)
    expect_maps_equal(maps[geom], golden$era[[nm]][geom], info = nm)
  }
})

test_that("labelled confound RDMs are aligned by key in the cached and lm paths", {
  fx <- perf_era_fixture()
  maps_ref <- perf_era_maps(fx)

  # Same confound RDMs, rows and columns stored in reverse key order and
  # labelled by key: aligned by label, they must give identical maps.
  ord <- rev(seq_along(fx$keys))
  fx_shuf <- fx
  fx_shuf$conf$time <- fx$conf$time[ord, ord, drop = FALSE]
  fx_shuf$conf$block <- fx$conf$block[ord, ord, drop = FALSE]
  maps_shuf <- perf_era_maps(fx_shuf)
  geom <- grep("^(geom_|beta_|sp_)", names(maps_ref), value = TRUE)
  expect_maps_equal(maps_shuf[geom], maps_ref[geom], info = "shuffled labels")
})

# ---------------------------------------------------------------------------
# Contrast RSA
# ---------------------------------------------------------------------------

test_that("contrast intra-block mask matches the pairwise intersect() definition", {
  blocks <- list(
    a = c(1, 2, 3),
    b = c(4, 5),
    c = c(3, 6),
    d = numeric(0),
    e = c(NA, 7),
    f = c(7, 8),
    g = c(4, 9)
  )
  K <- length(blocks)
  intra <- matrix(FALSE, K, K)
  for (r in 1:(K - 1L)) {
    for (cc in (r + 1L):K) {
      if (length(intersect(blocks[[r]], blocks[[cc]])) > 0) {
        intra[r, cc] <- TRUE
        intra[cc, r] <- TRUE
      }
    }
  }
  expect_identical(rMVPA:::.contrast_intra_block_lower(blocks), intra[lower.tri(intra)])
})

test_that("contrast RSA searchlight outputs equal the pre-change implementation", {
  golden <- readRDS(perf_golden_path())
  maps <- perf_contrast_maps(perf_contrast_fixture())
  expect_maps_equal(maps, golden$contrast$blocks_nuisance, info = "blocks_nuisance")

  maps2 <- perf_contrast_maps(perf_contrast_fixture(seed = 77L, with_nuisance = FALSE))
  expect_maps_equal(maps2, golden$contrast$blocks_only, info = "blocks_only")
})

test_that("missing condition_block_list warning is emitted once per run, not per sphere", {
  # keep_intra_run = FALSE and no block variable: the condition-to-block
  # mapping is unavailable, so intra-block cells cannot be excluded.
  fx <- perf_contrast_fixture(with_blocks = FALSE, with_nuisance = FALSE, seed = 3L)
  n_warn <- 0L
  count_handler <- function(w) {
    if (grepl("condition_block_list", conditionMessage(w), fixed = TRUE)) {
      n_warn <<- n_warn + 1L
    }
    invokeRestart("muffleWarning")
  }
  spec <- withCallingHandlers(
    contrast_rsa_model(
      dataset = fx$dataset,
      design = fx$design,
      output_metric = c("beta_delta", "beta_only")
    ),
    warning = count_handler
  )
  withCallingHandlers(
    run_searchlight(spec, radius = 1, method = "standard",
                    engine = "legacy", backend = "default", preflight = "off"),
    warning = count_handler
  )
  expect_identical(n_warn, 1L)
})

# ---------------------------------------------------------------------------
# Vector RSA
# ---------------------------------------------------------------------------

test_that("reindexed permutation distances equal recomputed distances", {
  vfx <- perf_vector_fixture()
  set.seed(11)
  perm <- sample(nrow(vfx$X))
  for (distfun in list(cordist(method = "pearson"), cordist(method = "spearman"), euclidean())) {
    D <- as.matrix(pairwise_dist(distfun, vfx$X))
    expect_equal(
      unname(D[perm, perm]),
      unname(as.matrix(pairwise_dist(distfun, vfx$X[perm, , drop = FALSE]))),
      tolerance = 1e-12
    )
  }
})

test_that("vector RSA permutation nulls equal recomputed within-block nulls", {
  # Oracle: draw each permutation with the same RNG stream, permute the rows
  # of X and recompute the scores from scratch with second_order_similarity().
  vfx <- perf_vector_fixture()
  for (cfg in list(list(cordist(), "pearson", 40L, 2024L),
                   list(euclidean(), "spearman", 25L, 7L))) {
    got <- perf_vector_run(vfx, cfg[[1]], cfg[[2]], cfg[[3]], cfg[[4]])
    set.seed(cfg[[4]])
    scores <- second_order_similarity(cfg[[1]], vfx$X, vfx$D, vfx$block, cfg[[2]])
    null <- vapply(seq_len(cfg[[3]]), function(i) {
      perm <- rMVPA:::.vector_rsa_block_perm(vfx$block, nrow(vfx$X))
      expect_identical(vfx$block[perm], vfx$block)
      mean(second_order_similarity(cfg[[1]], vfx$X[perm, , drop = FALSE],
                                   vfx$D, vfx$block, cfg[[2]]), na.rm = TRUE)
    }, numeric(1))
    expect_equal(got$scores, scores, tolerance = 1e-12)
    expect_equal(unname(got$perm$permutation_distributions$rsa_score),
                 null, tolerance = 1e-12)
    expect_equal(unname(got$perm$p_values[["rsa_score"]]),
                 (sum(null >= mean(scores)) + 1) / (cfg[[3]] + 1))
    expect_identical(got$rng_after, .Random.seed)
  }
})

test_that("second_order_similarity output equals the pre-change implementation", {
  golden <- readRDS(perf_golden_path())
  vfx <- perf_vector_fixture()
  got <- second_order_similarity(
    distfun = cordist(), X = vfx$X, Dref = vfx$D, block_var = vfx$block,
    method = "spearman"
  )
  expect_equal(got, golden$vector$second_order, tolerance = 1e-12)
})
