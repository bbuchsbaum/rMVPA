# Permutation searchlights use the exact engines (prepared once, rescored per
# permutation). The null distribution and p-values must equal those of the
# per-ROI path for the same seed.
testthat::skip_if_not_installed("neuroim2")

perm_pair <- function(model, strategy) {
  set.seed(1301)
  ds <- gen_sample_dataset(c(4, 4, 4), 40, nlevels = 3, blocks = 4)
  ms <- mvpa_model(load_model(model), ds$dataset, ds$design, "classification",
                   crossval = blocked_cross_validation(ds$design$block_var))
  pc <- permutation_control(n_perm = 3, seed = 7, perm_strategy = strategy,
                            subsample = 0.5, diagnose = FALSE)
  obs <- run_searchlight(ms, radius = 2, engine = "legacy", backend = "default")
  list(
    fast = suppressWarnings(run_permutation_searchlight(ms, observed = obs, radius = 2,
                                                        perm_ctrl = pc, metric = "Accuracy")),
    ref = suppressWarnings(run_permutation_searchlight(ms, observed = obs, radius = 2,
                                                       perm_ctrl = pc, metric = "Accuracy",
                                                       engine = "legacy"))
  )
}

for (model in c("corclass", "naive_bayes", "sda_notune")) {
  for (strategy in c("iterate", "searchlight")) {
    local({
      m <- model; st <- strategy
      test_that(sprintf("engine permutations equal the per-ROI path (%s, %s)", m, st), {
        pp <- perm_pair(m, st)
        expect_equal(sort(unlist(pp$fast$adj_null$bin_nulls)),
                     sort(unlist(pp$ref$adj_null$bin_nulls)), tolerance = 1e-12)
        expect_equal(pp$fast$p_values, pp$ref$p_values, tolerance = 1e-12)
        expect_equal(pp$fast$adj_null$bin_nulls, pp$ref$adj_null$bin_nulls,
                     tolerance = 1e-12)
        expect_identical(pp$fast$n_perm_used, pp$ref$n_perm_used)
      })
    })
  }
}

test_that("models without an exact engine keep the per-ROI permutation path", {
  set.seed(1302)
  ds <- gen_sample_dataset(c(4, 4, 4), 40, nlevels = 3, blocks = 4)
  ms <- mvpa_model(load_model("dual_lda"), ds$dataset, ds$design, "classification",
                   crossval = blocked_cross_validation(ds$design$block_var))
  expect_null(rMVPA:::.permutation_engine(ms, 2, get_center_ids(ds$dataset), "Accuracy", list()))
})

test_that("prepared permutations require fixed folds and leave RNG untouched when declined", {
  set.seed(1303)
  ds <- gen_sample_dataset(c(3, 3, 3), 48, nlevels = 3, blocks = 4)
  block <- ds$design$block_var
  fixed <- list(
    blocked = blocked_cross_validation(block),
    custom = custom_cross_validation(lapply(unique(block), function(b) {
      list(train = which(block != b), test = which(block == b))
    }))
  )
  random <- list(
    kfold = kfold_cross_validation(48, 4),
    twofold = twofold_blocked_cross_validation(block, nreps = 2),
    bootstrap = bootstrap_blocked_cross_validation(block, nreps = 2),
    sequential = sequential_blocked_cross_validation(block, nfolds = 2, nreps = 2)
  )
  for (model in c("corclass", "naive_bayes", "sda_notune")) {
    for (cv in fixed) {
      ms <- mvpa_model(load_model(model), ds$dataset, ds$design, "classification",
                       crossval = cv)
      expect_type(rMVPA:::.permutation_engine(ms, 2, get_center_ids(ds$dataset),
                                             "Accuracy", list()), "list")
    }
    for (cv in random) {
      ms <- mvpa_model(load_model(model), ds$dataset, ds$design, "classification",
                       crossval = cv)
      before <- .Random.seed
      expect_null(rMVPA:::.permutation_engine(ms, 2, get_center_ids(ds$dataset),
                                             "Accuracy", list()))
      expect_identical(.Random.seed, before)
    }
  }
})

test_that("randomized CV retains the per-permutation path and RNG for both strategies", {
  set.seed(1304)
  ds <- gen_sample_dataset(c(3, 3, 3), 48, nlevels = 3, blocks = 4)
  block <- ds$design$block_var
  cvs <- list(
    kfold = kfold_cross_validation(48, 4),
    twofold = twofold_blocked_cross_validation(block, nreps = 2),
    bootstrap = bootstrap_blocked_cross_validation(block, nreps = 2),
    sequential = sequential_blocked_cross_validation(block, nfolds = 2, nreps = 2)
  )
  for (cv in cvs) {
    ms <- mvpa_model(load_model("corclass"), ds$dataset, ds$design, "classification",
                     crossval = cv)
    for (strategy in c("iterate", "searchlight")) {
      pc <- permutation_control(n_perm = 3, seed = 7, perm_strategy = strategy,
                                subsample = 0.5, diagnose = FALSE)
      run <- function(...) suppressWarnings(run_permutation_searchlight(
        ms, observed = rep(0.5, 27), radius = 2, perm_ctrl = pc,
        metric = "Accuracy", ...
      ))
      set.seed(913)
      auto <- run()
      auto_seed <- .Random.seed
      set.seed(913)
      # Explicit engines bypass preparation reuse. Searchlight keeps its
      # existing engine dispatch; iterate keeps its per-ROI iterator.
      ref <- run(engine = if (strategy == "iterate") "legacy" else "aggregate_fast")
      expect_identical(auto$adj_null$bin_nulls, ref$adj_null$bin_nulls)
      expect_identical(auto$p_values, ref$p_values)
      expect_identical(.Random.seed, auto_seed)
    }
  }
})

test_that("null collection retains valid draws around skipped draws", {
  set.seed(1305)
  ds <- gen_sample_dataset(c(3, 3, 3), 20, nlevels = 2, blocks = 2)
  cv <- custom_cross_validation(list(list(train = 1:18, test = 19:20)))
  ms <- mvpa_model(load_model("corclass"), ds$dataset, ds$design, "classification",
                   crossval = cv)
  # Some global draws give the two assessment rows the same class: AUC is
  # missing everywhere, so these draws leave gaps in the collection list.
  valid <- vapply(1:10, function(i) {
    y <- permute_labels(ds$design, "global", seed = 30 + i)$y_train
    length(unique(y[19:20])) == 2L
  }, logical(1))
  expect_true(any(valid) && any(!valid))
  for (strategy in c("iterate", "searchlight")) {
    run <- function(seed, n_perm) suppressWarnings(run_permutation_searchlight(
      ms, observed = rep(0, 27), radius = 2, metric = "AUC",
      perm_ctrl = permutation_control(n_perm = n_perm, seed = seed,
        shuffle = "global", perm_strategy = strategy, subsample = 1,
        null_method = "global", diagnose = FALSE)
    ))
    combined <- run(30, 10)
    separate <- lapply(which(valid), function(i) run(30 + i - 1L, 1))
    expected <- sort(unlist(lapply(separate, function(x) x$adj_null$bin_nulls[[1L]])))
    expect_identical(combined$adj_null$bin_nulls[[1L]], expected)
    expect_equal(combined$n_null_vals, 27 * sum(valid))
    expect_error(run(30 + which(!valid)[1L] - 1L, 1), "No valid null values collected")
  }
})

test_that("RSA permutations through rsa_fast equal the per-ROI path", {
  set.seed(1303)
  ds <- gen_sample_dataset(c(4, 4, 4), 20, blocks = 4)
  D1 <- dist(matrix(rnorm(20 * 4), 20)); D2 <- dist(matrix(rnorm(20 * 4), 20))
  rdes <- rsa_design(~ D1 + D2, list(D1 = D1, D2 = D2, block = ds$design$block_var),
                     block_var = "block")
  for (measure in c("distance", "similarity"))
    for (rt in c("pearson", "lm")) for (st in c("iterate", "searchlight")) {
    ms <- rsa_model(ds$dataset, rdes, distmethod = "pearson", regtype = rt,
                     measure = measure, check_collinearity = FALSE)
    pc <- permutation_control(n_perm = 3, seed = 5, perm_strategy = st, subsample = 0.5,
                              diagnose = FALSE, rsa_null = if (rt == "lm") "joint" else "individual")
    obs <- suppressWarnings(run_searchlight(ms, radius = 2, engine = "legacy", backend = "default"))
    fast <- suppressWarnings(run_permutation_searchlight(ms, observed = obs, radius = 2,
                                                         perm_ctrl = pc, metric = "D1"))
    ref <- suppressWarnings(run_permutation_searchlight(ms, observed = obs, radius = 2,
                                                        perm_ctrl = pc, metric = "D1", engine = "legacy"))
    # rsa_fast reindexes each sphere's cached RDM per permutation; values
    # equal recomputation up to BLAS summation order.
    expect_equal(sort(unlist(fast$adj_null$bin_nulls)), sort(unlist(ref$adj_null$bin_nulls)),
                 tolerance = 1e-12, info = paste(measure, rt, st))
    expect_equal(fast$p_values, ref$p_values, info = paste(measure, rt, st))
    # Without the cache every permutation recomputes its RDMs: identical.
    old_opt <- options(rMVPA.rsa_perm_cache_bytes = 0)
    nocache <- suppressWarnings(run_permutation_searchlight(ms, observed = obs, radius = 2,
                                                            perm_ctrl = pc, metric = "D1"))
    options(old_opt)
    expect_identical(sort(unlist(nocache$adj_null$bin_nulls)), sort(unlist(ref$adj_null$bin_nulls)),
                     info = paste(measure, rt, st))
    expect_identical(nocache$p_values, ref$p_values, info = paste(measure, rt, st))
  }
})
