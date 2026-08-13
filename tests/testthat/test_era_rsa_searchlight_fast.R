context("ERA-RSA dedicated searchlight engine")

.make_era_fast_fixture <- function(K = 12L, dims = c(4L, 4L, 4L),
                                   seed = 731L,
                                   retrieval_na = FALSE,
                                   training_na = FALSE) {
  set.seed(seed)
  keys <- sprintf("item_%02d", seq_len(K))
  p <- prod(dims)
  E <- matrix(stats::rnorm(K * p), nrow = K)
  R <- E + matrix(stats::rnorm(K * p, sd = 0.35), nrow = K)

  # Exercise the filter_roi() center-preservation contract: feature 1 is
  # constant during encoding but informative during retrieval.
  E[, 1L] <- 0
  R[, 1L] <- seq_len(K)
  if (isTRUE(training_na)) {
    E[3L, 2L] <- NA_real_
  }
  if (isTRUE(retrieval_na)) {
    # Training remains finite, so this feature survives the train-only filter;
    # fit_roi() must then use its pairwise-complete similarity fallback.
    R[3L, 2L] <- NA_real_
  }

  enc_vec <- neuroim2::NeuroVec(
    array(as.numeric(t(E)), dim = c(dims, K)),
    neuroim2::NeuroSpace(c(dims, K))
  )
  ret_vec <- neuroim2::NeuroVec(
    array(as.numeric(t(R)), dim = c(dims, K)),
    neuroim2::NeuroSpace(c(dims, K))
  )
  mask <- neuroim2::NeuroVol(
    array(1, dim = dims),
    neuroim2::NeuroSpace(dims)
  )

  vividness <- sample(seq_len(K))
  vividness[c(4L, K - 1L)] <- NA_real_
  train_tab <- data.frame(
    item = factor(keys, levels = keys),
    vividness = NA_real_,
    retrieval_run = factor(rep(c("run_1", "run_2"), length.out = K)),
    trial_order = seq_len(K)
  )
  test_tab <- train_tab
  test_tab$vividness <- vividness

  dataset <- mvpa_dataset(enc_vec, test_data = ret_vec, mask = mask)
  design <- mvpa_design(
    train_tab,
    test_design = test_tab,
    y_train = ~ item,
    y_test = ~ item
  )
  list(dataset = dataset, design = design)
}

.make_era_fast_model <- function(fixture, simfun = "pearson",
                                 components = "item") {
  suppressWarnings(era_rsa_model(
    dataset = fixture$dataset,
    design = fixture$design,
    key_var = ~ item,
    pairing = "one_to_one",
    era_simfun = simfun,
    era_components = components,
    era_association = ~ vividness + retrieval_run + trial_order,
    era_effects = ~ vividness,
    era_min_complete = 4L
  ))
}

.expect_searchlight_maps_equal <- function(actual, expected, tolerance = 1e-12) {
  expect_identical(actual$metrics, expected$metrics)
  expect_identical(names(actual$results), names(expected$results))
  for (nm in expected$metrics) {
    expect_equal(
      as.numeric(neuroim2::values(actual$results[[nm]])),
      as.numeric(neuroim2::values(expected$results[[nm]])),
      tolerance = tolerance,
      info = nm
    )
  }
}

test_that("auto selects the dedicated ERA-RSA engine only for eligible runs", {
  fx <- .make_era_fast_fixture()
  model <- .make_era_fast_model(fx)

  expect_true(rMVPA:::.is_era_rsa_fast_path(model, "standard"))
  expect_identical(
    rMVPA:::.resolve_searchlight_engine(model, "standard", "auto"),
    "era_rsa_fast"
  )
  expect_identical(
    rMVPA:::.resolve_searchlight_engine(model, "randomized", "auto"),
    "legacy"
  )

  repeated_item_model <- .make_era_fast_model(fx)
  repeated_item_model$pairing <- "average"
  expect_false(rMVPA:::.is_era_rsa_fast_path(repeated_item_model, "standard"))

  audit <- explain_searchlight_engine(model, method = "standard", engine = "auto")
  expect_identical(audit$engine[audit$selected], "era_rsa_fast")
  expect_true(audit$eligible[audit$engine == "era_rsa_fast"])

  expect_error(
    run_searchlight(
      model,
      radius = 1,
      method = "randomized",
      engine = "era_rsa_fast",
      backend = "default",
      preflight = "off"
    ),
    "not eligible"
  )
})

test_that("direct Pearson engine matches the legacy searchlight exactly", {
  fx <- .make_era_fast_fixture()
  model <- .make_era_fast_model(fx, simfun = "pearson")

  legacy <- run_searchlight(
    model, radius = 1, method = "standard",
    engine = "legacy", backend = "default", preflight = "off"
  )
  fast <- run_searchlight(
    model, radius = 1, method = "standard",
    engine = "era_rsa_fast", backend = "default", preflight = "off"
  )

  expect_identical(attr(legacy, "searchlight_engine"), "legacy")
  expect_identical(attr(fast, "searchlight_engine"), "era_rsa_fast")
  .expect_searchlight_maps_equal(fast, legacy)
  expect_true(any(is.finite(neuroim2::values(
    fast$results$era_assoc_part_r_vividness
  ))))
})

test_that("direct Spearman engine preserves retrieval missing-data semantics", {
  fx <- .make_era_fast_fixture(retrieval_na = TRUE)
  model <- .make_era_fast_model(fx, simfun = "spearman")

  legacy <- run_searchlight(
    model, radius = 1, method = "standard",
    engine = "legacy", backend = "default", preflight = "off"
  )
  fast <- run_searchlight(
    model, radius = 1, method = "standard",
    engine = "era_rsa_fast", backend = "default", preflight = "off"
  )

  .expect_searchlight_maps_equal(fast, legacy)
})

test_that("dedicated engine has exact differential parity across seeded missing-data fixtures", {
  cases <- data.frame(
    seed = c(101L, 211L, 307L, 419L),
    method = c("pearson", "spearman", "pearson", "spearman"),
    training_na = c(FALSE, TRUE, TRUE, FALSE),
    retrieval_na = c(TRUE, TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(cases))) {
    case <- cases[i, ]
    fx <- .make_era_fast_fixture(
      K = 9L,
      dims = c(3L, 3L, 3L),
      seed = case$seed,
      training_na = case$training_na,
      retrieval_na = case$retrieval_na
    )
    model <- .make_era_fast_model(fx, simfun = case$method)

    legacy <- run_searchlight(
      model, radius = 1, method = "standard",
      engine = "legacy", backend = "default", preflight = "off"
    )
    fast <- run_searchlight(
      model, radius = 1, method = "standard",
      engine = "era_rsa_fast", backend = "default", preflight = "off"
    )

    .expect_searchlight_maps_equal(fast, legacy)
  }
})

test_that("dedicated engine remains exact for identification and geometry", {
  fx <- .make_era_fast_fixture(K = 10L)
  model <- .make_era_fast_model(
    fx,
    simfun = "pearson",
    components = c("item", "identification", "geometry")
  )

  legacy <- run_searchlight(
    model, radius = 1, method = "standard",
    engine = "legacy", backend = "default", preflight = "off"
  )
  fast <- run_searchlight(
    model, radius = 1, method = "standard",
    engine = "era_rsa_fast", backend = "default", preflight = "off"
  )

  .expect_searchlight_maps_equal(fast, legacy)
})

test_that("shard-backed dedicated engine matches the legacy path", {
  skip_if_not_installed("shard")
  fx <- .make_era_fast_fixture(K = 10L, dims = c(3L, 3L, 3L))
  model <- .make_era_fast_model(fx, simfun = "pearson")

  legacy <- run_searchlight(
    model, radius = 1, method = "standard",
    engine = "legacy", backend = "default", preflight = "off"
  )
  fast <- run_searchlight(
    model, radius = 1, method = "standard",
    engine = "era_rsa_fast", backend = "shard", preflight = "off"
  )

  .expect_searchlight_maps_equal(fast, legacy)
})

test_that("parallel shard execution is numerically identical", {
  skip_if_not_installed("shard")
  skip_on_cran()
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  old_options <- options(parallelly.maxWorkers.localhost = 2)
  on.exit(options(old_options), add = TRUE)

  fx <- .make_era_fast_fixture(K = 10L, dims = c(3L, 3L, 3L))
  model <- .make_era_fast_model(fx, simfun = "spearman")
  reference <- run_searchlight(
    model, radius = 1, method = "standard",
    engine = "era_rsa_fast", backend = "default", preflight = "off"
  )

  suppressWarnings({
    future::plan(future::multisession, workers = 2L)
    parent_engine_signature <- paste(deparse(body(get(
      ".era_fast_fit_task",
      envir = asNamespace("rMVPA"),
      inherits = FALSE
    ))), collapse = "\n")
    worker_engine_signature <- future::value(future::future({
      ns <- asNamespace("rMVPA")
      if (!exists(".era_fast_fit_task", envir = ns, inherits = FALSE)) {
        return(NULL)
      }
      paste(deparse(body(get(
        ".era_fast_fit_task",
        envir = ns,
        inherits = FALSE
      ))), collapse = "\n")
    }))
    if (!identical(worker_engine_signature, parent_engine_signature)) {
      skip(paste0(
        "multisession parity requires workers to load the current installed ",
        "rMVPA package; source-loaded tests cannot replace a stale worker namespace"
      ))
    }
    parallel <- run_searchlight(
      model, radius = 1, method = "standard",
      engine = "era_rsa_fast", backend = "shard", preflight = "off"
    )
  })

  .expect_searchlight_maps_equal(parallel, reference)
})
