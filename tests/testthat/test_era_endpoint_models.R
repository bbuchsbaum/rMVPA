context("Endpoint model and spatial-output integration")

.endpoint_model_fixture <- function(seed = 908L) {
  set.seed(seed)
  n <- 12L; dims <- c(2L, 2L, 2L); p <- prod(dims)
  keys <- sprintf("item%02d", seq_len(n))
  E <- matrix(rnorm(n * p), n)
  R <- E + matrix(rnorm(n * p, sd = 0.4), n)
  space <- neuroim2::NeuroSpace(dims)
  mask <- neuroim2::NeuroVol(array(1, dims), space)
  vec <- function(X) neuroim2::NeuroVec(array(as.numeric(t(X)), c(dims, n)),
                                       neuroim2::NeuroSpace(c(dims, n)))
  tab <- data.frame(item = keys, vividness = rnorm(n))
  list(E = E, R = R, keys = keys, block = rep(1:3, each = 4), tab = tab,
    mask = mask, dataset = mvpa_dataset(vec(E), test_data = vec(R), mask = mask),
    feature_dataset = mvpa_dataset(vec(R), mask = mask),
    design = mvpa_design(tab, test_design = tab, y_train = ~ item, y_test = ~ item),
    F = matrix(rnorm(n * 3), n, 3))
}

.endpoint_era_model <- function(fx, endpoints = NULL, score = "item_specificity") {
  era_rsa_model(fx$dataset, fx$design, key_var = ~ item, pairing = "one_to_one",
    era_components = "item", item_block = setNames(fx$block, fx$keys),
    era_correlates = ~ vividness, era_cor_method = "pearson",
    era_association_score = score, endpoints = endpoints)
}

.endpoint_feature_model <- function(fx, endpoints = NULL, return_predictions = FALSE) {
  design <- feature_rsa_design(F = fx$F, labels = fx$keys, max_comps = 3)
  feature_rsa_model(fx$feature_dataset, design, method = "ridge",
    crossval = blocked_cross_validation(fx$block), lambda = 0.1,
    lambda_selection = "fixed", endpoints = endpoints, return_predictions = return_predictions)
}

test_that("ERA opt-in preserves legacy scores and exposes joint excess", {
  fx <- .endpoint_model_fixture()
  spec <- era_endpoint_spec(fx$keys, fx$block, behavior = fx$tab$vividness,
                            return_items = TRUE)
  off <- .endpoint_era_model(fx)
  on <- .endpoint_era_model(fx, spec)
  joint <- .endpoint_era_model(fx, spec, "joint_excess")
  roi <- list(train_data = fx$E, test_data = fx$R, indices = 1:8)
  a <- fit_roi(off, roi, list(id = 1L))
  b <- fit_roi(on, roi, list(id = 1L))
  c <- fit_roi(joint, roi, list(id = 1L))
  expect_equal(b$metrics[names(a$metrics)], a$metrics)
  expect_setequal(names(b$metrics), names(output_schema(on)))
  raw <- era_endpoint_decomposition(fx$E, fx$R, item_key = fx$keys, competitor_block = fx$block)
  expect_equal(b$result$predictor$endpoints$decomposition$items, raw$items)
  expect_equal(unname(c$metrics["era_vividness_cor"]), cor(raw$items$excess, fx$tab$vividness))
  expect_true(on$return_fits)
  expect_null(a$result)
  expect_error(.endpoint_era_model(fx, score = "joint_excess"), "requires endpoints")
  expect_error(era_rsa_model(fx$dataset, fx$design, key_var = ~ item,
    endpoints = spec, pairing = "average"), "one_to_one")
  expect_error(era_rsa_model(fx$dataset, fx$design, key_var = ~ item,
    endpoints = spec, pairing = "one_to_one", era_simfun = "spearman"), "pearson")
  expect_error(.endpoint_era_model(fx, list()), "era_endpoint_spec")
})

test_that("ERA regional and both searchlight engines agree and export structural NA", {
  fx <- .endpoint_model_fixture()
  spec <- era_endpoint_spec(fx$keys, fx$block, behavior = fx$tab$vividness, return_items = TRUE)
  model <- .endpoint_era_model(fx, spec)
  regional <- run_regional(model, fx$mask, preflight = "off", backend = "default")
  expect_equal(length(regional$fits), 1L)
  expect_s3_class(regional$fits[[1]]$endpoints$decomposition, "era_endpoints")
  scalar_names <- grep("^era_endpoint_", names(output_schema(model)), value = TRUE)
  for (engine in c("legacy", "era_rsa_fast")) {
    sl <- run_searchlight(model, radius = 10, engine = engine,
                         backend = "default", preflight = "off")
    for (nm in scalar_names) {
      actual <- as.numeric(neuroim2::values(sl$results[[nm]]))
      expect_equal(actual, rep(regional$performance_table[[nm]][1], 8L), tolerance = 1e-10, info = nm)
    }
  }
  expect_true(all(is.na(neuroim2::values(sl$results$era_endpoint_cor_shared))))
  directory <- tempfile("endpoint-export-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  save_results(sl, directory, quiet = TRUE)
  files <- list.files(directory, pattern = "\\.nii(\\.gz)?$", full.names = TRUE, recursive = TRUE)
  mean_map <- files[grepl("era_endpoint_mean_excess", files)]
  na_map <- files[grepl("era_endpoint_cor_shared", files)]
  count_map <- files[grepl("era_endpoint_assoc_n\\.", files)]
  expect_length(mean_map, 1L)
  expect_length(na_map, 1L)
  expect_length(count_map, 1L)
  expect_equal(as.numeric(neuroim2::values(neuroim2::read_vol(mean_map))),
               rep(regional$performance_table$era_endpoint_mean_excess, 8), tolerance = 1e-6)
  expect_true(all(is.na(neuroim2::values(neuroim2::read_vol(na_map)))))
  expect_equal(as.numeric(neuroim2::values(neuroim2::read_vol(count_map))), rep(12, 8))
})

test_that("feature endpoints preserve OOF predictions, folds, tuning and legacy metrics", {
  fx <- .endpoint_model_fixture()
  spec <- era_endpoint_spec(fx$keys, fx$block, behavior = fx$tab$vividness, return_items = TRUE)
  off <- .endpoint_feature_model(fx, return_predictions = TRUE)
  on <- .endpoint_feature_model(fx, spec, return_predictions = TRUE)
  a <- run_regional(off, fx$mask, preflight = "off", backend = "default")
  b <- run_regional(on, fx$mask, preflight = "off", backend = "default")
  expect_equal(b$performance_table[, names(a$performance_table)], a$performance_table)
  ap <- feature_rsa_predictions(a)
  bp <- feature_rsa_predictions(b)
  expect_identical(ap$predicted, bp$predicted)
  expect_identical(ap$observed, bp$observed)
  expect_identical(ap$fold_id, bp$fold_id)
  expect_identical(ap$observation_index, bp$observation_index)
  endpoints <- b$fits[[1]]$endpoints
  expect_true(all(endpoints$decomposition$items$template_source == "held_out_prediction"))
  expect_identical(endpoints$decomposition$items$fold_id, bp$fold_id[[1]])
  ref <- era_endpoint_decomposition(bp$predicted[[1]], bp$observed[[1]],
    item_key = fx$keys[bp$observation_index[[1]]], competitor_block = fx$block,
    template_source = "held_out_prediction")
  expect_equal(endpoints$decomposition$items[, names(ref$items)], ref$items)
  expect_lt(endpoints$association$identity_error, 1e-12)
  # Endpoint scoring itself does not require retaining prediction matrices.
  scalars <- .endpoint_feature_model(fx, spec)
  reg <- run_regional(scalars, fx$mask, preflight = "off", backend = "default")
  expect_equal(reg$performance_table, b$performance_table)
  sl <- run_searchlight(scalars, radius = 10, backend = "default", preflight = "off")
  for (nm in grep("^feature_endpoint_", names(output_schema(scalars)), value = TRUE)) {
    expect_equal(as.numeric(neuroim2::values(sl$results[[nm]])),
      rep(reg$performance_table[[nm]][1], 8L), tolerance = 1e-10, info = nm)
  }
})

test_that("feature endpoint scoring refuses cross-fold pools and repeated assessments", {
  fx <- .endpoint_model_fixture()
  spec <- era_endpoint_spec(fx$keys, rep("pooled", 12))
  model <- .endpoint_feature_model(fx, spec)
  result <- tibble::tibble(observed = list(fx$R[1:6, ], fx$R[7:12, ]),
    predicted = list(fx$E[1:6, ], fx$E[7:12, ]), test_index = list(1:6, 7:12),
    result = list(list(lambda = 0.1), list(lambda = 0.1)),
    error = FALSE, error_message = "~")
  expect_error(merge_results(model, result, indices = 1:8, id = 1L), "nested within outer")
  result$test_index <- list(1:6, 1:6)
  expect_error(merge_results(model, result, indices = 1:8, id = 1L), "exactly one out-of-fold")
  bad <- spec; bad$item_key[1] <- "unmatched"
  expect_error(.endpoint_feature_model(fx, bad), "match model item keys")
})
