# Golden characterisation scenarios for the core method set.
#
# Each scenario is a deterministic function returning a compact list of plain
# R values: numeric vectors, character vectors and data frames, never S4 or
# NeuroVol objects. data-raw/golden/make_golden.R evaluates them once and saves
# the result to tests/testthat/fixtures/golden/; test_golden_characterisation.R
# re-evaluates them and compares. A fixture records current behaviour, not
# correctness. When a deliberate fix changes a scenario, regenerate the
# fixture in the same commit and explain the change in NEWS.md.

golden_fixture_dir <- function() {
  testthat::test_path("fixtures", "golden")
}

golden_quiet <- function(code) {
  old <- futile.logger::flog.threshold()
  futile.logger::flog.threshold(futile.logger::ERROR)
  on.exit(futile.logger::flog.threshold(old), add = TRUE)
  suppressMessages(suppressWarnings(force(code)))
}

golden_find_extdata <- function(...) {
  installed <- system.file("extdata", ..., package = "rMVPA")
  if (nzchar(installed)) {
    return(installed)
  }
  candidates <- c(
    file.path("inst", "extdata", ...),
    file.path("..", "..", "inst", "extdata", ...)
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) {
    stop("Cannot locate rMVPA extdata: ", file.path(...), call. = FALSE)
  }
  hit[[1]]
}

# Haxby (2001) subject 1, VT mask: 96 block-mean patterns
# (8 categories x 12 runs) over 577 voxels.
golden_haxby <- function() {
  bundle <- readRDS(golden_find_extdata("haxby2001_subj1", "patterns.rds"))
  mask_arr <- array(0L, bundle$mask_dim)
  mask_arr[bundle$mask_idx] <- 1L
  mask <- neuroim2::LogicalNeuroVol(mask_arr, bundle$mask_space)

  vec_data <- t(bundle$patterns)
  storage.mode(vec_data) <- "double"
  vec_space <- neuroim2::add_dim(bundle$mask_space, ncol(vec_data))
  bold_vec <- neuroim2::SparseNeuroVec(vec_data, space = vec_space, mask = as.logical(mask_arr))
  dataset <- mvpa_dataset(bold_vec, mask = mask)

  # Four deterministic ROIs: contiguous quarters of the VT mask in index order.
  roi_arr <- array(0L, bundle$mask_dim)
  roi_arr[bundle$mask_idx] <- as.integer(cut(seq_along(bundle$mask_idx), 4L, labels = FALSE))
  rois <- neuroim2::NeuroVol(roi_arr, bundle$mask_space)

  design_df <- data.frame(category = bundle$category, run = bundle$run)
  design <- mvpa_design(design_df, y_train = ~ category, block_var = ~ run)

  list(bundle = bundle, dataset = dataset, design = design, rois = rois,
       crossval = blocked_cross_validation(bundle$run))
}

golden_synthetic <- function() {
  set.seed(20260930)
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 48, nlevels = 3, blocks = 4)
  list(dataset = ds$dataset, design = ds$design,
       crossval = blocked_cross_validation(ds$design$block_var))
}

golden_searchlight_maps <- function(res) {
  lapply(res$results, function(m) as.numeric(neuroim2::values(m)))
}

golden_regional_tables <- function(res) {
  perf <- as.data.frame(res$performance_table)
  pred <- as.data.frame(res$prediction_table)
  keep <- c("roinum", "observed", "predicted", grep("^prob_", names(pred), value = TRUE))
  pred <- pred[, intersect(keep, names(pred)), drop = FALSE]
  pred$observed <- as.character(pred$observed)
  pred$predicted <- as.character(pred$predicted)
  pred <- pred[order(pred$roinum, seq_len(nrow(pred))), , drop = FALSE]
  rownames(pred) <- NULL
  list(performance = perf, predictions = pred)
}

golden_searchlight_scenario <- function(model_name, engine) {
  force(model_name)
  force(engine)
  function() {
    s <- golden_synthetic()
    mspec <- mvpa_model(load_model(model_name), s$dataset, s$design, "classification",
                        crossval = s$crossval)
    set.seed(11)
    res <- run_searchlight(mspec, radius = 2, method = "standard",
                           engine = engine, backend = "default")
    list(engine = attr(res, "searchlight_engine"), maps = golden_searchlight_maps(res))
  }
}

golden_regional_scenario <- function(model_name, model = NULL) {
  force(model_name)
  force(model)
  function() {
    h <- golden_haxby()
    spec <- if (is.null(model)) load_model(model_name) else model
    mspec <- mvpa_model(spec, h$dataset, h$design, "classification",
                        crossval = h$crossval, return_predictions = TRUE)
    set.seed(12)
    golden_regional_tables(run_regional(mspec, h$rois, verbose = FALSE))
  }
}

# The original SDA estimator is an independent reference for the native fit.
# Pass a local model specification through the public runner; do not alter the
# model registry or namespace used by the implementation under test.
golden_external_sda <- function() {
  model <- load_model("sda_notune")
  model$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    fit <- sda::sda(as.matrix(x), y, verbose = FALSE, ...)
    fit$obsLevels <- lev
    fit
  }
  model$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    predict(modelFit, as.matrix(newdata), verbose = FALSE)$class
  }
  model$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    predict(modelFit, as.matrix(newdata), verbose = FALSE)$posterior
  }
  golden_quiet(golden_regional_scenario("sda_notune", model)())
}

# Independent Mann-Whitney pair counting: exact score ties receive half credit.
# Use the probabilities as returned, without rounding or tolerance-based ties.
golden_pairwise_auc <- function(tables) {
  vapply(tables$performance$roinum, function(roi) {
    pred <- tables$predictions[tables$predictions$roinum == roi, , drop = FALSE]
    columns <- grep("^prob_", names(pred), value = TRUE)
    auc <- vapply(columns, function(column) {
      positive <- pred$observed == sub("^prob_", "", column)
      difference <- outer(pred[[column]][positive], pred[[column]][!positive], "-")
      mean((difference > 0) + 0.5 * (difference == 0))
    }, numeric(1))
    mean(2 * auc - 1)
  }, numeric(1))
}

golden_rsa_scenario <- function(distmethod, regtype) {
  force(distmethod)
  force(regtype)
  function() {
    h <- golden_haxby()
    cat_ <- h$bundle$category
    same_category <- stats::as.dist(1 - outer(cat_, cat_, "=="))
    run_distance <- stats::dist(h$bundle$run)
    rdes <- rsa_design(~ same_category + run_distance,
                       list(same_category = same_category, run_distance = run_distance,
                            run = h$bundle$run),
                       block_var = "run")
    mspec <- rsa_model(h$dataset, rdes, distmethod = distmethod, regtype = regtype,
                       check_collinearity = FALSE)
    set.seed(13)
    res <- run_regional(mspec, h$rois, verbose = FALSE)
    list(performance = as.data.frame(res$performance_table))
  }
}

golden_vector_rsa_scenario <- function() {
  h <- golden_haxby()
  cats <- levels(h$bundle$category)
  # Category-level model RDM from the category order (an arbitrary but fixed
  # geometry); vector RSA compares it to each trial's neural similarity profile.
  D <- as.matrix(stats::dist(seq_along(cats)))
  dimnames(D) <- list(cats, cats)
  vdes <- vector_rsa_design(D = D, labels = as.character(h$bundle$category),
                            block_var = h$bundle$run)
  mspec <- vector_rsa_model(h$dataset, vdes, distfun = cordist(), rsa_simfun = "pearson")
  set.seed(14)
  # No verbose argument: run_regional.vector_rsa_model passes it twice to
  # mvpa_iterate (fixed separately).
  res <- run_regional(mspec, h$rois)
  list(performance = as.data.frame(res$performance_table))
}

golden_crossnobis_scenario <- function() {
  h <- golden_haxby()
  b <- h$bundle
  cats <- levels(b$category)
  runs <- sort(unique(b$run))
  U <- array(NA_real_, dim = c(length(cats), ncol(b$patterns), length(runs)),
             dimnames = list(cats, NULL, NULL))
  for (m in seq_along(runs)) {
    for (k in seq_along(cats)) {
      rows <- which(b$run == runs[m] & b$category == cats[k])
      U[k, , m] <- colMeans(b$patterns[rows, , drop = FALSE])
    }
  }
  d <- rMVPA:::compute_crossnobis_distances_sl(U)
  list(distances = as.numeric(d), names = names(d))
}

golden_metric_scenario <- function() {
  # Canonical ASCII name order, independent of R CMD check's C collation.
  ordered_metrics <- function(x) {
    unclass(x)[order(tolower(names(x)), method = "radix")]
  }
  set.seed(15)
  obs <- seq(-2, 2, length.out = 20)
  regression <- list(
    exact = performance(regression_result(obs, obs)),
    offset = performance(regression_result(obs, obs + 3)),
    rescaled = performance(regression_result(obs, 2 * obs)),
    noisy = performance(regression_result(obs, obs + stats::rnorm(20, sd = 0.5)))
  )

  lv <- c("a", "b", "c")
  observed <- factor(rep(lv, length.out = 12), levels = lv)
  probs <- matrix(c(
    0.6, 0.2, 0.2,  0.1, 0.8, 0.1,  0.3, 0.3, 0.4,  0.5, 0.5, 0.0,
    0.2, 0.2, 0.6,  0.1, 0.1, 0.8,  0.4, 0.4, 0.2,  0.3, 0.3, 0.4,
    0.7, 0.2, 0.1,  1/3, 1/3, 1/3,  0.2, 0.6, 0.2,  0.1, 0.45, 0.45
  ), ncol = 3, byrow = TRUE, dimnames = list(NULL, lv))
  predicted <- factor(lv[max.col(probs, ties.method = "first")], levels = lv)
  multiclass <- performance(classification_result(observed, predicted, probs))

  bin_obs <- factor(rep(c("a", "b"), 6), levels = c("a", "b"))
  bin_probs <- cbind(a = c(0.9, 0.4, 0.5, 0.5, 0.2, 0.3, 0.8, 0.6, 0.5, 0.1, 0.7, 0.45),
                     b = NA_real_)
  bin_probs[, "b"] <- 1 - bin_probs[, "a"]
  bin_pred <- factor(ifelse(bin_probs[, "a"] >= 0.5, "a", "b"), levels = c("a", "b"))
  binary <- performance(classification_result(bin_obs, bin_pred, bin_probs))

  list(
    regression = lapply(regression, ordered_metrics),
    multiclass = ordered_metrics(multiclass),
    binary = ordered_metrics(binary)
  )
}

golden_scenarios <- function() {
  list(
    sl_corclass_legacy = golden_searchlight_scenario("corclass", "legacy"),
    sl_sda_notune_legacy = golden_searchlight_scenario("sda_notune", "legacy"),
    sl_dual_lda_legacy = golden_searchlight_scenario("dual_lda", "legacy"),
    sl_dual_lda_fast = golden_searchlight_scenario("dual_lda", "dual_lda_fast"),
    regional_corclass = golden_regional_scenario("corclass"),
    regional_sda_notune = golden_regional_scenario("sda_notune"),
    regional_naive_bayes = golden_regional_scenario("naive_bayes"),
    regional_dual_lda = golden_regional_scenario("dual_lda"),
    rsa_spearman_pearson = golden_rsa_scenario("spearman", "pearson"),
    rsa_pearson_lm = golden_rsa_scenario("pearson", "lm"),
    vector_rsa = golden_vector_rsa_scenario,
    crossnobis_haxby = golden_crossnobis_scenario,
    metrics = golden_metric_scenario
  )
}

golden_run_scenario <- function(name) {
  golden_quiet(golden_scenarios()[[name]]())
}
