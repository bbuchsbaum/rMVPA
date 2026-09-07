library(testthat)

# helper-test-logging.R already pins the log threshold; this only hides the
# iteration banners that the runners print to stdout.
quiet_run <- function(expr) {
  utils::capture.output(res <- suppressMessages(expr))
  res
}

test_that("pattern_model constructor validates inputs and fixes the penalty API", {
  sim <- sim_pattern_data(n = 45, dims = c(5, 5, 3), K = 3, seed = 11)
  spec <- pattern_model(sim$dataset, sim$design, max_rank = 2)
  expect_s3_class(spec, "pattern_model")
  expect_s3_class(spec, "model_spec")
  expect_equal(spec$target_type, "categorical")
  expect_s3_class(spec$crossval, "blocked_cross_validation")
  expect_true(spec$compute_performance)
  expect_false(spec$return_fits)
  expect_output(print(spec), "pattern_model specification")

  expect_error(pattern_model(sim$dataset, sim$design, penalty = list(sparse = 0.1)), "not\\s+implemented")
  expect_error(pattern_model(sim$dataset, sim$design, penalty = list(support_smooth = "auto")), "not\\s+implemented")
  expect_error(pattern_model(sim$dataset, sim$design, penalty = list(banana = 1)), "unknown penalty")
  expect_silent(pattern_model(sim$dataset, sim$design, penalty = list(sparse = 0, signed_smooth = NULL)))
  expect_error(pattern_model(sim$dataset, sim$design, rank = 0), "positive integer")
  expect_error(pattern_model(sim$dataset, sim$design, control = list()), "pattern_control")
  expect_error(pattern_model("nope", sim$design))

  # retention flags: predictions imply result retention; fits imply both
  s1 <- pattern_model(sim$dataset, sim$design, return_predictions = TRUE)
  expect_true(s1$return_predictions); expect_true(s1$return_fits); expect_false(s1$keep_fold_fits)
  s2 <- pattern_model(sim$dataset, sim$design, return_fits = TRUE)
  expect_true(s2$return_fits); expect_true(s2$keep_fold_fits); expect_false(s2$return_predictions)
})

test_that("pattern_model passes the plugin contract and schema validation", {
  sim <- sim_pattern_data(n = 45, dims = c(5, 5, 3), K = 3, seed = 12)
  spec <- pattern_model(sim$dataset, sim$design, rank = 2)
  rd <- mock_roi_data(train_data = sim$X[, 1:30], indices = 1:30)
  v <- validate_plugin_model(spec, roi_data = rd)
  expect_true(v$valid)
  expect_equal(v$metric_names, c("Accuracy", "AUC", "logloss", "rank_mean"))
  expect_true(validate_model_spec(spec, require_schema = TRUE, roi_data = rd)$valid)

  res <- fit_roi(spec, rd, mock_context(design = spec$design, cv_spec = spec$crossval, id = 7L))
  expect_s3_class(res, "roi_result")
  expect_false(res$error)
  expect_equal(res$id, 7L)
  expect_equal(res$indices, 1:30)
  expect_equal(names(res$metrics), c("Accuracy", "AUC", "logloss", "rank_mean"))
  expect_equal(unname(res$metrics["rank_mean"]), 2)
  expect_s3_class(res$result, "classification_result")
  expect_s3_class(res$result$predictor, "pattern_roi_fits")
  expect_equal(length(res$result$observed), 45L)
  expect_null(res$result$predictor$fold_fits)

  # continuous schema
  simc <- sim_pattern_data(n = 45, dims = c(5, 5, 3), q = 3, r = 2, seed = 13)
  specc <- pattern_model(simc$dataset, simc$design, rank = 2)
  expect_equal(names(output_schema(specc)), c("R2", "RMSE", "cor", "rank_mean"))
  resc <- fit_roi(specc, mock_roi_data(train_data = simc$X[, 1:30], indices = 1:30),
                  mock_context(design = specc$design, cv_spec = specc$crossval))
  expect_false(resc$error)
  expect_equal(names(resc$metrics), c("R2", "RMSE", "cor", "rank_mean"))
  expect_s3_class(resc$result$predictor$ledger, "pattern_ledger")

  # an unusable ROI yields an error result, not an exception
  bad <- mock_roi_data(train_data = matrix(1, 45, 3), indices = 1:3)
  rb <- fit_roi(spec, bad, mock_context(design = spec$design, cv_spec = spec$crossval))
  expect_true(rb$error)
  expect_match(rb$error_message, "no usable")
})

test_that("run_global recovers rank and patterns on known-truth data and retains fits", {
  # snr chosen so the second task dimension carries real, non-saturating
  # information: at high snr rank 1 already separates the classes perfectly and
  # the parsimony tie-break correctly returns the smaller rank (covered below).
  sim <- sim_pattern_data(n = 90, dims = c(8, 8, 4), K = 3, snr = 0.5, seed = 14)
  spec <- pattern_model(sim$dataset, sim$design, max_rank = 3)
  res <- quiet_run(run_global(spec, return_fits = TRUE, refit = TRUE))
  expect_s3_class(res, "pattern_global_result")
  expect_output(print(res), "pattern_global_result")
  pt <- performance(res)
  expect_equal(names(pt), c("Accuracy", "AUC", "logloss", "rank_mean"))
  expect_gt(pt$Accuracy, 0.9)
  expect_equal(unname(pt$rank_mean), 2)          # true rank recovered in every fold
  expect_equal(res$ranks, c(2L, 2L, 2L))
  expect_length(res$fold_fits, 3L)
  expect_s3_class(res$fold_fits[[1]], "pattern_fit")
  expect_s3_class(res$refit, "pattern_fit")
  expect_equal(res$refit$rank, 2L)
  expect_equal(res$feature_ids, which(sim$dataset$mask > 0))
  # pattern subspace recovery on the refit
  cs <- svd(crossprod(qr.Q(qr(sim$A)), qr.Q(qr(res$refit$A))))$d
  expect_true(all(cs > 0.9))
  # ledger covers every observation exactly once
  expect_equal(sort(res$ledger$observation), 1:90)
  expect_equal(res$ledger$partition, "cv")
  # default keeps no fits
  res0 <- quiet_run(run_global(spec))
  expect_null(res0$fold_fits); expect_null(res0$refit)
})

test_that("rank selection prefers the smaller rank when extra dimensions add nothing", {
  # A separable problem: every rank reaches zero held-out loss, so the
  # deliberate tie-break returns the most parsimonious rank.
  sim <- sim_pattern_data(n = 90, dims = c(8, 8, 4), K = 3, snr = 3, seed = 23)
  sel <- rMVPA:::.pattern_select_rank(sim$X, sim$targets, sim$design$block_var,
                                      pattern_control(max_rank = 3))
  expect_equal(sel$losses, c(0, 0))
  expect_equal(sel$rank, 1L)
  res <- quiet_run(run_global(pattern_model(sim$dataset, sim$design, max_rank = 3)))
  expect_equal(performance(res)$Accuracy, 1)
  expect_true(all(res$ranks == 1L))
})

test_that("run_global handles continuous targets, external test sets, and fixed ranks", {
  simc <- sim_pattern_data(n = 90, dims = c(8, 8, 4), q = 4, r = 2, snr = 2, seed = 15)
  specc <- pattern_model(simc$dataset, simc$design, max_rank = 4)
  resc <- quiet_run(run_global(specc))
  pt <- performance(resc)
  expect_equal(names(pt), c("R2", "RMSE", "cor", "rank_mean"))
  expect_gt(pt$R2, 0.2)
  expect_gt(pt$cor, 0.4)
  expect_equal(dim(resc$ledger$prediction), c(90L, 4L))
  expect_equal(dim(resc$ledger$truth), c(90L, 4L))

  sime <- sim_pattern_data(n = 60, dims = c(6, 6, 3), K = 2, snr = 1.5, external_test = TRUE, n_test = 40, seed = 16)
  spece <- pattern_model(sime$dataset, sime$design, rank = 5)   # capped at K - 1 = 1
  expect_true(spece$has_test_set)
  rese <- quiet_run(run_global(spece))
  expect_equal(rese$ledger$partition, "external")
  expect_equal(length(rese$ledger$observation), 40L)
  expect_equal(rese$ranks, 1L)
  expect_gt(performance(rese)$Accuracy, 0.9)

  # continuous external targets via targets_test
  simce <- sim_pattern_data(n = 60, dims = c(6, 6, 3), q = 2, r = 2, snr = 2, external_test = TRUE, n_test = 30, seed = 17)
  specce <- pattern_model(simce$dataset, simce$design, rank = 2)
  resce <- quiet_run(run_global(specce))
  expect_equal(dim(resce$ledger$truth), c(30L, 2L))
  expect_gt(performance(resce)$cor, 0.5)
})

test_that("run_global works with feature_sets_design targets and clustered datasets", {
  set.seed(18)
  ds <- gen_clustered_sample_dataset(c(8, 8, 8), nobs = 60, K = 6, nlevels = 2, blocks = 3)
  spec <- pattern_model(ds$dataset, ds$design, max_rank = 1)
  res <- quiet_run(run_global(spec))
  expect_equal(res$n_features, 6L)
  expect_equal(names(performance(res)), c("Accuracy", "AUC", "logloss", "rank_mean"))

  simc <- sim_pattern_data(n = 60, dims = c(6, 6, 3), q = 5, r = 2, snr = 2, seed = 19)
  fs <- feature_sets(simc$targets, blocks(low = 2, high = 3))
  fdes <- feature_sets_design(fs, block_var_train = simc$design$block_var)
  specf <- pattern_model(simc$dataset, fdes, max_rank = 3)
  expect_equal(specf$target_type, "continuous")
  # validate_analysis() reads design$block_var, which feature_sets_design does
  # not set, so it warns about blocking even though the CV is blocked.
  resf <- suppressWarnings(quiet_run(run_global(specf)))
  expect_gt(performance(resf)$cor, 0.3)
})

test_that("run_regional and run_searchlight dispatch pattern_model through fit_roi", {
  sim <- sim_pattern_data(n = 60, dims = c(6, 6, 3), K = 3, snr = 1.5, seed = 20)
  spec <- pattern_model(sim$dataset, sim$design, rank = 2, return_predictions = TRUE, return_fits = TRUE)
  region <- neuroim2::NeuroVol(array(rep(1:3, each = 36), sim$dims), neuroim2::space(sim$dataset$mask))

  reg <- quiet_run(run_regional(spec, region))
  expect_s3_class(reg, "regional_mvpa_result")
  expect_equal(nrow(reg$performance_table), 3L)
  expect_equal(names(reg$performance_table), c("roinum", "Accuracy", "AUC", "logloss", "rank_mean"))
  expect_true(all(reg$performance_table$rank_mean == 2))
  expect_true(!is.null(reg$prediction_table) && nrow(reg$prediction_table) == 180L)
  expect_length(reg$fits, 3L)
  expect_s3_class(reg$fits[[1]], "pattern_roi_fits")
  expect_length(reg$fits[[1]]$fold_fits, 3L)
  expect_equal(reg$fits[[1]]$fold_fits[[1]]$p_input, 36L)
  expect_true(all(c("Accuracy", "logloss") %in% names(reg$vol_results)))

  # continuous targets: no prediction table, ledger retained under $fits
  simc <- sim_pattern_data(n = 60, dims = c(6, 6, 3), q = 2, r = 1, snr = 2, seed = 21)
  specc <- pattern_model(simc$dataset, simc$design, rank = 1, return_predictions = TRUE)
  regc <- quiet_run(run_regional(specc, region))
  expect_equal(names(regc$performance_table), c("roinum", "R2", "RMSE", "cor", "rank_mean"))
  expect_s3_class(regc$fits[[1]]$ledger, "pattern_ledger")
  expect_null(regc$fits[[1]]$fold_fits)

  sl <- quiet_run(run_searchlight(pattern_model(sim$dataset, sim$design, rank = 1), radius = 2, method = "standard"))
  expect_equal(names(sl$results), c("Accuracy", "AUC", "logloss", "rank_mean"))
  acc <- neuroim2::values(sl$results$Accuracy)
  expect_true(all(is.finite(acc)))
  # informative voxels decode better than the rest
  expect_gt(mean(acc[sim$informative]), mean(acc[-sim$informative]))
})

test_that("blocked CV is taken from feature_sets_design block_var_train", {
  simc <- sim_pattern_data(n = 60, dims = c(5, 5, 3), q = 4, r = 2, seed = 30)
  fs <- feature_sets(simc$targets, blocks(a = 2, b = 2))
  fdes <- feature_sets_design(fs, block_var_train = simc$design$block_var)
  expect_null(fdes$block_var)                       # blocks live in block_var_train
  spec <- pattern_model(simc$dataset, fdes, rank = 1)
  expect_s3_class(spec$crossval, "blocked_cross_validation")
  expect_equal(spec$crossval$block_var, simc$design$block_var)
  expect_equal(spec$crossval$nfolds, 3L)

  # a design with no blocks at all warns rather than silently using random folds
  nb <- mvpa_design(data.frame(id = seq_len(60)), cv_labels = seq_len(60),
                    targets = simc$targets)
  expect_warning(pattern_model(simc$dataset, nb, rank = 1), "no block variable")
})

test_that("a fold missing a class is handled rather than crashing", {
  set.seed(31)
  sim <- sim_pattern_data(n = 60, dims = c(5, 5, 3), K = 3, snr = 1.5, seed = 31)
  # confine class "c" to block 3, so leave-one-block-out trains without it
  y <- as.character(sim$targets)
  blk <- sim$design$block_var
  y[y == "c"] <- rep(c("a", "b"), length.out = sum(y == "c"))
  y[blk == 3] <- "c"
  des <- mvpa_design(data.frame(y = factor(y), block = blk), y_train = ~ y, block_var = ~ block)
  spec <- pattern_model(sim$dataset, des, rank = 1)

  # this design trips several legitimate preflight warnings by construction
  res <- suppressWarnings(quiet_run(run_global(spec)))
  expect_s3_class(res, "pattern_global_result")
  pt <- performance(res)
  expect_true(all(is.finite(c(pt$Accuracy, pt$logloss))))
  # the ledger keeps every class column, including the one some folds never saw
  expect_equal(colnames(res$ledger$prediction), c("a", "b", "c"))
  expect_equal(levels(res$ledger$truth), c("a", "b", "c"))
  expect_equal(sort(res$ledger$observation), 1:60)
  # leave-one-block-out: fold 3 trains on blocks 1 and 2, which contain no "c",
  # and its test rows are exactly the "c" rows. That fold must give "c" zero
  # probability (padded column) and be scored as wrong rather than crash.
  f3 <- which(res$ledger$fold == 3)
  expect_true(all(as.character(res$ledger$truth[f3]) == "c"))
  expect_true(all(res$ledger$prediction[f3, "c"] == 0))
  expect_equal(unname(rowSums(res$ledger$prediction[f3, , drop = FALSE])), rep(1, length(f3)))
  expect_true(all(res$ledger$prediction[-f3, "c"] > 0))
  # rank selection survives inner folds that cannot score every class
  expect_true(all(res$ranks >= 1L))
})

test_that("repeated-CV predictions are pooled the way wrap_result pools them", {
  sim <- sim_pattern_data(n = 60, dims = c(5, 5, 3), K = 3, snr = 1.5, seed = 50)
  bv <- sim$design$block_var
  schemes <- list(
    blocked = blocked_cross_validation(bv),
    twofold = twofold_blocked_cross_validation(bv, nreps = 4),
    sequential = sequential_blocked_cross_validation(bv, nfolds = 2, nreps = 3),
    bootstrap = bootstrap_blocked_cross_validation(bv, nreps = 4)
  )
  for (nm in names(schemes)) {
    spec <- pattern_model(sim$dataset, sim$design, rank = 2, crossval = schemes[[nm]])
    res <- suppressWarnings(quiet_run(run_global(spec)))
    pooled <- res$ledger
    folded <- res$fold_ledger

    # one record per tested observation, in sorted order, like wrap_result()
    info <- paste("scheme", nm)
    expect_true(isTRUE(pooled$pooled), info = info)
    expect_equal(pooled$observation, sort(unique(folded$observation)), info = info)
    expect_false(anyDuplicated(pooled$observation) > 0, info = info)
    expect_equal(nrow(pooled$prediction), length(pooled$observation), info = info)
    expect_equal(unname(rowSums(pooled$prediction)), rep(1, length(pooled$observation)),
                 info = info)
    expect_equal(sum(pooled$n_repeats), length(folded$observation), info = info)

    # the pooled probability of each observation is the mean over its repeats
    i <- pooled$observation[1]
    reps <- which(folded$observation == i)
    expect_equal(unname(pooled$prediction[1, ]),
                 unname(colMeans(folded$prediction[reps, , drop = FALSE])),
                 tolerance = 1e-10, info = info)
    # truth survives pooling
    expect_equal(as.character(pooled$truth[1]), as.character(folded$truth[reps[1]]), info = info)
  }
  # only the repeated schemes actually repeat
  expect_true(all(res$ledger$n_repeats == 4L))          # bootstrap, nreps = 4
})

test_that("pooling averages continuous predictions and their baseline", {
  simc <- sim_pattern_data(n = 60, dims = c(5, 5, 3), q = 3, r = 2, snr = 2, seed = 51)
  spec <- pattern_model(simc$dataset, simc$design, rank = 2,
                        crossval = bootstrap_blocked_cross_validation(simc$design$block_var, nreps = 3))
  res <- suppressWarnings(quiet_run(run_global(spec)))
  pooled <- res$ledger; folded <- res$fold_ledger
  expect_equal(length(pooled$observation), 60L)
  expect_equal(dim(pooled$prediction), c(60L, 3L))
  expect_equal(dim(pooled$baseline), c(60L, 3L))
  reps <- which(folded$observation == pooled$observation[1])
  expect_gt(length(reps), 1L)
  expect_equal(unname(pooled$prediction[1, ]),
               unname(colMeans(folded$prediction[reps, , drop = FALSE])), tolerance = 1e-10)
  expect_equal(unname(pooled$baseline[1, ]),
               unname(colMeans(folded$baseline[reps, , drop = FALSE])), tolerance = 1e-10)
  expect_true(is.finite(performance(res)$R2))
})

test_that("a repeatedly tested observation is not double counted in the metrics", {
  # Under bootstrap CV a naive fold-order ledger would weight repeatedly tested
  # rows more heavily. Scoring the pooled ledger gives each observation one vote.
  sim <- sim_pattern_data(n = 60, dims = c(5, 5, 3), K = 3, snr = 0.4, seed = 52)
  spec <- pattern_model(sim$dataset, sim$design, rank = 2,
                        crossval = bootstrap_blocked_cross_validation(sim$design$block_var, nreps = 4))
  res <- suppressWarnings(quiet_run(run_global(spec)))
  pooled <- res$ledger
  pred <- factor(colnames(pooled$prediction)[max.col(pooled$prediction, ties.method = "first")],
                 levels = levels(pooled$truth))
  expect_equal(performance(res)$Accuracy, mean(pred == pooled$truth), tolerance = 1e-12)
  # and that differs from the fold-order accuracy when repeats disagree
  f <- res$fold_ledger
  fpred <- factor(colnames(f$prediction)[max.col(f$prediction, ties.method = "first")],
                  levels = levels(f$truth))
  expect_equal(length(f$observation), 240L)
  expect_false(isTRUE(all.equal(mean(fpred == f$truth), performance(res)$Accuracy)))
})

test_that("the regional prediction table has one row per observation per ROI under repeats", {
  sim <- sim_pattern_data(n = 60, dims = c(5, 5, 3), K = 3, snr = 1.5, seed = 53)
  region <- neuroim2::NeuroVol(array(rep(1:3, length.out = 75), sim$dims),
                               neuroim2::space(sim$dataset$mask))
  spec <- pattern_model(sim$dataset, sim$design, rank = 2, return_predictions = TRUE,
                        crossval = bootstrap_blocked_cross_validation(sim$design$block_var, nreps = 4))
  reg <- suppressWarnings(quiet_run(run_regional(spec, region)))
  pt <- reg$prediction_table
  expect_equal(nrow(pt), 180L)                       # 60 observations x 3 ROIs
  expect_false(any(duplicated(pt[, c("roinum", ".rownum")])))
  expect_equal(sort(unique(pt$.rownum)), 1:60)
  expect_s3_class(reg$fits[[1]]$pooled_ledger, "pattern_ledger")
  expect_equal(length(reg$fits[[1]]$ledger$observation), 240L)   # fold-resolved kept too
})

test_that("run_global preflight and run_analysis wiring", {
  sim <- sim_pattern_data(n = 45, dims = c(5, 5, 3), K = 2, seed = 22)
  spec <- pattern_model(sim$dataset, sim$design, rank = 1)
  expect_error(quiet_run(run_global(spec, preflight = "error")), NA)
  run <- quiet_run(run_analysis(mvpa_config(mode = "global", model_spec = spec)))
  expect_s3_class(run, "rmvpa_analysis_run")
})
