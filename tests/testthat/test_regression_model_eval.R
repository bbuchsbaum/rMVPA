# Regression tests for model-evaluation correctness bugs.
#
# Each block targets one confirmed bug. Expected values are computed in the
# test from independent oracles (hand-written arithmetic, aov(), Mann-Whitney
# ranks, or a repeated-prediction invariant), never by re-running the code
# path under test and comparing it to itself.

mk_matrix <- function(rows, cols) {
  matrix(unlist(rows), nrow = length(rows), byrow = TRUE,
         dimnames = list(NULL, cols))
}

# One-way ANOVA F statistic for a single feature column (independent oracle).
oneway_F <- function(x, g) {
  summary(aov(x ~ g))[[1]][["F value"]][1]
}

# Mann-Whitney AUC: P(score_pos > score_neg) + 0.5 * P(tie). Independent of yardstick.
auc_mann_whitney <- function(score, is_pos) {
  r <- rank(score)
  n1 <- sum(is_pos)
  n0 <- sum(!is_pos)
  (sum(r[is_pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

abc_design <- function(y = factor(c("a", "b", "c", "a", "b", "c"), levels = c("a", "b", "c"))) {
  mvpa_design(train_design = data.frame(x = seq_along(y)), y_train = y)
}

# ---------------------------------------------------------------------------
# Bug 1: predict.class_model_fit must align probability columns by class name
# and give classes absent from the fitted model probability 0.
# ---------------------------------------------------------------------------

test_that("align_class_probs maps named subset/reordered columns onto all levels", {
  probs <- mk_matrix(list(c(0.8, 0.2), c(0.3, 0.7)), c("c", "b"))  # subset, reordered
  full <- .align_class_probs(probs, c("a", "b", "c"))
  expect_equal(colnames(full), c("a", "b", "c"))
  expect_equal(unname(full[, "a"]), c(0, 0))   # absent class -> 0
  expect_equal(unname(full[, "b"]), c(0.2, 0.7))
  expect_equal(unname(full[, "c"]), c(0.8, 0.3))
})

test_that("align_class_probs keeps positional fallback for unnamed full-width output", {
  probs <- matrix(c(0.1, 0.2, 0.7, 0.5, 0.25, 0.25), nrow = 2, byrow = TRUE)
  full <- .align_class_probs(probs, c("a", "b", "c"))
  expect_equal(colnames(full), c("a", "b", "c"))
  expect_equal(unname(full), probs)
})

test_that("align_class_probs errors on names it cannot place", {
  probs <- mk_matrix(list(c(0.5, 0.5)), c("x", "y"))
  expect_error(.align_class_probs(probs, c("a", "b", "c")), "probability columns")
})

test_that("predict.class_model_fit labels and argmaxes a model that reports a class subset", {
  # Stub model: returns probability columns only for the classes seen in training
  # (in reversed order), as naive_bayes / sda_notune do when a fold lacks a class.
  P <- mk_matrix(list(c(0.2, 0.8), c(0.6, 0.4)), c("c", "b"))
  obj <- structure(
    list(fit = list(), model = list(prob = function(fit, mat) P),
         y = factor(c("a", "b", "c"), levels = c("a", "b", "c")),
         feature_mask = NULL),
    class = c("class_model_fit", "model_fit"))
  pred <- predict.class_model_fit(obj, matrix(0, 2, 3))
  expect_equal(colnames(pred$probs), c("a", "b", "c"))
  expect_equal(unname(pred$probs[1, ]), c(0, 0.8, 0.2))
  expect_equal(unname(pred$probs[2, ]), c(0, 0.4, 0.6))
  expect_equal(as.character(pred$class), c("b", "c"))
})

test_that("naive_bayes fold lacking a class gives that class zero probability (public path)", {
  set.seed(2)
  dset <- gen_sample_dataset(c(5, 5, 5), nobs = 24, nlevels = 3, data_mode = "image",
                             response_type = "categorical")
  blk <- rep(1:2, each = 12)
  # Class "c" occurs only in block 1, so the fold trained on block 2 lacks it.
  y <- factor(c(rep(c("a", "b", "c"), 4), rep(c("a", "b"), 6)), levels = c("a", "b", "c"))
  des <- mvpa_design(train_design = data.frame(block_var = blk), y_train = y)
  cval <- twofold_blocked_cross_validation(blk)
  region_mask <- neuroim2::NeuroVol(rep(1L, length(dset$dataset$mask)),
                                    neuroim2::space(dset$dataset$mask))
  mspec <- mvpa_model(load_model("naive_bayes"), dset$dataset, des,
                      model_type = "classification", crossval = cval)
  pt <- run_regional(mspec, region_mask)$prediction_table
  expect_true(all(c("prob_a", "prob_b", "prob_c") %in% names(pt)))
  # Rows 1:12 are tested by the fold whose training block has no class "c".
  block1 <- pt[pt$.rownum %in% 1:12, ]
  expect_equal(nrow(block1), 12)
  expect_true(all(block1$prob_c == 0))
  expect_equal(unname(rowSums(as.matrix(pt[, c("prob_a", "prob_b", "prob_c")]))),
               rep(1, nrow(pt)), tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# Bug 2: group_means must divide by the sizes of the groups rowsum() returns.
# ---------------------------------------------------------------------------

test_that("group_means uses per-group counts when a middle level is unused", {
  X <- matrix(1:12, nrow = 6, ncol = 2)  # rows: (1,7) (2,8) (3,9) (4,10) (5,11) (6,12)
  g <- factor(c("a", "a", "c", "c", "c", "e"), levels = c("a", "b", "c", "d", "e"))
  m <- group_means(X, margin = 1, group = g)
  expect_equal(rownames(m), c("a", "c", "e"))
  # Hand-computed group means.
  expected <- rbind(c(1.5, 7.5), c(4, 10), c(6, 12))
  expect_equal(unname(m), expected, tolerance = 1e-12)
  expect_true(all(is.finite(m)))
})

test_that("group_means margin = 2 averages columns by group with the same counts", {
  X <- rbind(1:6, 7:12)  # 2 features x 6 observations
  g <- factor(c("a", "a", "c", "c", "c", "e"), levels = c("a", "b", "c", "d", "e"))
  m <- group_means(X, margin = 2, group = g)
  expect_equal(dim(m), c(2L, 3L))
  expect_equal(unname(m[1, ]), c(1.5, 4, 6))     # feature 1: 1,2 | 3,4,5 | 6
  expect_equal(unname(m[2, ]), c(7.5, 10, 12))   # feature 2: 7,8 | 9,10,11 | 12
})

# ---------------------------------------------------------------------------
# Bug 3: classification test labels must use the training class levels.
# ---------------------------------------------------------------------------

test_that("design aligns y_test levels to training levels (vector and formula)", {
  ytr <- factor(rep(c("a", "b", "c"), each = 2), levels = c("a", "b", "c"))
  yte_raw <- c("a", "b", "a", "b")
  des_vec <- mvpa_design(train_design = data.frame(block = 1:6),
                         test_design = data.frame(dummy = 1:4),
                         y_train = ytr, y_test = yte_raw)
  expect_equal(levels(des_vec$y_test), c("a", "b", "c"))
  expect_equal(as.character(des_vec$y_test), yte_raw)

  te_df <- data.frame(cond = factor(yte_raw, levels = c("a", "b")))
  des_f <- mvpa_design(train_design = data.frame(block = 1:6), test_design = te_df,
                       y_train = ytr, y_test = ~ cond)
  expect_equal(levels(des_f$y_test), c("a", "b", "c"))
  expect_equal(as.character(des_f$y_test), yte_raw)
})

test_that("design errors when test labels contain a class absent from training", {
  ytr <- factor(rep(c("a", "b", "c"), each = 2), levels = c("a", "b", "c"))
  expect_error(
    mvpa_design(train_design = data.frame(block = 1:6),
                test_design = data.frame(dummy = 1:3),
                y_train = ytr, y_test = c("a", "z", "b")),
    "not present in the training labels.*z")
})

test_that("design still drops unused training levels (zero-count check unchanged)", {
  ytr <- factor(c("a", "a", "b", "b", "c", "c"), levels = c("a", "b", "c", "d"))
  des <- mvpa_design(train_design = data.frame(block = 1:6), y_train = ytr)
  expect_equal(levels(des$y_train), c("a", "b", "c"))
})

test_that("external-test run with a test level subset labels rows by test truth (public path)", {
  set.seed(3)
  dset <- gen_sample_dataset(c(5, 5, 5), nobs = 24, nlevels = 3, data_mode = "image",
                             response_type = "categorical", external_test = TRUE, ntest_obs = 12)
  blk <- rep(1:2, each = 12)
  y <- factor(rep(c("a", "b", "c"), length.out = 24), levels = c("a", "b", "c"))
  yt_chr <- rep(c("a", "b"), length.out = 12)
  des <- mvpa_design(train_design = data.frame(block_var = blk), y_train = y,
                     test_design = data.frame(yt = factor(yt_chr, levels = c("a", "b"))),
                     y_test = ~ yt)
  region_mask <- neuroim2::NeuroVol(rep(1L, length(dset$dataset$mask)),
                                    neuroim2::space(dset$dataset$mask))
  mspec <- mvpa_model(load_model("sda_notune"), dset$dataset, des,
                      model_type = "classification", crossval = blocked_cross_validation(blk))
  pt <- run_regional(mspec, region_mask)$prediction_table
  expect_equal(nrow(pt), 12)
  # Observed labels are the test truth for each test row (.rownum indexes the test set).
  expect_equal(as.character(pt$observed), yt_chr[pt$.rownum])
  prob_cols <- as.matrix(pt[, c("prob_a", "prob_b", "prob_c")])
  expect_equal(unname(rowSums(prob_cols)), rep(1, nrow(pt)), tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# Bug 4: wrap_result must align probs by name, sum duplicated test indices,
# average regression predictions by true counts, and subset observed.
# ---------------------------------------------------------------------------

test_that("wrap_result aligns class probabilities by name across folds", {
  des <- abc_design()
  # Fold 1 reports only b and a, in that order (class c absent in its training fold).
  p1 <- mk_matrix(list(c(0.9, 0.1), c(0.3, 0.7), c(0.5, 0.5)), c("b", "a"))
  # Fold 2 reports all classes in a different order.
  p2 <- mk_matrix(list(c(0.2, 0.3, 0.5), c(0.6, 0.3, 0.1), c(0.1, 0.1, 0.8)), c("c", "a", "b"))
  rt <- list(test_ind = list(1:3, 4:6), probs = list(p1, p2))
  res <- wrap_result(rt, des)
  expect_equal(res$testind, 1:6)
  P <- as.matrix(res$probs)
  expect_equal(colnames(P), c("a", "b", "c"))
  # Hand-computed rows in (a, b, c) order.
  expect_equal(unname(P[1, ]), c(0.1, 0.9, 0))   # fold 1: b=.9, a=.1, c absent
  expect_equal(unname(P[2, ]), c(0.7, 0.3, 0))   # fold 1: b=.3, a=.7
  expect_equal(unname(P[3, ]), c(0.5, 0.5, 0))   # fold 1: b=.5, a=.5
  expect_equal(unname(P[4, ]), c(0.3, 0.5, 0.2)) # fold 2: c=.2, a=.3, b=.5
  expect_equal(unname(P[5, ]), c(0.3, 0.1, 0.6)) # fold 2: c=.6, a=.3, b=.1
  expect_equal(unname(P[6, ]), c(0.1, 0.8, 0.1)) # fold 2: c=.1, a=.1, b=.8
  expect_equal(as.character(res$predicted), c("b", "a", "a", "b", "c", "b"))
})

test_that("wrap_result sums repeated test indices within a fold (oversampled test set)", {
  des <- abc_design()
  # Observation 1 appears twice in one fold with different probabilities.
  p <- mk_matrix(list(c(0.2, 0.3, 0.5), c(0.6, 0.2, 0.2), c(0.1, 0.1, 0.8)), c("a", "b", "c"))
  rt <- list(test_ind = list(c(1L, 1L, 2L)), probs = list(p))
  res <- wrap_result(rt, des)
  expect_equal(res$testind, c(1L, 2L))
  P <- as.matrix(res$probs)
  # Summed row for obs 1 is (0.8, 0.5, 0.7), normalised to sum one.
  expect_equal(unname(P[1, ]), c(0.8, 0.5, 0.7) / 2, tolerance = 1e-12)
  expect_equal(unname(P[2, ]), c(0.1, 0.1, 0.8))
})

test_that("wrap_result regression averages duplicates by true counts and subsets observed", {
  y <- c(1.1, 2.2, 3.3, 4.4, 5.5, 6.6)
  des <- mvpa_design(train_design = data.frame(x = 1:6), y_train = y)
  # Obs 1 predicted twice (1 and 3), obs 2 once (5); tested obs are 1:4 only.
  rt <- list(test_ind = list(c(1L, 1L, 2L), c(3L, 4L)), preds = list(c(1, 3, 5), c(10, 20)))
  res <- wrap_result(rt, des)
  expect_equal(res$testind, 1:4)
  expect_equal(res$predicted, c((1 + 3) / 2, 5, 10, 20), tolerance = 1e-12)
  # observed must be subset to the tested rows, not the full design.
  expect_equal(res$observed, y[1:4])
  expect_length(res$observed, length(res$testind))
})

# Minimal deterministic least-squares regression spec, passed directly to
# mvpa_model() (no registry changes). Installed regression backends are not
# available in every environment, so the public regression path uses this.
make_ols_regression_model <- function() {
  list(
    type = "Regression",
    label = "ols_test",
    library = NULL,
    loop = NULL,
    parameters = data.frame(parameter = "none", class = "character", label = "none"),
    grid = function(x, y, len = NULL) data.frame(none = 1),
    fit = function(x, y, wts = NULL, param = NULL, lev = NULL, last = FALSE,
                   classProbs = FALSE, ...) {
      X <- cbind(1, as.matrix(x))
      list(coef = solve(crossprod(X) + diag(1e-8, ncol(X)), crossprod(X, y)))
    },
    predict = function(modelFit, newdata, ...) {
      as.numeric(cbind(1, as.matrix(newdata)) %*% modelFit$coef)
    }
  )
}

test_that("regression run leaves untested rows out and observed matches tested rows (public path)", {
  set.seed(4)
  dset <- gen_sample_dataset(c(5, 5, 5), 48, response_type = "continuous")
  y <- dset$design$y_train
  region_mask <- neuroim2::NeuroVol(rep(1L, length(dset$dataset$mask)),
                                    neuroim2::space(dset$dataset$mask))
  cv <- custom_cross_validation(list(
    list(train = 1:24, test = 25:36),
    list(train = 25:48, test = 1:12)))   # rows 13:24 and 37:48 are never tested
  mspec <- mvpa_model(make_ols_regression_model(), dset$dataset, dset$design,
                      model_type = "regression", crossval = cv)
  pt <- run_regional(mspec, region_mask)$prediction_table
  pt <- pt[order(pt$.rownum), ]
  expect_equal(pt$.rownum, c(1:12, 25:36))
  expect_equal(pt$observed, unname(y[pt$.rownum]), tolerance = 1e-12)
})

test_that("regression run with a repeated test observation keeps its prediction (public path)", {
  set.seed(4)
  dset <- gen_sample_dataset(c(5, 5, 5), 48, response_type = "continuous")
  region_mask <- neuroim2::NeuroVol(rep(1L, length(dset$dataset$mask)),
                                    neuroim2::space(dset$dataset$mask))
  run_with <- function(fold2_test) {
    cv <- custom_cross_validation(list(
      list(train = 1:24, test = 25:48),
      list(train = 25:48, test = fold2_test)))
    mspec <- mvpa_model(make_ols_regression_model(), dset$dataset, dset$design,
                        model_type = "regression", crossval = cv)
    run_regional(mspec, region_mask)$prediction_table
  }
  # All rows tested. The repeated run tests observations 1:4 twice in fold 2.
  # A repeated observation gets the same model output, so its averaged
  # prediction must equal the single-test run.
  pu <- run_with(1:24)
  pd <- run_with(c(1:24, 1:4))
  pu <- pu[order(pu$.rownum), ]
  pd <- pd[order(pd$.rownum), ]
  expect_equal(pd$.rownum, 1:48)
  expect_equal(pd$predicted, pu$predicted, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# Bug 5: matrixAnova with unused levels; cutoffs via validate_cutoff.
# ---------------------------------------------------------------------------

test_that("matrixAnova F statistics match aov() when a middle level is unused", {
  set.seed(5)
  X <- matrix(rnorm(30 * 4), 30, 4)
  g <- factor(rep(c("a", "c", "d"), length.out = 30), levels = c("a", "b", "c", "d"))
  tab <- matrixAnova(g, X)
  expect_false(any(is.nan(tab)))
  oracle <- apply(X, 2, function(col) oneway_F(col, droplevels(g)))
  expect_equal(unname(tab[, "Ftest"]), unname(oracle), tolerance = 1e-8)
})

test_that("FTest selects the top-k features by ANOVA p-value when a level is unused", {
  set.seed(6)
  X <- matrix(rnorm(30 * 8), 30, 8)
  g <- factor(rep(c("a", "c", "d"), length.out = 30), levels = c("a", "b", "c", "d"))
  X[g == "d", 2] <- X[g == "d", 2] + 3   # strong signal in feature 2
  X[g == "c", 5] <- X[g == "c", 5] + 2   # moderate signal in feature 5
  pvals <- apply(X, 2, function(col) summary(aov(col ~ droplevels(g)))[[1]][["Pr(>F)"]][1])
  keep <- select_features.FTest(feature_selector("FTest", "topk", 3), X, g)
  expect_equal(which(keep), sort(order(pvals)[1:3]))
})

test_that("FTest top_p rounds up like validate_cutoff (ceiling)", {
  set.seed(7)
  X <- matrix(rnorm(30 * 10), 30, 10)
  g <- factor(rep(c("a", "b"), 15))
  keep <- select_features.FTest(feature_selector("FTest", "top_p", 0.25), X, g)
  expect_equal(sum(keep), 3L)  # ceiling(0.25 * 10) = 3
  expect_equal(sum(keep), validate_cutoff("top_p", 0.25, 10))
})

test_that("catscore accepts topk and top_p with the shared cutoff rounding", {
  skip_if_not_installed("sda")
  set.seed(8)
  X <- matrix(rnorm(30 * 10), 30, 10)
  g <- factor(rep(c("a", "b", "c"), length.out = 30))
  ranking <- sda::sda.ranking(X, g, ranking.score = "entropy", fdr = FALSE, verbose = FALSE)[, "idx"]
  keep_topk <- select_features.catscore(feature_selector("catscore", "topk", 3), X, g)
  expect_equal(which(keep_topk), sort(ranking[1:3]))
  keep_topp <- select_features.catscore(feature_selector("catscore", "top_p", 0.25), X, g)
  expect_equal(sum(keep_topp), 3L)  # ceiling(0.25 * 10), not truncated to 2
  expect_equal(which(keep_topp), sort(ranking[1:3]))
})

test_that("feature selection emits no per-fold messages", {
  set.seed(9)
  X <- matrix(rnorm(20 * 5), 20, 5)
  g <- factor(rep(c("a", "b"), 10))
  expect_silent(select_features.FTest(feature_selector("FTest", "topk", 2), X, g))
})

# ---------------------------------------------------------------------------
# Bug 6: combine_regional_results must skip failed ROIs and index by name.
# ---------------------------------------------------------------------------

test_that("combine_regional_results skips a failed first ROI and uses name-indexed pobserved", {
  ok_res <- list(
    observed = factor(c("a", "b", "a"), levels = c("a", "b")),
    predicted = c("a", "b", "b"),
    # Columns deliberately reversed relative to the factor levels.
    # Row 1: b=.7, a=.3 ; row 2: b=.2, a=.8 ; row 3: b=.6, a=.4
    probs = mk_matrix(list(c(0.7, 0.3), c(0.2, 0.8), c(0.6, 0.4)), c("b", "a")),
    testind = c(4L, 5L, 6L))
  results <- tibble::tibble(id = c(1L, 2L), result = list(NULL, ok_res),
                            error = c(TRUE, FALSE), error_message = c("boom", "~"))
  out <- as.data.frame(combine_regional_results(results))
  expect_equal(nrow(out), 3)
  expect_equal(unique(out$roinum), 2L)
  expect_equal(out$.rownum, c(4L, 5L, 6L))
  # Probability of each observation's true class (a, b, a): .3, .2, .4.
  expect_equal(out$pobserved, c(0.3, 0.2, 0.4), tolerance = 1e-12)
  expect_equal(out$prob_a, c(0.3, 0.8, 0.4), tolerance = 1e-12)
  expect_equal(out$predicted, c("a", "b", "b"))
  expect_equal(out$correct, c(TRUE, TRUE, FALSE))
})

test_that("combine_regional_results picks the regression branch from the first non-NULL result", {
  reg_res <- list(observed = c(1.5, 2.5), predicted = c(1, 3), testind = c(2L, 3L))
  results <- tibble::tibble(id = c(1L, 2L), result = list(NULL, reg_res),
                            error = c(TRUE, FALSE), error_message = c("boom", "~"))
  out <- as.data.frame(combine_regional_results(results))
  expect_equal(nrow(out), 2)
  expect_equal(out$.rownum, c(2L, 3L))
  expect_equal(out$observed, c(1.5, 2.5))
  expect_equal(out$predicted, c(1, 3))
})

test_that("combine_regional_results returns an empty table when every ROI failed", {
  results <- tibble::tibble(id = 1:2, result = list(NULL, NULL),
                            error = c(TRUE, TRUE), error_message = c("a", "b"))
  expect_equal(nrow(as.data.frame(combine_regional_results(results))), 0)
})

# ---------------------------------------------------------------------------
# Bug 7: performance helpers index probability columns by class name.
# ---------------------------------------------------------------------------

test_that("prob_observed looks up the true-class column by name", {
  obs <- factor(c("neg", "pos", "neg"), levels = c("neg", "pos"))
  # Columns reversed (pos, neg). Row 1: neg=.9 ; row 2: pos=.6 ; row 3: neg=.25.
  probs <- mk_matrix(list(c(0.1, 0.9), c(0.6, 0.4), c(0.75, 0.25)), c("pos", "neg"))
  r <- binary_classification_result(obs, c("neg", "pos", "pos"), probs)
  expect_equal(unname(prob_observed(r)), c(0.9, 0.6, 0.25), tolerance = 1e-12)
})

test_that("multiclass AUC scores each class by its own probability column", {
  set.seed(10)
  n <- 30
  obs <- factor(rep(c("a", "b", "c"), length.out = n), levels = c("a", "b", "c"))
  P <- matrix(runif(n * 3), n, 3)
  P <- P / rowSums(P)
  colnames(P) <- c("c", "a", "b")          # deliberately not in level order
  pred <- factor(c("a", "b", "c")[max.col(P[, c("a", "b", "c")], ties.method = "first")],
                 levels = c("a", "b", "c"))
  perf <- multiclass_perf(obs, pred, P, class_metrics = TRUE)
  auc_k <- sapply(c("a", "b", "c"), function(k) auc_mann_whitney(P[, k], obs == k))
  expect_equal(unname(perf[["AUC_a"]]), 2 * auc_k[["a"]] - 1, tolerance = 1e-8)
  expect_equal(unname(perf[["AUC_b"]]), 2 * auc_k[["b"]] - 1, tolerance = 1e-8)
  expect_equal(unname(perf[["AUC_c"]]), 2 * auc_k[["c"]] - 1, tolerance = 1e-8)
  expect_equal(unname(perf[["AUC"]]), mean(2 * auc_k - 1), tolerance = 1e-8)
  expect_equal(unname(perf[["Accuracy"]]), mean(obs == pred))
})

test_that("dead combinedACC helper is removed", {
  expect_false(exists("combinedACC", envir = asNamespace("rMVPA"), inherits = FALSE))
})

# ---------------------------------------------------------------------------
# Bug 8: return_fits must be read by exact name (no partial matching).
# ---------------------------------------------------------------------------

test_that("merge_results.mvpa_model reads return_fits exactly, not via partial match", {
  des <- abc_design()
  rs <- tibble::tibble(error = FALSE, error_message = "~",
                       test_ind = list(1:3),
                       probs = list(mk_matrix(list(c(0.2, 0.3, 0.5), c(0.6, 0.2, 0.2),
                                                   c(0.1, 0.1, 0.8)), c("a", "b", "c"))),
                       fit = list(NULL))
  # Two fields start with "return_fit": `obj$return_fit` is then ambiguous and
  # returns NULL, so the old code errors. The explicit field is FALSE.
  obj <- list(design = des, return_fits = FALSE, return_fit_debug = TRUE,
              compute_performance = FALSE)
  out <- merge_results.mvpa_model(obj, rs, indices = 1:3, id = 1L)
  expect_false(out$error)
  # Row 1 of the probs is (a, b, c) = (.2, .3, .5); it already sums to one.
  expect_equal(unname(as.matrix(out$result[[1]]$probs)[1, ]), c(0.2, 0.3, 0.5), tolerance = 1e-12)
})
