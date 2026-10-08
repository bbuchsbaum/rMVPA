# Multiclass AUC is one-vs-rest on each class's own probability, as in
# scikit-learn's roc_auc_score(multi_class = "ovr"). Reference values were
# computed with scikit-learn 1.9.1 (tools/bench/.venv) on these inputs.

test_that("per-class one-vs-rest AUC equals scikit-learn's", {
  lv <- c("a", "b", "c")
  observed <- factor(rep(lv, length.out = 12), levels = lv)
  probs <- matrix(c(
    0.6, 0.2, 0.2,  0.1, 0.8, 0.1,  0.3, 0.3, 0.4,  0.5, 0.5, 0.0,
    0.2, 0.2, 0.6,  0.1, 0.1, 0.8,  0.4, 0.4, 0.2,  0.3, 0.3, 0.4,
    0.7, 0.2, 0.1,  1/3, 1/3, 1/3,  0.2, 0.6, 0.2,  0.1, 0.45, 0.45
  ), ncol = 3, byrow = TRUE, dimnames = list(NULL, lv))
  pred <- lv[max.col(probs, ties.method = "first")]
  m <- rMVPA:::multiclass_perf(observed, pred, probs, class_metrics = TRUE)
  sklearn_auc <- c(0.875, 0.671875, 0.71875)
  expect_equal(unname((m[paste0("AUC_", lv)] + 1) / 2), sklearn_auc, tolerance = 1e-12)
  expect_equal(unname(m["AUC"]), mean(2 * sklearn_auc - 1), tolerance = 1e-12)
})

test_that("tiny but distinct class probabilities are ranked, not tied", {
  # p_k - mean(p_-k) rounds these to the same -1/2; ranking p_k keeps them apart.
  observed <- factor(c("a", "b", "b", "c"), levels = c("a", "b", "c"))
  probs <- rbind(c(1e-30, 0.5, 0.5), c(1e-40, 0.5, 0.5), c(0, 1, 0), c(0, 0, 1))
  colnames(probs) <- c("a", "b", "c")
  m <- rMVPA:::multiclass_perf(observed, c("b", "b", "b", "c"), probs, class_metrics = TRUE)
  expect_equal(unname(m["AUC_a"]), 1)  # the positive has the largest p_a
})
