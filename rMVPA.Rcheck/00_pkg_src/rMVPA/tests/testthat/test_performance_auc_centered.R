context("chance-centered AUC metrics")

test_that("binary AUC is chance-centered", {
  obs <- factor(c("A", "B", "A", "B"), levels = c("A", "B"))

  # Perfect separation -> AUCc = 1
  probs_perfect <- matrix(c(0.9, 0.1,
                             0.1, 0.9,
                             0.8, 0.2,
                             0.2, 0.8), byrow = TRUE, ncol = 2,
                           dimnames = list(NULL, c("A", "B")))
  res_perfect <- classification_result(obs, obs, probs_perfect)
  perf_perfect <- performance(res_perfect)
  expect_equal(unname(perf_perfect["AUC"]), 1, tolerance = 1e-6)
  expect_equal(unname(perf_perfect["Accuracy"]), 1)

  # Reversed separation -> AUCc = -1
  # Construct scores that flip the ranking: positive class gets lowest scores
  prob_pos <- c(0.9, 0.1, 0.8, 0.2)  # high for negatives (A), low for positives (B)
  probs_rev <- cbind(A = 1 - prob_pos, B = prob_pos)
  pred_rev <- factor(ifelse(prob_pos > 0.5, "B", "A"), levels = c("A", "B"))
  res_rev <- classification_result(obs, pred_rev, probs_rev)
  perf_rev <- performance(res_rev)
  expect_equal(unname(perf_rev["AUC"]), -1, tolerance = 1e-6)
})


test_that("multiclass AUC is chance-centered and class metrics are opt-in", {
  obs <- factor(c("A", "B", "C", "A", "B", "C"), levels = c("A", "B", "C"))

  # Perfect classifier: diag probabilities
  probs_perfect <- rbind(
    c(0.9, 0.05, 0.05),
    c(0.05, 0.9, 0.05),
    c(0.05, 0.05, 0.9),
    c(0.9, 0.05, 0.05),
    c(0.05, 0.9, 0.05),
    c(0.05, 0.05, 0.9)
  )
  colnames(probs_perfect) <- c("A", "B", "C")
  res_mc <- classification_result(obs, obs, probs_perfect)

  # No class metrics by default
  perf_mc <- performance(res_mc)
  expect_equal(unname(perf_mc["Accuracy"]), 1)
  expect_equal(unname(perf_mc["AUC"]), 1, tolerance = 1e-6)
  expect_false(any(grepl("AUC_A", names(perf_mc))))

  # With per-class metrics
  perf_mc_class <- performance(res_mc, class_metrics = TRUE)
  expect_true(all(c("AUC_A", "AUC_B", "AUC_C") %in% names(perf_mc_class)))
  expect_true(all(abs(perf_mc_class[c("AUC_A", "AUC_B", "AUC_C")] - 1) < 1e-6))
})
