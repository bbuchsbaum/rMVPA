library(testthat)
library(tibble)

context("hrfdecoder merge_results zero-probability filtering")

test_that("merge_results filters zero-probability duplicate events across folds", {
  # Minimal obj with required class and flags
  obj <- structure(list(compute_performance = FALSE, performance = NULL), class = "hrfdecoder_model")

  # Two-class problem with three events
  classes <- c("A", "B")

  # Fold 1: e1, e2 have predictions; e3 is zero-row (belongs to other fold)
  probs1 <- matrix(c(0.8, 0.2,
                     0.1, 0.9,
                     0.0, 0.0), nrow = 3, byrow = TRUE)
  colnames(probs1) <- classes
  y_true1 <- factor(c("A", "B", "A"), levels = classes)
  class1 <- factor(c("A", "B", "A"), levels = classes)

  # Fold 2: e1,e2 zeros; e3 has prediction
  probs2 <- matrix(c(0.0, 0.0,
                     0.0, 0.0,
                     0.3, 0.7), nrow = 3, byrow = TRUE)
  colnames(probs2) <- classes
  y_true2 <- factor(c("A", "B", "B"), levels = classes)
  class2 <- factor(c("A", "B", "B"), levels = classes)

  result_set <- tibble(
    probs = list(probs1, probs2),
    y_true = list(y_true1, y_true2),
    class = list(class1, class2),
    error = c(FALSE, FALSE),
    error_message = c("~", "~")
  )

  out <- merge_results(obj, result_set, indices = list(1:5), id = 1L)
  expect_false(isTRUE(out$error))
  cres <- out$result[[1]]
  # Should keep only the three non-zero rows total
  expect_equal(length(cres$observed), 3L)
  expect_equal(nrow(cres$probs), 3L)
  expect_identical(colnames(cres$probs), classes)
})

