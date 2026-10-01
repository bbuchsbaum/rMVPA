# Predicted classes are the exact argmax, with exact ties going to the first
# class (as in numpy/scikit-learn argmax and CoSMoMVPA). max.col()'s default
# broke ties at random *and* treated values within a relative 1e-5 as tied, so
# it could return a class that was not the maximum.

test_that("corclass prediction is the exact argmax and ignores the RNG", {
  fit <- list(conditionMeans = rbind(a = c(1, 0, 0), b = c(0, 1, 0)), levs = c("a", "b"),
              method = "pearson", robust = FALSE)
  newdata <- rbind(c(0.5, 0.5, -1),   # exact tie between a and b -> first ("a")
                   c(0.2, 1, 0))      # clear b
  model <- rMVPA:::MVPAModels$corclass
  preds <- lapply(1:20, function(s) { set.seed(s); as.character(model$predict(fit, newdata)) })
  expect_true(all(vapply(preds, identical, logical(1), preds[[1]])))
  expect_identical(preds[[1]], c("a", "b"))
})

test_that("predicted_class picks the exact maximum, not a near-tie, regardless of the RNG", {
  probs <- rbind(c(0.1250849, 0.1250855, 0.7498296),   # clear third class
                 c(0.5000049, 0.4999951, 0),            # x by a relative 2e-5
                 c(0.50000001, 0.49999999, 0),          # x by < 1e-5: was a coin flip
                 c(0.4, 0.4, 0.2))                      # exact tie: first
  colnames(probs) <- c("x", "y", "z")
  out <- lapply(1:25, function(s) { set.seed(s); predicted_class(probs) })
  expect_true(all(vapply(out, identical, logical(1), out[[1]])))
  expect_identical(out[[1]], c("z", "x", "x", "x"))
})
