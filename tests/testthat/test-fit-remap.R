test_that("frozen scores match spatial Pearson for single trials and batches", {
  set.seed(141)
  source <- matrix(rnorm(36), 9, 4, dimnames = list(NULL, letters[1:4]))
  keys <- factor(rep(c("red", "green", "blue"), 3), levels = c("blue", "green", "red", "unused"))
  test <- matrix(rnorm(20), 5, 4, dimnames = list(NULL, colnames(source)))
  fit <- fit_remap(source, keys, shrink_whiten = FALSE)
  prototypes <- t(vapply(levels(droplevels(keys)), function(k) colMeans(source[keys == k, ]), numeric(4)))
  expect_equal(unname(fit$mu), unname(colMeans(prototypes)))
  # Derive the target and template scores from the source sample covariance.
  eig <- eigen(cov(prototypes), symmetric = TRUE)
  W <- eig$vectors %*% diag(1 / sqrt(pmax(eig$values, 1e-7))) %*% t(eig$vectors)
  X <- sweep(prototypes, 2, colMeans(prototypes)) %*% W
  Y <- sweep(test, 2, colMeans(prototypes)) %*% W
  oracle <- t(vapply(seq_len(nrow(Y)), function(i) apply(X, 1, cor, y = Y[i, ]), numeric(nrow(X))))
  colnames(oracle) <- rownames(X)
  scores <- predict(fit, test)
  expect_equal(scores, oracle, tolerance = 1e-8)
  expect_equal(predict(fit, test[2, , drop = FALSE]), scores[2, , drop = FALSE])
  expect_equal(predict(unserialize(serialize(fit, NULL)), test), scores)
  expect_equal(rowSums(predict(fit, test, "prob")), rep(1, 5))
  expect_equal(as.character(predict(fit, test, "class")), colnames(scores)[max.col(scores)])
  expect_error(predict(fit, test[, 4:1]), "Features must match")
  expect_error(predict(fit, unname(test)), "Features must match")
  expect_error(fit_remap(source, keys, source[, 4:1], keys), "Features must match")
  expect_error(fit_remap(source, keys[-1]), "one nonmissing key")
  expect_error(fit_remap(source, keys, rank = 1), "Without calibration")
  expect_error(fit_remap(source, keys, lambda_grid = NA_real_), "lambda_grid")
})

test_that("two-key source-only whitening has a well-defined default", {
  x <- rbind(left = c(-1, 1, 0), right = c(1, -1, 0))
  fit <- fit_remap(x)
  expect_equal(unname(predict(fit, x)), matrix(c(1, -1, -1, 1), 2))
  expect_equal(unname(predict(fit, matrix(0, 1, 3))), matrix(0, 1, 2))
})

test_that("positive rank recovers a known reflection on fresh target observations", {
  skip_if_not_installed("rrpack")
  set.seed(822)
  source <- rbind(diag(4), -diag(4))
  rownames(source) <- paste0("class", seq_len(nrow(source)))
  calibration <- source
  calibration[, 1] <- -calibration[, 1]
  held <- calibration[rep(seq_len(8), 3), ] + matrix(rnorm(96, sd = .03), 24, 4)
  truth <- rownames(held)
  source_fit <- fit_remap(source, shrink_whiten = FALSE)
  for (ridge in list(NULL, .01)) {
    fit <- fit_remap(source, calibration = calibration, rank = 1, lambda_grid = 1,
                     ridge_rrr_lambda = ridge, shrink_whiten = FALSE)
    expect_equal(fit$rank, 1)
    expect_gt(sum(fit$Delta^2), 0)
    expect_equal(as.character(predict(fit, held, "class")), truth)
    expect_lt(mean(as.character(predict(source_fit, held, "class")) == truth), 1)
  }
  partial <- fit_remap(source, calibration = calibration[-c(2, 6), ], rank = 1,
                       lambda_grid = 1, shrink_whiten = FALSE)
  expect_setequal(colnames(predict(partial, held)), rownames(source))
  expect_equal(predict(partial, held[9, , drop = FALSE]), predict(partial, held)[9, , drop = FALSE])
  unknown <- calibration
  rownames(unknown)[1] <- "unknown"
  expect_error(fit_remap(source, calibration = unknown), "present in source")
})

test_that("the ROI split freezes preprocessing before observing assessment rows or labels", {
  skip_if_not_installed("rrpack")
  set.seed(44)
  sim <- gen_sample_dataset(D = c(4, 4, 4), nobs = 60, nlevels = 3, blocks = 3, external_test = TRUE)
  design <- sim$design
  keys <- y_test(design)
  # Two observations of each class in calibration, every other row in assessment.
  cal <- sort(unname(unlist(lapply(split(seq_along(keys), keys), head, 2L))))
  test <- setdiff(seq_along(keys), cal)
  make <- function(d, ...) remap_rrr_model(sim$dataset, d, rank = 1, lambda_grid = 1,
      return_adapter = TRUE, shrink_whiten = FALSE, compute_naive_baseline = TRUE, ...)
  roi <- extract_roi(which(sim$dataset$mask > 0)[1:8], sim$dataset)
  model <- make(design, calibration_rows = cal)
  inputs <- list(train_data = as.matrix(neuroim2::values(roi$train_roi)),
                 test_data = as.matrix(neuroim2::values(roi$test_roi)),
                 indices = neuroim2::indices(roi$train_roi))
  run <- function(m, x) fit_roi(m, x, list(id = 1L))
  outer <- run(model, inputs)
  a <- process_roi(model, roi, 1)$result[[1]]
  expect_equal(a$probs, outer$result$probs)
  baseline <- fit_remap(inputs$train_data, y_train(design), shrink_whiten = FALSE)
  expected <- predict(baseline, inputs$test_data[test, , drop = FALSE], "class")
  expect_equal(unname(outer$metrics["source_only_Accuracy"]),
               mean(as.character(expected) == as.character(keys[test])))
  expect_true(model$.preserve_roi_features)
  expect_s3_class(a$predictor, "remap_fit")
  expect_equal(as.integer(a$testind), test)
  expect_equal(unname(a$probs), unname(predict(a$predictor, inputs$test_data[test, , drop = FALSE], "prob")))
  # Reassign only the withheld labels, preserving the class set.
  changed <- design
  changed$y_test[test] <- rev(keys[test])
  b <- run(make(changed, calibration_rows = cal), inputs)$result
  expect_identical(a$predictor, b$predictor)
  expect_equal(a$probs, b$probs)
  changed_roi <- inputs
  changed_roi$test_data[test, ] <- 4 + 3 * inputs$test_data[test, ]
  mutated <- run(model, changed_roi)$result
  expect_identical(mutated$predictor, a$predictor)
  expect_error(make(design, calibration_rows = c(cal, cal)), "distinct")
  expect_error(make(design, calibration_rows = seq_along(keys)), "proper subset")
  expect_error(make(design, calibration_rows = cal, leave_one_key_out = TRUE), "leave_one_key_out")
})
