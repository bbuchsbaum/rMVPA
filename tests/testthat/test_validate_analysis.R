test_that("validate_analysis passes for well-formed blocked design", {
  des_df <- data.frame(
    condition = factor(rep(c("A","B"), each = 50)),
    run = rep(1:5, each = 20)
  )
  des <- mvpa_design(des_df, y_train = ~ condition, block_var = ~ run)
  cv  <- blocked_cross_validation(des$block_var)

  result <- validate_analysis(des, crossval = cv, verbose = FALSE)

  expect_s3_class(result, "validation_result")
  expect_equal(result$n_fail, 0)
  expect_true(result$n_pass > 0)
})


test_that("validate_analysis detects kfold CV with blocking variable", {
  des_df <- data.frame(
    condition = factor(rep(c("A","B"), each = 50)),
    run = rep(1:5, each = 20)
  )
  des <- mvpa_design(des_df, y_train = ~ condition, block_var = ~ run)
  cv_bad <- kfold_cross_validation(len = 100, nfolds = 5)

  result <- validate_analysis(des, crossval = cv_bad, verbose = FALSE)

  expect_s3_class(result, "validation_result")
  expect_true(result$n_fail >= 1)
  # Should flag cv_block_alignment as FAIL
  check_names <- vapply(result$checks, `[[`, character(1), "name")
  statuses <- vapply(result$checks, `[[`, character(1), "status")
  expect_true("cv_block_alignment" %in% check_names)
  expect_equal(statuses[check_names == "cv_block_alignment"], "fail")
})


test_that("validate_analysis warns about too few blocks", {
  des_df <- data.frame(
    condition = factor(rep(c("A","B"), each = 20)),
    run = rep(1:2, each = 20)
  )
  des <- mvpa_design(des_df, y_train = ~ condition, block_var = ~ run)
  cv  <- blocked_cross_validation(des$block_var)

  result <- validate_analysis(des, crossval = cv, verbose = FALSE)

  check_names <- vapply(result$checks, `[[`, character(1), "name")
  statuses <- vapply(result$checks, `[[`, character(1), "status")
  expect_true("block_count" %in% check_names)
  expect_equal(statuses[check_names == "block_count"], "warn")
})


test_that("validate_analysis detects single-class test folds", {
  # Create a design where one class only appears in one block
  des_df <- data.frame(
    condition = factor(c(rep("A", 30), rep("B", 10))),
    run = c(rep(1, 10), rep(2, 10), rep(3, 10), rep(3, 10))
  )
  des <- mvpa_design(des_df, y_train = ~ condition, block_var = ~ run)
  cv  <- blocked_cross_validation(des$block_var)

  result <- validate_analysis(des, crossval = cv, verbose = FALSE)

  # Block 1 and 2 have only class A — when they're held out as test,
  # the test fold has only one class
  check_names <- vapply(result$checks, `[[`, character(1), "name")
  expect_true("single_class_test" %in% check_names ||
              "block_class_coverage" %in% check_names)
})


test_that("validate_analysis warns about very few observations per class", {
  des_df <- data.frame(
    condition = factor(c("A","A","A","B","B","B")),
    run = rep(1:3, each = 2)
  )
  des <- mvpa_design(des_df, y_train = ~ condition, block_var = ~ run)
  cv  <- blocked_cross_validation(des$block_var)

  result <- validate_analysis(des, crossval = cv, verbose = FALSE)

  check_names <- vapply(result$checks, `[[`, character(1), "name")
  statuses <- vapply(result$checks, `[[`, character(1), "status")
  expect_true("min_observations" %in% check_names)
  expect_equal(statuses[check_names == "min_observations"], "warn")
})


test_that("validate_analysis works with mvpa_model object", {
  dat <- gen_sample_dataset(c(5,5,5), 100, blocks = 5)
  model <- mvpa_model(
    load_model("corclass"),
    dataset = dat$dataset,
    design = dat$design,
    crossval = blocked_cross_validation(dat$design$block_var)
  )

  result <- validate_analysis(model, verbose = FALSE)
  expect_s3_class(result, "validation_result")
  expect_equal(result$n_fail, 0)
})


test_that("validate_analysis print method works", {
  des_df <- data.frame(
    condition = factor(rep(c("A","B"), each = 50)),
    run = rep(1:5, each = 20)
  )
  des <- mvpa_design(des_df, y_train = ~ condition, block_var = ~ run)
  cv  <- blocked_cross_validation(des$block_var)

  result <- validate_analysis(des, crossval = cv, verbose = FALSE)
  expect_output(print(result), "MVPA Analysis Validation")
  expect_output(print(result), "Summary:")
})


test_that("validate_analysis warns about no block_var", {
  des_df <- data.frame(
    condition = factor(rep(c("A","B"), each = 50))
  )
  des <- mvpa_design(des_df, y_train = ~ condition)
  cv  <- kfold_cross_validation(len = 100, nfolds = 5)

  result <- validate_analysis(des, crossval = cv, verbose = FALSE)

  check_names <- vapply(result$checks, `[[`, character(1), "name")
  # Should NOT flag cv_block_alignment as fail (no block_var means kfold is OK)
  # But should warn about no blocking for fMRI
  expect_false("cv_block_alignment" %in% check_names[
    vapply(result$checks, `[[`, character(1), "status") == "fail"])
})
