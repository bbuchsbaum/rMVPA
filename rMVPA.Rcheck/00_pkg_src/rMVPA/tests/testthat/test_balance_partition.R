library(testthat)
library(rMVPA)
library(tibble)
library(dplyr)

context("balance_partitions")

# Helper function to check if counts are balanced within a vector
is_balanced <- function(targets) {
  counts <- table(targets)
  if (length(counts) <= 1) return(TRUE) # Cannot balance 0 or 1 class
  all(counts == counts[1])
}

# --- Test Setup --- 

# Design 1: Imbalanced Binary (more 'b')
design_df_bin <- data.frame(
  condition = factor(rep(c("a", "b", "b", "b"), 25)), # 25 'a', 75 'b'
  block = rep(1:5, each = 20) # 5 blocks, 20 samples each
)
des_bin <- mvpa_design(design_df_bin, y_train = ~ condition, block_var = ~ block)
y_bin <- des_bin$y_train

# Design 2: Imbalanced Multi-class (more 'c')
design_df_multi <- data.frame(
  condition = factor(rep(c("a", "b", "c", "c", "c", "c"), 20)), # 20 'a', 20 'b', 80 'c'
  block = rep(1:6, each = 20) # 6 blocks, 20 samples each
)
des_multi <- mvpa_design(design_df_multi, y_train = ~ condition, block_var = ~ block)
y_multi <- des_multi$y_train

# Create various CV objects
cval_blocked_bin <- blocked_cross_validation(des_bin$block_var)
cval_kfold_bin <- kfold_cross_validation(length(y_bin), nfolds = 5)
cval_twofold_bin <- twofold_blocked_cross_validation(des_bin$block_var, nreps = 4) # Less reps for speed

cval_blocked_multi <- blocked_cross_validation(des_multi$block_var)

# --- Test Subsampling --- 

test_that("balance_partitions subsample works for blocked_cv (binary)", {
  cval_bal <- balance_partitions(cval_blocked_bin, des_bin, method = "subsample", seed = 1)
  expect_s3_class(cval_bal, "custom_cross_validation")
  
  samples_bal <- crossval_samples(cval_bal, design_df_bin, y_bin)
  
  for (i in 1:nrow(samples_bal)) {
    # Check training set balance
    ytrain_bal <- samples_bal$ytrain[[i]]
    expect_true(is_balanced(ytrain_bal))
    # Check test set balance (default balance_test_set = TRUE)
    ytest_bal <- samples_bal$ytest[[i]]
    expect_true(is_balanced(ytest_bal))
    # Check total size (should be smaller)
    samples_unbal <- crossval_samples(cval_blocked_bin, design_df_bin, y_bin)
    expect_lt(length(ytrain_bal), length(samples_unbal$ytrain[[i]]))
    # Cannot guarantee test set is smaller if a block was already balanced or had only one class
  }
})

test_that("balance_partitions subsample works for kfold_cv (binary)", {
  cval_bal <- balance_partitions(cval_kfold_bin, des_bin, method = "subsample", seed = 2)
  expect_s3_class(cval_bal, "custom_cross_validation")
  samples_bal <- crossval_samples(cval_bal, design_df_bin, y_bin)
  for (i in 1:nrow(samples_bal)) {
    expect_true(is_balanced(samples_bal$ytrain[[i]]))
    expect_true(is_balanced(samples_bal$ytest[[i]]))
  }
})

test_that("balance_partitions subsample works for twofold_cv (binary)", {
  cval_bal <- balance_partitions(cval_twofold_bin, des_bin, method = "subsample", seed = 3)
  expect_s3_class(cval_bal, "custom_cross_validation")
  samples_bal <- crossval_samples(cval_bal, design_df_bin, y_bin)
  expect_equal(nrow(samples_bal), cval_twofold_bin$nreps)
  for (i in 1:nrow(samples_bal)) {
    expect_true(is_balanced(samples_bal$ytrain[[i]]))
    expect_true(is_balanced(samples_bal$ytest[[i]]))
  }
})

test_that("balance_partitions subsample works for multi-class", {
  cval_bal <- balance_partitions(cval_blocked_multi, des_multi, method = "subsample", seed = 4)
  expect_s3_class(cval_bal, "custom_cross_validation")
  samples_bal <- crossval_samples(cval_bal, design_df_multi, y_multi)
  for (i in 1:nrow(samples_bal)) {
    expect_true(is_balanced(samples_bal$ytrain[[i]]))
    expect_true(is_balanced(samples_bal$ytest[[i]]))
  }
})

test_that("balance_partitions subsample respects balance_test_set = FALSE", {
  cval_bal <- balance_partitions(cval_blocked_bin, des_bin, method = "subsample", balance_test_set = FALSE, seed = 5)
  expect_s3_class(cval_bal, "custom_cross_validation")
  samples_bal <- crossval_samples(cval_bal, design_df_bin, y_bin)
  samples_unbal <- crossval_samples(cval_blocked_bin, design_df_bin, y_bin)
  
  test_sets_balanced = TRUE
  for (i in 1:nrow(samples_bal)) {
    expect_true(is_balanced(samples_bal$ytrain[[i]])) # Train should be balanced
    # Test should be *identical* to original unbalanced test set
    expect_equal(sort(samples_bal$test[[i]]$idx), sort(samples_unbal$test[[i]]$idx))
    if (!is_balanced(samples_bal$ytest[[i]])) { 
      test_sets_balanced = FALSE
    }
  }
  expect_false(test_sets_balanced) # Ensure at least one test set was not balanced
})

# --- Test Oversampling --- 

test_that("balance_partitions oversample works for blocked_cv (binary)", {
  expect_warning(
    cval_bal <- balance_partitions(cval_blocked_bin, des_bin, method = "oversample", seed = 6),
    "Oversampling the test set"
  )
  expect_s3_class(cval_bal, "custom_cross_validation")
  
  samples_bal <- crossval_samples(cval_bal, design_df_bin, y_bin)
  
  for (i in 1:nrow(samples_bal)) {
    ytrain_bal <- samples_bal$ytrain[[i]]
    expect_true(is_balanced(ytrain_bal))
    ytest_bal <- samples_bal$ytest[[i]]
    expect_true(is_balanced(ytest_bal))
    
    # Check total size (should be larger)
    samples_unbal <- crossval_samples(cval_blocked_bin, design_df_bin, y_bin)
    expect_gt(length(ytrain_bal), length(samples_unbal$ytrain[[i]]))
    # Cannot guarantee test set is larger if a block was already balanced 
  }
})

test_that("balance_partitions oversample works for multi-class", {
   expect_warning(
      cval_bal <- balance_partitions(cval_blocked_multi, des_multi, method = "oversample", seed = 7),
     "Oversampling the test set"
   )
  expect_s3_class(cval_bal, "custom_cross_validation")
  samples_bal <- crossval_samples(cval_bal, design_df_multi, y_multi)
  for (i in 1:nrow(samples_bal)) {
    expect_true(is_balanced(samples_bal$ytrain[[i]]))
    expect_true(is_balanced(samples_bal$ytest[[i]]))
  }
})

test_that("balance_partitions oversample respects balance_test_set = FALSE", {
  # No warning expected here
  cval_bal <- balance_partitions(cval_blocked_bin, des_bin, method = "oversample", balance_test_set = FALSE, seed = 8)
  expect_s3_class(cval_bal, "custom_cross_validation")
  samples_bal <- crossval_samples(cval_bal, design_df_bin, y_bin)
  samples_unbal <- crossval_samples(cval_blocked_bin, design_df_bin, y_bin)
  
  test_sets_balanced = TRUE
  for (i in 1:nrow(samples_bal)) {
    expect_true(is_balanced(samples_bal$ytrain[[i]])) # Train should be balanced
    # Test should be *identical* to original unbalanced test set
    expect_equal(sort(samples_bal$test[[i]]$idx), sort(samples_unbal$test[[i]]$idx))
     if (!is_balanced(samples_bal$ytest[[i]])) { 
      test_sets_balanced = FALSE
    }
  }
   expect_false(test_sets_balanced) # Ensure at least one test set was not balanced
})

# --- Test Edge Cases & Options --- 

test_that("balance_partitions handles folds with only one class", {
  # Create a scenario where one block (fold) has only one class ('b')
  design_df_single <- data.frame(
      condition = factor(c(rep("a", 5), rep("b", 25))), # 5 'a' in Blk 1, 25 'b' in Blk 2+3
      block = c(rep(1, 5), rep(2, 15), rep(3, 10)) # Block 3 has only 'b'
  )
  des_single <- mvpa_design(design_df_single, y_train = ~ condition, block_var = ~ block)
  y_single <- des_single$y_train
  cval_blocked_single <- blocked_cross_validation(des_single$block_var)
  samples_unbal <- crossval_samples(cval_blocked_single, design_df_single, y_single)
  
  # Run balancing - expect warnings but don't test for specific text
  suppressWarnings({
    cval_bal <- balance_partitions(cval_blocked_single, des_single, method = "subsample", seed = 9)
  })
  
  samples_bal <- crossval_samples(cval_bal, design_df_single, y_single)
  
  # Fold 1: Test=Blk1('a' only), Train=Blk2+3('b' only)
  expect_false(is_balanced(samples_bal$ytrain[[1]])) # Train cannot be balanced (only 'b')
  expect_false(is_balanced(samples_bal$ytest[[1]]))  # Test cannot be balanced (only 'a')
  expect_equal(sort(samples_bal$test[[1]]$idx), sort(samples_unbal$test[[1]]$idx)) # Test unchanged
  
  # Fold 2: Test=Blk2('b' only), Train=Blk1('a')+Blk3('b')
  expect_true(is_balanced(samples_bal$ytrain[[2]]))  # Train *can* be balanced (has 'a' and 'b')
  expect_false(is_balanced(samples_bal$ytest[[2]])) # Test cannot be balanced (only 'b')
  expect_equal(sort(samples_bal$test[[2]]$idx), sort(samples_unbal$test[[2]]$idx)) # Test unchanged

  # Fold 3: Test=Blk3('b' only), Train=Blk1('a')+Blk2('b')
  expect_true(is_balanced(samples_bal$ytrain[[3]])) # Train *can* be balanced (has 'a' and 'b')
  expect_false(is_balanced(samples_bal$ytest[[3]])) # Test cannot be balanced (only 'b')
  expect_equal(sort(samples_bal$test[[3]]$idx), sort(samples_unbal$test[[3]]$idx)) # Test unchanged
  
})

test_that("balance_partitions gives reproducible results with seed", {
  cval_bal1 <- balance_partitions(cval_blocked_bin, des_bin, method = "subsample", seed = 123)
  cval_bal2 <- balance_partitions(cval_blocked_bin, des_bin, method = "subsample", seed = 123)
  cval_bal3 <- balance_partitions(cval_blocked_bin, des_bin, method = "subsample", seed = 456)
  
  expect_identical(cval_bal1$sample_set, cval_bal2$sample_set)
  expect_false(identical(cval_bal1$sample_set, cval_bal3$sample_set))
  
  # Also test oversampling seed
   expect_warning(cval_bal_over1 <- balance_partitions(cval_blocked_bin, des_bin, method = "oversample", seed = 789))
   expect_warning(cval_bal_over2 <- balance_partitions(cval_blocked_bin, des_bin, method = "oversample", seed = 789))
   expect_warning(cval_bal_over3 <- balance_partitions(cval_blocked_bin, des_bin, method = "oversample", seed = 101))
   
  expect_identical(cval_bal_over1$sample_set, cval_bal_over2$sample_set)
  expect_false(identical(cval_bal_over1$sample_set, cval_bal_over3$sample_set))
})

# --- Test Error Handling --- 

test_that("balance_partitions throws error for invalid method", {
  expect_error(
    balance_partitions(cval_blocked_bin, des_bin, method = "invalid_method"),
    "`method` must be one of: subsample, oversample"
  )
})

test_that("balance_partitions throws error for non-mvpa_design object", {
  expect_error(
    balance_partitions(cval_blocked_bin, design_df_bin), # Pass dataframe instead of design
    "`design` argument must be an 'mvpa_design' object."
  )
})

# Test balancing an already custom_cross_validation object
test_that("balance_partitions can re-balance a custom_cross_validation object", {
  # Create an unbalanced custom CV
  samples_unbal <- crossval_samples(cval_blocked_bin, design_df_bin, y_bin)
  custom_unbal_set <- lapply(1:nrow(samples_unbal), function(i) {
      list(train = samples_unbal$train[[i]]$idx, test = samples_unbal$test[[i]]$idx)
  })
  cval_custom_unbal <- custom_cross_validation(custom_unbal_set)

  # Balance it
  cval_custom_bal <- balance_partitions(cval_custom_unbal, des_bin, method = "subsample", seed = 10)
  expect_s3_class(cval_custom_bal, "custom_cross_validation")
  samples_custom_bal <- crossval_samples(cval_custom_bal, design_df_bin, y_bin)

  # Check if balanced
  for (i in 1:nrow(samples_custom_bal)) {
      expect_true(is_balanced(samples_custom_bal$ytrain[[i]]))
      expect_true(is_balanced(samples_custom_bal$ytest[[i]]))
  }
}) 