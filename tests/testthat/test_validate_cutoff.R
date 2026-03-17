test_that("validate_cutoff is exported and normalizes top_k/top_p cutoffs", {
  expect_equal(validate_cutoff("top_k", 5, 20), 5)
  expect_equal(validate_cutoff("topk", 50, 20), 20)
  expect_equal(validate_cutoff("top_p", 0.1, 100), 10)
  expect_equal(validate_cutoff("topp", 0.01, 10), 1)
})

test_that("validate_cutoff errors on unsupported type and invalid percentage", {
  expect_error(validate_cutoff("bad_type", 1, 10), "Cutoff type")
  expect_error(validate_cutoff("top_p", 0, 10), "percentage cutoff")
  expect_error(validate_cutoff("top_p", 2, 10), "percentage cutoff")
})
