test_that("rsa_model detects dataset/design mismatch", {
  ds <- gen_sample_dataset(c(4,4,4), nobs = 10, blocks = 2)
  D <- dist(matrix(rnorm(8 * 8), 8, 8))
  rdes <- rsa_design(~ D, list(D = D))
  expect_error(
    rsa_model(ds$dataset, design = rdes),
    "Mismatch between dataset observations"
  )
})
