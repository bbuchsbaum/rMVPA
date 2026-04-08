library(rMVPA)

test_that("rmvpa_api_lifecycle exposes expected lifecycle tiers", {
  reg <- rmvpa_api_lifecycle()

  expect_true(is.data.frame(reg))
  expect_true(all(c("symbol", "lifecycle", "notes") %in% names(reg)))
  expect_true(all(c("stable", "experimental", "developer") %in% reg$lifecycle))
})

test_that("rmvpa_stable_api is a subset of exports", {
  stable <- rmvpa_stable_api()
  exports <- getNamespaceExports("rMVPA")

  expect_true(length(stable) > 0L)
  expect_true(all(stable %in% exports))
  expect_true(all(c("mvpa_config", "build_analysis", "run_analysis") %in% stable))
})
