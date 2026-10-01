# Direct counts are an independent oracle for the sorted lookup. Include ties,
# endpoints, missing observations, infinities, empty bins and covariate edges.
test_that("sorted p-values preserve exact upper-tail counts and ties", {
  set.seed(1401)
  for (method in c("global", "adjusted")) {
    for (discrete in c(FALSE, TRUE)) {
      vals <- c(-Inf, if (discrete) sample(0:20, 500, TRUE) / 20 else rnorm(500),
                Inf, NA_real_, NaN)
      null_cov <- data.frame(nfeatures = sample(10:100, length(vals), TRUE))
      adj <- rMVPA:::build_adjusted_null(vals, null_cov, n_bins = 5, method = method)
      obs <- c(-Inf, -2, 0, 0.5, 1, 2, Inf, NA_real_, NaN, vals[2:30])
      cov <- data.frame(nfeatures = rep(c(0, adj$breaks, 200), length.out = length(obs)))
      bins <- if (method == "global") rep(1L, length(obs)) else {
        pmax(1L, pmin(findInterval(cov$nfeatures, adj$breaks, rightmost.closed = TRUE),
                     adj$n_bins))
      }
      expected <- vapply(seq_along(obs), function(i) {
        if (is.na(obs[i])) return(NA_real_)
        null <- adj$bin_nulls[[bins[i]]]
        (1 + sum(null >= obs[i])) / (1 + length(null))
      }, numeric(1))
      expect_identical(rMVPA:::score_observed(obs, adj, cov), expected)
    }
  }
})

test_that("sorted p-values handle empty nulls, empty bins and no observations", {
  empty <- rMVPA:::build_adjusted_null(numeric(0), data.frame(nfeatures = numeric(0)),
                                       method = "global")
  obs <- c(-Inf, 0, Inf, NA_real_, NaN)
  expect_identical(rMVPA:::score_observed(obs, empty, data.frame(nfeatures = rep(10, 5))),
                   c(1, 1, 1, NA_real_, NA_real_))
  expect_identical(rMVPA:::score_observed(numeric(0), empty,
                                        data.frame(nfeatures = numeric(0))), numeric(0))

  adj <- rMVPA:::build_adjusted_null(c(0, 1, 2, 3), data.frame(nfeatures = c(1, 2, 100, 101)),
                                    n_bins = 8, method = "adjusted")
  midpoints <- (head(adj$breaks, -1) + tail(adj$breaks, -1)) / 2
  expected <- vapply(adj$bin_nulls, function(null) (1 + sum(null >= 1)) / (1 + length(null)),
                     numeric(1))
  expect_true(any(lengths(adj$bin_nulls) == 0L))
  expect_identical(rMVPA:::score_observed(rep(1, adj$n_bins), adj,
                                        data.frame(nfeatures = midpoints)), expected)
})
