# The vectorised column screens must agree exactly with the apply()-based
# reference they replaced.
ref_nonzero_var <- function(M) {
  ret <- apply(M, 2, sd, na.rm = TRUE) > 0
  ret[is.na(ret)] <- FALSE
  ret
}
ref_na_cols <- function(M) apply(M, 2, function(x) any(is.na(x)))

test_that("vectorised column screens match the apply() reference", {
  set.seed(4401)
  for (rep in 1:200) {
    n <- sample(1:12, 1)
    p <- sample(1:15, 1)
    M <- matrix(rnorm(n * p), n, p)
    for (j in seq_len(p)) {
      switch(sample(8, 1),
        M[, j] <- 3.25,                                    # constant
        M[sample(n, max(1, n %/% 2)), j] <- NA,             # partial NA
        M[, j] <- NA,                                        # all NA
        M[sample(n, 1), j] <- NaN,                           # NaN
        M[sample(n, 1), j] <- sample(c(Inf, -Inf), 1),       # infinite
        { M[, j] <- 1; M[sample(n, 1), j] <- 1 + 1e-12 },    # near-constant
        { M[, j] <- NA; M[sample(n, 1), j] <- 2 },           # single value
        NULL                                                 # untouched
      )
    }
    expect_identical(rMVPA:::nonzeroVarianceColumns2(M), ref_nonzero_var(M), info = rep)
    expect_identical(rMVPA:::na_cols(M), ref_na_cols(M), info = rep)
  }
})

test_that("vectorised column screens keep column names", {
  M <- matrix(c(1, 2, 3, 5, 5, 5, NA, 1, 2), 3, dimnames = list(NULL, c("a", "b", "c")))
  expect_identical(rMVPA:::nonzeroVarianceColumns2(M), ref_nonzero_var(M))
  expect_identical(rMVPA:::na_cols(M), ref_na_cols(M))
})
