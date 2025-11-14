library(testthat)
library(rMVPA)

context("hrfdecoder defensive handling")

make_mock_dataset <- function(Tlen = 120, V = 10) {
  train_mat <- matrix(rnorm(Tlen * V), Tlen, V)
  class(train_mat) <- "NeuroVec"
  mask <- array(rep(1, V), dim = c(V, 1, 1))
  class(mask) <- "NeuroVol"
  mvpa_dataset(train_mat, mask = mask)
}

make_mock_design <- function(Tlen = 120, n_runs = 2) {
  # Minimal mock event_model with sampling_frame attrs
  TR <- 2
  blocklens <- rep(Tlen / n_runs, n_runs)
  sframe <- list(blocklens = blocklens, TR = TR)
  attr(sframe, "blocklens") <- sframe$blocklens
  attr(sframe, "TR") <- sframe$TR
  class(sframe) <- "sampling_frame"

  evmod <- list()
  attr(evmod, "sampling_frame") <- sframe
  class(evmod) <- "event_model"

  # Simple events within total time
  events <- data.frame(
    onset = c(5, 15, 25),
    condition = factor(c("A", "B", "A")),
    run = c(1, 1, 2)
  )

  block_var <- rep(seq_len(n_runs), each = Tlen / n_runs)
  hrfdecoder_design(evmod, events, block_var)
}

test_that("invalid window resets to default with warning", {
  dset <- make_mock_dataset()
  design <- make_mock_design()
  expect_warning({
    mspec <- hrfdecoder_model(dataset = dset, design = design, window = c(NA_real_, NA_real_))
  }, "invalid 'window'; using default")
  expect_equal(mspec$window, c(4, 8))
})

test_that("single-block design warns about CV", {
  dset <- make_mock_dataset(Tlen = 60, V = 8)
  design <- make_mock_design(Tlen = 60, n_runs = 1)
  expect_warning({
    mspec <- hrfdecoder_model(dataset = dset, design = design)
  }, "only one run/block detected")
})

test_that("basis handling depends on fmrihrf availability", {
  dset <- make_mock_dataset()
  design <- make_mock_design()
  fake_basis <- list(name = "dummy_hrf")
  if (!requireNamespace("fmrihrf", quietly = TRUE)) {
    expect_warning({
      mspec <- hrfdecoder_model(dataset = dset, design = design, basis = fake_basis)
    }, "'fmrihrf' not available; ignoring provided basis")
    expect_null(mspec$basis)
  } else {
    mspec <- hrfdecoder_model(dataset = dset, design = design, basis = fake_basis)
    expect_identical(mspec$basis, fake_basis)
  }
})

test_that("train_model requires hrfdecoder when not installed", {
  skip_if(requireNamespace("hrfdecoder", quietly = TRUE), "hrfdecoder installed; skipping negative path test")
  dset <- make_mock_dataset()
  design <- make_mock_design()
  mspec <- hrfdecoder_model(dataset = dset, design = design)

  # Use a small ROI matrix as train data
  train_dat <- matrix(rnorm(60 * 5), 60, 5)
  expect_error(
    train_model(mspec, train_dat, y = 1:60, sl_info = NULL, cv_spec = NULL, indices = 1:5),
    regexp = "hrfdecoder_model: The 'hrfdecoder' package is required",
    fixed = FALSE
  )
})

test_that("format_result returns error tibble when hrfdecoder missing", {
  skip_if(requireNamespace("hrfdecoder", quietly = TRUE), "hrfdecoder installed; skipping negative path test")
  dset <- make_mock_dataset()
  design <- make_mock_design()
  mspec <- hrfdecoder_model(dataset = dset, design = design)

  # Build minimal context resembling internal_crossval context
  context <- list(
    test = 1:10,
    ytest = 1:10
  )
  out <- format_result(mspec, result = list(fit = NULL), error_message = NULL, context = context)
  expect_true(isTRUE(out$error))
  expect_match(out$error_message[[1]], "Missing required package 'hrfdecoder'")
})
