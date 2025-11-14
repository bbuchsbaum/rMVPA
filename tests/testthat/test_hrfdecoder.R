library(testthat)
library(rMVPA)

context("hrfdecoder integration")

# Skip all tests if hrfdecoder not available
skip_if_not_installed("hrfdecoder")
skip_if_not_installed("fmridesign")
skip_if_not_installed("fmrihrf")

test_that("hrfdecoder_design requires correct event table columns", {
  # Mock sampling frame
  sframe <- list(blocklens = c(60, 60), TR = 2)
  attr(sframe, "blocklens") <- sframe$blocklens
  attr(sframe, "TR") <- sframe$TR
  class(sframe) <- "sampling_frame"

  # Mock event model
  evmod <- list()
  attr(evmod, "sampling_frame") <- sframe
  class(evmod) <- "event_model"

  # Missing 'onset' column
  events_bad <- data.frame(
    condition = c("A", "B"),
    run = c(1, 1)
  )

  block_var <- rep(1:2, each = 60)

  expect_error(
    hrfdecoder_design(evmod, events_bad, block_var),
    "events table must contain columns.*onset"
  )

  # Missing 'condition' column
  events_bad2 <- data.frame(
    onset = c(5, 15),
    run = c(1, 1)
  )

  expect_error(
    hrfdecoder_design(evmod, events_bad2, block_var),
    "events table must contain columns.*condition"
  )
})

test_that("hrfdecoder_design validates block_var length", {
  sframe <- list(blocklens = c(60, 60), TR = 2)
  attr(sframe, "blocklens") <- sframe$blocklens
  attr(sframe, "TR") <- sframe$TR
  class(sframe) <- "sampling_frame"

  evmod <- list()
  attr(evmod, "sampling_frame") <- sframe
  class(evmod) <- "event_model"

  events <- data.frame(
    onset = c(5, 15),
    condition = c("A", "B"),
    run = c(1, 1)
  )

  # Wrong length block_var
  block_var_wrong <- rep(1:2, each = 50)  # 100 instead of 120

  expect_warning(
    hrfdecoder_design(evmod, events, block_var_wrong),
    "block_var length.*does not match sampling_frame"
  )
})

test_that("hrfdecoder_design warns about events outside acquisition time", {
  sframe <- list(blocklens = c(60, 60), TR = 2)
  attr(sframe, "blocklens") <- sframe$blocklens
  attr(sframe, "TR") <- sframe$TR
  class(sframe) <- "sampling_frame"

  evmod <- list()
  attr(evmod, "sampling_frame") <- sframe
  class(evmod) <- "event_model"

  # Events beyond 120 TRs * 2 sec = 240 seconds
  events <- data.frame(
    onset = c(5, 15, 250),  # Last event is outside
    condition = c("A", "B", "A"),
    run = c(1, 1, 2)
  )

  block_var <- rep(1:2, each = 60)

  expect_warning(
    hrfdecoder_design(evmod, events, block_var),
    "event\\(s\\) have onsets beyond total acquisition time"
  )
})

test_that("hrfdecoder_design creates proper class hierarchy", {
  sframe <- list(blocklens = c(60, 60), TR = 2)
  attr(sframe, "blocklens") <- sframe$blocklens
  attr(sframe, "TR") <- sframe$TR
  class(sframe) <- "sampling_frame"

  evmod <- list()
  attr(evmod, "sampling_frame") <- sframe
  class(evmod) <- "event_model"

  events <- data.frame(
    onset = c(5, 15, 25),
    condition = c("A", "B", "A"),
    run = c(1, 1, 2)
  )

  block_var <- rep(1:2, each = 60)

  design <- hrfdecoder_design(evmod, events, block_var)

  expect_true(inherits(design, "hrfdecoder_design"))
  expect_true(inherits(design, "mvpa_design"))
  expect_true(inherits(design, "list"))

  # Check fields
  expect_false(is.null(design$event_model))
  expect_false(is.null(design$events))
  expect_equal(nrow(design$events), 3)
})

test_that("print.hrfdecoder_design works without error", {
  sframe <- list(blocklens = c(60, 60), TR = 2)
  attr(sframe, "blocklens") <- sframe$blocklens
  attr(sframe, "TR") <- sframe$TR
  class(sframe) <- "sampling_frame"

  evmod <- list()
  attr(evmod, "sampling_frame") <- sframe
  class(evmod) <- "event_model"

  events <- data.frame(
    onset = c(5, 15, 25),
    condition = factor(c("A", "B", "A")),
    run = c(1, 1, 2)
  )

  block_var <- rep(1:2, each = 60)
  design <- hrfdecoder_design(evmod, events, block_var)

  # Should not error
  expect_output(print(design), "hrfdecoder_design")
  expect_output(print(design), "TRs: 120")
  expect_output(print(design), "Events: 3")
})

test_that("hrfdecoder_model requires hrfdecoder_design", {
  # Mock dataset
  train_mat <- matrix(rnorm(120 * 20), 120, 20)
  class(train_mat) <- "NeuroVec"
  mask <- array(rep(1, 20), dim = c(4, 5, 1))
  class(mask) <- "NeuroVol"

  dset <- mvpa_dataset(train_mat, mask = mask)

  # Regular mvpa_design (not hrfdecoder_design)
  regular_design <- mvpa_design(
    train_design = data.frame(y = 1:3, block = c(1,1,2)),
    y_train = ~ y,
    block_var = ~ block
  )

  expect_error(
    hrfdecoder_model(dset, regular_design),
    "design must be created with hrfdecoder_design"
  )
})

test_that("hrfdecoder_model creates proper model spec", {
  skip_if_not(requireNamespace("hrfdecoder", quietly = TRUE))

  # Mock dataset
  train_mat <- matrix(rnorm(120 * 20), 120, 20)
  class(train_mat) <- "NeuroVec"
  mask <- array(rep(1, 20), dim = c(4, 5, 1))
  class(mask) <- "NeuroVol"

  dset <- mvpa_dataset(train_mat, mask = mask)

  # Mock design
  sframe <- list(blocklens = c(60, 60), TR = 2)
  attr(sframe, "blocklens") <- sframe$blocklens
  attr(sframe, "TR") <- sframe$TR
  class(sframe) <- "sampling_frame"

  evmod <- list()
  attr(evmod, "sampling_frame") <- sframe
  class(evmod) <- "event_model"

  events <- data.frame(
    onset = c(5, 15, 25),
    condition = factor(c("A", "B", "A")),
    run = c(1, 1, 2)
  )

  block_var <- rep(1:2, each = 60)
  design <- hrfdecoder_design(evmod, events, block_var)

  mspec <- hrfdecoder_model(
    dataset = dset,
    design = design,
    lambda_W = 10,
    lambda_HRF = 1,
    lambda_smooth = 5
  )

  expect_true(inherits(mspec, "hrfdecoder_model"))
  expect_true(inherits(mspec, "model_spec"))
  expect_equal(mspec$lambda_W, 10)
  expect_equal(mspec$lambda_HRF, 1)
  expect_equal(mspec$lambda_smooth, 5)
  expect_equal(mspec$window, c(4, 8))  # default
})

test_that("y_train.hrfdecoder_model returns TR sequence", {
  skip_if_not(requireNamespace("hrfdecoder", quietly = TRUE))

  # Mock dataset with 120 TRs
  train_mat <- matrix(rnorm(120 * 20), 120, 20)
  class(train_mat) <- "NeuroVec"
  mask <- array(rep(1, 20), dim = c(4, 5, 1))
  class(mask) <- "NeuroVol"

  dset <- mvpa_dataset(train_mat, mask = mask)

  # Mock design
  sframe <- list(blocklens = c(60, 60), TR = 2)
  attr(sframe, "blocklens") <- sframe$blocklens
  attr(sframe, "TR") <- sframe$TR
  class(sframe) <- "sampling_frame"

  evmod <- list()
  attr(evmod, "sampling_frame") <- sframe
  class(evmod) <- "event_model"

  events <- data.frame(
    onset = c(5, 15, 25),
    condition = factor(c("A", "B", "A")),
    run = c(1, 1, 2)
  )

  block_var <- rep(1:2, each = 60)
  design <- hrfdecoder_design(evmod, events, block_var)

  mspec <- hrfdecoder_model(dataset = dset, design = design)

  ytrain <- y_train(mspec)

  expect_equal(length(ytrain), 120)  # Number of TRs
  expect_equal(ytrain, 1:120)        # Sequence 1:T
})
