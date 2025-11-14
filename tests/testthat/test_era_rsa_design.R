context("era_rsa_design helper")

test_that("era_rsa_design builds item summaries and confound RDMs from a two-phase table", {
  # Synthetic two-phase (enc/ret) table in a single train_design
  set.seed(1)
  K <- 6
  items <- letters[1:K]
  Phase <- rep(c("enc", "ret"), each = K)
  Item  <- rep(items, times = 2)
  # Blocks/runs differ between items
  Block_enc <- sample(1:2, K, replace = TRUE)
  Block_ret <- sample(1:2, K, replace = TRUE)
  Block <- c(Block_enc, Block_ret)
  # Times: enc earlier, ret later with constant lag 10
  Time_enc <- seq_len(K)
  Time_ret <- Time_enc + 10
  Time <- c(Time_enc, Time_ret)

  train_df <- data.frame(Phase = factor(Phase, levels = c("enc", "ret")),
                         Item  = Item,
                         Block = Block,
                         Time  = Time,
                         stringsAsFactors = FALSE)

  # Minimal mvpa_design using only train_design
  des <- mvpa_design(train_design = train_df, y_train = ~ Phase, block_var = ~ Block)

  out <- era_rsa_design(design = des,
                        key_var = ~ Item,
                        phase_var = ~ Phase,
                        encoding_level = "enc",
                        retrieval_level = "ret",
                        block_var = ~ Block,
                        time_var  = ~ Time)

  expect_true(is.list(out))
  expect_true(all(c("items", "confound_rdms") %in% names(out)))
  expect_equal(length(out$items), K)

  # Check item-level summaries
  expect_true(is.factor(out$item_block))
  expect_equal(names(out$item_block), items)
  expect_equal(names(out$item_time_enc), items)
  expect_equal(names(out$item_time_ret), items)
  expect_equal(unname(out$item_time_enc), Time_enc)
  expect_equal(unname(out$item_time_ret), Time_ret)
  expect_equal(unname(out$item_lag), rep(10, K))

  # Confounds: block (matrix) and time_enc (dist)
  cr <- out$confound_rdms
  expect_true(all(c("block", "time_enc", "run_enc", "run_ret") %in% names(cr)))
  expect_true(is.matrix(cr$block))
  expect_equal(dim(cr$block), c(K, K))
  expect_true(inherits(cr$time_enc, "dist"))
  expect_equal(attr(cr$time_enc, "Size"), K)
  expect_true(is.matrix(cr$run_enc) && is.matrix(cr$run_ret))
  expect_equal(dim(cr$run_enc), c(K, K))
  expect_equal(dim(cr$run_ret), c(K, K))
})


test_that("era_rsa_design outputs can be passed into era_rsa_model", {
  skip_on_cran()
  set.seed(2)

  # External-test toy dataset for model-level smoke test
  toy <- gen_sample_dataset(D = c(3,3,3), nobs = 30, nlevels = 3, blocks = 3, external_test = TRUE)
  # Inject a shared key into both designs
  toy$design$train_design$Item <- toy$design$train_design$Y
  toy$design$test_design$Item  <- toy$design$test_design$Ytest

  # Build helper outputs from a two-phase-like train table (reuse train_design; fake Phase to satisfy helper)
  td <- toy$design$train_design
  # Duplicate enc/ret rows by reusing train_design for both and stacking (small K ensures speed)
  # For this smoke test, simply assign a dummy Phase and Time; the helper will still return a list
  td$Phase <- factor(rep("enc", nrow(td)))
  td$Time  <- seq_len(nrow(td))
  des2 <- toy$design
  des2$train_design <- td

  helper <- era_rsa_design(design = des2,
                           key_var = ~ Item,
                           phase_var = ~ Phase,
                           encoding_level = "enc",
                           retrieval_level = "enc",
                           block_var = ~ block_var,
                           time_var  = ~ Time)

  ms <- era_rsa_model(dataset = toy$dataset,
                      design  = toy$design,
                      key_var = ~ Item,
                      phase_var = ~ block_var,
                      confound_rdms = helper$confound_rdms,
                      item_block    = helper$item_block,
                      item_lag      = helper$item_lag)

  # Simple two-region mask
  mask <- neuroim2::NeuroVol(sample(1:2, size = length(toy$dataset$mask), replace = TRUE),
                             neuroim2::space(toy$dataset$mask))
  res <- run_regional(ms, mask)
  expect_s3_class(res, "regional_mvpa_result")
})
