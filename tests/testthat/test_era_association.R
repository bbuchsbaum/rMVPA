context("ERA-RSA item associations")

.make_era_association_fixture <- function(K = 10L, dims = c(2L, 2L, 2L), seed = 41L) {
  set.seed(seed)
  keys <- sprintf("item_%02d", seq_len(K))
  p <- prod(dims)

  E <- matrix(stats::rnorm(K * p), nrow = K, dimnames = list(keys, NULL))
  noise <- seq(0.05, 0.8, length.out = K)
  R <- E + matrix(stats::rnorm(K * p), nrow = K) * noise

  combined <- matrix(NA_real_, nrow = 2L * K, ncol = p)
  combined[seq.int(1L, 2L * K, by = 2L), ] <- E
  combined[seq.int(2L, 2L * K, by = 2L), ] <- R

  arr <- array(as.numeric(t(combined)), dim = c(dims, 2L * K))
  vec <- neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(c(dims, 2L * K)))
  mask <- neuroim2::NeuroVol(
    array(1, dim = dims),
    neuroim2::NeuroSpace(dims)
  )

  vividness <- sample(seq_len(K))
  vividness[intersect(c(3L, 8L), seq_len(K))] <- NA_real_
  retrieval_run <- factor(rep(c("run_1", "run_2"), length.out = K))
  trial_order <- seq_len(K)

  tab <- data.frame(
    item = factor(rep(keys, each = 2L), levels = keys),
    phase = factor(rep(c("E", "R"), times = K), levels = c("E", "R")),
    vividness = as.numeric(rbind(rep(NA_real_, K), vividness)),
    retrieval_run = factor(as.vector(rbind(rep(NA_character_, K), as.character(retrieval_run))),
                           levels = levels(retrieval_run)),
    trial_order = as.numeric(rbind(rep(NA_real_, K), trial_order)),
    stringsAsFactors = FALSE
  )

  design <- mvpa_design(tab, y_train = ~ item)
  dataset <- mvpa_dataset(vec, mask = mask)

  enc_idx <- which(tab$phase == "E")
  ret_idx <- which(tab$phase == "R")
  enc_tab <- tab[enc_idx, , drop = FALSE]
  ret_tab <- tab[ret_idx, , drop = FALSE]
  split_design <- mvpa_design(
    enc_tab,
    test_design = ret_tab,
    y_train = ~ item,
    y_test = ~ item
  )
  split_dataset <- mvpa_dataset(
    neuroim2::sub_vector(vec, enc_idx),
    test_data = neuroim2::sub_vector(vec, ret_idx),
    mask = mask
  )

  list(
    keys = keys,
    E = E,
    R = R,
    table = tab,
    dataset = dataset,
    design = design,
    split_dataset = split_dataset,
    split_design = split_design,
    vividness = vividness,
    retrieval_run = retrieval_run,
    trial_order = trial_order
  )
}

test_that("combined alternating phases normalize to the external-test representation", {
  fx <- .make_era_association_fixture()

  combined <- suppressWarnings(era_rsa_model(
    dataset = fx$dataset,
    design = fx$design,
    key_var = ~ item,
    phase_var = ~ phase,
    encoding_level = "E",
    retrieval_level = "R",
    pairing = "one_to_one"
  ))
  external <- suppressWarnings(era_rsa_model(
    dataset = fx$split_dataset,
    design = fx$split_design,
    key_var = ~ item,
    pairing = "one_to_one"
  ))

  expect_true(has_test_set(combined$dataset))
  expect_equal(as.matrix(combined$dataset$train_data), as.matrix(external$dataset$train_data))
  expect_equal(as.matrix(combined$dataset$test_data), as.matrix(external$dataset$test_data))
  expect_identical(as.character(combined$design$train_design$item), fx$keys)
  expect_identical(as.character(combined$design$test_design$item), fx$keys)
})

test_that("strict one-to-one pairing rejects duplicate and unmatched keys", {
  fx <- .make_era_association_fixture(K = 6L)
  bad <- fx$design
  bad$train_design$item[1L] <- bad$train_design$item[3L]

  expect_error(
    suppressWarnings(era_rsa_model(
      fx$dataset, bad,
      key_var = ~ item,
      phase_var = ~ phase,
      encoding_level = "E",
      retrieval_level = "R",
      pairing = "one_to_one"
    )),
    "one-to-one"
  )
})

test_that("zero-order ERA correlation uses its own complete cases", {
  fx <- .make_era_association_fixture()
  model <- suppressWarnings(era_rsa_model(
    dataset = fx$split_dataset,
    design = fx$split_design,
    key_var = ~ item,
    pairing = "one_to_one",
    era_simfun = "spearman",
    era_correlates = ~ vividness,
    era_cor_method = "spearman",
    era_min_complete = 4L
  ))

  out <- fit_roi(
    model,
    roi_data = list(train_data = fx$E, test_data = fx$R, indices = seq_len(ncol(fx$E))),
    context = list(id = 1L)
  )
  S <- stats::cor(t(fx$E), t(fx$R), method = "spearman", use = "pairwise.complete.obs")
  expected <- stats::cor(diag(S), fx$vividness, method = "spearman", use = "complete.obs")

  expect_false(out$error)
  expect_equal(out$metrics[["era_diag_mean"]], mean(diag(S)), tolerance = 1e-12)
  expect_equal(out$metrics[["era_vividness_cor"]], expected, tolerance = 1e-12)
  expect_equal(out$metrics[["era_vividness_n"]], sum(is.finite(diag(S)) & is.finite(fx$vividness)))
})

test_that("item-key alignment makes associations invariant to phase row order", {
  fx <- .make_era_association_fixture(K = 10L)
  base_model <- suppressWarnings(era_rsa_model(
    fx$split_dataset,
    fx$split_design,
    key_var = ~ item,
    pairing = "one_to_one",
    era_correlates = ~ vividness,
    era_association = ~ vividness + retrieval_run + trial_order,
    era_effects = ~ vividness
  ))
  base <- fit_roi(
    base_model,
    list(train_data = fx$E, test_data = fx$R, indices = seq_len(ncol(fx$E))),
    list(id = 1L)
  )

  perm <- c(4L, 1L, 9L, 3L, 10L, 2L, 8L, 5L, 7L, 6L)
  perm_design <- fx$split_design
  perm_design$test_design <- perm_design$test_design[perm, , drop = FALSE]
  perm_design$y_test <- perm_design$y_test[perm]
  perm_dataset <- mvpa_dataset(
    fx$split_dataset$train_data,
    test_data = neuroim2::sub_vector(fx$split_dataset$test_data, perm),
    mask = fx$split_dataset$mask
  )
  perm_model <- suppressWarnings(era_rsa_model(
    perm_dataset,
    perm_design,
    key_var = ~ item,
    pairing = "one_to_one",
    era_correlates = ~ vividness,
    era_association = ~ vividness + retrieval_run + trial_order,
    era_effects = ~ vividness
  ))
  permuted <- fit_roi(
    perm_model,
    list(train_data = fx$E, test_data = fx$R[perm, , drop = FALSE], indices = seq_len(ncol(fx$E))),
    list(id = 1L)
  )

  metrics <- c("era_diag_mean", "era_vividness_cor", "era_assoc_part_r_vividness")
  expect_equal(base$metrics[metrics], permuted$metrics[metrics], tolerance = 1e-12)
})

test_that("adjusted ERA effect is signed part-r and invariant to predictor scale", {
  fx <- .make_era_association_fixture(K = 12L)
  make_model <- function(design) suppressWarnings(era_rsa_model(
    dataset = fx$split_dataset,
    design = design,
    key_var = ~ item,
    pairing = "one_to_one",
    era_simfun = "pearson",
    era_association = ~ vividness + retrieval_run + trial_order,
    era_effects = ~ vividness,
    era_min_complete = 4L
  ))

  model <- make_model(fx$split_design)
  out <- fit_roi(
    model,
    roi_data = list(train_data = fx$E, test_data = fx$R, indices = seq_len(ncol(fx$E))),
    context = list(id = 1L)
  )

  sim <- diag(stats::cor(t(fx$E), t(fx$R), use = "pairwise.complete.obs"))
  dat <- data.frame(
    sim = sim,
    vividness = fx$vividness,
    retrieval_run = fx$retrieval_run,
    trial_order = fx$trial_order
  )
  keep <- stats::complete.cases(dat)
  x_resid <- stats::resid(stats::lm(
    vividness ~ retrieval_run + trial_order,
    data = dat,
    subset = keep
  ))
  expected <- stats::cor(dat$sim[keep], x_resid)

  expect_equal(out$metrics[["era_assoc_part_r_vividness"]], expected, tolerance = 1e-12)
  expect_equal(out$metrics[["era_assoc_n"]], sum(keep))

  scaled_design <- fx$split_design
  scaled_design$test_design$vividness <- scaled_design$test_design$vividness * 100
  scaled <- make_model(scaled_design)
  scaled_out <- fit_roi(
    scaled,
    roi_data = list(train_data = fx$E, test_data = fx$R, indices = seq_len(ncol(fx$E))),
    context = list(id = 1L)
  )
  expect_equal(
    scaled_out$metrics[["era_assoc_part_r_vividness"]],
    out$metrics[["era_assoc_part_r_vividness"]],
    tolerance = 1e-12
  )
})

test_that("association schema reports only focal directional effects", {
  fx <- .make_era_association_fixture(K = 8L)
  model <- suppressWarnings(era_rsa_model(
    dataset = fx$split_dataset,
    design = fx$split_design,
    key_var = ~ item,
    pairing = "one_to_one",
    era_correlates = ~ vividness,
    era_association = ~ vividness + retrieval_run + trial_order,
    era_effects = ~ vividness
  ))
  nms <- names(output_schema(model))

  expect_true(all(c(
    "era_vividness_cor", "era_vividness_n",
    "era_assoc_part_r_vividness", "era_assoc_n", "era_assoc_df_resid"
  ) %in% nms))
  expect_false(any(grepl("beta|r2|retrieval_run|trial_order", nms)))
})

test_that("association metrics flow through regional and searchlight results", {
  skip_on_cran()
  fx <- .make_era_association_fixture(K = 8L)
  model <- suppressWarnings(era_rsa_model(
    dataset = fx$dataset,
    design = fx$design,
    key_var = ~ item,
    phase_var = ~ phase,
    encoding_level = "E",
    retrieval_level = "R",
    pairing = "one_to_one",
    era_correlates = ~ vividness,
    era_association = ~ vividness + retrieval_run + trial_order,
    era_effects = ~ vividness
  ))

  region_mask <- neuroim2::NeuroVol(
    array(1, dim = c(2L, 2L, 2L)),
    neuroim2::NeuroSpace(c(2L, 2L, 2L))
  )
  regional <- run_regional(model, region_mask)
  searchlight <- run_searchlight(model, radius = 3, method = "standard")

  expected <- c("era_vividness_cor", "era_vividness_n", "era_assoc_part_r_vividness")
  expect_true(all(expected %in% names(regional$performance_table)))
  expect_true(all(expected %in% searchlight$metrics))
})

test_that("association degeneracies return stable NA metrics", {
  fx <- .make_era_association_fixture(K = 8L)
  constant <- fx$split_design
  constant$test_design$vividness <- 1
  model <- suppressWarnings(era_rsa_model(
    dataset = fx$split_dataset,
    design = constant,
    key_var = ~ item,
    pairing = "one_to_one",
    era_correlates = ~ vividness,
    era_association = ~ vividness + retrieval_run + trial_order,
    era_effects = ~ vividness
  ))

  out <- fit_roi(
    model,
    roi_data = list(train_data = fx$E, test_data = fx$R, indices = seq_len(ncol(fx$E))),
    context = list(id = 1L)
  )
  expect_true(is.na(out$metrics[["era_vividness_cor"]]))
  expect_true(is.na(out$metrics[["era_assoc_part_r_vividness"]]))
  expect_true(is.finite(out$metrics[["era_assoc_n"]]))
})
