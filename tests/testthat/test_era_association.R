context("ERA-RSA item associations")

.reference_era_item_scores <- function(E,
                                       R,
                                       method = "pearson",
                                       min_voxels = 2L,
                                       item_block = NULL) {
  S <- suppressWarnings(stats::cor(
    t(E), t(R), method = method, use = "pairwise.complete.obs"
  ))
  finite_counts <- (is.finite(E) * 1L) %*% t(is.finite(R) * 1L)
  S[finite_counts < min_voxels] <- NA_real_
  matched <- diag(S)
  background <- vapply(seq_len(nrow(S)), function(i) {
    values <- S[-i, i]
    values <- values[is.finite(values)]
    if (length(values)) mean(values) else NA_real_
  }, numeric(1L))
  same_specificity <- diff_specificity <- rep(NA_real_, nrow(S))
  if (!is.null(item_block)) {
    block <- as.character(item_block)
    for (i in seq_len(nrow(S))) {
      if (is.na(block[[i]]) || !is.finite(matched[[i]])) next
      same_idx <- which(!is.na(block) & block == block[[i]] & seq_len(nrow(S)) != i)
      diff_idx <- which(!is.na(block) & block != block[[i]])
      same_values <- S[same_idx, i]
      diff_values <- S[diff_idx, i]
      same_values <- same_values[is.finite(same_values)]
      diff_values <- diff_values[is.finite(diff_values)]
      if (length(same_values)) {
        same_specificity[[i]] <- matched[[i]] - mean(same_values)
      }
      if (length(diff_values)) {
        diff_specificity[[i]] <- matched[[i]] - mean(diff_values)
      }
    }
  }
  list(
    matrix = S,
    matched = matched,
    background = background,
    specificity = matched - background,
    same_block_specificity = same_specificity,
    diff_block_specificity = diff_specificity
  )
}

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

test_that("minimum count controls require finite whole-number scalars", {
  fx <- .make_era_association_fixture(K = 6L)
  base_args <- list(
    dataset = fx$split_dataset,
    design = fx$split_design,
    key_var = ~ item,
    pairing = "one_to_one"
  )
  invalid <- list(2.5, Inf, -Inf, NaN, NA_real_, c(2, 3), "3", TRUE)

  for (value in invalid) {
    expect_error(
      suppressWarnings(do.call(
        era_rsa_model,
        c(base_args, list(era_min_voxels = value))
      )),
      "era_min_voxels.*single integer"
    )
    expect_error(
      suppressWarnings(do.call(
        era_rsa_model,
        c(base_args, list(era_min_complete = value))
      )),
      "era_min_complete.*single integer"
    )
  }
  expect_error(
    suppressWarnings(do.call(era_rsa_model, c(base_args, list(era_min_voxels = 1)))),
    "era_min_voxels.*single integer"
  )
  expect_error(
    suppressWarnings(do.call(era_rsa_model, c(base_args, list(era_min_complete = 2)))),
    "era_min_complete.*single integer"
  )

  expect_s3_class(suppressWarnings(do.call(
    era_rsa_model,
    c(base_args, list(era_min_voxels = 2, era_min_complete = 3))
  )), "era_rsa_model")
})

test_that("association formulas reject ambiguous or unsupported specifications", {
  fx <- .make_era_association_fixture(K = 9L)
  make_model <- function(...) suppressWarnings(era_rsa_model(
    dataset = fx$split_dataset,
    design = fx$split_design,
    key_var = ~ item,
    pairing = "one_to_one",
    ...
  ))

  expect_error(make_model(era_correlates = vividness ~ trial_order),
               "era_correlates.*right-hand-side")
  expect_error(make_model(era_association = vividness ~ trial_order,
                          era_effects = ~ trial_order),
               "era_association.*right-hand-side")
  expect_error(make_model(era_effects = ~ vividness),
               "era_effects.*requires.*era_association")
  expect_error(make_model(era_effects_block = list(ratings = ~ vividness)),
               "era_effects_block.*requires.*era_association")
  expect_error(make_model(era_association = ~ vividness),
               "era_association.*requires.*era_effects")
  expect_error(make_model(
    era_association = ~ vividness,
    era_effects_block = list(~ vividness)
  ), "named, non-empty list")
  expect_error(make_model(
    era_association = ~ vividness,
    era_effects_block = list(other = ~ trial_order)
  ), "terms must occur exactly.*trial_order")
  expect_error(make_model(era_association = ~ 0 + vividness,
                          era_effects = ~ vividness),
               "must include an intercept")
  expect_error(make_model(era_association = ~ vividness,
                          era_effects = ~ trial_order),
               "terms must occur exactly")
  expect_error(make_model(era_correlates = ~ retrieval_run),
               "must be numeric")
  expect_error(make_model(era_correlates = ~ absent_rating),
               "not found")

  three_level <- fx$split_design
  three_level$test_design$retrieval_run <- factor(
    rep(c("run_1", "run_2", "run_3"), length.out = 9L)
  )
  expect_error(
    suppressWarnings(era_rsa_model(
      dataset = fx$split_dataset,
      design = three_level,
      key_var = ~ item,
      pairing = "one_to_one",
      era_association = ~ retrieval_run,
      era_effects = ~ retrieval_run
    )),
    "one-degree-of-freedom"
  )
})

test_that("prepared repeated-item prototypes preserve item means and key order", {
  train_design <- data.frame(item = c("b", "a", "b", "c"))
  test_design <- data.frame(item = c("c", "a", "c", "b"))
  design <- list(train_design = train_design, test_design = test_design)
  pairing <- rMVPA:::.era_pairing_info(design, ~ item, pairing = "average")
  E_trials <- matrix(seq_len(12L), nrow = 4L)
  R_trials <- matrix(seq_len(12L) + 20, nrow = 4L)

  got <- rMVPA:::.era_item_prototypes(E_trials, R_trials, pairing)
  expected_E <- rbind(
    a = E_trials[2L, ],
    b = colMeans(E_trials[c(1L, 3L), , drop = FALSE]),
    c = E_trials[4L, ]
  )
  expected_R <- rbind(
    a = R_trials[2L, ],
    b = R_trials[4L, ],
    c = colMeans(R_trials[c(1L, 3L), , drop = FALSE])
  )

  expect_identical(got$keys, c("a", "b", "c"))
  expect_equal(got$E, expected_E)
  expect_equal(got$R, expected_R)
})

test_that("repeated-item associations use prototypes and require constant metadata", {
  set.seed(3400)
  dims <- c(2L, 2L, 2L)
  p <- prod(dims)
  enc_keys <- c("b", "a", "b", "c", "d", "d")
  ret_keys <- c("c", "a", "c", "b", "d", "d")
  E_trials <- matrix(stats::rnorm(length(enc_keys) * p), nrow = length(enc_keys))
  R_trials <- matrix(stats::rnorm(length(ret_keys) * p), nrow = length(ret_keys))
  ratings <- c(7, 2, 7, 5, 9, 9)
  mask <- neuroim2::NeuroVol(array(1, dims), neuroim2::NeuroSpace(dims))
  enc <- neuroim2::NeuroVec(
    array(as.numeric(t(E_trials)), c(dims, length(enc_keys))),
    neuroim2::NeuroSpace(c(dims, length(enc_keys)))
  )
  ret <- neuroim2::NeuroVec(
    array(as.numeric(t(R_trials)), c(dims, length(ret_keys))),
    neuroim2::NeuroSpace(c(dims, length(ret_keys)))
  )
  train_tab <- data.frame(item = enc_keys)
  test_tab <- data.frame(item = ret_keys, vividness = ratings)
  dataset <- mvpa_dataset(enc, test_data = ret, mask = mask)
  design <- mvpa_design(
    train_tab,
    test_design = test_tab,
    y_train = ~ item,
    y_test = ~ item
  )
  model <- suppressWarnings(era_rsa_model(
    dataset,
    design,
    key_var = ~ item,
    pairing = "average",
    era_components = "item",
    era_correlates = ~ vividness,
    era_association = ~ vividness,
    era_effects = ~ vividness,
    era_min_complete = 3L
  ))
  out <- fit_roi(
    model,
    list(train_data = E_trials, test_data = R_trials, indices = seq_len(p)),
    list(id = 1L)
  )

  pairing <- rMVPA:::.era_pairing_info(design, ~ item, pairing = "average")
  prototypes <- rMVPA:::.era_item_prototypes(E_trials, R_trials, pairing)
  ref <- .reference_era_item_scores(prototypes$E, prototypes$R)
  item_ratings <- ratings[match(prototypes$keys, ret_keys)]
  expected_cor <- stats::cor(ref$specificity, item_ratings, method = "spearman")
  expected_part_r <- stats::cor(ref$specificity, item_ratings)

  expect_false(out$error)
  expect_equal(out$metrics[["era_vividness_cor"]], expected_cor, tolerance = 1e-12)
  expect_equal(
    out$metrics[["era_assoc_part_r_vividness"]],
    expected_part_r,
    tolerance = 1e-12
  )
  expect_equal(out$metrics[["era_assoc_n"]], 4)

  nonconstant <- design
  nonconstant$test_design$vividness[3L] <- 8
  expect_error(
    suppressWarnings(era_rsa_model(
      dataset,
      nonconstant,
      key_var = ~ item,
      pairing = "average",
      era_correlates = ~ vividness
    )),
    "not constant within retrieval item"
  )
})

test_that("zero-order ERA correlation uses trial-specific nonmatch background", {
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
  ref <- .reference_era_item_scores(fx$E, fx$R, method = "spearman")
  expected <- stats::cor(
    ref$specificity, fx$vividness, method = "spearman", use = "complete.obs"
  )

  expect_false(out$error)
  expect_equal(out$metrics[["era_diag_mean"]], mean(ref$matched), tolerance = 1e-12)
  expect_equal(
    out$metrics[["era_diag_minus_off"]],
    mean(ref$specificity),
    tolerance = 1e-12
  )
  expect_equal(out$metrics[["era_vividness_cor"]], expected, tolerance = 1e-12)
  expect_equal(
    out$metrics[["era_vividness_n"]],
    sum(is.finite(ref$specificity) & is.finite(fx$vividness))
  )
})

test_that("raw matched similarity remains an explicit association estimand", {
  fx <- .make_era_association_fixture(K = 10L)
  model <- suppressWarnings(era_rsa_model(
    dataset = fx$split_dataset,
    design = fx$split_design,
    key_var = ~ item,
    pairing = "one_to_one",
    era_correlates = ~ vividness,
    era_association_score = "matched"
  ))
  out <- fit_roi(
    model,
    roi_data = list(train_data = fx$E, test_data = fx$R, indices = seq_len(ncol(fx$E))),
    context = list(id = 1L)
  )
  ref <- .reference_era_item_scores(fx$E, fx$R)
  expected <- stats::cor(
    ref$matched, fx$vividness, method = "spearman", use = "complete.obs"
  )

  expect_equal(out$metrics[["era_vividness_cor"]], expected, tolerance = 1e-12)
})

test_that("optimized item scores match the full matrix with finite and missing data", {
  set.seed(3401)
  E <- matrix(stats::rnorm(12L * 17L), nrow = 12L)
  R <- E + matrix(stats::rnorm(12L * 17L, sd = 0.6), nrow = 12L)

  for (method in c("pearson", "spearman")) {
    ref <- .reference_era_item_scores(E, R, method = method, min_voxels = 3L)
    got <- rMVPA:::.era_item_similarity_scores(
      E, R, method = method, min_voxels = 3L, need_matrix = FALSE
    )
    expect_null(got$matrix)
    expect_equal(got$matched, ref$matched, tolerance = 1e-12)
    expect_equal(got$background, ref$background, tolerance = 1e-12)
    expect_equal(got$specificity, ref$specificity, tolerance = 1e-12)

    got_matrix <- rMVPA:::.era_item_similarity_scores(
      E, R, method = method, min_voxels = 3L, need_matrix = TRUE
    )
    expect_equal(got_matrix$matrix, ref$matrix, tolerance = 1e-12)
  }

  block <- rep(c("run_1", "run_2", "run_3"), each = 4L)
  ref <- .reference_era_item_scores(E, R, method = "pearson", min_voxels = 3L)
  same_ref <- diff_ref <- rep(NA_real_, nrow(E))
  for (i in seq_len(nrow(E))) {
    same_idx <- which(block == block[[i]] & seq_len(nrow(E)) != i)
    diff_idx <- which(block != block[[i]])
    same_ref[[i]] <- ref$matched[[i]] - mean(ref$matrix[same_idx, i])
    diff_ref[[i]] <- ref$matched[[i]] - mean(ref$matrix[diff_idx, i])
  }
  got_block <- rMVPA:::.era_item_similarity_scores(
    E,
    R,
    method = "pearson",
    min_voxels = 3L,
    need_matrix = FALSE,
    item_block = block
  )
  expect_equal(got_block$same_block_specificity, same_ref, tolerance = 1e-12)
  expect_equal(got_block$diff_block_specificity, diff_ref, tolerance = 1e-12)

  E[cbind(c(1L, 1L, 4L, 7L), c(1L, 3L, 8L, 12L))] <- NA_real_
  R[cbind(c(2L, 4L, 4L, 9L), c(2L, 6L, 15L, 17L))] <- NA_real_
  ref <- .reference_era_item_scores(E, R, method = "pearson", min_voxels = 3L)
  got <- rMVPA:::.era_item_similarity_scores(
    E, R, method = "pearson", min_voxels = 3L, need_matrix = FALSE
  )
  expect_equal(got$matched, ref$matched, tolerance = 1e-12)
  expect_equal(got$background, ref$background, tolerance = 1e-12)
  expect_equal(got$specificity, ref$specificity, tolerance = 1e-12)
})

test_that("missing block labels are excluded from block-specific backgrounds", {
  set.seed(3402)
  E <- matrix(stats::rnorm(10L * 9L), nrow = 10L)
  R <- 0.6 * E + matrix(stats::rnorm(10L * 9L), nrow = 10L)
  block <- c("run_1", "run_1", "run_1", "run_2", NA,
             "run_1", "run_2", "run_2", "run_2", "run_2")

  assert_matches_reference <- function(E, R, method) {
    ref <- .reference_era_item_scores(
      E, R, method = method, min_voxels = 3L, item_block = block
    )
    got <- rMVPA:::.era_item_similarity_scores(
      E,
      R,
      method = method,
      min_voxels = 3L,
      item_block = block
    )
    expect_equal(got$matched, ref$matched, tolerance = 1e-12)
    expect_equal(got$background, ref$background, tolerance = 1e-12)
    expect_equal(
      got$same_block_specificity,
      ref$same_block_specificity,
      tolerance = 1e-12
    )
    expect_equal(
      got$diff_block_specificity,
      ref$diff_block_specificity,
      tolerance = 1e-12
    )
  }

  for (method in c("pearson", "spearman")) {
    assert_matches_reference(E, R, method)
  }

  E[2L, 3L] <- NA_real_
  R[8L, 6L] <- NA_real_
  for (method in c("pearson", "spearman")) {
    assert_matches_reference(E, R, method)
  }
})

test_that("optimized item scores satisfy randomized full-matrix oracles", {
  set.seed(3403)
  for (case in seq_len(24L)) {
    K <- sample(4:10, 1L)
    p <- sample(4:14, 1L)
    E <- matrix(stats::rnorm(K * p), nrow = K)
    R <- 0.4 * E + matrix(stats::rnorm(K * p), nrow = K)
    if (case %% 3L == 0L) {
      E <- round(E, 1L)
      R <- round(R, 1L)
    }
    if (case %% 2L == 0L) {
      E[sample(length(E), 1L)] <- NA_real_
      R[sample(length(R), 1L)] <- NA_real_
    }
    block <- rep(c("run_1", "run_2", "run_3"), length.out = K)
    if (case %% 4L == 0L) block[sample(K, 1L)] <- NA_character_

    for (method in c("pearson", "spearman")) {
      ref <- .reference_era_item_scores(
        E, R, method = method, min_voxels = 3L, item_block = block
      )
      got <- rMVPA:::.era_item_similarity_scores(
        E,
        R,
        method = method,
        min_voxels = 3L,
        need_matrix = TRUE,
        item_block = block
      )
      expect_equal(got$matrix, ref$matrix, tolerance = 5e-12,
                   info = paste("matrix case", case, method))
      expect_equal(got$matched, ref$matched, tolerance = 5e-12,
                   info = paste("matched case", case, method))
      expect_equal(got$background, ref$background, tolerance = 5e-12,
                   info = paste("background case", case, method))
      expect_equal(got$specificity, ref$specificity, tolerance = 5e-12,
                   info = paste("specificity case", case, method))
      expect_equal(got$same_block_specificity, ref$same_block_specificity,
                   tolerance = 5e-12, info = paste("same-block case", case, method))
      expect_equal(got$diff_block_specificity, ref$diff_block_specificity,
                   tolerance = 5e-12, info = paste("different-block case", case, method))
    }
  }
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

test_that("adjusted ERA effect is signed part-r with correct affine behavior", {
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

  sim <- .reference_era_item_scores(fx$E, fx$R)$specificity
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

  shifted_design <- fx$split_design
  shifted_design$test_design$vividness <- shifted_design$test_design$vividness + 100
  shifted_out <- fit_roi(
    make_model(shifted_design),
    roi_data = list(train_data = fx$E, test_data = fx$R, indices = seq_len(ncol(fx$E))),
    context = list(id = 1L)
  )
  expect_equal(
    shifted_out$metrics[["era_assoc_part_r_vividness"]],
    out$metrics[["era_assoc_part_r_vividness"]],
    tolerance = 1e-12
  )

  reversed_design <- fx$split_design
  reversed_design$test_design$vividness <- -reversed_design$test_design$vividness
  reversed_out <- fit_roi(
    make_model(reversed_design),
    roi_data = list(train_data = fx$E, test_data = fx$R, indices = seq_len(ncol(fx$E))),
    context = list(id = 1L)
  )
  expect_equal(
    reversed_out$metrics[["era_assoc_part_r_vividness"]],
    -out$metrics[["era_assoc_part_r_vividness"]],
    tolerance = 1e-12
  )
})

test_that("signed part-r agrees with incremental R-squared across randomized designs", {
  set.seed(3404)
  for (case in seq_len(24L)) {
    n <- sample(10:25, 1L)
    keys <- sprintf("item_%02d", seq_len(n))
    item_data <- data.frame(
      .era_key = keys,
      focal_1 = stats::rnorm(n),
      focal_2 = stats::rnorm(n),
      nuisance = stats::rnorm(n),
      run = factor(rep(c("a", "b"), length.out = n))
    )
    similarity <- 0.45 * item_data$focal_1 - 0.3 * item_data$focal_2 +
      0.2 * item_data$nuisance + stats::rnorm(n)
    if (case %% 2L == 0L) item_data$focal_1[sample(n, 1L)] <- NA_real_
    if (case %% 3L == 0L) item_data$nuisance[sample(n, 1L)] <- NA_real_
    if (case %% 4L == 0L) similarity[sample(n, 1L)] <- NA_real_

    spec <- rMVPA:::.era_prepare_association(
      ~ focal_1 + focal_2 + nuisance + run,
      ~ focal_1 + focal_2,
      item_data
    )
    got <- rMVPA:::.era_association_metrics(
      similarity, keys, spec, min_complete = 4L
    )
    X <- spec$matrix
    keep <- is.finite(similarity) & apply(X, 1L, function(row) all(is.finite(row)))
    Xk <- X[keep, , drop = FALSE]
    yk <- similarity[keep]
    full <- stats::lm.fit(Xk, yk)
    sst <- sum((yk - mean(yk))^2)

    for (i in seq_along(spec$effect_columns)) {
      col <- spec$effect_columns[[i]]
      reduced <- stats::lm.fit(Xk[, -col, drop = FALSE], yk)
      delta_r2 <- max(0, sum(reduced$residuals^2) - sum(full$residuals^2)) / sst
      expected <- sign(full$coefficients[[col]]) * sqrt(delta_r2)
      metric <- paste0("era_assoc_part_r_", spec$effect_metric_labels[[i]])
      expect_equal(got[[metric]], expected, tolerance = 2e-12,
                   info = paste("case", case, spec$effect_labels[[i]]))
    }
    expect_equal(got[["era_assoc_n"]], sum(keep), info = paste("n case", case))
    expect_equal(got[["era_assoc_df_resid"]], sum(keep) - qr(Xk)$rank,
                 info = paste("df case", case))
  }
})

test_that("association blocks match nested LM tests on a shared complete-case set", {
  set.seed(6901)
  n <- 36L
  keys <- sprintf("item_%02d", seq_len(n))
  item_data <- data.frame(
    .era_key = keys,
    focal_1 = stats::rnorm(n),
    focal_2 = stats::rnorm(n),
    nuisance = stats::rnorm(n),
    group = factor(rep(c("a", "b", "c"), each = n / 3L))
  )
  similarity <- with(
    item_data,
    0.7 * focal_1 - 0.4 * focal_2 + 0.2 * nuisance +
      c(0, 0.5, -0.3)[group] + stats::rnorm(n, sd = 0.6)
  )
  item_data$focal_2[[4L]] <- NA_real_
  item_data$nuisance[[11L]] <- NA_real_
  similarity[[19L]] <- NA_real_

  blocks <- rMVPA:::.era_effect_blocks(list(
    signal = ~ focal_1 + focal_2,
    category = ~ group
  ))
  spec <- rMVPA:::.era_prepare_association(
    ~ focal_1 + focal_2 + nuisance + group,
    ~ focal_1,
    item_data,
    blocks
  )
  got <- rMVPA:::.era_association_metrics(
    similarity, keys, spec, min_complete = 4L
  )

  X <- spec$matrix
  keep <- is.finite(similarity) & apply(X, 1L, function(row) all(is.finite(row)))
  Xk <- X[keep, , drop = FALSE]
  yk <- similarity[keep]
  full <- stats::lm.fit(Xk, yk)
  tss <- sum((yk - mean(yk))^2)

  for (i in seq_along(spec$block_columns)) {
    label <- spec$block_metric_labels[[i]]
    reduced <- stats::lm.fit(Xk[, -spec$block_columns[[i]], drop = FALSE], yk)
    delta_rss <- sum(reduced$residuals^2) - sum(full$residuals^2)
    expected_df1 <- full$rank - reduced$rank
    expected_dr2 <- delta_rss / tss
    expected_F <- (delta_rss / expected_df1) /
      (sum(full$residuals^2) / (sum(keep) - full$rank))

    expect_equal(got[[paste0("era_assoc_dr2_", label)]], expected_dr2,
                 tolerance = 2e-12)
    expect_equal(got[[paste0("era_assoc_F_", label)]], expected_F,
                 tolerance = 2e-12)
    expect_equal(got[[paste0("era_assoc_df1_", label)]], expected_df1)
  }
  expect_equal(got[["era_assoc_n"]], n - 3L)
  expect_equal(got[["era_assoc_df1_category"]], 2)
})

test_that("association blocks report aliased rank contributions without failing", {
  set.seed(6902)
  n <- 24L
  x <- stats::rnorm(n)
  item_data <- data.frame(
    .era_key = sprintf("item_%02d", seq_len(n)),
    x = x,
    duplicate = x,
    nuisance = stats::rnorm(n)
  )
  similarity <- 0.5 * x + stats::rnorm(n)
  blocks <- rMVPA:::.era_effect_blocks(list(aliased = ~ duplicate))
  spec <- rMVPA:::.era_prepare_association(
    ~ x + duplicate + nuisance, NULL, item_data, blocks
  )

  got <- rMVPA:::.era_association_metrics(
    similarity, item_data$.era_key, spec, min_complete = 4L
  )

  expect_equal(got[["era_assoc_df1_aliased"]], 0)
  expect_equal(got[["era_assoc_dr2_aliased"]], 0, tolerance = 1e-14)
  expect_true(is.na(got[["era_assoc_F_aliased"]]))
})

test_that("association block metrics appear in the stable output schema", {
  fx <- .make_era_association_fixture(K = 9L)
  model <- suppressWarnings(era_rsa_model(
    dataset = fx$split_dataset,
    design = fx$split_design,
    key_var = ~ item,
    pairing = "one_to_one",
    era_components = "item",
    era_association = ~ vividness + retrieval_run + trial_order,
    era_effects_block = list(context = ~ retrieval_run + trial_order)
  ))

  nms <- names(output_schema(model))
  expect_true(all(c(
    "era_assoc_dr2_context",
    "era_assoc_F_context",
    "era_assoc_df1_context",
    "era_assoc_n",
    "era_assoc_df_resid"
  ) %in% nms))
  expect_false(any(grepl("era_assoc_part_r", nms)))

  out <- fit_roi(
    model,
    roi_data = list(
      train_data = fx$E,
      test_data = fx$R,
      indices = seq_len(ncol(fx$E))
    ),
    context = list(id = 1L)
  )
  expect_false(out$error)
  expect_true(all(c(
    "era_assoc_dr2_context",
    "era_assoc_F_context",
    "era_assoc_df1_context"
  ) %in% names(out$metrics)))
  expect_true(is.finite(out$metrics[["era_assoc_dr2_context"]]))
  expect_true(is.finite(out$metrics[["era_assoc_F_context"]]))
  expect_equal(out$metrics[["era_assoc_df1_context"]], 2)
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

test_that("item-only components omit identification and geometry metrics", {
  fx <- .make_era_association_fixture(K = 8L)
  model <- suppressWarnings(era_rsa_model(
    dataset = fx$split_dataset,
    design = fx$split_design,
    key_var = ~ item,
    pairing = "one_to_one",
    era_components = "item",
    era_correlates = ~ vividness,
    era_association = ~ vividness + retrieval_run + trial_order,
    era_effects = ~ vividness
  ))
  nms <- names(output_schema(model))

  expect_true(all(c(
    "n_items", "era_diag_mean", "era_diag_minus_off",
    "era_vividness_cor", "era_assoc_part_r_vividness"
  ) %in% nms))
  expect_false(any(c("era_top1_acc", "geom_cor", "geom_cor_partial") %in% nms))
  expect_error(
    suppressWarnings(era_rsa_model(
      fx$split_dataset,
      fx$split_design,
      key_var = ~ item,
      era_components = "geometry",
      era_correlates = ~ vividness
    )),
    "item.*component"
  )
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
    era_components = "item",
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

test_that("item-only ERA associations run through the shard searchlight backend", {
  skip_on_cran()
  skip_if_not_installed("shard")

  fx <- .make_era_association_fixture(K = 8L)
  model <- suppressWarnings(era_rsa_model(
    dataset = fx$split_dataset,
    design = fx$split_design,
    key_var = ~ item,
    pairing = "one_to_one",
    era_components = "item",
    era_correlates = ~ vividness,
    era_association = ~ vividness + retrieval_run + trial_order,
    era_effects = ~ vividness
  ))

  result <- run_searchlight(
    model,
    radius = 2,
    method = "standard",
    backend = "shard"
  )

  expect_s3_class(result, "searchlight_result")
  expect_true(all(c(
    "era_diag_minus_off",
    "era_vividness_cor",
    "era_assoc_part_r_vividness"
  ) %in% result$metrics))
  expect_false(any(c("era_top1_acc", "geom_cor") %in% result$metrics))
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
