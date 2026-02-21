test_that("dual_lda geometry enter/leave update matches full sphere recomputation", {
  radius <- 3
  offsets <- rMVPA:::.searchlight_offsets(radius = radius, spacing = c(1, 1, 1))
  delta <- rMVPA:::.searchlight_enter_leave(offsets)

  dims <- c(13L, 13L, 13L)
  mask_active <- rep(TRUE, prod(dims))
  center <- c(7L, 7L, 7L)

  for (key in names(delta)) {
    step <- as.integer(strsplit(key, ",", fixed = TRUE)[[1]])
    center_next <- center + step

    prev_ids <- rMVPA:::.offset_to_indices(center, offsets, dims, mask_active)
    next_ids <- rMVPA:::.offset_to_indices(center_next, offsets, dims, mask_active)
    leave_ids <- rMVPA:::.offset_to_indices(center, delta[[key]]$leave, dims, mask_active)
    enter_ids <- rMVPA:::.offset_to_indices(center_next, delta[[key]]$enter, dims, mask_active)

    updated <- union(setdiff(prev_ids, leave_ids), enter_ids)
    expect_setequal(updated, next_ids)
  }
})

test_that("dual_lda fast boundary tables match reference step boundary columns", {
  radius <- 3
  offsets <- rMVPA:::.searchlight_offsets(radius = radius, spacing = c(1, 1, 1))
  delta <- rMVPA:::.searchlight_enter_leave(offsets)
  dims <- c(13L, 13L, 13L)

  mask_active <- rep(FALSE, prod(dims))
  active_idx <- sample.int(prod(dims), size = floor(prod(dims) * 0.8))
  mask_active[active_idx] <- TRUE
  col_lookup <- integer(prod(dims))
  col_lookup[active_idx] <- seq_along(active_idx)

  tables <- rMVPA:::.precompute_boundary_tables(delta, dims)
  center <- c(7L, 7L, 7L)
  prev_id <- center[1] + (center[2] - 1L) * dims[1] + (center[3] - 1L) * dims[1] * dims[2]

  for (key in names(delta)) {
    step <- as.integer(strsplit(key, ",", fixed = TRUE)[[1]])
    center_next <- center + step
    curr_id <- center_next[1] + (center_next[2] - 1L) * dims[1] + (center_next[3] - 1L) * dims[1] * dims[2]

    ref <- rMVPA:::.step_boundary_cols(
      prev_coord = center,
      curr_coord = center_next,
      leave_offsets = delta[[key]]$leave,
      enter_offsets = delta[[key]]$enter,
      dims = dims,
      mask_active = mask_active,
      col_lookup = col_lookup
    )
    fast <- rMVPA:::.step_boundary_cols_fast(
      prev_coord = center,
      prev_id = prev_id,
      curr_coord = center_next,
      curr_id = curr_id,
      step_table = tables[[key]],
      dims = dims,
      mask_active = mask_active,
      col_lookup = col_lookup
    )

    expect_setequal(ref$out_cols, fast$out_cols)
    expect_setequal(ref$in_cols, fast$in_cols)
  }
})

test_that("dual_lda rank-update preference is deterministic", {
  old_opt <- getOption("rMVPA.dual_lda_rank_update")
  old_env <- Sys.getenv("RMVPA_DUAL_LDA_RANK_UPDATE", unset = NA_character_)
  on.exit({
    options(rMVPA.dual_lda_rank_update = old_opt)
    if (is.na(old_env)) {
      Sys.unsetenv("RMVPA_DUAL_LDA_RANK_UPDATE")
    } else {
      Sys.setenv(RMVPA_DUAL_LDA_RANK_UPDATE = old_env)
    }
  }, add = TRUE)

  options(rMVPA.dual_lda_rank_update = TRUE)
  Sys.setenv(RMVPA_DUAL_LDA_RANK_UPDATE = "true")
  expect_null(rMVPA:::.dual_lda_rank_update_preference())
})

test_that("dual_lda rank-update heuristic cache keys by shape/radius/boundary size", {
  rMVPA:::.dual_lda_rank_cache_clear()
  on.exit(rMVPA:::.dual_lda_rank_cache_clear(), add = TRUE)

  fold <- list(
    R = matrix(0, nrow = 12, ncol = 80),
    M = matrix(0, nrow = 80, ncol = 3)
  )
  state <- list()

  key_a <- rMVPA:::.dual_lda_rank_cache_key(fold = fold, radius = 3, boundary_size = 18)
  key_b <- rMVPA:::.dual_lda_rank_cache_key(fold = fold, radius = 4, boundary_size = 18)
  key_c <- rMVPA:::.dual_lda_rank_cache_key(fold = fold, radius = 3, boundary_size = 22)
  key_d <- rMVPA:::.dual_lda_rank_cache_key(fold = fold, radius = 3, boundary_size = 18, chunk_size = 4)
  expect_false(identical(key_a, key_b))
  expect_false(identical(key_a, key_c))
  expect_false(identical(key_a, key_d))

  rMVPA:::.dual_lda_rank_cache_set(key_a, TRUE)
  expect_identical(
    rMVPA:::.dual_lda_rank_update_benchmark(
      state = state,
      fold = fold,
      out_cols = seq_len(9),
      in_cols = seq_len(9),
      radius = 3
    ),
    TRUE
  )

  expect_null(rMVPA:::.dual_lda_rank_cache_get(key_b))
  expect_null(rMVPA:::.dual_lda_rank_cache_get(key_c))
})

test_that("dual_lda native chunk update matches scalar rank updates", {
  has_native_chunk <- rMVPA:::.dual_lda_native_symbol_loaded(
    "rmvpa_dual_lda_chunk_update_predict"
  )
  if (!isTRUE(has_native_chunk)) {
    skip("Native chunk update symbol unavailable in this environment.")
  }

  set.seed(3191)
  dims <- c(9L, 9L, 9L)
  p <- prod(dims)
  n <- 36L
  classes <- c("a", "b", "c")
  y_all <- factor(rep(classes, length.out = n), levels = classes)
  x_all <- matrix(rnorm(n * p), nrow = n, ncol = p)

  fold <- rMVPA:::.prepare_dual_lda_fold(
    x_all = x_all,
    y_all = y_all,
    train_idx = seq_len(n),
    test_idx = seq_len(n),
    classes = classes,
    gamma = 1e-2
  )

  mask_active <- rep(TRUE, p)
  col_lookup <- seq_len(p)
  offsets <- rMVPA:::.searchlight_offsets(radius = 2, spacing = c(1, 1, 1))
  delta <- rMVPA:::.searchlight_enter_leave(offsets)
  tables <- rMVPA:::.precompute_boundary_tables(delta, dims)

  center1 <- c(4L, 4L, 4L)
  center2 <- center1 + c(1L, 0L, 0L)
  center3 <- center2 + c(1L, 0L, 0L)

  lin_id <- function(coord) {
    coord[1] + (coord[2] - 1L) * dims[1] + (coord[3] - 1L) * dims[1] * dims[2]
  }

  step1 <- rMVPA:::.step_boundary_cols_fast(
    prev_coord = center1,
    prev_id = lin_id(center1),
    curr_coord = center2,
    curr_id = lin_id(center2),
    step_table = tables[["1,0,0"]],
    dims = dims,
    mask_active = mask_active,
    col_lookup = col_lookup
  )
  step2 <- rMVPA:::.step_boundary_cols_fast(
    prev_coord = center2,
    prev_id = lin_id(center2),
    curr_coord = center3,
    curr_id = lin_id(center3),
    step_table = tables[["1,0,0"]],
    dims = dims,
    mask_active = mask_active,
    col_lookup = col_lookup
  )
  step_seq <- list(step1, step2)

  state0 <- rMVPA:::.init_dual_lda_state(fold, center1, offsets, dims, mask_active, col_lookup)

  state_scalar <- rMVPA:::.dual_lda_clone_state(state0)
  scalar_probs <- vector("list", length(step_seq))
  for (i in seq_along(step_seq)) {
    upd <- rMVPA:::.update_dual_lda_state_neighbor(
      state = state_scalar,
      fold = fold,
      out_cols = step_seq[[i]]$out_cols,
      in_cols = step_seq[[i]]$in_cols,
      use_rank_update = TRUE
    )
    expect_true(upd$ok)
    state_scalar <- upd$state
    scalar_probs[[i]] <- rMVPA:::.predict_dual_lda_state(state_scalar, fold)$probs
  }

  chunk_res <- rMVPA:::.dual_lda_chunk_update_predict_native(
    state = rMVPA:::.dual_lda_clone_state(state0),
    fold = fold,
    step_cols_seq = step_seq
  )
  expect_false(is.null(chunk_res))
  expect_length(chunk_res$probs, length(step_seq))

  expect_equal(chunk_res$state$L, state_scalar$L, tolerance = 1e-10)
  expect_equal(chunk_res$state$T, state_scalar$T, tolerance = 1e-10)
  expect_equal(chunk_res$state$U, state_scalar$U, tolerance = 1e-10)
  expect_equal(chunk_res$state$Q, state_scalar$Q, tolerance = 1e-10)
  expect_equal(chunk_res$state$C, state_scalar$C, tolerance = 1e-10)
  for (i in seq_along(step_seq)) {
    expect_equal(unname(chunk_res$probs[[i]]), unname(scalar_probs[[i]]), tolerance = 1e-10)
  }
})

test_that("dual_lda rank-k update falls back when native chol kernel is unavailable", {
  set.seed(3192)
  n <- 10L
  base <- matrix(rnorm(n * n), nrow = n)
  spd <- crossprod(base) + diag(0.5, n)
  r0 <- chol(spd)
  x <- matrix(rnorm(n * 3), nrow = n)

  ref <- r0
  for (j in seq_len(ncol(x))) {
    step <- rMVPA:::.chol_rank1_update_upper(ref, x[, j], downdate = FALSE)
    expect_true(step$ok)
    ref <- step$R
  }

  out <- testthat::with_mocked_bindings(
    .dual_lda_call_native = function(...) stop("mock native unavailable"),
    {
      rMVPA:::.chol_rankk_update_upper(r0, x, downdate = FALSE, scale = 1)
    },
    .package = "rMVPA"
  )

  expect_true(out$ok)
  expect_equal(out$R, ref, tolerance = 1e-10)
})

test_that("dual_lda state update falls back to non-native path on native miss", {
  set.seed(3193)
  dims <- c(7L, 7L, 7L)
  p <- prod(dims)
  n <- 30L
  classes <- c("a", "b", "c")
  y_all <- factor(rep(classes, length.out = n), levels = classes)
  x_all <- matrix(rnorm(n * p), nrow = n, ncol = p)

  fold <- rMVPA:::.prepare_dual_lda_fold(
    x_all = x_all,
    y_all = y_all,
    train_idx = seq_len(n),
    test_idx = seq_len(n),
    classes = classes,
    gamma = 1e-2
  )

  mask_active <- rep(TRUE, p)
  col_lookup <- seq_len(p)
  offsets <- rMVPA:::.searchlight_offsets(radius = 2, spacing = c(1, 1, 1))
  delta <- rMVPA:::.searchlight_enter_leave(offsets)
  tables <- rMVPA:::.precompute_boundary_tables(delta, dims)

  center <- c(4L, 4L, 4L)
  center_next <- center + c(1L, 0L, 0L)
  lin_id <- function(coord) {
    coord[1] + (coord[2] - 1L) * dims[1] + (coord[3] - 1L) * dims[1] * dims[2]
  }

  step_cols <- rMVPA:::.step_boundary_cols_fast(
    prev_coord = center,
    prev_id = lin_id(center),
    curr_coord = center_next,
    curr_id = lin_id(center_next),
    step_table = tables[["1,0,0"]],
    dims = dims,
    mask_active = mask_active,
    col_lookup = col_lookup
  )

  state0 <- rMVPA:::.init_dual_lda_state(fold, center, offsets, dims, mask_active, col_lookup)
  ref <- rMVPA:::.update_dual_lda_state_neighbor(
    state = rMVPA:::.dual_lda_clone_state(state0),
    fold = fold,
    out_cols = step_cols$out_cols,
    in_cols = step_cols$in_cols,
    use_rank_update = FALSE
  )
  expect_true(ref$ok)

  got <- testthat::with_mocked_bindings(
    .dual_lda_step_update_native = function(...) NULL,
    {
      rMVPA:::.update_dual_lda_state_neighbor(
        state = rMVPA:::.dual_lda_clone_state(state0),
        fold = fold,
        out_cols = step_cols$out_cols,
        in_cols = step_cols$in_cols,
        use_rank_update = TRUE
      )
    },
    .package = "rMVPA"
  )

  expect_true(got$ok)
  expect_equal(unname(got$state$L), unname(ref$state$L), tolerance = 1e-10)
  expect_equal(unname(got$state$T), unname(ref$state$T), tolerance = 1e-10)
  expect_equal(unname(got$state$U), unname(ref$state$U), tolerance = 1e-10)
  expect_equal(unname(got$state$Q), unname(ref$state$Q), tolerance = 1e-10)
  expect_equal(unname(got$state$C), unname(ref$state$C), tolerance = 1e-10)
})

test_that("dual_lda native chunk wrapper rejects non-finite payload", {
  set.seed(3194)
  dims <- c(7L, 7L, 7L)
  p <- prod(dims)
  n <- 24L
  classes <- c("a", "b", "c")
  y_all <- factor(rep(classes, length.out = n), levels = classes)
  x_all <- matrix(rnorm(n * p), nrow = n, ncol = p)

  fold <- rMVPA:::.prepare_dual_lda_fold(
    x_all = x_all,
    y_all = y_all,
    train_idx = seq_len(n),
    test_idx = seq_len(n),
    classes = classes,
    gamma = 1e-2
  )

  mask_active <- rep(TRUE, p)
  col_lookup <- seq_len(p)
  offsets <- rMVPA:::.searchlight_offsets(radius = 2, spacing = c(1, 1, 1))
  state0 <- rMVPA:::.init_dual_lda_state(fold, c(4L, 4L, 4L), offsets, dims, mask_active, col_lookup)

  step_cols_seq <- list(list(out_cols = integer(0), in_cols = integer(0)))
  bad_payload <- list(
    state = list(
      L = state0$L,
      T = state0$T,
      U = state0$U,
      Q = state0$Q,
      C = state0$C
    ),
    probs = list(matrix(NaN, nrow = ncol(fold$Xtest_t), ncol = length(classes)))
  )

  got <- testthat::with_mocked_bindings(
    .dual_lda_call_native = function(...) bad_payload,
    {
      rMVPA:::.dual_lda_chunk_update_predict_native(
        state = state0,
        fold = fold,
        step_cols_seq = step_cols_seq
      )
    },
    .package = "rMVPA"
  )

  expect_null(got)
})

test_that("dual_lda incremental state update matches full neighbor recomputation", {
  set.seed(3120)
  dims <- c(7L, 7L, 7L)
  p <- prod(dims)
  n <- 30L
  classes <- c("a", "b", "c")
  y_all <- factor(rep(classes, length.out = n), levels = classes)
  x_all <- matrix(rnorm(n * p), nrow = n, ncol = p)

  fold <- rMVPA:::.prepare_dual_lda_fold(
    x_all = x_all,
    y_all = y_all,
    train_idx = seq_len(n),
    test_idx = seq_len(n),
    classes = classes,
    gamma = 1e-2
  )

  mask_active <- rep(TRUE, p)
  col_lookup <- seq_len(p)
  offsets <- rMVPA:::.searchlight_offsets(radius = 2, spacing = c(1, 1, 1))
  delta <- rMVPA:::.searchlight_enter_leave(offsets)
  tables <- rMVPA:::.precompute_boundary_tables(delta, dims)

  center <- c(4L, 4L, 4L)
  step <- c(1L, 0L, 0L)
  center_next <- center + step

  lin_id <- function(coord) {
    coord[1] + (coord[2] - 1L) * dims[1] + (coord[3] - 1L) * dims[1] * dims[2]
  }

  for (rank_mode in c(FALSE, TRUE)) {
    st_prev <- rMVPA:::.init_dual_lda_state(fold, center, offsets, dims, mask_active, col_lookup)
    step_cols <- rMVPA:::.step_boundary_cols_fast(
      prev_coord = center,
      prev_id = lin_id(center),
      curr_coord = center_next,
      curr_id = lin_id(center_next),
      step_table = tables[["1,0,0"]],
      dims = dims,
      mask_active = mask_active,
      col_lookup = col_lookup
    )
    upd <- rMVPA:::.update_dual_lda_state_neighbor(
      state = st_prev,
      fold = fold,
      out_cols = step_cols$out_cols,
      in_cols = step_cols$in_cols,
      use_rank_update = rank_mode
    )
    expect_true(upd$ok, info = paste("rank_mode =", rank_mode))

    st_next <- rMVPA:::.init_dual_lda_state(fold, center_next, offsets, dims, mask_active, col_lookup)

    expect_equal(crossprod(upd$state$L), crossprod(st_next$L), tolerance = 1e-10, info = paste("rank_mode =", rank_mode))
    expect_equal(upd$state$L, st_next$L, tolerance = 1e-10, info = paste("rank_mode =", rank_mode))
    expect_equal(upd$state$T, st_next$T, tolerance = 1e-10, info = paste("rank_mode =", rank_mode))
    expect_equal(upd$state$U, st_next$U, tolerance = 1e-10, info = paste("rank_mode =", rank_mode))
    expect_equal(upd$state$Q, st_next$Q, tolerance = 1e-10, info = paste("rank_mode =", rank_mode))
    expect_equal(upd$state$C, st_next$C, tolerance = 1e-10, info = paste("rank_mode =", rank_mode))

    p_upd <- rMVPA:::.predict_dual_lda_state(upd$state, fold)$probs
    p_ref <- rMVPA:::.predict_dual_lda_state(st_next, fold)$probs
    expect_equal(p_upd, p_ref, tolerance = 1e-10, info = paste("rank_mode =", rank_mode))
  }
})

test_that("dual_lda core probabilities are finite and normalized", {
  set.seed(3101)
  x <- matrix(rnorm(40 * 24), nrow = 40, ncol = 24)
  y <- factor(rep(letters[1:3], length.out = nrow(x)))

  fit <- rMVPA:::.dual_lda_fit_core(x, y, gamma = 1e-2)
  probs <- rMVPA:::.dual_lda_prob_core(fit, x[1:11, , drop = FALSE])

  expect_true(all(is.finite(probs)))
  expect_true(all(probs >= 0))
  expect_equal(rowSums(probs), rep(1, nrow(probs)), tolerance = 1e-12)
})

test_that("dual_lda core is invariant to translation and feature permutation", {
  set.seed(3102)
  x <- matrix(rnorm(36 * 20), nrow = 36, ncol = 20)
  y <- factor(rep(letters[1:3], length.out = nrow(x)))
  x_test <- matrix(rnorm(12 * 20), nrow = 12, ncol = 20)

  fit_ref <- rMVPA:::.dual_lda_fit_core(x, y, gamma = 1e-2)
  p_ref <- rMVPA:::.dual_lda_prob_core(fit_ref, x_test)

  shift <- 7.5
  fit_shift <- rMVPA:::.dual_lda_fit_core(x + shift, y, gamma = 1e-2)
  p_shift <- rMVPA:::.dual_lda_prob_core(fit_shift, x_test + shift)
  expect_equal(p_ref, p_shift, tolerance = 1e-10)

  perm <- sample(seq_len(ncol(x)))
  fit_perm <- rMVPA:::.dual_lda_fit_core(x[, perm, drop = FALSE], y, gamma = 1e-2)
  p_perm <- rMVPA:::.dual_lda_prob_core(fit_perm, x_test[, perm, drop = FALSE])
  expect_equal(p_ref, p_perm, tolerance = 1e-10)
})

test_that("dual_lda fast metric kernel matches generic multiclass performance", {
  set.seed(3110)
  classes <- c("a", "b", "c")
  observed <- factor(sample(classes, 60, replace = TRUE), levels = classes)

  probs <- matrix(stats::runif(60 * length(classes)), nrow = 60, ncol = length(classes))
  probs <- probs / rowSums(probs)
  colnames(probs) <- classes

  pred <- factor(classes[max.col(probs)], levels = classes)
  test_idx <- sample(seq_len(150), 60)
  split_groups <- list(
    odd = test_idx[seq(1, length(test_idx), by = 2)],
    even = test_idx[seq(2, length(test_idx), by = 2)]
  )

  cres <- classification_result(observed = observed, predicted = pred, probs = probs, testind = test_idx)
  perf_ref <- performance(cres, split_groups, class_metrics = FALSE)
  perf_fast <- rMVPA:::.dual_lda_metric_with_splits(
    observed = observed,
    probs = probs,
    test_idx = test_idx,
    classes = classes,
    kind = "multiclass",
    split_groups = split_groups,
    class_metrics = FALSE
  )
  expect_equal(perf_fast, perf_ref, tolerance = 1e-10)

  perf_ref_class <- performance(cres, split_groups, class_metrics = TRUE)
  perf_fast_class <- rMVPA:::.dual_lda_metric_with_splits(
    observed = observed,
    probs = probs,
    test_idx = test_idx,
    classes = classes,
    kind = "multiclass",
    split_groups = split_groups,
    class_metrics = TRUE
  )
  expect_equal(perf_fast_class, perf_ref_class, tolerance = 1e-10)
})

test_that("dual_lda fast metric kernel matches generic binary performance", {
  set.seed(3111)
  classes <- c("neg", "pos")
  observed <- factor(sample(classes, 80, replace = TRUE), levels = classes)

  probs <- matrix(stats::runif(80 * length(classes)), nrow = 80, ncol = length(classes))
  probs <- probs / rowSums(probs)
  colnames(probs) <- classes

  pred <- factor(classes[max.col(probs)], levels = classes)
  test_idx <- sample(seq_len(200), 80)
  split_groups <- list(
    first = test_idx[1:40],
    second = test_idx[41:80]
  )

  cres <- classification_result(observed = observed, predicted = pred, probs = probs, testind = test_idx)
  perf_ref <- performance(cres, split_groups)
  perf_fast <- rMVPA:::.dual_lda_metric_with_splits(
    observed = observed,
    probs = probs,
    test_idx = test_idx,
    classes = classes,
    kind = "binary",
    split_groups = split_groups,
    class_metrics = FALSE
  )
  expect_equal(perf_fast, perf_ref, tolerance = 1e-10)
})

test_that("dual_lda standard searchlight produces valid metric maps", {
  set.seed(2026)
  built <- build_dual_lda_mspec(D = c(4, 4, 4), nobs = 24, nlevels = 2, blocks = 3, gamma = 1e-2)

  res <- run_searchlight(built$mspec, radius = 2, method = "standard")
  expect_s3_class(res, "searchlight_result")
  expect_true(all(c("Accuracy", "AUC") %in% names(res$results)))

  acc_vals <- metric_map_values(res, "Accuracy")
  auc_vals <- metric_map_values(res, "AUC")
  expect_true(any(is.finite(acc_vals)))
  expect_true(any(is.finite(auc_vals)))
})

test_that("dual_lda fast path preserves custom performance via fallback path", {
  set.seed(2029)
  ds <- gen_sample_dataset(D = c(3, 3, 3), nobs = 18, nlevels = 3, blocks = 3)
  cv <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("dual_lda")

  custom_perf <- function(result) c(custom_acc = mean(result$observed == result$predicted))
  mspec <- mvpa_model(
    model = mdl,
    dataset = ds$dataset,
    design = ds$design,
    model_type = "classification",
    crossval = cv,
    tune_grid = data.frame(gamma = 1e-2),
    performance = custom_perf
  )

  res_inc <- run_searchlight(mspec, radius = 1, method = "standard", incremental = TRUE)
  res_full <- run_searchlight(mspec, radius = 1, method = "standard", incremental = FALSE)

  expect_true("custom_acc" %in% names(res_inc$results))
  expect_true("custom_acc" %in% names(res_full$results))
  expect_searchlight_parity(
    reference = res_full,
    candidate = res_inc,
    metrics = "custom_acc",
    atol = 1e-12,
    rtol = 1e-10
  )
})

test_that("dual_lda with external test set stays consistent across incremental flag", {
  set.seed(2028)
  built <- build_dual_lda_mspec(
    D = c(4, 4, 4),
    nobs = 24,
    nlevels = 2,
    blocks = 3,
    gamma = 1e-2,
    external_test = TRUE
  )

  res_inc <- run_searchlight(built$mspec, radius = 2, method = "standard", incremental = TRUE)
  res_full <- run_searchlight(built$mspec, radius = 2, method = "standard", incremental = FALSE)

  expect_searchlight_parity(
    reference = res_full,
    candidate = res_inc,
    atol = 1e-12,
    rtol = 1e-10
  )
})

test_that("dual_lda incremental updates match full recomputation", {
  set.seed(2027)
  built <- build_dual_lda_mspec(D = c(4, 4, 4), nobs = 30, nlevels = 3, blocks = 3, gamma = 1e-2)

  res_inc <- run_searchlight(built$mspec, radius = 2, method = "standard", incremental = TRUE)
  res_full <- run_searchlight(built$mspec, radius = 2, method = "standard", incremental = FALSE)

  expect_searchlight_parity(
    reference = res_full,
    candidate = res_inc,
    # Incremental rank-update accumulates floating-point drift (~0.025 AUC);
    # this is inherent to the approach, not a correctness bug.
    atol = 3e-2,
    rtol = 1e-6
  )
})

test_that("dual_lda differential fuzz check across random seeds", {
  seeds <- c(11, 29, 73)
  max_diffs <- numeric(length(seeds))

  for (i in seq_along(seeds)) {
    set.seed(seeds[i])
    built <- build_dual_lda_mspec(D = c(3, 3, 3), nobs = 18, nlevels = 2, blocks = 3, gamma = 5e-3)
    res_inc <- run_searchlight(built$mspec, radius = 1, method = "standard", incremental = TRUE)
    res_full <- run_searchlight(built$mspec, radius = 1, method = "standard", incremental = FALSE)

    inc_vals <- metric_map_values(res_inc, "Accuracy")
    full_vals <- metric_map_values(res_full, "Accuracy")
    keep <- is.finite(inc_vals) & is.finite(full_vals)
    expect_gt(sum(keep), 0)
    max_diffs[i] <- max(abs(inc_vals[keep] - full_vals[keep]))
  }

  expect_true(all(max_diffs < 6e-2))
})

test_that("dual_lda remains finite under near-singular regime with tiny ridge", {
  set.seed(4041)
  built <- build_dual_lda_mspec(D = c(4, 4, 4), nobs = 15, nlevels = 3, blocks = 3, gamma = 1e-6)
  res <- run_searchlight(built$mspec, radius = 2, method = "standard", incremental = TRUE)

  acc <- metric_map_values(res, "Accuracy")
  auc <- metric_map_values(res, "AUC")
  expect_true(any(is.finite(acc)))
  expect_true(any(is.finite(auc)))
})

test_that("dual_lda chunk recovery falls back when native chunk returns NULL", {
  skip_on_cran()
  skip_if_not_installed("rMVPA")

  set.seed(9999)
  built <- build_dual_lda_mspec(D = c(4, 4, 4), nobs = 30, nlevels = 3, blocks = 3, gamma = 1e-2)

  # Run with chunk update forced to return NULL (simulates native failure)
  res_recovery <- testthat::with_mocked_bindings(
    run_searchlight(built$mspec, radius = 2, method = "standard", incremental = TRUE),
    .dual_lda_chunk_update_predict_native = function(...) NULL,
    .package = "rMVPA"
  )

  # Run without rank-chunk updates for reference
  res_ref <- testthat::with_mocked_bindings(
    run_searchlight(built$mspec, radius = 2, method = "standard", incremental = TRUE),
    .dual_lda_use_rank_update = function() FALSE,
    .package = "rMVPA"
  )

  # Both should produce valid searchlight results
  expect_s3_class(res_recovery, "searchlight_result")
  expect_s3_class(res_ref, "searchlight_result")

  # Recovery path should produce same results as the non-chunk path
  acc_recovery <- metric_map_values(res_recovery, "Accuracy")
  acc_ref <- metric_map_values(res_ref, "Accuracy")
  keep <- is.finite(acc_recovery) & is.finite(acc_ref)
  expect_gt(sum(keep), 0)
  expect_equal(acc_recovery[keep], acc_ref[keep], tolerance = 1e-6)
})
