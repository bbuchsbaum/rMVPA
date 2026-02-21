#' @keywords internal
#' @noRd
.is_dual_lda_fast_path <- function(model_spec, method) {
  inherits(model_spec, "mvpa_model") &&
    method %in% c("standard", "randomized", "resampled") &&
    isTRUE(has_crossval(model_spec)) &&
    !isTRUE(has_test_set(model_spec)) &&
    !is.null(model_spec$model) &&
    identical(model_spec$model$label, "dual_lda")
}

#' @keywords internal
#' @noRd
.dual_lda_gamma <- function(model_spec, gamma = NULL) {
  if (!is.null(gamma)) {
    return(as.numeric(gamma)[1])
  }
  tg <- model_spec$tune_grid
  if (is.data.frame(tg) && "gamma" %in% names(tg) && nrow(tg) > 0L) {
    return(as.numeric(tg$gamma[[1]]))
  }
  1e-2
}

#' @keywords internal
#' @noRd
.dual_lda_rank_update_preference <- function() {
  NULL
}

#' @keywords internal
#' @noRd
.dual_lda_use_rank_update <- function() {
  isTRUE(.dual_lda_rank_update_preference())
}

#' @keywords internal
#' @noRd
.dual_lda_rank_chunk_size <- function() {
  4L
}

#' @keywords internal
#' @noRd
.dual_lda_rank_heuristic_cache <- local({
  new.env(hash = TRUE, parent = emptyenv())
})

#' @keywords internal
#' @noRd
.dual_lda_rank_cache_key <- function(fold, radius, boundary_size, chunk_size = 1L) {
  paste(
    as.integer(nrow(fold$R)),
    as.integer(ncol(fold$M)),
    format(as.numeric(radius)[1], scientific = FALSE, trim = TRUE),
    as.integer(boundary_size),
    as.integer(chunk_size),
    sep = ":"
  )
}

#' @keywords internal
#' @noRd
.dual_lda_rank_cache_get <- function(key) {
  if (exists(key, envir = .dual_lda_rank_heuristic_cache, inherits = FALSE)) {
    get(key, envir = .dual_lda_rank_heuristic_cache, inherits = FALSE)
  } else {
    NULL
  }
}

#' @keywords internal
#' @noRd
.dual_lda_rank_cache_set <- function(key, value) {
  assign(key, isTRUE(value), envir = .dual_lda_rank_heuristic_cache)
  invisible(value)
}

#' @keywords internal
#' @noRd
.dual_lda_rank_cache_clear <- function() {
  rm(list = ls(envir = .dual_lda_rank_heuristic_cache, all.names = TRUE), envir = .dual_lda_rank_heuristic_cache)
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.dual_lda_clone_state <- function(state) {
  lapply(state, function(x) {
    if (is.matrix(x) && is.numeric(x)) {
      y <- matrix(0, nrow = nrow(x), ncol = ncol(x))
      y[] <- x
      y
    } else {
      x
    }
  })
}

#' @keywords internal
#' @noRd
.dual_lda_first_step_cols <- function(ord_coords, ord_ids, boundary_tables,
                                      dims, mask_active, col_lookup) {
  if (length(ord_ids) < 2L) {
    return(NULL)
  }

  for (pos in 2:length(ord_ids)) {
    prev_coord <- ord_coords[pos - 1L, ]
    curr_coord <- ord_coords[pos, ]
    dx <- as.integer(curr_coord - prev_coord)
    if (!.is_unit_axis_step(dx)) {
      next
    }
    key <- .dir_key(dx)
    return(.step_boundary_cols_fast(
      prev_coord = prev_coord,
      prev_id = ord_ids[pos - 1L],
      curr_coord = curr_coord,
      curr_id = ord_ids[pos],
      step_table = boundary_tables[[key]],
      dims = dims,
      mask_active = mask_active,
      col_lookup = col_lookup
    ))
  }

  NULL
}

#' @keywords internal
#' @noRd
.dual_lda_rank_update_benchmark <- function(state, fold, out_cols, in_cols,
                                            radius = NA_real_, chunk_size = 1L,
                                            trials = 8L, verbose = FALSE) {
  pref <- .dual_lda_rank_update_preference()
  if (!is.null(pref)) {
    return(isTRUE(pref))
  }

  boundary_size <- length(out_cols) + length(in_cols)
  if (boundary_size == 0L) {
    return(FALSE)
  }

  cache_key <- .dual_lda_rank_cache_key(
    fold = fold,
    radius = radius,
    boundary_size = boundary_size,
    chunk_size = chunk_size
  )
  cache_val <- .dual_lda_rank_cache_get(cache_key)
  if (!is.null(cache_val)) {
    if (isTRUE(verbose)) {
      futile.logger::flog.info(
        "dual_lda rank-update heuristic cache hit: key=%s use_rank=%s",
        cache_key,
        as.character(cache_val)
      )
    }
    return(isTRUE(cache_val))
  }

  has_chunk_native <- .dual_lda_native_symbol_loaded("rmvpa_dual_lda_chunk_update_predict")
  can_chunk <- isTRUE(has_chunk_native) && chunk_size > 1L

  bench_mode <- function(mode) {
    st <- .dual_lda_clone_state(state)

    if (isTRUE(mode) && isTRUE(can_chunk)) {
      seq_steps <- rep(list(list(out_cols = out_cols, in_cols = in_cols)), trials)
      t0 <- proc.time()[3]
      upd <- .dual_lda_chunk_update_predict_native(
        state = st,
        fold = fold,
        step_cols_seq = seq_steps
      )
      if (is.null(upd)) {
        return(Inf)
      }
      return(proc.time()[3] - t0)
    }

    t0 <- proc.time()[3]
    for (i in seq_len(trials)) {
      upd <- .update_dual_lda_state_neighbor(
        state = st,
        fold = fold,
        out_cols = out_cols,
        in_cols = in_cols,
        use_rank_update = mode
      )
      if (!isTRUE(upd$ok)) {
        return(Inf)
      }
      st <- upd$state
      .predict_dual_lda_state(st, fold)
    }
    proc.time()[3] - t0
  }

  t_chol <- bench_mode(FALSE)
  t_rank <- bench_mode(TRUE)

  use_rank <- is.finite(t_chol) && is.finite(t_rank) && t_rank < t_chol * 0.95
  if (isTRUE(verbose)) {
    futile.logger::flog.info(
      "dual_lda rank-update heuristic: chol=%.4fs rank=%.4fs use_rank=%s",
      t_chol, t_rank, as.character(use_rank)
    )
  }

  .dual_lda_rank_cache_set(cache_key, use_rank)
  use_rank
}

#' @keywords internal
#' @noRd
.dir_key <- function(dx) paste(as.integer(dx), collapse = ",")

#' @keywords internal
#' @noRd
.is_unit_axis_step <- function(dx) {
  dx <- as.integer(dx)
  sum(abs(dx)) == 1L && all(abs(dx) <= 1L)
}

#' @keywords internal
#' @noRd
.searchlight_offsets <- function(radius, spacing) {
  spacing <- as.numeric(spacing)
  if (length(spacing) < 3L) spacing <- c(spacing, rep(1, 3L - length(spacing)))
  spacing <- spacing[1:3]
  max_off <- floor(radius / pmax(spacing, .Machine$double.eps))

  grid <- expand.grid(
    dx = seq.int(-max_off[1], max_off[1]),
    dy = seq.int(-max_off[2], max_off[2]),
    dz = seq.int(-max_off[3], max_off[3])
  )

  dist2 <- (grid$dx * spacing[1])^2 +
    (grid$dy * spacing[2])^2 +
    (grid$dz * spacing[3])^2
  keep <- dist2 <= (radius^2 + 1e-12)

  as.matrix(grid[keep, c("dx", "dy", "dz"), drop = FALSE])
}

#' @keywords internal
#' @noRd
.searchlight_enter_leave <- function(offsets) {
  keys <- paste(offsets[, 1], offsets[, 2], offsets[, 3], sep = ",")
  key_env <- new.env(hash = TRUE, parent = emptyenv())
  for (k in keys) key_env[[k]] <- TRUE

  has_offset <- function(mat) {
    if (is.null(dim(mat))) {
      mat <- matrix(mat, ncol = 3)
    }
    qk <- paste(mat[, 1], mat[, 2], mat[, 3], sep = ",")
    vapply(qk, exists, logical(1), envir = key_env, inherits = FALSE)
  }

  dirs <- rbind(
    c(1L, 0L, 0L),
    c(-1L, 0L, 0L),
    c(0L, 1L, 0L),
    c(0L, -1L, 0L),
    c(0L, 0L, 1L),
    c(0L, 0L, -1L)
  )

  out <- vector("list", nrow(dirs))
  names(out) <- apply(dirs, 1, .dir_key)
  for (i in seq_len(nrow(dirs))) {
    d <- dirs[i, ]
    leave_mask <- !has_offset(sweep(offsets, 2, d, "-"))
    enter_mask <- !has_offset(sweep(offsets, 2, d, "+"))
    out[[i]] <- list(
      leave = offsets[leave_mask, , drop = FALSE],
      enter = offsets[enter_mask, , drop = FALSE]
    )
  }
  out
}

#' @keywords internal
#' @noRd
.snake_center_order <- function(coords) {
  df <- data.frame(
    idx = seq_len(nrow(coords)),
    x = coords[, 1],
    y = coords[, 2],
    z = coords[, 3]
  )

  z_vals <- sort(unique(df$z))
  ord <- integer(0)
  flip_z <- FALSE

  for (z in z_vals) {
    yz <- df[df$z == z, , drop = FALSE]
    y_vals <- sort(unique(yz$y))
    if (flip_z) y_vals <- rev(y_vals)

    flip_row <- FALSE
    for (y in y_vals) {
      row <- yz[yz$y == y, , drop = FALSE]
      if (nrow(row) == 0L) next
      row <- row[order(row$x, decreasing = flip_row), , drop = FALSE]
      ord <- c(ord, row$idx)
      flip_row <- !flip_row
    }
    flip_z <- !flip_z
  }

  ord
}

#' @keywords internal
#' @noRd
.offset_to_indices <- function(center_coord, offsets, dims, mask_active) {
  if (nrow(offsets) == 0L) {
    return(integer(0))
  }

  coords <- sweep(offsets, 2, center_coord, "+")
  in_bounds <- coords[, 1] >= 1L & coords[, 1] <= dims[1] &
    coords[, 2] >= 1L & coords[, 2] <= dims[2] &
    coords[, 3] >= 1L & coords[, 3] <= dims[3]

  if (!any(in_bounds)) {
    return(integer(0))
  }
  coords <- coords[in_bounds, , drop = FALSE]

  idx <- coords[, 1] +
    (coords[, 2] - 1L) * dims[1] +
    (coords[, 3] - 1L) * dims[1] * dims[2]

  idx <- idx[mask_active[idx]]
  unique(as.integer(idx))
}

#' @keywords internal
#' @noRd
.chol_rank1_update_upper <- function(R, x, downdate = FALSE) {
  p <- nrow(R)
  x <- as.numeric(x)

  for (k in seq_len(p)) {
    rkk <- R[k, k]
    xk <- x[k]
    r2 <- if (downdate) rkk * rkk - xk * xk else rkk * rkk + xk * xk

    if (!is.finite(r2) || r2 <= .Machine$double.eps) {
      return(list(R = R, ok = FALSE))
    }

    r <- sqrt(r2)
    c <- r / rkk
    s <- xk / rkk
    R[k, k] <- r

    if (k < p) {
      idx <- (k + 1L):p
      old_row <- R[k, idx]
      new_row <- if (downdate) {
        (old_row - s * x[idx]) / c
      } else {
        (old_row + s * x[idx]) / c
      }
      R[k, idx] <- new_row
      x[idx] <- c * x[idx] - s * new_row
    }
  }

  list(R = R, ok = TRUE)
}

#' @keywords internal
#' @noRd
.dual_lda_call_native <- function(symbol, ...) {
  ccall <- .Call
  ccall(symbol, ..., PACKAGE = "rMVPA")
}

#' @keywords internal
#' @noRd
.dual_lda_native_symbol_loaded <- function(symbol) {
  symbol <- as.character(symbol)[1]
  isTRUE(tryCatch(
    is.loaded(symbol, PACKAGE = "rMVPA"),
    error = function(...) FALSE
  ))
}

#' @keywords internal
#' @noRd
.chol_rankk_update_upper <- function(R, X, downdate = FALSE, scale = 1) {
  scale <- as.numeric(scale)[1]
  if (length(X) == 0L) {
    return(list(R = R, ok = TRUE))
  }

  X <- as.matrix(X)
  if (ncol(X) == 0L) {
    return(list(R = R, ok = TRUE))
  }
  if (!all(is.finite(X))) {
    return(list(R = R, ok = FALSE))
  }

  if (.dual_lda_native_symbol_loaded("rmvpa_chol_rankk_update")) {
    native_out <- try(
      .dual_lda_call_native(
        "rmvpa_chol_rankk_update",
        R,
        X,
        as.logical(downdate),
        as.numeric(scale)
      ),
      silent = TRUE
    )
    if (!inherits(native_out, "try-error") && !is.null(native_out)) {
      return(list(R = native_out, ok = TRUE))
    }
  }

  # Fallback for environments where native symbols are unavailable.
  Rout <- R
  for (j in seq_len(ncol(X))) {
    step <- .chol_rank1_update_upper(Rout, X[, j] * scale, downdate = downdate)
    if (!isTRUE(step$ok)) {
      return(list(R = R, ok = FALSE))
    }
    Rout <- step$R
  }

  list(R = Rout, ok = TRUE)
}

#' @keywords internal
#' @noRd
.dual_lda_finite_matrix <- function(x, nrow = NULL, ncol = NULL) {
  if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x))) {
    return(FALSE)
  }

  if (!is.null(nrow) && nrow(x) != as.integer(nrow)[1]) {
    return(FALSE)
  }
  if (!is.null(ncol) && ncol(x) != as.integer(ncol)[1]) {
    return(FALSE)
  }

  TRUE
}

#' @keywords internal
#' @noRd
.dual_lda_valid_native_state <- function(state, n, k, m) {
  is.list(state) &&
    .dual_lda_finite_matrix(state$L, nrow = n, ncol = n) &&
    .dual_lda_finite_matrix(state$T, nrow = n, ncol = k) &&
    .dual_lda_finite_matrix(state$U, nrow = k, ncol = k) &&
    .dual_lda_finite_matrix(state$Q, nrow = n, ncol = m) &&
    .dual_lda_finite_matrix(state$C, nrow = m, ncol = k)
}

#' @keywords internal
#' @noRd
.dual_lda_valid_prob_matrix <- function(p, classes) {
  if (!.dual_lda_finite_matrix(p, ncol = length(classes))) {
    return(FALSE)
  }
  if (any(p < -1e-12)) {
    return(FALSE)
  }

  rs <- rowSums(p)
  all(is.finite(rs)) && all(abs(rs - 1) <= 1e-6)
}

#' @keywords internal
#' @noRd
.dual_lda_step_update_native <- function(state, fold, out_cols, in_cols) {
  if (length(out_cols) == 0L && length(in_cols) == 0L) {
    return(state)
  }
  if (!.dual_lda_native_symbol_loaded("rmvpa_dual_lda_step_update")) {
    return(NULL)
  }

  out <- try(
    .dual_lda_call_native(
      "rmvpa_dual_lda_step_update",
      state,
      fold$R,
      fold$M,
      fold$Xtest_t,
      as.integer(out_cols),
      as.integer(in_cols),
      as.numeric(fold$gamma_inv)
    ),
    silent = TRUE
  )

  if (inherits(out, "try-error") || is.null(out)) {
    NULL
  } else if (!.dual_lda_valid_native_state(
    out,
    n = nrow(fold$R),
    k = ncol(fold$M),
    m = ncol(fold$Xtest_t)
  )) {
    NULL
  } else {
    out
  }
}

#' @keywords internal
#' @noRd
.dual_lda_pack_boundary_steps <- function(step_cols_seq) {
  n_steps <- length(step_cols_seq)
  if (n_steps == 0L) {
    return(list(
      out_ptr = integer(1L),
      out_idx = integer(0),
      in_ptr = integer(1L),
      in_idx = integer(0)
    ))
  }

  out_len <- sum(vapply(step_cols_seq, function(x) length(x$out_cols), integer(1)))
  in_len <- sum(vapply(step_cols_seq, function(x) length(x$in_cols), integer(1)))

  out_ptr <- integer(n_steps + 1L)
  in_ptr <- integer(n_steps + 1L)
  out_idx <- integer(out_len)
  in_idx <- integer(in_len)

  out_cursor <- 0L
  in_cursor <- 0L
  for (i in seq_len(n_steps)) {
    out_ptr[i] <- out_cursor
    in_ptr[i] <- in_cursor

    ocols <- as.integer(step_cols_seq[[i]]$out_cols)
    icols <- as.integer(step_cols_seq[[i]]$in_cols)

    if (length(ocols) > 0L) {
      rng <- seq.int(out_cursor + 1L, out_cursor + length(ocols))
      out_idx[rng] <- ocols
      out_cursor <- out_cursor + length(ocols)
    }
    if (length(icols) > 0L) {
      rng <- seq.int(in_cursor + 1L, in_cursor + length(icols))
      in_idx[rng] <- icols
      in_cursor <- in_cursor + length(icols)
    }
  }

  out_ptr[n_steps + 1L] <- out_cursor
  in_ptr[n_steps + 1L] <- in_cursor

  list(
    out_ptr = out_ptr,
    out_idx = out_idx,
    in_ptr = in_ptr,
    in_idx = in_idx
  )
}

#' @keywords internal
#' @noRd
.dual_lda_chunk_update_predict_native <- function(state, fold, step_cols_seq) {
  if (length(step_cols_seq) == 0L) {
    return(list(state = state, probs = list()))
  }
  if (!.dual_lda_native_symbol_loaded("rmvpa_dual_lda_chunk_update_predict")) {
    return(NULL)
  }

  packed <- .dual_lda_pack_boundary_steps(step_cols_seq)
  out <- try(
    .dual_lda_call_native(
      "rmvpa_dual_lda_chunk_update_predict",
      state,
      fold$R,
      fold$M,
      fold$Xtest_t,
      packed$out_ptr,
      packed$out_idx,
      packed$in_ptr,
      packed$in_idx,
      as.numeric(fold$gamma_inv),
      as.numeric(fold$log_priors)
    ),
    silent = TRUE
  )

  if (inherits(out, "try-error") || is.null(out) || !is.list(out)) {
    return(NULL)
  }

  if (!.dual_lda_valid_native_state(
    out$state,
    n = nrow(fold$R),
    k = ncol(fold$M),
    m = ncol(fold$Xtest_t)
  )) {
    return(NULL)
  }

  if (!is.list(out$probs) || length(out$probs) != length(step_cols_seq)) {
    return(NULL)
  }
  if (!all(vapply(out$probs, .dual_lda_valid_prob_matrix, logical(1), classes = fold$classes))) {
    NULL
  } else {
    out
  }
}

#' @keywords internal
#' @noRd
.solve_from_chol_upper <- function(U, B) {
  # Solve U' U X = B without materializing t(U) each call.
  y <- backsolve(U, B, transpose = TRUE)
  backsolve(U, y)
}

#' @keywords internal
#' @noRd
.dual_lda_fast_metric_mode <- function(model_spec, classes) {
  perf_fun <- model_spec$performance
  perf_kind <- attr(perf_fun, "rmvpa_perf_kind", exact = TRUE)

  if (is.null(perf_kind)) {
    return(list(enabled = FALSE))
  }

  if (identical(perf_kind, "binary") && length(classes) == 2L) {
    return(list(
      enabled = TRUE,
      kind = "binary",
      class_metrics = FALSE,
      split_groups = model_spec$design$split_groups
    ))
  }

  if (identical(perf_kind, "multiclass") && length(classes) > 2L) {
    return(list(
      enabled = TRUE,
      kind = "multiclass",
      class_metrics = isTRUE(attr(perf_fun, "rmvpa_class_metrics", exact = TRUE)),
      split_groups = model_spec$design$split_groups
    ))
  }

  list(enabled = FALSE)
}

#' @keywords internal
#' @noRd
.dual_lda_auc_from_scores <- function(score, positive) {
  positive <- as.logical(positive)
  n_pos <- sum(positive)
  n_neg <- length(positive) - n_pos

  if (n_pos == 0L || n_neg == 0L || any(!is.finite(score))) {
    return(NA_real_)
  }

  r <- rank(score, ties.method = "average")
  (sum(r[positive]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

#' @keywords internal
#' @noRd
.dual_lda_metric_core <- function(observed, probs, classes, kind, class_metrics = FALSE) {
  pred <- factor(classes[max.col(probs)], levels = classes)
  acc <- mean(pred == observed)

  if (identical(kind, "binary")) {
    score <- if (ncol(probs) >= 2L) probs[, 2L] else probs[, 1L]
    auc <- .dual_lda_auc_from_scores(score, observed == classes[2L])
    auc_centered <- if (is.na(auc)) NA_real_ else 2 * auc - 1
    return(c(Accuracy = acc, AUC = auc_centered))
  }

  k <- length(classes)
  aucres <- rep(NA_real_, k)
  for (i in seq_len(k)) {
    pclass <- probs[, i]
    pother <- rowMeans(probs[, -i, drop = FALSE])
    score <- pclass - pother
    aucres[i] <- .dual_lda_auc_from_scores(score, observed == classes[i])
  }

  auc_centered <- 2 * aucres - 1
  mean_auc <- mean(auc_centered, na.rm = TRUE)
  if (is.nan(mean_auc)) {
    mean_auc <- NA_real_
  }

  metrics <- c(Accuracy = acc, AUC = mean_auc)
  if (isTRUE(class_metrics)) {
    names(auc_centered) <- paste0("AUC_", classes)
    c(metrics, auc_centered)
  } else {
    metrics
  }
}

#' @keywords internal
#' @noRd
.dual_lda_metric_with_splits <- function(observed, probs, test_idx, classes, kind,
                                         split_groups = NULL, class_metrics = FALSE) {
  base_vals <- .dual_lda_metric_core(
    observed = observed,
    probs = probs,
    classes = classes,
    kind = kind,
    class_metrics = class_metrics
  )

  if (is.null(split_groups) || length(split_groups) == 0L) {
    return(base_vals)
  }

  tags <- names(split_groups)
  if (is.null(tags) || length(tags) == 0L) {
    return(base_vals)
  }

  split_vecs <- lapply(tags, function(tag) {
    design_idx <- split_groups[[tag]]
    res_idx <- which(test_idx %in% design_idx)

    if (!length(res_idx)) {
      vals <- rep(NA_real_, length(base_vals))
      names(vals) <- paste0(names(base_vals), "_", tag)
      return(vals)
    }

    vals <- .dual_lda_metric_core(
      observed = observed[res_idx],
      probs = probs[res_idx, , drop = FALSE],
      classes = classes,
      kind = kind,
      class_metrics = class_metrics
    )
    names(vals) <- paste0(names(vals), "_", tag)
    vals
  })

  c(base_vals, unlist(split_vecs))
}

#' @keywords internal
#' @noRd
.combine_dual_lda_fold_probs <- function(fold_preds, fold_row_index, n_rows, n_classes) {
  prob <- matrix(0, nrow = n_rows, ncol = n_classes)
  for (i in seq_along(fold_preds)) {
    ridx <- fold_row_index[[i]]
    prob[ridx, ] <- prob[ridx, ] + fold_preds[[i]]$probs
  }

  row_sum <- rowSums(prob)
  nz <- row_sum > 0
  prob[nz, ] <- prob[nz, , drop = FALSE] / row_sum[nz]
  if (any(!nz)) {
    prob[!nz, ] <- matrix(1 / n_classes, nrow = sum(!nz), ncol = n_classes)
  }

  prob
}

#' @keywords internal
#' @noRd
.prepare_dual_lda_fold <- function(x_all, y_all, train_idx, test_idx, classes, gamma) {
  x_train <- x_all[train_idx, , drop = FALSE]
  x_test <- x_all[test_idx, , drop = FALSE]
  y_train <- factor(y_all[train_idx], levels = classes)

  if (any(!is.finite(x_train)) || any(!is.finite(x_test))) {
    x_train[!is.finite(x_train)] <- 0
    x_test[!is.finite(x_test)] <- 0
  }

  class_counts <- as.numeric(table(y_train)[classes])
  if (any(class_counts == 0L)) {
    stop("dual_lda fast path requires every class to appear in every training fold.")
  }

  k <- length(classes)
  p <- ncol(x_train)
  means <- matrix(0, nrow = k, ncol = p, dimnames = list(classes, NULL))
  for (i in seq_len(k)) {
    idx <- which(y_train == classes[i])
    means[i, ] <- colMeans(x_train[idx, , drop = FALSE])
  }

  class_id <- match(y_train, classes)
  residual <- x_train - means[class_id, , drop = FALSE]

  priors <- class_counts / sum(class_counts)
  list(
    train_idx = train_idx,
    test_idx = test_idx,
    y_test = factor(y_all[test_idx], levels = classes),
    classes = classes,
    gamma = gamma,
    gamma_inv = 1 / gamma,
    gamma_inv2 = 1 / (gamma * gamma),
    priors = priors,
    log_priors = log(pmax(priors, .Machine$double.eps)),
    R = residual,
    M = t(means),
    Xtest_t = t(x_test)
  )
}

#' @keywords internal
#' @noRd
.init_dual_lda_state <- function(fold, center_coord, offsets, dims, mask_active, col_lookup) {
  ids <- .offset_to_indices(center_coord, offsets, dims, mask_active)
  cols <- col_lookup[ids]
  cols <- cols[cols > 0L]
  cols <- unique(cols)

  n <- nrow(fold$R)
  k <- ncol(fold$M)
  m <- ncol(fold$Xtest_t)

  if (length(cols) == 0L) {
    return(list(
      A = diag(1, n),
      L = diag(1, n),
      T = matrix(0, n, k),
      U = matrix(0, k, k),
      Q = matrix(0, n, m),
      C = matrix(0, m, k)
    ))
  }

  r_sub <- fold$R[, cols, drop = FALSE]
  m_sub <- fold$M[cols, , drop = FALSE]
  x_sub_t <- fold$Xtest_t[cols, , drop = FALSE]

  gmat <- tcrossprod(r_sub)
  tmat <- r_sub %*% m_sub
  umat <- crossprod(m_sub)
  qmat <- r_sub %*% x_sub_t
  cmat <- crossprod(x_sub_t, m_sub)

  a <- diag(1, n) + gmat * fold$gamma_inv
  lmat <- try(chol(a), silent = TRUE)
  if (inherits(lmat, "try-error")) {
    lmat <- chol(a + diag(1e-8, n))
  }

  list(
    A = a,
    L = lmat,
    T = tmat,
    U = umat,
    Q = qmat,
    C = cmat
  )
}

#' @keywords internal
#' @noRd
.step_boundary_cols <- function(prev_coord, curr_coord, leave_offsets, enter_offsets,
                                dims, mask_active, col_lookup) {
  out_ids <- .offset_to_indices(prev_coord, leave_offsets, dims, mask_active)
  in_ids <- .offset_to_indices(curr_coord, enter_offsets, dims, mask_active)

  out_cols <- unique(col_lookup[out_ids])
  out_cols <- out_cols[out_cols > 0L]
  in_cols <- unique(col_lookup[in_ids])
  in_cols <- in_cols[in_cols > 0L]
  list(out_cols = out_cols, in_cols = in_cols)
}

#' @keywords internal
#' @noRd
.precompute_boundary_tables <- function(delta_map, dims) {
  d1 <- as.integer(dims[1])
  d12 <- as.integer(dims[1] * dims[2])

  build_table <- function(offsets) {
    if (is.null(offsets) || nrow(offsets) == 0L) {
      return(list(
        dx = integer(0),
        dy = integer(0),
        dz = integer(0),
        dlin = integer(0)
      ))
    }

    dx <- as.integer(offsets[, 1])
    dy <- as.integer(offsets[, 2])
    dz <- as.integer(offsets[, 3])
    list(
      dx = dx,
      dy = dy,
      dz = dz,
      dlin = as.integer(dx + dy * d1 + dz * d12)
    )
  }

  out <- vector("list", length(delta_map))
  names(out) <- names(delta_map)
  for (key in names(delta_map)) {
    out[[key]] <- list(
      leave = build_table(delta_map[[key]]$leave),
      enter = build_table(delta_map[[key]]$enter)
    )
  }
  out
}

#' @keywords internal
#' @noRd
.boundary_cols_from_table <- function(center_coord, center_id, table, dims, mask_active, col_lookup) {
  if (length(table$dlin) == 0L) {
    return(integer(0))
  }

  x <- center_coord[1] + table$dx
  y <- center_coord[2] + table$dy
  z <- center_coord[3] + table$dz
  in_bounds <- x >= 1L & x <= dims[1] &
    y >= 1L & y <= dims[2] &
    z >= 1L & z <= dims[3]

  if (!any(in_bounds)) {
    return(integer(0))
  }

  idx <- as.integer(center_id + table$dlin[in_bounds])
  keep <- mask_active[idx]
  if (!any(keep)) {
    return(integer(0))
  }

  cols <- col_lookup[idx[keep]]
  cols[cols > 0L]
}

#' @keywords internal
#' @noRd
.step_boundary_cols_fast <- function(prev_coord, prev_id, curr_coord, curr_id, step_table,
                                     dims, mask_active, col_lookup) {
  list(
    out_cols = .boundary_cols_from_table(
      center_coord = prev_coord,
      center_id = prev_id,
      table = step_table$leave,
      dims = dims,
      mask_active = mask_active,
      col_lookup = col_lookup
    ),
    in_cols = .boundary_cols_from_table(
      center_coord = curr_coord,
      center_id = curr_id,
      table = step_table$enter,
      dims = dims,
      mask_active = mask_active,
      col_lookup = col_lookup
    )
  )
}

#' @keywords internal
#' @noRd
.update_dual_lda_state_neighbor <- function(state, fold, out_cols, in_cols,
                                            use_rank_update = .dual_lda_use_rank_update()) {
  changed <- FALSE
  if (isTRUE(use_rank_update)) {
    native_state <- .dual_lda_step_update_native(state, fold, out_cols, in_cols)
    if (!is.null(native_state)) {
      return(list(state = native_state, ok = TRUE))
    }
  }

  lnew <- state$L
  chol_scale <- sqrt(fold$gamma_inv)

  if (length(out_cols) > 0L) {
    changed <- TRUE
    r_out <- fold$R[, out_cols, drop = FALSE]
    m_out <- fold$M[out_cols, , drop = FALSE]
    x_out_t <- fold$Xtest_t[out_cols, , drop = FALSE]

    if (isTRUE(use_rank_update)) {
      upd <- .chol_rankk_update_upper(lnew, r_out, downdate = TRUE, scale = chol_scale)
      if (!isTRUE(upd$ok)) {
        return(list(state = state, ok = FALSE))
      }
      lnew <- upd$R
    } else {
      state$A <- state$A - tcrossprod(r_out) * fold$gamma_inv
    }
    state$T <- state$T - r_out %*% m_out
    state$U <- state$U - crossprod(m_out)
    state$Q <- state$Q - r_out %*% x_out_t
    state$C <- state$C - crossprod(x_out_t, m_out)
  }

  if (length(in_cols) > 0L) {
    changed <- TRUE
    r_in <- fold$R[, in_cols, drop = FALSE]
    m_in <- fold$M[in_cols, , drop = FALSE]
    x_in_t <- fold$Xtest_t[in_cols, , drop = FALSE]

    if (isTRUE(use_rank_update)) {
      upd <- .chol_rankk_update_upper(lnew, r_in, downdate = FALSE, scale = chol_scale)
      if (!isTRUE(upd$ok)) {
        return(list(state = state, ok = FALSE))
      }
      lnew <- upd$R
    } else {
      state$A <- state$A + tcrossprod(r_in) * fold$gamma_inv
    }
    state$T <- state$T + r_in %*% m_in
    state$U <- state$U + crossprod(m_in)
    state$Q <- state$Q + r_in %*% x_in_t
    state$C <- state$C + crossprod(x_in_t, m_in)
  }

  if (changed) {
    if (isTRUE(use_rank_update)) {
      state$L <- lnew
    } else {
      lnew <- try(chol(state$A), silent = TRUE)
      if (inherits(lnew, "try-error")) {
        lnew <- try(chol(state$A + diag(1e-8, nrow(state$A))), silent = TRUE)
        if (inherits(lnew, "try-error")) {
          return(list(state = state, ok = FALSE))
        }
      }
      state$L <- lnew
    }
  }

  list(state = state, ok = TRUE)
}

#' @keywords internal
#' @noRd
.predict_dual_lda_state <- function(state, fold) {
  z <- .solve_from_chol_upper(state$L, state$T)

  lin <- state$C * fold$gamma_inv - crossprod(state$Q, z) * fold$gamma_inv2
  quad <- diag(state$U) * fold$gamma_inv - colSums(state$T * z) * fold$gamma_inv2

  scores <- sweep(lin, 2, 0.5 * quad, "-")
  scores <- sweep(scores, 2, fold$log_priors, "+")

  probs <- .dual_lda_row_softmax(scores)
  colnames(probs) <- fold$classes
  list(test_idx = fold$test_idx, probs = probs)
}

#' @keywords internal
#' @noRd
.combine_dual_lda_folds <- function(fold_preds, y_all, classes, design) {
  test_idx <- sort(unique(unlist(lapply(fold_preds, `[[`, "test_idx"))))
  fold_row_index <- lapply(fold_preds, function(fp) match(fp$test_idx, test_idx))
  prob <- .combine_dual_lda_fold_probs(
    fold_preds = fold_preds,
    fold_row_index = fold_row_index,
    n_rows = length(test_idx),
    n_classes = length(classes)
  )
  colnames(prob) <- classes

  obs <- factor(y_all[test_idx], levels = classes)
  pred <- factor(classes[max.col(prob)], levels = classes)
  classification_result(
    observed = obs,
    predicted = pred,
    probs = prob,
    testind = test_idx,
    test_design = design$test_design,
    predictor = NULL
  )
}

#' @keywords internal
#' @noRd
.validate_dual_lda_inputs <- function(model_spec, gamma) {
  ds <- model_spec$dataset
  if (!inherits(ds, "mvpa_image_dataset")) {
    stop("dual_lda fast path currently supports mvpa_image_dataset only.")
  }
  if (inherits(ds, "mvpa_multibasis_image_dataset")) {
    stop("dual_lda fast path does not currently support multibasis datasets.")
  }

  y_all <- y_train(model_spec)
  if (!is.factor(y_all)) {
    stop("dual_lda fast path requires factor responses.")
  }
  classes <- levels(y_all)
  if (length(classes) < 2L) {
    stop("dual_lda fast path requires at least two classes.")
  }

  gamma <- .dual_lda_gamma(model_spec, gamma = gamma)
  if (!is.finite(gamma) || gamma <= 0) {
    stop("dual_lda fast path requires `gamma > 0`.")
  }

  list(ds = ds, y_all = y_all, classes = classes, gamma = gamma)
}

#' @keywords internal
#' @noRd
run_searchlight_dual_lda_fast <- function(model_spec, radius, incremental = TRUE, gamma = NULL,
                                          verbose = FALSE, scan_chunk = NULL) {
  validated <- .validate_dual_lda_inputs(model_spec, gamma)
  ds <- validated$ds
  y_all <- validated$y_all
  classes <- validated$classes
  gamma <- validated$gamma

  space_obj <- resolve_volume_space(ds)
  dims <- spatial_dim_shape(space_obj)
  spacing <- neuroim2::spacing(space_obj)[1:3]

  mask_indices <- ds$mask_indices
  if (is.null(mask_indices)) {
    mask_indices <- compute_mask_indices(ds$mask)
  }
  if (length(mask_indices) == 0L) {
    return(empty_searchlight_result(ds))
  }

  centers <- get_center_ids(ds)
  centers <- intersect(centers, mask_indices)
  if (length(centers) == 0L) {
    return(empty_searchlight_result(ds))
  }

  n_space <- prod(dims)
  mask_active <- logical(n_space)
  mask_active[mask_indices] <- TRUE

  col_lookup <- integer(n_space)
  col_lookup[mask_indices] <- seq_along(mask_indices)

  center_coords <- neuroim2::index_to_grid(space_obj, centers)
  center_order <- .snake_center_order(center_coords)
  ord_ids <- centers[center_order]
  ord_coords <- center_coords[center_order, , drop = FALSE]

  offsets <- .searchlight_offsets(radius, spacing = spacing)
  if (nrow(offsets) == 0L) {
    return(empty_searchlight_result(ds))
  }
  delta_map <- .searchlight_enter_leave(offsets)
  boundary_tables <- .precompute_boundary_tables(delta_map, dims)

  x_all <- as.matrix(neuroim2::series(ds$train_data, mask_indices))
  if (nrow(x_all) != length(y_all)) {
    stop("dual_lda fast path: mismatch between train rows and y_train length.")
  }

  folds <- generate_folds(
    model_spec$crossval,
    tibble::tibble(.row = seq_len(nrow(x_all))),
    y_all
  )
  if (nrow(folds) == 0L) {
    stop("dual_lda fast path: cross-validation produced zero folds.")
  }

  fold_list <- vector("list", nrow(folds))
  for (i in seq_len(nrow(folds))) {
    tr_idx <- as.integer(folds$train[[i]])
    te_idx <- as.integer(folds$test[[i]])
    fold_list[[i]] <- .prepare_dual_lda_fold(
      x_all = x_all,
      y_all = y_all,
      train_idx = tr_idx,
      test_idx = te_idx,
      classes = classes,
      gamma = gamma
    )
  }

  states <- lapply(fold_list, function(fold) {
    .init_dual_lda_state(fold, ord_coords[1, ], offsets, dims, mask_active, col_lookup)
  })

  rank_chunk <- if (is.null(scan_chunk)) .dual_lda_rank_chunk_size() else as.integer(scan_chunk)[1]
  if (!is.finite(rank_chunk) || rank_chunk < 1L) {
    rank_chunk <- 1L
  }
  has_chunk_native <- .dual_lda_native_symbol_loaded("rmvpa_dual_lda_chunk_update_predict")

  rank_pref <- .dual_lda_rank_update_preference()
  if (!is.null(rank_pref)) {
    use_rank_update <- isTRUE(rank_pref)
  } else {
    probe_step <- .dual_lda_first_step_cols(
      ord_coords = ord_coords,
      ord_ids = ord_ids,
      boundary_tables = boundary_tables,
      dims = dims,
      mask_active = mask_active,
      col_lookup = col_lookup
    )

    if (is.null(probe_step)) {
      use_rank_update <- FALSE
    } else {
      use_rank_update <- .dual_lda_rank_update_benchmark(
        state = states[[1]],
        fold = fold_list[[1]],
        out_cols = probe_step$out_cols,
        in_cols = probe_step$in_cols,
        radius = radius,
        chunk_size = rank_chunk,
        trials = 8L,
        verbose = verbose
      )
    }
  }

  use_rank_chunk <- isTRUE(incremental) &&
    isTRUE(use_rank_update) &&
    rank_chunk > 1L &&
    isTRUE(has_chunk_native)
  if (isTRUE(verbose) && isTRUE(use_rank_chunk)) {
    futile.logger::flog.info("dual_lda rank-update chunk path enabled (chunk=%d)", rank_chunk)
  }

  metric_mode <- .dual_lda_fast_metric_mode(model_spec, classes)
  if (isTRUE(metric_mode$enabled)) {
    test_idx_union <- sort(unique(unlist(lapply(fold_list, `[[`, "test_idx"))))
    fold_row_index <- lapply(fold_list, function(fold) match(fold$test_idx, test_idx_union))
    obs_union <- factor(y_all[test_idx_union], levels = classes)
  }

  perf_rows <- vector("list", length(centers))
  bad_rows <- list()

  score_center <- function(pos, fold_preds) {
    out_idx <- center_order[pos]
    center_id <- centers[out_idx]

    perf <- try({
      if (isTRUE(metric_mode$enabled)) {
        prob <- .combine_dual_lda_fold_probs(
          fold_preds = fold_preds,
          fold_row_index = fold_row_index,
          n_rows = length(test_idx_union),
          n_classes = length(classes)
        )
        colnames(prob) <- classes
        .dual_lda_metric_with_splits(
          observed = obs_union,
          probs = prob,
          test_idx = test_idx_union,
          classes = classes,
          kind = metric_mode$kind,
          split_groups = metric_mode$split_groups,
          class_metrics = metric_mode$class_metrics
        )
      } else {
        cres <- .combine_dual_lda_folds(fold_preds, y_all, classes, model_spec$design)
        compute_performance(model_spec, cres)
      }
    }, silent = TRUE)

    if (inherits(perf, "try-error")) {
      return(list(
        success = FALSE,
        bad_row = tibble::tibble(
          id = center_id,
          error = TRUE,
          error_message = conditionMessage(attr(perf, "condition"))
        )
      ))
    }

    list(success = TRUE, out_idx = out_idx, perf = perf)
  }

  n_ord <- length(ord_ids)
  pos <- 1L
  while (pos <= n_ord) {
    if (isTRUE(verbose) && (pos %% 250L == 0L || pos == n_ord)) {
      futile.logger::flog.info("dual_lda fast path: processed %d/%d centers", pos, n_ord)
    }

    if (pos == 1L) {
      fold_preds <- lapply(seq_along(fold_list), function(f) {
        .predict_dual_lda_state(states[[f]], fold_list[[f]])
      })
      sc <- score_center(pos, fold_preds)
      if (sc$success) {
        perf_rows[[sc$out_idx]] <- sc$perf
      } else {
        bad_rows[[length(bad_rows) + 1L]] <- sc$bad_row
      }
      pos <- pos + 1L
      next
    }

    chunk_done <- FALSE
    if (isTRUE(use_rank_chunk)) {
      pos_last <- min(n_ord, pos + rank_chunk - 1L)
      step_cols_seq <- list()
      step_positions <- integer(0L)

      for (pidx in pos:pos_last) {
        prev_coord <- ord_coords[pidx - 1L, ]
        curr_coord <- ord_coords[pidx, ]
        dx <- as.integer(curr_coord - prev_coord)
        if (!.is_unit_axis_step(dx)) {
          break
        }

        key <- .dir_key(dx)
        step_cols_seq[[length(step_cols_seq) + 1L]] <- .step_boundary_cols_fast(
          prev_coord = prev_coord,
          prev_id = ord_ids[pidx - 1L],
          curr_coord = curr_coord,
          curr_id = ord_ids[pidx],
          step_table = boundary_tables[[key]],
          dims = dims,
          mask_active = mask_active,
          col_lookup = col_lookup
        )
        step_positions <- c(step_positions, pidx)
      }

      if (length(step_cols_seq) > 0L) {
        chunk_prob <- vector("list", length(fold_list))
        chunk_ok <- TRUE

        for (f in seq_along(fold_list)) {
          chunk_res <- .dual_lda_chunk_update_predict_native(
            state = states[[f]],
            fold = fold_list[[f]],
            step_cols_seq = step_cols_seq
          )
          if (is.null(chunk_res)) {
            chunk_ok <- FALSE
            break
          }
          states[[f]] <- chunk_res$state
          chunk_prob[[f]] <- chunk_res$probs
        }

        if (!isTRUE(chunk_ok)) {
          # Recover to the pre-chunk center and continue via the scalar fallback.
          for (f in seq_along(fold_list)) {
            states[[f]] <- .init_dual_lda_state(
              fold_list[[f]], ord_coords[pos - 1L, ], offsets, dims, mask_active, col_lookup
            )
          }
        } else {
          for (s in seq_along(step_positions)) {
            pidx <- step_positions[s]
            fold_preds <- lapply(seq_along(fold_list), function(f) {
              list(
                test_idx = fold_list[[f]]$test_idx,
                probs = chunk_prob[[f]][[s]]
              )
            })
            sc <- score_center(pidx, fold_preds)
            if (sc$success) {
              perf_rows[[sc$out_idx]] <- sc$perf
            } else {
              bad_rows[[length(bad_rows) + 1L]] <- sc$bad_row
            }
          }
          pos <- step_positions[length(step_positions)] + 1L
          chunk_done <- TRUE
        }
      }
    }

    if (isTRUE(chunk_done)) {
      next
    }

    prev_coord <- ord_coords[pos - 1L, ]
    curr_coord <- ord_coords[pos, ]
    dx <- as.integer(curr_coord - prev_coord)
    can_update <- isTRUE(incremental) && .is_unit_axis_step(dx)

    if (can_update) {
      key <- .dir_key(dx)
      step_cols <- .step_boundary_cols_fast(
        prev_coord = prev_coord,
        prev_id = ord_ids[pos - 1L],
        curr_coord = curr_coord,
        curr_id = ord_ids[pos],
        step_table = boundary_tables[[key]],
        dims = dims,
        mask_active = mask_active,
        col_lookup = col_lookup
      )
    }

    for (f in seq_along(fold_list)) {
      if (can_update) {
        upd <- .update_dual_lda_state_neighbor(
          state = states[[f]],
          fold = fold_list[[f]],
          out_cols = step_cols$out_cols,
          in_cols = step_cols$in_cols,
          use_rank_update = use_rank_update
        )
        if (isTRUE(upd$ok)) {
          states[[f]] <- upd$state
        } else {
          states[[f]] <- .init_dual_lda_state(
            fold_list[[f]], curr_coord, offsets, dims, mask_active, col_lookup
          )
        }
      } else {
        states[[f]] <- .init_dual_lda_state(
          fold_list[[f]], curr_coord, offsets, dims, mask_active, col_lookup
        )
      }
    }

    fold_preds <- lapply(seq_along(fold_list), function(f) {
      .predict_dual_lda_state(states[[f]], fold_list[[f]])
    })
    sc <- score_center(pos, fold_preds)
    if (sc$success) {
      perf_rows[[sc$out_idx]] <- sc$perf
    } else {
      bad_rows[[length(bad_rows) + 1L]] <- sc$bad_row
    }
    pos <- pos + 1L
  }

  good_idx <- which(!vapply(perf_rows, is.null, logical(1)))
  if (length(good_idx) == 0L) {
    out <- empty_searchlight_result(ds)
    if (length(bad_rows) > 0L) {
      attr(out, "bad_results") <- dplyr::bind_rows(bad_rows)
    }
    return(out)
  }

  metric_names <- names(perf_rows[[good_idx[1]]])
  perf_mat <- matrix(NA_real_, nrow = length(good_idx), ncol = length(metric_names))
  colnames(perf_mat) <- metric_names
  for (i in seq_along(good_idx)) {
    vals <- perf_rows[[good_idx[i]]]
    if (!identical(names(vals), metric_names)) {
      stop("dual_lda fast path: inconsistent performance metric names across centers.")
    }
    perf_mat[i, ] <- as.numeric(vals)
  }

  out <- wrap_out(perf_mat, ds, ids = centers[good_idx])
  attr(out, "bad_results") <- if (length(bad_rows) > 0L) dplyr::bind_rows(bad_rows) else tibble::tibble()
  out
}

#' @keywords internal
#' @noRd
.dual_lda_sampled_combiner <- function(combiner) {
  .engine_sampled_combiner(combiner, "dual_lda_fast")
}

#' @keywords internal
#' @noRd
.dual_lda_extract_roi_indices <- .engine_extract_roi_indices

#' @keywords internal
#' @noRd
.dual_lda_extract_roi_ids <- .engine_extract_roi_ids

#' @keywords internal
#' @noRd
.init_dual_lda_state_cols <- function(fold, cols) {
  n <- nrow(fold$R)
  k <- ncol(fold$M)
  m <- ncol(fold$Xtest_t)

  if (length(cols) == 0L) {
    return(list(
      L = diag(1, n),
      T = matrix(0, n, k),
      U = matrix(0, k, k),
      Q = matrix(0, n, m),
      C = matrix(0, m, k)
    ))
  }

  r_sub <- fold$R[, cols, drop = FALSE]
  m_sub <- fold$M[cols, , drop = FALSE]
  x_sub_t <- fold$Xtest_t[cols, , drop = FALSE]

  gmat <- tcrossprod(r_sub)
  tmat <- r_sub %*% m_sub
  umat <- crossprod(m_sub)
  qmat <- r_sub %*% x_sub_t
  cmat <- crossprod(x_sub_t, m_sub)

  a <- diag(1, n) + gmat * fold$gamma_inv
  lmat <- try(chol(a), silent = TRUE)
  if (inherits(lmat, "try-error")) {
    lmat <- try(chol(a + diag(1e-8, n)), silent = TRUE)
    if (inherits(lmat, "try-error")) {
      return(NULL)
    }
  }

  list(
    L = lmat,
    T = tmat,
    U = umat,
    Q = qmat,
    C = cmat
  )
}

#' @keywords internal
#' @noRd
run_searchlight_dual_lda_sampled_fast <- function(model_spec,
                                                  radius,
                                                  method = c("randomized", "resampled"),
                                                  niter = 4L,
                                                  combiner = "average",
                                                  drop_probs = FALSE,
                                                  fail_fast = FALSE,
                                                  backend = c("default", "shard", "auto"),
                                                  gamma = NULL,
                                                  verbose = FALSE,
                                                  ...) {
  method <- match.arg(method)
  combiner_fun <- .dual_lda_sampled_combiner(combiner)
  niter <- as.integer(niter)[1]
  if (!is.finite(niter) || is.na(niter) || niter < 1L) {
    stop("dual_lda sampled fast path requires niter >= 1.", call. = FALSE)
  }
  if (isTRUE(drop_probs)) {
    warning("drop_probs is ignored in dual_lda sampled fast path.", call. = FALSE)
  }
  if (isTRUE(fail_fast)) {
    warning("fail_fast is ignored in dual_lda sampled fast path.", call. = FALSE)
  }
  backend <- match.arg(backend)

  validated <- .validate_dual_lda_inputs(model_spec, gamma)
  ds <- validated$ds
  y_all <- validated$y_all
  classes <- validated$classes
  gamma <- validated$gamma

  mask_indices <- ds$mask_indices
  if (is.null(mask_indices)) {
    mask_indices <- compute_mask_indices(ds$mask)
  }
  if (length(mask_indices) == 0L) {
    return(empty_searchlight_result(ds))
  }

  space_obj <- resolve_volume_space(ds)
  dims <- spatial_dim_shape(space_obj)
  n_space <- prod(dims)
  col_lookup <- integer(n_space)
  col_lookup[mask_indices] <- seq_along(mask_indices)

  x_all <- as.matrix(neuroim2::series(ds$train_data, mask_indices))
  if (nrow(x_all) != length(y_all)) {
    stop("dual_lda sampled fast path: mismatch between training rows and y_train length.")
  }

  folds <- generate_folds(
    model_spec$crossval,
    tibble::tibble(.row = seq_len(nrow(x_all))),
    y_all
  )
  if (nrow(folds) == 0L) {
    stop("dual_lda sampled fast path: cross-validation produced zero folds.")
  }

  fold_list <- vector("list", nrow(folds))
  n_valid_folds <- 0L
  bad_rows <- list()

  for (i in seq_len(nrow(folds))) {
    tr_idx <- as.integer(folds$train[[i]])
    te_idx <- as.integer(folds$test[[i]])

    fold <- try(
      .prepare_dual_lda_fold(
        x_all = x_all,
        y_all = y_all,
        train_idx = tr_idx,
        test_idx = te_idx,
        classes = classes,
        gamma = gamma
      ),
      silent = TRUE
    )
    if (inherits(fold, "try-error")) {
      bad_rows[[length(bad_rows) + 1L]] <- tibble::tibble(
        id = NA_integer_,
        error = TRUE,
        error_message = conditionMessage(attr(fold, "condition"))
      )
      next
    }

    n_valid_folds <- n_valid_folds + 1L
    fold_list[[n_valid_folds]] <- fold
  }

  if (n_valid_folds == 0L) {
    out <- empty_searchlight_result(ds)
    attr(out, "bad_results") <- if (length(bad_rows) > 0L) dplyr::bind_rows(bad_rows) else tibble::tibble()
    return(out)
  }
  fold_list <- fold_list[seq_len(n_valid_folds)]

  metric_mode <- .dual_lda_fast_metric_mode(model_spec, classes)
  test_idx_union <- NULL
  fold_row_index <- NULL
  obs_union <- NULL
  if (isTRUE(metric_mode$enabled)) {
    test_idx_union <- sort(unique(unlist(lapply(fold_list, `[[`, "test_idx"))))
    fold_row_index <- lapply(fold_list, function(fold) match(fold$test_idx, test_idx_union))
    obs_union <- factor(y_all[test_idx_union], levels = classes)
  }

  collect_good <- list()
  collect_n <- 0L
  iter_seq <- if (identical(method, "randomized")) seq_len(niter) else 1L

  for (iter in iter_seq) {
    slight <- if (identical(method, "randomized")) {
      get_searchlight(ds, type = "randomized", radius = radius)
    } else {
      get_searchlight(ds, type = "resampled", radius = radius, iter = niter)
    }
    if (length(slight) == 0L) {
      next
    }

    roi_indices <- .dual_lda_extract_roi_indices(slight)
    roi_ids <- .dual_lda_extract_roi_ids(slight)
    roi_cols <- lapply(roi_indices, function(ids) {
      cols <- col_lookup[ids]
      unique(cols[cols > 0L])
    })

    perf_rows <- vector("list", length(roi_cols))
    for (r in seq_along(roi_cols)) {
      fold_preds <- vector("list", length(fold_list))
      pred_ok <- TRUE

      for (f in seq_along(fold_list)) {
        state <- .init_dual_lda_state_cols(fold_list[[f]], roi_cols[[r]])
        if (is.null(state)) {
          pred_ok <- FALSE
          break
        }
        pred <- try(.predict_dual_lda_state(state, fold_list[[f]]), silent = TRUE)
        if (inherits(pred, "try-error")) {
          pred_ok <- FALSE
          break
        }
        fold_preds[[f]] <- pred
      }

      if (!isTRUE(pred_ok)) {
        bad_rows[[length(bad_rows) + 1L]] <- tibble::tibble(
          id = roi_ids[r],
          error = TRUE,
          error_message = "dual_lda sampled fast path: failed to compute fold predictions."
        )
        next
      }

      perf <- try({
        if (isTRUE(metric_mode$enabled)) {
          prob <- .combine_dual_lda_fold_probs(
            fold_preds = fold_preds,
            fold_row_index = fold_row_index,
            n_rows = length(test_idx_union),
            n_classes = length(classes)
          )
          colnames(prob) <- classes
          .dual_lda_metric_with_splits(
            observed = obs_union,
            probs = prob,
            test_idx = test_idx_union,
            classes = classes,
            kind = metric_mode$kind,
            split_groups = metric_mode$split_groups,
            class_metrics = metric_mode$class_metrics
          )
        } else {
          cres <- .combine_dual_lda_folds(fold_preds, y_all, classes, model_spec$design)
          compute_performance(model_spec, cres)
        }
      }, silent = TRUE)

      if (inherits(perf, "try-error")) {
        bad_rows[[length(bad_rows) + 1L]] <- tibble::tibble(
          id = roi_ids[r],
          error = TRUE,
          error_message = conditionMessage(attr(perf, "condition"))
        )
        next
      }
      perf_rows[[r]] <- perf
    }

    good_idx <- which(!vapply(perf_rows, is.null, logical(1)))
    if (length(good_idx) > 0L) {
      perf_list <- lapply(good_idx, function(i) {
        vals <- as.numeric(perf_rows[[i]])
        names(vals) <- names(perf_rows[[i]])
        vals
      })
      collect_n <- collect_n + 1L
      collect_good[[collect_n]] <- tibble::tibble(
        result = rep(list(NULL), length(good_idx)),
        indices = roi_indices[good_idx],
        performance = perf_list,
        id = roi_ids[good_idx],
        error = FALSE,
        error_message = "~"
      )
    }

    if (identical(method, "resampled")) {
      break
    }
  }

  if (collect_n == 0L) {
    out <- empty_searchlight_result(ds)
    attr(out, "bad_results") <- if (length(bad_rows) > 0L) dplyr::bind_rows(bad_rows) else tibble::tibble()
    return(out)
  }

  good_results <- dplyr::bind_rows(collect_good)
  bad_results <- if (length(bad_rows) > 0L) dplyr::bind_rows(bad_rows) else tibble::tibble()
  out <- combiner_fun(model_spec, good_results, bad_results)
  attr(out, "bad_results") <- bad_results
  attr(out, "searchlight_engine") <- "dual_lda_fast"
  out
}
