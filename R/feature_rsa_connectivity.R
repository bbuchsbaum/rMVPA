#' Extract Per-ROI Predicted and Observed RDM Vectors from Feature RSA Results
#'
#' Convenience helper to pull the compact lower-triangle predicted (and
#' optionally observed) RDM vectors stored by
#' \code{feature_rsa_model(..., return_rdm_vectors = TRUE)} from a
#' \code{regional_mvpa_result}. Results can come either from in-memory
#' \code{$fits} or from file-backed batches written via
#' \code{run_regional(..., save_rdm_vectors_dir = ...)}.
#'
#' @param x A \code{regional_mvpa_result} returned by \code{run_regional()} for a
#'   \code{feature_rsa_model}, or a tibble/data frame with columns
#'   \code{roinum} and \code{rdm_vec}. A \code{regional_mvpa_result} may store
#'   vectors either in-memory or on disk in \code{$rdm_batch_dir}.
#'
#' @return A tibble with one row per ROI and columns:
#'   \describe{
#'     \item{roinum}{ROI id.}
#'     \item{n_obs}{Number of observations contributing to the vector.}
#'     \item{observation_index}{List-column of observation ordering used for the
#'       predicted RDM.}
#'     \item{rdm_vec}{List-column containing the lower-triangle predicted RDM
#'       vector for that ROI.}
#'     \item{observed_rdm_vec}{List-column containing the lower-triangle
#'       observed RDM vector for that ROI (if available).}
#'   }
#'
#' @examples
#' \dontrun{
#' res <- run_regional(
#'   feature_rsa_model(dataset, design, method = "pls", return_rdm_vectors = TRUE),
#'   region_mask
#' )
#' vecs <- feature_rsa_rdm_vectors(res)
#' }
#' @export
feature_rsa_rdm_vectors <- function(x) {
  if (is.data.frame(x) && all(c("roinum", "rdm_vec") %in% names(x))) {
    return(tibble::as_tibble(x))
  }

  if (!inherits(x, "regional_mvpa_result")) {
    stop("feature_rsa_rdm_vectors: pass a regional_mvpa_result or a tibble with columns `roinum` and `rdm_vec`.")
  }

  if (!is.null(x$rdm_batch_dir)) {
    files <- .feature_rsa_batch_files(x$rdm_batch_dir)
    if (!length(files)) {
      stop("feature_rsa_rdm_vectors: no RDM batch files found in `rdm_batch_dir`.")
    }
    rows <- lapply(files, readRDS)
    rows <- Filter(function(tbl) is.data.frame(tbl) && nrow(tbl) > 0L, rows)
    if (!length(rows)) {
      stop("feature_rsa_rdm_vectors: no predicted RDM vectors found in `rdm_batch_dir`.")
    }
    return(dplyr::bind_rows(rows))
  }

  fits <- x$fits
  if (is.null(fits) || !length(fits)) {
    stop("feature_rsa_rdm_vectors: no ROI diagnostics found; re-run feature_rsa_model(..., return_rdm_vectors=TRUE).")
  }

  roi_ids <- seq_along(fits)
  if (!is.null(x$performance_table) &&
      is.data.frame(x$performance_table) &&
      "roinum" %in% names(x$performance_table) &&
      nrow(x$performance_table) == length(fits)) {
    roi_ids <- x$performance_table$roinum
  }

  rows <- lapply(seq_along(fits), function(i) {
    pred <- fits[[i]]
    vec <- tryCatch(pred$predicted_rdm_vec, error = function(...) NULL)
    if (is.null(vec)) {
      return(NULL)
    }
    obs_vec <- tryCatch(pred$observed_rdm_vec, error = function(...) NULL)
    tibble::tibble(
      roinum = roi_ids[[i]],
      n_obs = as.integer(if (is.null(pred$n_obs)) NA_integer_ else pred$n_obs),
      observation_index = list(pred$observation_index),
      rdm_vec = list(as.numeric(vec)),
      observed_rdm_vec = list(if (!is.null(obs_vec)) as.numeric(obs_vec) else NULL)
    )
  })

  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    stop("feature_rsa_rdm_vectors: no predicted RDM vectors found; re-run feature_rsa_model(..., return_rdm_vectors=TRUE).")
  }

  dplyr::bind_rows(rows)
}

.feature_rsa_sparsify_connectivity <- function(mat, keep = 1, absolute = FALSE) {
  if (!is.numeric(keep) || length(keep) != 1L || !is.finite(keep) || keep <= 0 || keep > 1) {
    stop("feature_rsa_connectivity: `keep` must be a scalar in (0, 1].")
  }
  if (keep >= 1 || nrow(mat) < 2L) {
    return(mat)
  }

  out <- mat
  tri <- upper.tri(out)
  vals <- out[tri]
  ok <- is.finite(vals)

  if (!any(ok)) {
    return(out)
  }

  scores <- if (isTRUE(absolute)) abs(vals[ok]) else vals[ok]
  n_keep <- max(1L, ceiling(keep * length(scores)))
  cutoff <- sort(scores, decreasing = TRUE)[n_keep]

  drop_mask <- rep(FALSE, length(vals))
  drop_mask[ok] <- scores < cutoff
  vals[drop_mask] <- 0
  out[tri] <- vals
  out[lower.tri(out)] <- t(out)[lower.tri(out)]
  diag(out) <- 1
  out
}

.feature_rsa_validate_obs_order <- function(obs_idx, fn_label) {
  non_null_obs_idx <- Filter(Negate(is.null), obs_idx)
  if (length(non_null_obs_idx) <= 1L) {
    return(invisible(NULL))
  }

  ref_idx <- non_null_obs_idx[[1]]
  same_idx <- vapply(non_null_obs_idx[-1], identical, logical(1), y = ref_idx)
  if (!all(same_idx)) {
    stop(fn_label, ": ROI RDM vectors do not share the same observation ordering.")
  }

  invisible(NULL)
}

.feature_rsa_double_center <- function(mat) {
  if (!is.matrix(mat)) {
    mat <- as.matrix(mat)
  }

  row_mean <- apply(mat, 1L, function(x) {
    if (any(is.finite(x))) mean(x, na.rm = TRUE) else NA_real_
  })
  col_mean <- apply(mat, 2L, function(x) {
    if (any(is.finite(x))) mean(x, na.rm = TRUE) else NA_real_
  })
  grand_mean <- if (any(is.finite(mat))) mean(mat, na.rm = TRUE) else NA_real_

  centered <- sweep(sweep(mat, 1L, row_mean, "-"), 2L, col_mean, "-") + grand_mean

  list(
    matrix = centered,
    source_offset = row_mean - grand_mean,
    target_offset = col_mean - grand_mean,
    grand_mean = grand_mean
  )
}

.feature_rsa_residualize_vector <- function(y, x) {
  y <- as.numeric(y)
  x <- as.numeric(x)
  out <- rep(NA_real_, length(y))
  ok <- is.finite(y) & is.finite(x)

  if (sum(ok) < 2L) {
    return(out)
  }

  x_ok <- x[ok]
  y_ok <- y[ok]
  x_sd <- stats::sd(x_ok)

  if (!is.finite(x_sd) || x_sd <= 1e-12) {
    out[ok] <- y_ok - mean(y_ok)
    return(out)
  }

  fit <- stats::lm.fit(x = cbind(1, x_ok), y = y_ok)
  out[ok] <- fit$residuals
  out
}

.feature_rsa_residualize_columns <- function(mat, ref_vec) {
  if (!is.matrix(mat)) {
    mat <- as.matrix(mat)
  }
  if (ncol(mat) < 1L) {
    return(mat)
  }

  out <- vapply(
    seq_len(ncol(mat)),
    function(i) .feature_rsa_residualize_vector(mat[, i], ref_vec),
    numeric(nrow(mat))
  )

  if (!is.matrix(out)) {
    out <- matrix(out, nrow = nrow(mat), ncol = ncol(mat))
  }

  dimnames(out) <- dimnames(mat)
  out
}

.feature_rsa_get_pred_vec <- function(x) {
  tryCatch(x$predicted_rdm_vec, error = function(...) NULL)
}

.feature_rsa_get_obs_vec <- function(x) {
  tryCatch(x$observed_rdm_vec, error = function(...) NULL)
}

.feature_rsa_get_obs_idx <- function(x) {
  tryCatch(x$observation_index, error = function(...) NULL)
}

.feature_rsa_bind_numeric_columns <- function(vecs) {
  n_cols <- length(vecs)
  if (n_cols == 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = 0L))
  }

  mat <- do.call(cbind, vecs)
  if (!is.matrix(mat)) {
    mat <- matrix(mat, ncol = n_cols)
  }
  storage.mode(mat) <- "double"
  mat
}

.feature_rsa_batch_manifest <- function(dir) {
  manifest_file <- file.path(dir, "manifest.rds")
  if (!file.exists(manifest_file)) {
    return(NULL)
  }
  readRDS(manifest_file)
}

.feature_rsa_batch_files <- function(dir) {
  manifest <- .feature_rsa_batch_manifest(dir)
  if (!is.null(manifest) && length(manifest$batches) > 0L) {
    return(file.path(dir, vapply(manifest$batches, `[[`, character(1), "file")))
  }

  list.files(
    dir,
    pattern = "^batch_[0-9]+[.]rds$",
    full.names = TRUE
  ) |> sort()
}

.feature_rsa_batch_index <- function(dir) {
  manifest <- .feature_rsa_batch_manifest(dir)
  if (is.null(manifest) || !length(manifest$batches)) {
    files <- .feature_rsa_batch_files(dir)
    if (!length(files)) {
      return(list(
        files = character(0),
        starts = integer(0),
        ends = integer(0),
        roi_ids = character(0),
        pred_lengths = integer(0),
        obs_lengths = integer(0),
        has_obs = logical(0),
        observation_index = list()
      ))
    }

    batch_rows <- lapply(files, function(path) {
      tbl <- readRDS(path)
      list(
        n_rows = nrow(tbl),
        roi_ids = as.character(tbl$roinum),
        pred_lengths = vapply(tbl$rdm_vec, length, integer(1)),
        obs_lengths = vapply(tbl$observed_rdm_vec, function(v) if (is.null(v)) NA_integer_ else length(v), integer(1)),
        has_obs = vapply(tbl$observed_rdm_vec, function(v) !is.null(v), logical(1)),
        observation_index = tbl$observation_index
      )
    })
  } else {
    files <- file.path(dir, vapply(manifest$batches, `[[`, character(1), "file"))
    batch_rows <- lapply(manifest$batches, function(batch) {
      list(
        n_rows = as.integer(batch$n_rows),
        roi_ids = as.character(batch$roi_ids),
        pred_lengths = as.integer(batch$pred_lengths),
        obs_lengths = as.integer(batch$obs_lengths),
        has_obs = as.logical(batch$has_obs),
        observation_index = batch$observation_index
      )
    })
  }

  n_rows <- vapply(batch_rows, `[[`, integer(1), "n_rows")
  starts <- cumsum(c(1L, head(n_rows, -1L)))
  ends <- cumsum(n_rows)

  list(
    files = files,
    starts = starts,
    ends = ends,
    roi_ids = unlist(lapply(batch_rows, `[[`, "roi_ids"), use.names = FALSE),
    pred_lengths = unlist(lapply(batch_rows, `[[`, "pred_lengths"), use.names = FALSE),
    obs_lengths = unlist(lapply(batch_rows, `[[`, "obs_lengths"), use.names = FALSE),
    has_obs = unlist(lapply(batch_rows, `[[`, "has_obs"), use.names = FALSE),
    observation_index = unlist(lapply(batch_rows, `[[`, "observation_index"), recursive = FALSE, use.names = FALSE)
  )
}

.feature_rsa_batch_loader <- function(dir, column, vec_length, fill_missing = TRUE) {
  index <- .feature_rsa_batch_index(dir)
  cache <- new.env(parent = emptyenv())

  function(idx) {
    if (!length(idx)) {
      return(matrix(numeric(0), nrow = vec_length, ncol = 0L))
    }

    idx <- as.integer(idx)
    out <- vector("list", length(idx))
    for (i in seq_along(idx)) {
      batch_id <- which(index$starts <= idx[[i]] & index$ends >= idx[[i]])
      if (!length(batch_id)) {
        stop("feature_rsa batch loader: index out of bounds.")
      }

      batch_key <- as.character(batch_id[[1]])
      if (!exists(batch_key, envir = cache, inherits = FALSE)) {
        assign(batch_key, readRDS(index$files[[batch_id[[1]]]]), envir = cache)
      }
      batch_tbl <- get(batch_key, envir = cache, inherits = FALSE)
      local_idx <- idx[[i]] - index$starts[[batch_id[[1]]]] + 1L
      value <- batch_tbl[[column]][[local_idx]]
      if (is.null(value) && isTRUE(fill_missing)) {
        value <- rep(NA_real_, vec_length)
      }
      out[[i]] <- if (is.null(value)) NULL else as.numeric(value)
    }

    .feature_rsa_bind_numeric_columns(out)
  }
}

.feature_rsa_default_block_size <- function(vec_length, n_roi) {
  target_bytes <- getOption("rMVPA.feature_rsa_block_target_bytes", 128 * 1024^2)
  max_cols <- getOption("rMVPA.feature_rsa_block_max_cols", 64L)

  if (!is.numeric(target_bytes) || length(target_bytes) != 1L || !is.finite(target_bytes) || target_bytes <= 0) {
    target_bytes <- 128 * 1024^2
  }
  if (!is.numeric(max_cols) || length(max_cols) != 1L || !is.finite(max_cols) || max_cols <= 0) {
    max_cols <- 64L
  }

  bytes_per_col <- max(8, as.double(vec_length) * 8)
  block_size <- as.integer(floor(as.double(target_bytes) / bytes_per_col))
  block_size <- max(1L, block_size)
  block_size <- min(block_size, as.integer(max_cols), as.integer(n_roi))
  as.integer(block_size)
}

.feature_rsa_block_slices <- function(n_items, block_size) {
  split(seq_len(n_items), ceiling(seq_len(n_items) / block_size))
}

.feature_rsa_rank_columns <- function(mat) {
  if (!length(mat)) {
    return(mat)
  }

  ranked <- apply(mat, 2L, rank, ties.method = "average", na.last = "keep")
  if (!is.matrix(ranked)) {
    ranked <- matrix(ranked, nrow = nrow(mat), ncol = ncol(mat))
  }
  storage.mode(ranked) <- "double"
  ranked
}

.feature_rsa_log_message <- function(verbose, label, ..., .append_elapsed = NULL) {
  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }

  msg <- sprintf(...)
  if (!is.null(.append_elapsed)) {
    msg <- sprintf("%s [%0.1fs]", msg, .append_elapsed)
  }
  message(label, ": ", msg)
  invisible(NULL)
}

.feature_rsa_block_cor <- function(x_loader,
                                   y_loader = x_loader,
                                   n_x,
                                   n_y = n_x,
                                   vec_length,
                                   method,
                                   use,
                                   symmetric = FALSE,
                                   verbose = FALSE,
                                   label = "feature_rsa_block_cor") {
  block_size <- .feature_rsa_default_block_size(vec_length = vec_length, n_roi = max(n_x, n_y))
  x_slices <- .feature_rsa_block_slices(n_x, block_size)
  y_slices <- .feature_rsa_block_slices(n_y, block_size)
  out <- matrix(NA_real_, nrow = n_x, ncol = n_y)
  cor_method <- if (identical(method, "spearman")) "pearson" else method
  x_cache <- new.env(parent = emptyenv())
  y_cache <- new.env(parent = emptyenv())
  started <- proc.time()[3]
  total_pairs <- if (isTRUE(symmetric)) {
    length(x_slices) * (length(x_slices) + 1L) / 2L
  } else {
    length(x_slices) * length(y_slices)
  }
  progress_every <- max(1L, as.integer(ceiling(total_pairs / 10L)))
  completed_pairs <- 0L

  .feature_rsa_log_message(
    verbose,
    label,
    "starting %s correlation across %d ROI(s) with %d block pair(s) (block_size=%d, vec_length=%d)",
    method,
    max(n_x, n_y),
    total_pairs,
    block_size,
    vec_length
  )

  load_block <- function(loader, idx, cache_env) {
    cache_key <- sprintf("%d:%d", idx[[1]], idx[[length(idx)]])
    if (!exists(cache_key, envir = cache_env, inherits = FALSE)) {
      block <- loader(idx)
      if (identical(method, "spearman")) {
        block <- .feature_rsa_rank_columns(block)
      }
      assign(cache_key, block, envir = cache_env)
    }
    get(cache_key, envir = cache_env, inherits = FALSE)
  }

  for (i in seq_along(x_slices)) {
    idx_i <- x_slices[[i]]
    x_block <- load_block(x_loader, idx_i, x_cache)

    j_seq <- if (isTRUE(symmetric)) seq.int(i, length(y_slices)) else seq_along(y_slices)
    for (j in j_seq) {
      idx_j <- y_slices[[j]]
      y_block <- load_block(y_loader, idx_j, y_cache)

      block_cor <- suppressWarnings(stats::cor(x_block, y_block, method = cor_method, use = use))
      if (!is.matrix(block_cor)) {
        block_cor <- matrix(block_cor, nrow = length(idx_i), ncol = length(idx_j))
      }

      out[idx_i, idx_j] <- block_cor
      if (isTRUE(symmetric) && i != j) {
        out[idx_j, idx_i] <- t(block_cor)
      }

       completed_pairs <- completed_pairs + 1L
       if (isTRUE(verbose) &&
           (completed_pairs == 1L ||
            completed_pairs %% progress_every == 0L ||
            completed_pairs == total_pairs)) {
         .feature_rsa_log_message(
           TRUE,
           label,
           "completed %d/%d block pair(s) (%0.1f%%)",
           completed_pairs,
           total_pairs,
           100 * completed_pairs / total_pairs,
           .append_elapsed = proc.time()[3] - started
         )
       }
    }
  }

  .feature_rsa_log_message(
    verbose,
    label,
    "finished %s correlation across %d ROI(s)",
    method,
    max(n_x, n_y),
    .append_elapsed = proc.time()[3] - started
  )

  out
}

.feature_rsa_row_means <- function(loader,
                                   n_cols,
                                   vec_length,
                                   verbose = FALSE,
                                   label = "feature_rsa_row_means") {
  block_size <- .feature_rsa_default_block_size(vec_length = vec_length, n_roi = n_cols)
  sums <- numeric(vec_length)
  counts <- integer(vec_length)
  slices <- .feature_rsa_block_slices(n_cols, block_size)
  started <- proc.time()[3]

  .feature_rsa_log_message(
    verbose,
    label,
    "starting row-mean accumulation across %d ROI block(s) (block_size=%d, vec_length=%d)",
    length(slices),
    block_size,
    vec_length
  )

  for (i in seq_along(slices)) {
    idx <- slices[[i]]
    block <- loader(idx)
    ok <- is.finite(block)
    block[!ok] <- 0
    sums <- sums + rowSums(block)
    counts <- counts + rowSums(ok)

    if (isTRUE(verbose) &&
        (i == 1L || i %% max(1L, ceiling(length(slices) / 4L)) == 0L || i == length(slices))) {
      .feature_rsa_log_message(
        TRUE,
        label,
        "processed %d/%d ROI block(s) (%0.1f%%)",
        i,
        length(slices),
        100 * i / length(slices),
        .append_elapsed = proc.time()[3] - started
      )
    }
  }

  out <- sums / counts
  out[counts == 0L] <- NA_real_
  .feature_rsa_log_message(verbose, label, "finished row-mean accumulation", .append_elapsed = proc.time()[3] - started)
  out
}

.feature_rsa_prepare_vectors <- function(x, fn_label) {
  if (is.data.frame(x) && all(c("roinum", "rdm_vec") %in% names(x))) {
    vec_tbl <- tibble::as_tibble(x)
    if (!"observed_rdm_vec" %in% names(vec_tbl)) {
      vec_tbl$observed_rdm_vec <- rep(list(NULL), nrow(vec_tbl))
    }
    if (!"observation_index" %in% names(vec_tbl)) {
      vec_tbl$observation_index <- rep(list(NULL), nrow(vec_tbl))
    }

    pred_lengths <- vapply(vec_tbl$rdm_vec, length, integer(1))
    unique_pred_lengths <- unique(pred_lengths)
    if (length(unique_pred_lengths) != 1L) {
      stop(fn_label, ": ROI RDM vectors must all have the same length.")
    }

    obs_lengths <- vapply(
      vec_tbl$observed_rdm_vec,
      function(v) if (is.null(v)) NA_integer_ else length(v),
      integer(1)
    )

    return(list(
      roi_labels = as.character(vec_tbl$roinum),
      observation_index = vec_tbl$observation_index,
      pred_lengths = pred_lengths,
      obs_lengths = obs_lengths,
      vec_length = unique_pred_lengths[[1]],
      has_obs = vapply(vec_tbl$observed_rdm_vec, function(v) !is.null(v), logical(1)),
      get_pred_block = function(idx) {
        .feature_rsa_bind_numeric_columns(lapply(vec_tbl$rdm_vec[idx], as.numeric))
      },
      get_obs_block = function(idx, fill_missing = TRUE) {
        vecs <- lapply(vec_tbl$observed_rdm_vec[idx], function(v) {
          if (is.null(v)) {
            if (isTRUE(fill_missing)) {
              rep(NA_real_, unique_pred_lengths[[1]])
            } else {
              NULL
            }
          } else {
            as.numeric(v)
          }
        })
        .feature_rsa_bind_numeric_columns(vecs)
      }
    ))
  }

  if (!inherits(x, "regional_mvpa_result")) {
    stop(fn_label, ": pass a regional_mvpa_result or a tibble with columns `roinum` and `rdm_vec`.")
  }

  if (!is.null(x$rdm_batch_dir)) {
    batch_index <- .feature_rsa_batch_index(x$rdm_batch_dir)
    if (!length(batch_index$files)) {
      stop(fn_label, ": no RDM batch files found in `rdm_batch_dir`.")
    }
    if (!length(batch_index$pred_lengths)) {
      stop(fn_label, ": no predicted RDM vectors found in `rdm_batch_dir`.")
    }
    if (length(unique(batch_index$pred_lengths)) != 1L) {
      stop(fn_label, ": ROI RDM vectors must all have the same length.")
    }

    pred_loader <- .feature_rsa_batch_loader(
      x$rdm_batch_dir,
      column = "rdm_vec",
      vec_length = batch_index$pred_lengths[[1]],
      fill_missing = FALSE
    )
    obs_loader <- .feature_rsa_batch_loader(
      x$rdm_batch_dir,
      column = "observed_rdm_vec",
      vec_length = batch_index$pred_lengths[[1]],
      fill_missing = TRUE
    )

    return(list(
      roi_labels = batch_index$roi_ids,
      observation_index = batch_index$observation_index,
      pred_lengths = batch_index$pred_lengths,
      obs_lengths = batch_index$obs_lengths,
      vec_length = batch_index$pred_lengths[[1]],
      has_obs = batch_index$has_obs,
      get_pred_block = function(idx) {
        pred_loader(idx)
      },
      get_obs_block = function(idx, fill_missing = TRUE) {
        if (isTRUE(fill_missing)) {
          obs_loader(idx)
        } else {
          .feature_rsa_batch_loader(
            x$rdm_batch_dir,
            column = "observed_rdm_vec",
            vec_length = batch_index$pred_lengths[[1]],
            fill_missing = FALSE
          )(idx)
        }
      }
    ))
  }

  fits <- x$fits
  if (is.null(fits) || !length(fits)) {
    stop(fn_label, ": no ROI diagnostics found; re-run feature_rsa_model(..., return_rdm_vectors=TRUE).")
  }

  roi_ids <- seq_along(fits)
  if (!is.null(x$performance_table) &&
      is.data.frame(x$performance_table) &&
      "roinum" %in% names(x$performance_table) &&
      nrow(x$performance_table) == length(fits)) {
    roi_ids <- x$performance_table$roinum
  }

  has_pred <- vapply(fits, function(pred) !is.null(.feature_rsa_get_pred_vec(pred)), logical(1))
  keep <- which(has_pred)
  if (!length(keep)) {
    stop(fn_label, ": no predicted RDM vectors found; re-run feature_rsa_model(..., return_rdm_vectors=TRUE).")
  }

  kept_fits <- fits[keep]
  pred_lengths <- vapply(kept_fits, function(pred) length(.feature_rsa_get_pred_vec(pred)), integer(1))
  unique_pred_lengths <- unique(pred_lengths)
  if (length(unique_pred_lengths) != 1L) {
    stop(fn_label, ": ROI RDM vectors must all have the same length.")
  }

  obs_lengths <- vapply(
    kept_fits,
    function(pred) {
      obs_vec <- .feature_rsa_get_obs_vec(pred)
      if (is.null(obs_vec)) NA_integer_ else length(obs_vec)
    },
    integer(1)
  )

  list(
    roi_labels = as.character(roi_ids[keep]),
    observation_index = lapply(kept_fits, .feature_rsa_get_obs_idx),
    pred_lengths = pred_lengths,
    obs_lengths = obs_lengths,
    vec_length = unique_pred_lengths[[1]],
    has_obs = vapply(kept_fits, function(pred) !is.null(.feature_rsa_get_obs_vec(pred)), logical(1)),
    get_pred_block = function(idx) {
      fit_idx <- keep[idx]
      vecs <- lapply(fit_idx, function(i) as.numeric(.feature_rsa_get_pred_vec(fits[[i]])))
      .feature_rsa_bind_numeric_columns(vecs)
    },
    get_obs_block = function(idx, fill_missing = TRUE) {
      fit_idx <- keep[idx]
      vecs <- lapply(fit_idx, function(i) {
        obs_vec <- .feature_rsa_get_obs_vec(fits[[i]])
        if (is.null(obs_vec)) {
          if (isTRUE(fill_missing)) rep(NA_real_, unique_pred_lengths[[1]]) else NULL
        } else {
          as.numeric(obs_vec)
        }
      })
      .feature_rsa_bind_numeric_columns(vecs)
    }
  )
}

#' Compute ROI-by-ROI Representational Connectivity from Feature RSA Predictions
#'
#' Forms an ROI x ROI similarity matrix by correlating lower-triangle predicted
#' RDM vectors across ROIs. Sparsification, when requested, is applied only to
#' the final ROI x ROI matrix and never to the per-ROI RDM vectors themselves.
#'
#' @param x Either a \code{regional_mvpa_result} produced by
#'   \code{feature_rsa_model(..., return_rdm_vectors=TRUE)} or the tibble
#'   returned by \code{feature_rsa_rdm_vectors()}. Regional results may store
#'   RDM vectors either in-memory or in file-backed batches written by
#'   \code{run_regional(..., save_rdm_vectors_dir = ...)}.
#' @param method Correlation method used across ROI RDM vectors, one of
#'   \code{"spearman"} or \code{"pearson"}.
#' @param keep Proportion of ROI-ROI edges to retain after optional
#'   sparsification. \code{keep = 1} disables sparsification. For example,
#'   \code{keep = 0.1} retains the top 10\% of finite off-diagonal edges.
#' @param absolute Logical; when \code{TRUE}, rank edges by absolute magnitude
#'   during sparsification. Defaults to \code{FALSE}.
#' @param use Missing-value handling passed to \code{\link[stats]{cor}}.
#' @param verbose Logical; if \code{TRUE}, emit block-level progress messages
#'   while connectivity is being computed.
#'
#' @return A symmetric numeric matrix with ROIs in rows/columns.
#'
#' @examples
#' \dontrun{
#' vecs <- feature_rsa_rdm_vectors(res)
#' conn <- feature_rsa_connectivity(vecs, method = "spearman", keep = 0.1)
#' }
#' @export
feature_rsa_connectivity <- function(x,
                                     method = c("spearman", "pearson"),
                                     keep = 1,
                                     absolute = FALSE,
                                     use = "pairwise.complete.obs",
                                     verbose = FALSE) {
  method <- match.arg(method)
  vec_info <- .feature_rsa_prepare_vectors(x, fn_label = "feature_rsa_connectivity")

  .feature_rsa_validate_obs_order(
    vec_info$observation_index,
    fn_label = "feature_rsa_connectivity"
  )

  n_roi <- length(vec_info$roi_labels)
  conn <- .feature_rsa_block_cor(
    x_loader = vec_info$get_pred_block,
    n_x = n_roi,
    vec_length = vec_info$vec_length,
    method = method,
    use = use,
    symmetric = TRUE,
    verbose = verbose,
    label = "feature_rsa_connectivity"
  )

  dimnames(conn) <- list(vec_info$roi_labels, vec_info$roi_labels)
  diag(conn) <- 1

  .feature_rsa_sparsify_connectivity(conn, keep = keep, absolute = absolute)
}

#' Compute Cross-Connectivity: Predicted-Observed ROI x ROI Matrix
#'
#' Builds an asymmetric ROI x ROI matrix where entry (i, j) is the correlation
#' between the predicted RDM vector of ROI i and the observed RDM vector of
#' ROI j.  This captures how well the model-predicted representational geometry
#' in one ROI matches the data-driven geometry in another.
#'
#' @param x Either a \code{regional_mvpa_result} produced by
#'   \code{feature_rsa_model(..., return_rdm_vectors=TRUE)} or the tibble
#'   returned by \code{feature_rsa_rdm_vectors()}. Regional results may store
#'   RDM vectors either in-memory or in file-backed batches written by
#'   \code{run_regional(..., save_rdm_vectors_dir = ...)}.
#' @param method Correlation method, one of \code{"spearman"} or
#'   \code{"pearson"}.
#' @param adjust Optional adjustment for ROI-level source/target offsets. Use
#'   \code{"none"} (default) to return the raw ROI x ROI correlation matrix,
#'   \code{"double_center"} to subtract additive source and target main effects
#'   from that matrix, or \code{"residualize_mean"} to remove the grand-mean
#'   RDM component from predicted and observed ROI vectors before computing the
#'   cross-correlation.
#' @param return_components Logical; if \code{TRUE}, return a list containing
#'   the requested matrix, the raw matrix, the adjusted matrix, and the source
#'   and target offset terms estimated from the raw matrix.
#' @param use Missing-value handling passed to \code{\link[stats]{cor}}.
#' @param verbose Logical; if \code{TRUE}, emit block-level progress messages
#'   while cross-connectivity is being computed.
#'
#' @return By default, a numeric matrix of dimension n_ROI x n_ROI. Rows
#'   correspond to predicted RDM vectors and columns to observed RDM vectors.
#'   The matrix is \emph{not} necessarily symmetric. If
#'   \code{return_components = TRUE}, a list is returned with elements
#'   \code{matrix}, \code{raw_matrix}, \code{adjusted_matrix},
#'   \code{source_offset}, \code{target_offset}, \code{grand_mean},
#'   \code{method}, and \code{adjust}.
#'
#' @examples
#' \dontrun{
#' res <- run_regional(
#'   feature_rsa_model(dataset, design, method = "pls", return_rdm_vectors = TRUE),
#'   region_mask
#' )
#' cross_conn <- feature_rsa_cross_connectivity(res, method = "spearman")
#' cross_dc <- feature_rsa_cross_connectivity(
#'   res,
#'   method = "spearman",
#'   adjust = "double_center"
#' )
#' }
#' @export
feature_rsa_cross_connectivity <- function(x,
                                           method = c("spearman", "pearson"),
                                           adjust = c("none", "double_center", "residualize_mean"),
                                           return_components = FALSE,
                                           use = "pairwise.complete.obs",
                                           verbose = FALSE) {
  method <- match.arg(method)
  adjust <- match.arg(adjust)
  vec_info <- .feature_rsa_prepare_vectors(x, fn_label = "feature_rsa_cross_connectivity")

  ## Check that observed RDM vectors are available
  has_obs <- vec_info$has_obs
  if (!any(has_obs)) {
    stop("feature_rsa_cross_connectivity: no observed RDM vectors found. ",
         "Re-run with feature_rsa_model(..., return_rdm_vectors = TRUE) ",
         "using a version that stores observed RDM vectors.")
  }
  if (!all(has_obs)) {
    warning("feature_rsa_cross_connectivity: some ROIs lack observed RDM vectors; ",
            "their rows/columns will be NA.")
  }

  ## Validate vector lengths
  pred_lengths <- vec_info$pred_lengths
  obs_lengths  <- vec_info$obs_lengths
  all_lengths <- unique(c(pred_lengths, stats::na.omit(obs_lengths)))
  if (length(all_lengths) != 1L) {
    stop("feature_rsa_cross_connectivity: predicted and observed RDM vectors ",
         "must all have the same length.")
  }

  ## Validate observation ordering
  .feature_rsa_validate_obs_order(
    vec_info$observation_index,
    fn_label = "feature_rsa_cross_connectivity"
  )

  n_roi <- length(vec_info$roi_labels)

  ## Cross-correlation: cor(pred_col_i, obs_col_j)
  raw_cross <- .feature_rsa_block_cor(
    x_loader = vec_info$get_pred_block,
    y_loader = function(idx) vec_info$get_obs_block(idx, fill_missing = TRUE),
    n_x = n_roi,
    n_y = n_roi,
    vec_length = vec_info$vec_length,
    method = method,
    use = use,
    symmetric = FALSE,
    verbose = verbose,
    label = "feature_rsa_cross_connectivity [raw]"
  )

  roi_labels <- vec_info$roi_labels
  dimnames(raw_cross) <- list(predicted = roi_labels, observed = roi_labels)

  adjusted_cross <- switch(
    adjust,
    none = raw_cross,
    double_center = {
      centered <- .feature_rsa_double_center(raw_cross)
      mat <- centered$matrix
      dimnames(mat) <- dimnames(raw_cross)
      mat
    },
    residualize_mean = {
      pred_ref <- .feature_rsa_row_means(
        loader = vec_info$get_pred_block,
        n_cols = n_roi,
        vec_length = vec_info$vec_length,
        verbose = verbose,
        label = "feature_rsa_cross_connectivity [pred_ref]"
      )
      obs_ref <- .feature_rsa_row_means(
        loader = function(idx) vec_info$get_obs_block(idx, fill_missing = TRUE),
        n_cols = n_roi,
        vec_length = vec_info$vec_length,
        verbose = verbose,
        label = "feature_rsa_cross_connectivity [obs_ref]"
      )
      pred_ref[!is.finite(pred_ref)] <- NA_real_
      obs_ref[!is.finite(obs_ref)] <- NA_real_

      mat <- .feature_rsa_block_cor(
        x_loader = function(idx) {
          .feature_rsa_residualize_columns(vec_info$get_pred_block(idx), pred_ref)
        },
        y_loader = function(idx) {
          .feature_rsa_residualize_columns(
            vec_info$get_obs_block(idx, fill_missing = TRUE),
            obs_ref
          )
        },
        n_x = n_roi,
        n_y = n_roi,
        vec_length = vec_info$vec_length,
        method = method,
        use = use,
        symmetric = FALSE,
        verbose = verbose,
        label = "feature_rsa_cross_connectivity [residualized]"
      )
      dimnames(mat) <- dimnames(raw_cross)
      mat
    }
  )

  if (!isTRUE(return_components)) {
    return(adjusted_cross)
  }

  raw_components <- .feature_rsa_double_center(raw_cross)

  list(
    matrix = adjusted_cross,
    raw_matrix = raw_cross,
    adjusted_matrix = adjusted_cross,
    source_offset = stats::setNames(raw_components$source_offset, roi_labels),
    target_offset = stats::setNames(raw_components$target_offset, roi_labels),
    grand_mean = raw_components$grand_mean,
    method = method,
    adjust = adjust
  )
}
