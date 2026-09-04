# Public single-domain banded-ridge model ------------------------------------

.brm_active_ids <- function(dataset) {
  ids <- which(dataset$mask > 0)
  if (!length(ids)) stop("banded_ridge_model: dataset mask has no active responses.",
                         call. = FALSE)
  as.integer(ids)
}

.brm_n_observations <- function(dataset) {
  dims <- dim(dataset$train_data)
  if (is.null(dims) || length(dims) < 4L) {
    stop("banded_ridge_model: train_data must be a 4D NeuroVec.", call. = FALSE)
  }
  as.integer(utils::tail(dims, 1L))
}

# Share of alpha selections pinned to the ends of the grid above which the grid
# is reported as saturated. A selection at the largest available alpha means the
# inner objective was still improving when the grid ran out, so the refit is
# under-penalized; at the smallest it is over-penalized. The share is counted
# over both ends together as well as each end alone: under the default
# per-response alpha scope a heterogeneous mask splits its boundary mass across
# the two ends, so a grid that brackets nothing can leave neither end near the
# threshold while almost no response lands inside.
.brm_alpha_saturation_share <- 0.95

# Median outer out-of-fold R2 below which the full model is reported as
# predicting worse than the mean of the held-out data. .banded_ridge_score()
# forms ss_tot from the mean of the rows being scored, so that, not the
# training mean, is the baseline a negative R2 loses to.
.brm_r2_floor <- -0.05

# Scale at which the penalty starts to bite for this design.
#
# Every solver path standardizes each training column to unit standard
# deviation and then scales band b by sqrt(theta_b), so alpha is compared
# against the eigenvalues of Z'Z with Z = Xs diag(sqrt(theta_feature)). Those
# eigenvalues average trace(Z'Z) / rank, which under a uniform theta is
# (n - 1) * mean(D_b) / min(n - 1, p), where D_b counts the non-degenerate
# columns of band b. The anchor therefore grows with the number of training
# rows and with band width, which is why a fixed grid cannot be right for every
# design: a grid that brackets the anchor is in the range where the penalty
# changes the fit, and a grid far below it selects its own ceiling for every
# response.
.brm_alpha_anchor <- function(X, groups) {
  n <- nrow(X)
  if (is.null(n) || n < 2L) return(NA_real_)
  # colVars rather than .brc_column_sds: identical verdict, without
  # materialising a centred and a squared copy of the design
  vars <- matrixStats::colVars(X)
  usable <- is.finite(vars) & vars > 0
  p <- sum(usable)
  if (p < 1L) return(NA_real_)
  per_band <- vapply(split(usable, as.character(groups)),
                     function(z) as.numeric(sum(z)), numeric(1))
  ((n - 1) * mean(per_band)) / min(n - 1, p)
}

# Design-scaled default alpha grid: nine points spanning six decades centred on
# the anchor above.
.brm_default_alphas <- function(X, groups, n_points = 9L, decades = 3,
                                anchor = .brm_alpha_anchor(X, groups)) {
  if (!is.finite(anchor) || anchor <= 0) anchor <- 1
  anchor * 10^seq(-decades, decades, length.out = n_points)
}

.brm_resolve_alphas <- function(alphas, X, groups) {
  if (is.numeric(alphas)) {
    return(list(alphas = alphas, origin = "user", anchor = NA_real_))
  }
  named <- (is.character(alphas) || is.factor(alphas)) &&
    length(alphas) == 1L && !is.na(alphas) &&
    identical(as.character(alphas), "auto")
  if (!named) {
    stop(
      "alphas must be \"auto\" or a numeric vector of non-negative penalties.",
      call. = FALSE
    )
  }
  anchor <- .brm_alpha_anchor(X, groups)
  list(alphas = .brm_default_alphas(X, groups, anchor = anchor),
       origin = "auto", anchor = anchor)
}

# Alphas a response could actually have been given, after the alpha/theta scope
# constraints. Using the whole manifest here would report a fixed alpha scope as
# fully saturated.
.brm_selectable_alphas <- function(candidates,
                                   alpha_scope,
                                   theta_scope,
                                   fixed_alpha = NULL,
                                   fixed_theta = NULL) {
  available <- tryCatch(
    .brt_available_candidates(
      candidates, alpha_scope, theta_scope,
      fixed_alpha = fixed_alpha, fixed_theta = fixed_theta
    ),
    error = function(e) seq_along(candidates$alpha)
  )
  sort(unique(candidates$alpha[available]))
}

# One row summarising where a model's alpha selections landed in its grid.
.brm_alpha_diagnostics <- function(selections, alphas, label) {
  grid <- sort(unique(as.numeric(alphas)))
  selected <- as.numeric(selections$alpha)
  n_selections <- length(selected)
  lo <- grid[[1L]]
  hi <- grid[[length(grid)]]
  # counted against the grid itself, so modal_alpha is one of its values rather
  # than a value round-tripped through table()'s character names
  counts <- tabulate(match(selected, grid), nbins = length(grid))
  scored <- sum(counts) > 0L
  # a one-value grid offered no choice, so it has no boundary to be at: 1 in
  # both share columns would be arithmetically true and read as total saturation
  single <- length(grid) < 2L
  share_min <- if (single) NA_real_ else mean(selected == lo)
  share_max <- if (single) NA_real_ else mean(selected == hi)
  share_interior <- if (length(grid) < 3L) {
    NA_real_
  } else mean(selected > lo & selected < hi)
  data.frame(
    model = label,
    n_selections = n_selections,
    n_alpha_grid = length(grid),
    alpha_grid_min = lo,
    alpha_grid_max = hi,
    modal_alpha = if (scored) grid[[which.max(counts)]] else NA_real_,
    modal_share = if (scored) max(counts) / n_selections else NA_real_,
    share_at_grid_min = share_min,
    share_at_grid_max = share_max,
    share_interior = share_interior,
    # either end alone, or both ends together once there is an interior to
    # leave empty
    saturated = !single &&
      (max(share_min, share_max) >= .brm_alpha_saturation_share ||
         (length(grid) >= 3L &&
            (share_min + share_max) >= .brm_alpha_saturation_share)),
    stringsAsFactors = FALSE
  )
}

.brm_fit_diagnostics <- function(metrics, label) {
  # a constant response has no variance to explain and scores NaN; it is
  # counted, but it cannot move the summary
  r2 <- metrics$r2[is.finite(metrics$r2)]
  data.frame(
    model = label,
    n_responses = nrow(metrics),
    n_scored = length(r2),
    median_r2 = if (length(r2)) stats::median(r2) else NA_real_,
    mean_r2 = if (length(r2)) mean(r2) else NA_real_,
    share_r2_positive = if (length(r2)) mean(r2 > 0) else NA_real_,
    stringsAsFactors = FALSE
  )
}

.brm_one_sided <- function(row) {
  max(row$share_at_grid_min, row$share_at_grid_max) >=
    .brm_alpha_saturation_share
}

.brm_saturation_message <- function(row, subject) {
  if (.brm_one_sided(row)) {
    at_max <- row$share_at_grid_max >= row$share_at_grid_min
    return(paste0(
      subject, " took the ", if (at_max) "largest" else "smallest",
      " available alpha (",
      format(if (at_max) row$alpha_grid_max else row$alpha_grid_min,
             digits = 4L),
      ") for ",
      sprintf("%.1f%%",
              100 * max(row$share_at_grid_min, row$share_at_grid_max)),
      " of its ", row$n_selections, " selections"
    ))
  }
  interior <- if (is.na(row$share_interior)) {
    max(0, 1 - row$share_at_grid_min - row$share_at_grid_max)
  } else row$share_interior
  paste0(
    subject, " left only ", sprintf("%.1f%%", 100 * interior),
    " of its ", row$n_selections, " selections inside the grid (",
    sprintf("%.1f%%", 100 * row$share_at_grid_min), " took ",
    format(row$alpha_grid_min, digits = 4L), ", ",
    sprintf("%.1f%%", 100 * row$share_at_grid_max), " took ",
    format(row$alpha_grid_max, digits = 4L), ")"
  )
}

# Saturation is reported once per run: one warning naming the full model and any
# leave-one-band-out models whose selections pinned to a grid edge.
.brm_warn_alpha_saturation <- function(alpha_diagnostics,
                                       alpha_grid_origin = NA_character_) {
  hit <- alpha_diagnostics[which(alpha_diagnostics$saturated), , drop = FALSE]
  if (!nrow(hit)) return(invisible(NULL))
  full <- hit[hit$model == "full", , drop = FALSE]
  reduced <- hit[hit$model != "full", , drop = FALSE]
  parts <- character()
  if (nrow(full)) {
    parts <- c(parts, .brm_saturation_message(full, "The full model"))
  }
  if (nrow(reduced)) {
    worst <- reduced[which.max(
      reduced$share_at_grid_min + reduced$share_at_grid_max
    ), , drop = FALSE]
    bands <- sub("^without_", "", reduced$model)
    parts <- c(parts, .brm_saturation_message(worst, if (nrow(reduced) == 1L) {
      paste0("Its leave-one-band-out model (", bands, ")")
    } else {
      paste0(nrow(reduced), " leave-one-band-out models (",
             paste(bands, collapse = ", "), "), the worst of which")
    }))
  }
  directions <- vapply(seq_len(nrow(hit)), function(ii) {
    row <- hit[ii, , drop = FALSE]
    if (!.brm_one_sided(row)) return("both ends")
    if (row$share_at_grid_max >= row$share_at_grid_min) {
      "under-penalized"
    } else "over-penalized"
  }, character(1))
  direction <- if (length(unique(directions)) == 1L &&
                     directions[[1L]] != "both ends") {
    directions[[1L]]
  } else "mis-penalized"
  # `alphas` is ignored when a candidate manifest is supplied, so naming it
  # there would send the user to an argument with no effect
  remedy <- if (identical(alpha_grid_origin, "candidates")) {
    "Widen the alpha values in `candidates`"
  } else "Widen `alphas`"
  rlang::warn(paste0(
    "run_banded_ridge: the alpha grid is saturated. ",
    paste(parts, collapse = ". "), ". The inner optimum is very likely ",
    "outside the grid, so the refits are ", direction,
    " and any leave-one-band-out delta R2 compares two mis-tuned models. ",
    remedy, " until the modal selection is interior; see ",
    "`result$selection_diagnostics$alpha`."
  ))
}

.brm_warn_negative_r2 <- function(fit_diagnostics) {
  full <- fit_diagnostics[fit_diagnostics$model == "full", , drop = FALSE]
  if (!nrow(full) || !is.finite(full$median_r2) ||
      full$median_r2 >= .brm_r2_floor) {
    return(invisible(NULL))
  }
  rlang::warn(paste0(
    "run_banded_ridge: the median outer out-of-fold R2 across the ",
    full$n_scored, " scored responses is ", sprintf("%.4g", full$median_r2),
    " and only ", sprintf("%.1f%%", 100 * full$share_r2_positive),
    " of responses beat the mean of the data they were scored on. Nothing ",
    "downstream of this fit, ",
    "leave-one-band-out delta R2 included, describes explained variance. ",
    "Check the alpha grid, the design, and the fold structure before reading ",
    "the maps."
  ))
}

# Convert a cross_validation object into the internal outer-fold list. The
# object only expresses test partitions, so `purge` is applied to the training
# rows here, exactly as the integer path does.
.brm_folds_from_cross_validation <- function(spec, n, purge) {
  cls <- class(spec)[[1L]]
  if (!is.null(spec$block_var) && length(spec$block_var) != n) {
    stop(
      "banded_ridge_model: outer_crossval (", cls, ") has a block_var of ",
      "length ", length(spec$block_var), " but the dataset has ", n,
      " training observations.",
      call. = FALSE
    )
  }
  tryCatch(
    lapply(seq_len(get_nfolds(spec)), function(k) {
      test <- sort(as.integer(partition_indices(spec, k, n_samples = n)))
      train <- sort(as.integer(train_indices(spec, k)))
      list(
        id = paste0("fold_", k),
        train = .brt_purge_train(train, test, purge),
        test = test
      )
    }),
    error = function(e) {
      stop(
        "banded_ridge_model: outer_crossval of class '", cls, "' cannot be ",
        "converted to deterministic, non-overlapping folds: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
}

.brm_resolve_outer_folds <- function(spec, n, blocks, purge, seed) {
  if (is.null(spec)) {
    if (is.null(blocks) || length(unique(blocks)) < 2L) {
      stop(
        "banded_ridge_model: supply outer_crossval, or declare at least two ",
        "training blocks in feature_sets_design(block_var_train = ...).",
        call. = FALSE
      )
    }
    spec <- min(5L, length(unique(blocks)))
  }
  if (inherits(spec, "cross_validation")) {
    cls <- class(spec)[[1L]]
    folds <- .brm_folds_from_cross_validation(spec, n, purge)
    return(tryCatch(
      .brt_validate_folds(
        folds, seq_len(n), n, "outer_crossval", blocks = blocks,
        purge = purge, require_coverage = TRUE
      ),
      error = function(e) {
        stop(
          "banded_ridge_model: outer_crossval of class '", cls, "' cannot be ",
          "converted to outer folds; only schemes whose test partitions cover ",
          "every training row exactly once and respect declared blocks are ",
          "usable: ", conditionMessage(e),
          call. = FALSE
        )
      }
    ))
  } else if (is.data.frame(spec)) {
    stop(
      "outer_crossval must be an integer fold count, a cross_validation ",
      "object, or an explicit list of list(id =, train =, test =) folds; a ",
      "data.frame such as the output of crossval_samples() is not accepted.",
      call. = FALSE
    )
  } else if (is.list(spec)) {
    # Explicit fold lists are validated as supplied; see .brt_validate_folds.
  } else {
    if (!is.numeric(spec) || length(spec) != 1L || !is.finite(spec) ||
        spec != as.integer(spec) || spec < 2L) {
      stop(
        "outer_crossval must be an integer of at least two, a ",
        "cross_validation object, or an explicit fold list.",
        call. = FALSE
      )
    }
    return(.banded_ridge_make_folds(
      seq_len(n), as.integer(spec), blocks = blocks, purge = purge, seed = seed
    ))
  }
  .brt_validate_folds(
    spec, seq_len(n), n, "outer_crossval", blocks = blocks,
    purge = purge, require_coverage = TRUE
  )
}

# Default inner fold count: at most five, limited by the declared blocks (or
# rows, when no blocks are declared) available inside each outer-training set.
.brm_default_inner_v <- function(outer_folds, blocks) {
  units <- vapply(outer_folds, function(fold) {
    if (is.null(blocks)) length(fold$train) else length(unique(blocks[fold$train]))
  }, integer(1))
  inner_v <- min(5L, min(units))
  if (inner_v < 2L) {
    short <- outer_folds[[which.min(units)]]$id
    stop(
      "banded_ridge_model: tune_crossval has no usable default because outer ",
      "fold '", short, "' retains only ", min(units), " declared training ",
      "block(s); supply tune_crossval as an explicit nested fold list, or ",
      "declare more training blocks.",
      call. = FALSE
    )
  }
  inner_v
}

.brm_resolve_tuning <- function(spec, outer_folds, n, blocks, purge, seed,
                                n_available) {
  if (inherits(spec, "cross_validation")) {
    stop(
      "banded_ridge_model: tune_crossval must be an inner-fold count or one ",
      "explicit inner-fold list per outer fold; cross_validation objects are ",
      "accepted only for outer_crossval because inner folds are nested inside ",
      "each outer-training set.",
      call. = FALSE
    )
  }
  # A supplied value is type-checked even when tuning turns out to be
  # unnecessary, so a malformed argument is never silently ignored.
  if (!is.null(spec)) {
    if (is.list(spec)) {
      if (length(spec) != length(outer_folds)) {
        stop("tune_crossval must contain one inner-fold list per outer fold.",
             call. = FALSE)
      }
    } else if (!is.numeric(spec) || length(spec) != 1L || !is.finite(spec) ||
               spec != as.integer(spec) || spec < 2L) {
      stop("tune_crossval must be an inner-fold list or an integer of at least two.",
           call. = FALSE)
    }
  }
  if (n_available < 2L) {
    return(list(inner_folds = NULL, inner_v = NA_integer_, inner_tuning = "none"))
  }
  if (is.null(spec)) spec <- .brm_default_inner_v(outer_folds, blocks)
  tuning <- if (is.list(spec)) {
    list(inner_folds = spec, inner_v = 2L, inner_tuning = "nested")
  } else {
    list(inner_folds = NULL, inner_v = as.integer(spec), inner_tuning = "nested")
  }
  # Construct and validate every inner fold now, so an impossible request
  # fails at model creation rather than inside run_banded_ridge().
  for (oo in seq_along(outer_folds)) {
    tryCatch(
      .brt_inner_folds(tuning$inner_folds, oo, outer_folds[[oo]]$train, n,
                       tuning$inner_v, blocks, purge, seed),
      error = function(e) {
        hint <- if (is.null(blocks) && purge > 0L && is.null(tuning$inner_folds)) {
          paste0(
            " Without declared training blocks, inner folds are random ",
            "row-wise splits, so a temporal purge can remove most training ",
            "rows; declare contiguous blocks with feature_sets_design(",
            "block_var_train = ...) or supply explicit nested inner folds."
          )
        } else ""
        stop(
          "banded_ridge_model: tune_crossval cannot be applied inside outer ",
          "fold '", outer_folds[[oo]]$id, "': ", conditionMessage(e), hint,
          call. = FALSE
        )
      }
    )
  }
  if (!is.null(tuning$inner_folds)) {
    tuning$inner_v <- as.integer(max(lengths(tuning$inner_folds)))
  }
  tuning
}

.brm_retention_estimate <- function(n, p, v, n_outer,
                                    weight_retention, return_predictions,
                                    retain_diagnostics, n_candidates, inner_v,
                                    n_delta = 0L, n_groups = 1L) {
  if (is.na(inner_v)) inner_v <- 0L
  prediction <- if (return_predictions) 8 * n * v * (1L + n_delta) else 0
  weights <- switch(weight_retention,
    none = 0,
    primal = 8 * n_outer * p * v,
    dual = 8 * n_outer * ceiling(n * (n_outer - 1) / n_outer) * v
  )
  diagnostics <- if (retain_diagnostics) {
    8 * n_outer * v * (n_candidates * (inner_v + 1L) + n + p) *
      (1L + n_delta)
  } else 0
  attribution <- if (n_delta) {
    8 * n_delta * v * (6L + n_outer * (5L + n_groups))
  } else 0
  c(predictions = prediction, weights = weights, diagnostics = diagnostics,
    attribution = attribution,
    total = prediction + weights + diagnostics + attribution)
}

.brm_response_matrix <- function(dataset, ids, n_expected) {
  Y <- as.matrix(neuroim2::series(dataset$train_data, ids))
  if (nrow(Y) != n_expected && ncol(Y) == n_expected) Y <- t(Y)
  if (nrow(Y) != n_expected || ncol(Y) != length(ids)) {
    stop("banded_ridge_model: extracted response matrix does not match design rows or requested voxel IDs.",
         call. = FALSE)
  }
  colnames(Y) <- paste0("voxel_", ids)
  Y
}

.brm_validate_partitions <- function(partitions, active_ids, target_batch_size) {
  if (is.null(partitions)) partitions <- list(active_ids)
  if (!is.list(partitions) || !length(partitions)) {
    stop("response_partitions must be a non-empty list of active voxel IDs.",
         call. = FALSE)
  }
  partitions <- lapply(seq_along(partitions), function(ii) {
    x <- partitions[[ii]]
    if (!is.numeric(x) || !length(x) || any(!is.finite(x)) ||
        any(x != as.integer(x)) || anyDuplicated(x)) {
      stop("response_partitions[[", ii,
           "]] must contain unique finite integer voxel IDs.", call. = FALSE)
    }
    as.integer(x)
  })
  flattened <- unlist(partitions, use.names = FALSE)
  if (anyDuplicated(flattened) || !setequal(flattened, active_ids) ||
      length(flattened) != length(active_ids)) {
    stop("response_partitions must contain every active voxel exactly once with no overlap.",
         call. = FALSE)
  }
  batches <- unlist(lapply(partitions, function(ids) {
    split(ids, ceiling(seq_along(ids) / target_batch_size))
  }), recursive = FALSE)
  unname(batches)
}

.brm_map <- function(values, ids, mask) {
  full <- rep(NA_real_, prod(dim(mask)[1:3]))
  full[ids] <- as.numeric(values)
  neuroim2::NeuroVol(
    array(full, dim = dim(mask)[1:3]), neuroim2::space(mask)
  )
}

.brm_extract_primal_weights <- function(nested) {
  out <- list()
  for (outer in nested$outer_results) {
    for (entry in outer$outer_fits) {
      fit <- entry$fit
      if (is.null(fit$coefficients)) next
      for (jj in seq_along(entry$responses)) {
        response <- entry$responses[[jj]]
        out[[response]][[outer$fold_id]] <- list(
          outer_fold = outer$fold_id,
          candidate_id = entry$candidate_id,
          coefficients = fit$coefficients[, jj]
        )
      }
    }
  }
  out
}

.brm_extract_dual_weights <- function(nested,
                                      X,
                                      Y,
                                      groups,
                                      cache,
                                      cache_prefix) {
  out <- list()
  theta_columns <- paste0("theta_", nested$group_names)
  for (oo in seq_along(nested$outer_results)) {
    outer <- nested$outer_results[[oo]]
    selected <- outer$selected
    for (candidate_index in unique(selected$candidate_index)) {
      rows <- which(selected$candidate_index == candidate_index)
      alpha <- selected$alpha[rows[[1L]]]
      theta <- as.numeric(selected[rows[[1L]], theta_columns, drop = TRUE])
      names(theta) <- nested$group_names
      theta_key <- .brt_theta_keys(matrix(theta, nrow = 1L))
      path <- .banded_ridge_fit_dual_path(
        .brc_subset_rows(X, outer$train), Y[outer$train, rows, drop = FALSE],
        alphas = alpha, theta = theta, groups = groups,
        recover_primal = FALSE, cache = cache,
        cache_key = paste(cache_prefix, "outer", oo, "refit", "theta",
                          theta_key, sep = "::"),
        kernel_cache_key = paste(cache_prefix, "outer", oo, "refit", sep = "::")
      )
      fit <- path$fits[[1L]]
      for (jj in seq_along(rows)) {
        response <- selected$response[rows[[jj]]]
        out[[response]][[outer$fold_id]] <- list(
          outer_fold = outer$fold_id,
          candidate_id = selected$candidate_id[rows[[jj]]],
          train = outer$train,
          dual_weights = fit$dual_weights[, jj]
        )
      }
    }
  }
  out
}

#' First-class single-domain banded-ridge encoding model
#'
#' `banded_ridge_model()` defines leakage-safe nested tuning for grouped stimulus
#' predictors and independent brain responses. Feature-band membership and
#' training blocks come exclusively from `feature_sets_design()`.
#'
#' @param dataset An image `mvpa_dataset` whose training rows are observations
#'   and active mask locations are responses.
#' @param design A `feature_sets_design` containing `X_train`, band membership,
#'   and optional `block_var_train` / `time_series` metadata.
#' @param outer_crossval Outer fold specification: an integer fold count, a
#'   `cross_validation` object (for example
#'   \code{\link{blocked_cross_validation}},
#'   \code{\link{kfold_cross_validation}}, or
#'   \code{\link{custom_cross_validation}}; only schemes whose test
#'   partitions cover every training row exactly once convert, so resampled or
#'   non-deterministic schemes are refused), or an explicit
#'   list of `list(id =, train =, test =)` folds. If omitted, intact training
#'   blocks are required and define at most five folds. Integer counts and
#'   `cross_validation` objects have `purge` applied to their training rows;
#'   an explicit fold list is validated as supplied and must already respect
#'   `purge`.
#' @param tune_crossval Inner fold count, or one explicit inner-fold list per
#'   outer fold. `NULL` (the default) uses at most five inner folds, limited
#'   by the number of declared training blocks (or rows, when no blocks are
#'   declared) available inside each outer training set. A supplied value is
#'   still type-checked but otherwise ignored when only one candidate survives the alpha/theta
#'   scopes (for example `theta_method = "fixed"` with one `theta` and one
#'   `alphas` value, or fixed scopes with `fixed_alpha` and `fixed_theta`):
#'   nothing needs selecting, so inner tuning is skipped, the model's
#'   `inner_v` is `NA`, and reported `inner_score` values are `NA`. Inner folds
#'   are constructed and validated when the model is created, so an impossible
#'   request fails here rather than in `run_banded_ridge()`.
#' @param candidates Optional candidate manifest. By default one is created from
#'   `alphas`, `theta_method`, and the remaining theta arguments.
#' @param alphas Non-negative overall ridge penalties used when constructing
#'   candidates, or `"auto"` (the default) to scale the grid to the design.
#'   Every solver path standardizes each training column and then scales band
#'   `b` by `sqrt(theta_b)`, so `alpha` competes with the eigenvalues of the
#'   resulting cross-product, which average
#'   `(n - 1) * mean(D_b) / min(n - 1, p)` under a uniform `theta` for band
#'   widths `D_b`. `"auto"` places nine points spanning six decades centred on
#'   that anchor, so the grid brackets the range in which the penalty changes
#'   the fit however many rows and features the design has. A fixed grid cannot:
#'   a design with a few thousand timepoints needs penalties several orders of
#'   magnitude above the fixed `10^seq(-2, 2)` grid this default replaced, and
#'   every response then selects the largest available value. `run_banded_ridge()`
#'   warns when that happens; see `result$selection_diagnostics`. The resolved
#'   grid is on the returned model as `alpha_grid`.
#' @param theta_method One of `"grid"`, `"fixed"`, or `"random"`. The default
#'   constructs a small, deterministic simplex grid.
#' @param theta Fixed named simplex vector/matrix for `theta_method = "fixed"`.
#' @param theta_grid_points Grid resolution for simplex grid candidates.
#' @param n_theta Number of deterministic random simplex points.
#' @param metric Inner selection metric; MSE is the default estimand.
#' @param alpha_scope,theta_scope `"response"`, `"shared"`/`"roi"`, or
#'   `"fixed"` selection constraints.
#' @param fixed_alpha,fixed_theta Values required by corresponding fixed scopes.
#' @param delta_sets Optional feature-band names, or `"all"`, for independently
#'   retuned predictive leave-one-band-out outer-OOF delta R2. `NULL` (the
#'   default) performs no reduced-model work.
#' @param purge Non-negative temporal gap (in rows) around validation/test
#'   rows. Training rows within the gap are dropped when folds are constructed
#'   from an integer count or a `cross_validation` object; explicit fold lists
#'   (outer or inner) are checked against the gap and rejected if they violate
#'   it, never modified.
#' @param solver `"auto"`, `"direct"`, `"svd_primal"`, or `"dual_kernel"`.
#' @param target_batch_size Maximum responses evaluated together.
#' @param return_predictions Retain the complete outer out-of-fold prediction
#'   matrix.
#' @param weight_retention Retain no weights, outer-fold primal coefficients, or
#'   outer-fold dual weights.
#' @param retain_diagnostics Retain full candidate/fold diagnostics per chunk.
#' @param max_retained_mb Allocation contract for predictions, weights, and
#'   diagnostics.
#' @param weight_overflow Refuse an unsafe request or fall back to no weights.
#' @param seed Deterministic fold/candidate seed.
#' @param memory_limit_mb Solver intermediate-allocation contract in MiB. It
#'   bounds two separate pools, so peak solver memory can approach twice this
#'   value: per-fit intermediates are checked against it up front, and the
#'   retained decomposition cache of the optimized solvers is capped at it
#'   independently. Once cached band Grams and eigendecompositions would exceed
#'   the cap, the least recently used entries are evicted and recomputed on
#'   demand. Results are unchanged at any value; only reuse is lost, and
#'   provenance reports `solver_cache_evictions` and `solver_cache_peak_mb`.
#'
#' @return A `banded_ridge_model` specification for `run_banded_ridge()`.
#' @export
#' @examples
#' set.seed(71)
#' n <- 24L
#' dims <- c(3L, 3L, 3L)
#' X <- matrix(rnorm(n * 6), n, 6)
#' Y <- X %*% matrix(rnorm(6 * prod(dims), sd = 0.8), 6, prod(dims)) +
#'   matrix(rnorm(n * prod(dims), sd = 0.5), n, prod(dims))
#' dataset <- mvpa_dataset(
#'   neuroim2::NeuroVec(array(as.vector(t(Y)), c(dims, n)),
#'                      neuroim2::NeuroSpace(c(dims, n), c(1, 1, 1))),
#'   mask = neuroim2::NeuroVol(array(1, dims),
#'                             neuroim2::NeuroSpace(dims, c(1, 1, 1)))
#' )
#' fs <- feature_sets(X, blocks(low = 3, semantic = 3))
#' design <- feature_sets_design(
#'   fs, block_var_train = rep(1:4, each = 6), time_series = TRUE
#' )
#' model <- banded_ridge_model(
#'   dataset, design, outer_crossval = 4, tune_crossval = 2,
#'   theta_method = "fixed", theta = c(low = 0.5, semantic = 0.5),
#'   target_batch_size = 9, delta_sets = "all"
#' )
#' result <- run_banded_ridge(model)
#' head(result$metrics)
#' result$selection_diagnostics$alpha
#' head(result$predictive_leave_one_band_out$effects)
#'
#' @section Predictive leave-one-band-out delta R2:
#' When `delta_sets` is enabled, each requested reduced model is tuned again
#' inside every outer-training set, over a candidate manifest projected onto
#' the retained bands. The reported effect is `R2_full - R2_without_band` from
#' matched outer out-of-fold predictions. It is a predictive drop-out effect,
#' not an additive unique/shared variance partition: values are not clipped,
#' may be negative, and need not sum to the full-model R2.
banded_ridge_model <- function(dataset,
                               design,
                               outer_crossval = NULL,
                               tune_crossval = NULL,
                               candidates = NULL,
                               alphas = "auto",
                               theta_method = c("grid", "fixed", "random"),
                               theta = NULL,
                               theta_grid_points = 3L,
                               n_theta = 20L,
                               metric = c("mse", "correlation", "r2"),
                               alpha_scope = "response",
                               theta_scope = "response",
                               fixed_alpha = NULL,
                               fixed_theta = NULL,
                               delta_sets = NULL,
                               purge = 0L,
                               solver = c("auto", "direct", "svd_primal", "dual_kernel"),
                               target_batch_size = 256L,
                               return_predictions = FALSE,
                               weight_retention = c("none", "dual", "primal"),
                               retain_diagnostics = FALSE,
                               max_retained_mb = 1024,
                               weight_overflow = c("error", "none"),
                               seed = 1L,
                               memory_limit_mb = 1024) {
  if (!inherits(dataset, "mvpa_image_dataset")) {
    stop("banded_ridge_model currently supports image mvpa_dataset objects only.",
         call. = FALSE)
  }
  if (!inherits(design, "feature_sets_design")) {
    stop("banded_ridge_model: design must be created by feature_sets_design().",
         call. = FALSE)
  }
  metric <- match.arg(metric)
  theta_method <- match.arg(theta_method)
  solver <- match.arg(solver)
  weight_retention <- match.arg(weight_retention)
  weight_overflow <- match.arg(weight_overflow)
  if (!is.logical(return_predictions) || length(return_predictions) != 1L ||
      is.na(return_predictions) || !is.logical(retain_diagnostics) ||
      length(retain_diagnostics) != 1L || is.na(retain_diagnostics)) {
    stop("return_predictions and retain_diagnostics must be TRUE or FALSE.",
         call. = FALSE)
  }
  target_batch_size <- .brs_batch_size(target_batch_size,
                                       length(.brm_active_ids(dataset)),
                                       "target_batch_size")
  if (!is.numeric(purge) || length(purge) != 1L || !is.finite(purge) ||
      purge != as.integer(purge) || purge < 0L) {
    stop("purge must be one non-negative integer.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("seed must be one finite integer.", call. = FALSE)
  }
  n <- .brm_n_observations(dataset)
  if (nrow(design$X_train$X) != n) {
    stop("banded_ridge_model: design X_train rows must match dataset training observations.",
         call. = FALSE)
  }
  blocks <- design$block_var_train
  if (isTRUE(design$time_series) && is.null(blocks) && purge < 1L) {
    stop(
      "banded_ridge_model: time_series = TRUE requires intact block_var_train ",
      "or an explicit positive temporal purge; random row-wise CV is unsafe.",
      call. = FALSE
    )
  }
  outer_folds <- .brm_resolve_outer_folds(
    outer_crossval, n, blocks, as.integer(purge), as.integer(seed)
  )
  group_names <- design$X_train$set_order
  if (is.null(candidates)) {
    alpha_grid <- .brm_resolve_alphas(
      alphas, design$X_train$X, as.character(design$X_train$set)
    )
    candidates <- .banded_ridge_candidates(
      group_names, alpha_grid$alphas, method = theta_method, theta = theta,
      grid_points = theta_grid_points, n_theta = n_theta, seed = seed
    )
  } else {
    candidates <- .brt_validate_candidates(candidates, group_names)
    alpha_grid <- list(alphas = sort(unique(candidates$alpha)),
                       origin = "candidates", anchor = NA_real_)
  }
  available <- .brt_available_candidates(
    candidates, alpha_scope, theta_scope,
    fixed_alpha = fixed_alpha, fixed_theta = fixed_theta
  )
  tuning <- .brm_resolve_tuning(
    tune_crossval, outer_folds, n, blocks, as.integer(purge), as.integer(seed),
    n_available = length(available)
  )
  delta_sets <- .bra_validate_delta_sets(delta_sets, group_names)
  reduced_candidates <- setNames(lapply(delta_sets, function(band) {
    .bra_reduced_candidates(candidates, band)
  }), delta_sets)
  reduced_fixed_theta <- setNames(lapply(delta_sets, function(band) {
    .bra_reduced_fixed_theta(fixed_theta, group_names, band)
  }), delta_sets)
  active_ids <- .brm_active_ids(dataset)
  plan <- .banded_ridge_solver_plan(
    n = n, p = ncol(design$X_train$X), v = min(target_batch_size, length(active_ids)),
    n_bands = length(group_names), n_alphas = length(unique(candidates$alpha)),
    alpha_batch_size = 1L, response_batch_size = target_batch_size,
    memory_limit_mb = memory_limit_mb, solver = solver
  )
  resolved_solver <- plan$solver
  if (resolved_solver == "direct") resolved_solver <- "direct"
  retention <- .brm_retention_estimate(
    n, ncol(design$X_train$X), length(active_ids), length(outer_folds),
    weight_retention, return_predictions, retain_diagnostics,
    length(candidates$alpha), tuning$inner_v,
    n_delta = length(delta_sets), n_groups = length(group_names)
  )
  if (!is.numeric(max_retained_mb) || length(max_retained_mb) != 1L ||
      !is.finite(max_retained_mb) || max_retained_mb <= 0) {
    stop("max_retained_mb must be one positive finite value.", call. = FALSE)
  }
  retention_notice <- NULL
  if (retention[["total"]] > max_retained_mb * 1024^2) {
    if (weight_retention != "none" && weight_overflow == "none") {
      retention_notice <- paste0(
        "Requested weight retention was disabled because estimated retained ",
        "storage exceeded ", max_retained_mb, " MiB."
      )
      weight_retention <- "none"
      retention <- .brm_retention_estimate(
        n, ncol(design$X_train$X), length(active_ids), length(outer_folds),
        weight_retention, return_predictions, retain_diagnostics,
        length(candidates$alpha), tuning$inner_v,
        n_delta = length(delta_sets), n_groups = length(group_names)
      )
    }
    if (retention[["total"]] > max_retained_mb * 1024^2) {
      stop(
        "banded_ridge_model: requested retained outputs are estimated at ",
        sprintf("%.3f", retention[["total"]] / 1024^2),
        " MiB, above max_retained_mb = ", max_retained_mb,
        ". Reduce diagnostics/predictions/weights or raise the explicit limit.",
        call. = FALSE
      )
    }
  }
  create_model_spec(
    "banded_ridge_model", dataset, design,
    return_predictions = return_predictions, compute_performance = TRUE,
    outer_folds = outer_folds, inner_folds = tuning$inner_folds,
    inner_v = tuning$inner_v, inner_tuning = tuning$inner_tuning,
    candidates = candidates, metric = metric,
    alpha_grid = alpha_grid$alphas, alpha_grid_origin = alpha_grid$origin,
    alpha_anchor = alpha_grid$anchor,
    selectable_alphas = sort(unique(candidates$alpha[available])),
    alpha_scope = alpha_scope, theta_scope = theta_scope,
    fixed_alpha = fixed_alpha, fixed_theta = fixed_theta,
    delta_sets = delta_sets, reduced_candidates = reduced_candidates,
    reduced_fixed_theta = reduced_fixed_theta,
    purge = as.integer(purge), solver = resolved_solver, solver_plan = plan,
    target_batch_size = target_batch_size, weight_retention = weight_retention,
    retain_diagnostics = retain_diagnostics,
    retention_estimated_bytes = retention,
    retention_notice = retention_notice,
    max_retained_mb = max_retained_mb, memory_limit_mb = memory_limit_mb,
    seed = as.integer(seed),
    active_response_ids = active_ids, has_test_set = FALSE
  )
}

#' Run a chunked whole-brain banded-ridge encoding analysis
#'
#' `run_banded_ridge()` evaluates each active mask voxel exactly once. Optional
#' `response_partitions` may reorder non-overlapping response batches (the same
#' contract used by regional batching), but may not overlap or omit voxels.
#' When `model$delta_sets` is non-empty, the result additionally contains
#' `predictive_leave_one_band_out`: matched full/reduced OOF metrics and
#' predictions, independently selected reduced hyperparameters, spatial
#' `delta_cv_r2_<band>` maps, and B+1 model-cost provenance.
#'
#' @param model A `banded_ridge_model` specification.
#' @param target_batch_size Optional runtime override of the maximum response
#'   chunk size.
#' @param response_partitions Optional list of non-overlapping global voxel IDs
#'   that exactly cover the active mask.
#' @return A `banded_ridge_result` with spatial maps, metrics, exact outer-fold
#'   hyperparameters, selection diagnostics, optional predictions/weights, and
#'   allocation provenance.
#'
#' @section Selection diagnostics:
#' `result$selection_diagnostics$alpha` carries one row per fitted model (the
#' full model and each leave-one-band-out model) giving the alpha grid it could
#' select from, the modal selection and its share, the share pinned to each end
#' of the grid, and the share strictly interior. `$fit` gives the median and
#' mean outer out-of-fold R2 per model and the share of responses above zero.
#'
#' Two conditions are warned about rather than left to be discovered. The first
#' is a saturated grid: at least 95% of a model's alpha selections taking the
#' largest available alpha, or the smallest, or the two ends between them once
#' the grid has an interior to leave empty. Under the default per-response
#' alpha scope a heterogeneous mask splits its boundary mass across both ends,
#' so the combined share is what catches a grid that brackets nothing. In every
#' form the inner optimum lies outside the grid: the refits are mis-penalized
#' and leave-one-band-out delta R2 compares two mis-tuned models. The second is
#' a median outer out-of-fold R2 below -0.05, which means the fit predicts
#' worse than the mean of the data it was scored on, whatever the cause, and
#' nothing derived from it describes explained variance.
#' @export
run_banded_ridge <- function(model,
                             target_batch_size = model$target_batch_size,
                             response_partitions = NULL) {
  if (!inherits(model, "banded_ridge_model")) {
    stop("run_banded_ridge: model must be created by banded_ridge_model().",
         call. = FALSE)
  }
  active_ids <- model$active_response_ids
  batch_size <- .brs_batch_size(target_batch_size, length(active_ids),
                                "target_batch_size")
  batches <- .brm_validate_partitions(response_partitions, active_ids, batch_size)
  X <- model$design$X_train$X
  groups <- as.character(model$design$X_train$set)
  n <- nrow(X)
  v <- length(active_ids)
  runtime_plan <- .banded_ridge_solver_plan(
    n = n, p = ncol(X), v = min(batch_size, v),
    n_bands = length(model$candidates$group_names),
    n_alphas = length(unique(model$candidates$alpha)),
    alpha_batch_size = 1L, response_batch_size = batch_size,
    memory_limit_mb = model$memory_limit_mb, solver = model$solver
  )
  predictions <- if (model$return_predictions) {
    matrix(NA_real_, n, v,
           dimnames = list(as.character(seq_len(n)), paste0("voxel_", active_ids)))
  } else NULL
  delta_bands <- model$delta_sets
  reduced_metric_chunks <- setNames(lapply(delta_bands, function(x) {
    vector("list", length(batches))
  }), delta_bands)
  reduced_selection_chunks <- setNames(lapply(delta_bands, function(x) {
    vector("list", length(batches))
  }), delta_bands)
  reduced_predictions <- if (model$return_predictions && length(delta_bands)) {
    setNames(lapply(delta_bands, function(x) {
      matrix(NA_real_, n, v,
             dimnames = list(as.character(seq_len(n)),
                             paste0("voxel_", active_ids)))
    }), delta_bands)
  } else NULL
  reduced_diagnostics <- if (model$retain_diagnostics && length(delta_bands)) {
    setNames(lapply(delta_bands, function(x) vector("list", length(batches))),
             delta_bands)
  } else NULL
  metric_rows <- list()
  selection_rows <- list()
  diagnostics <- if (model$retain_diagnostics) list() else NULL
  weights <- if (model$weight_retention == "none") NULL else list()
  cache_bytes <- model$memory_limit_mb * 1024^2
  cache <- if (model$solver == "direct") {
    NULL
  } else .banded_ridge_solver_cache(max_bytes = cache_bytes)
  weight_cache <- if (model$weight_retention == "dual" && is.null(cache)) {
    .banded_ridge_solver_cache(max_bytes = cache_bytes)
  } else cache
  max_chunk_result_bytes <- 0
  peak_dimensions <- NULL
  chunk_manifest <- vector("list", length(batches))
  work_rows <- list()

  for (bb in seq_along(batches)) {
    ids <- batches[[bb]]
    Y <- .brm_response_matrix(model$dataset, ids, n)
    counters_before <- .bra_cache_counters(cache)
    nested <- .banded_ridge_nested_cv(
      X, Y, model$outer_folds, model$candidates,
      feature_groups = groups, inner_folds = model$inner_folds,
      inner_v = model$inner_v, blocks = model$design$block_var_train,
      purge = model$purge, metric = model$metric,
      alpha_scope = model$alpha_scope, theta_scope = model$theta_scope,
      fixed_alpha = model$fixed_alpha, fixed_theta = model$fixed_theta,
      response_batch_size = length(ids), seed = model$seed,
      solver = model$solver, solver_cache = cache,
      cache_prefix = "public"
    )
    counters_after <- .bra_cache_counters(cache)
    work_rows[[length(work_rows) + 1L]] <- .bra_work_row(
      "full", NULL, bb, ncol(X), length(ids),
      length(model$candidates$alpha), nested, counters_before, counters_after
    )
    max_chunk_result_bytes <- max(max_chunk_result_bytes,
                                  as.numeric(utils::object.size(nested)))
    peak_dimensions <- if (is.null(peak_dimensions)) {
      nested$peak_intermediate_dimensions
    } else pmax(peak_dimensions, nested$peak_intermediate_dimensions)
    ids_by_name <- setNames(ids, paste0("voxel_", ids))
    metrics <- nested$performance
    metrics$voxel_index <- unname(ids_by_name[metrics$response])
    metric_rows[[bb]] <- metrics
    selections <- nested$selections
    selections$voxel_index <- unname(ids_by_name[selections$response])
    selection_rows[[bb]] <- selections
    if (model$return_predictions) {
      predictions[, match(ids, active_ids)] <- nested$predictions
    }
    if (model$retain_diagnostics) diagnostics[[bb]] <- nested$outer_results

    for (band in delta_bands) {
      keep <- groups != band
      reduced_x <- X[, keep, drop = FALSE]
      reduced_groups <- groups[keep]
      reduced_manifest <- model$reduced_candidates[[band]]
      counters_before <- .bra_cache_counters(cache)
      reduced <- .banded_ridge_nested_cv(
        reduced_x, Y, model$outer_folds, reduced_manifest,
        feature_groups = reduced_groups, inner_folds = model$inner_folds,
        inner_v = model$inner_v, blocks = model$design$block_var_train,
        purge = model$purge, metric = model$metric,
        alpha_scope = model$alpha_scope, theta_scope = model$theta_scope,
        fixed_alpha = model$fixed_alpha,
        fixed_theta = model$reduced_fixed_theta[[band]],
        response_batch_size = length(ids), seed = model$seed,
        solver = model$solver, solver_cache = cache,
        cache_prefix = paste("public", "without", band, sep = "::"),
        kernel_cache_prefix = "public"
      )
      counters_after <- .bra_cache_counters(cache)
      work_rows[[length(work_rows) + 1L]] <- .bra_work_row(
        paste0("without_", band), band, bb, ncol(reduced_x), length(ids),
        length(reduced_manifest$alpha), reduced, counters_before, counters_after
      )
      max_chunk_result_bytes <- max(max_chunk_result_bytes,
                                    as.numeric(utils::object.size(reduced)))
      peak_dimensions <- pmax(peak_dimensions,
                              reduced$peak_intermediate_dimensions)
      reduced_metrics <- reduced$performance
      reduced_metrics$voxel_index <- unname(ids_by_name[reduced_metrics$response])
      reduced_metric_chunks[[band]][[bb]] <- reduced_metrics
      reduced_selections <- reduced$selections
      reduced_selections$voxel_index <-
        unname(ids_by_name[reduced_selections$response])
      reduced_selection_chunks[[band]][[bb]] <- reduced_selections
      if (model$return_predictions) {
        reduced_predictions[[band]][, match(ids, active_ids)] <-
          reduced$predictions
      }
      if (model$retain_diagnostics) {
        reduced_diagnostics[[band]][[bb]] <- reduced$outer_results
      }
    }

    if (model$weight_retention == "primal") {
      weights <- c(weights, .brm_extract_primal_weights(nested))
    } else if (model$weight_retention == "dual") {
      weights <- c(weights, .brm_extract_dual_weights(
        nested, X, Y, groups, weight_cache, "public"
      ))
    }
    chunk_manifest[[bb]] <- list(
      chunk = bb, voxel_ids = ids, n_responses = length(ids),
      solver = model$solver
    )
  }
  metrics <- do.call(rbind, metric_rows)
  metrics <- metrics[match(active_ids, metrics$voxel_index),
                     c("voxel_index", "response", "n_obs", "mse",
                       "correlation", "r2")]
  rownames(metrics) <- NULL
  selections <- do.call(rbind, selection_rows)
  selections <- selections[order(match(selections$voxel_index, active_ids),
                                 selections$outer_fold), ]
  rownames(selections) <- NULL
  # model specs built before selectable_alphas existed still describe their own
  # scope constraints, so the grid is recoverable rather than fatal
  full_alphas <- model$selectable_alphas %||% .brm_selectable_alphas(
    model$candidates, model$alpha_scope, model$theta_scope,
    fixed_alpha = model$fixed_alpha, fixed_theta = model$fixed_theta
  )
  alpha_diagnostic_rows <- list(.brm_alpha_diagnostics(
    selections, full_alphas, "full"
  ))
  fit_diagnostic_rows <- list(.brm_fit_diagnostics(metrics, "full"))
  maps <- list(
    mse = .brm_map(metrics$mse, active_ids, model$dataset$mask),
    correlation = .brm_map(metrics$correlation, active_ids, model$dataset$mask),
    r2 = .brm_map(metrics$r2, active_ids, model$dataset$mask)
  )
  selection_columns <- c("alpha", paste0("theta_", model$candidates$group_names))
  for (column in selection_columns) {
    averaged <- vapply(active_ids, function(id) {
      mean(selections[[column]][selections$voxel_index == id])
    }, numeric(1))
    maps[[paste0(column, "_mean")]] <- .brm_map(
      averaged, active_ids, model$dataset$mask
    )
  }
  predictive_leave_one_band_out <- NULL
  if (length(delta_bands)) {
    effects <- metrics[, c("voxel_index", "response"), drop = FALSE]
    reduced_performance <- list()
    reduced_hyperparameters <- list()
    effect_maps <- list()
    for (band in delta_bands) {
      reduced_metric <- do.call(rbind, reduced_metric_chunks[[band]])
      reduced_metric <- reduced_metric[
        match(active_ids, reduced_metric$voxel_index),
        c("voxel_index", "response", "n_obs", "mse", "correlation", "r2")
      ]
      rownames(reduced_metric) <- NULL
      reduced_selection <- do.call(rbind, reduced_selection_chunks[[band]])
      reduced_selection <- reduced_selection[order(
        match(reduced_selection$voxel_index, active_ids),
        reduced_selection$outer_fold
      ), ]
      rownames(reduced_selection) <- NULL
      effect_name <- paste0("delta_cv_r2_", band)
      effect <- metrics$r2 - reduced_metric$r2
      effects[[effect_name]] <- effect
      effect_maps[[effect_name]] <- .brm_map(
        effect, active_ids, model$dataset$mask
      )
      reduced_performance[[band]] <- reduced_metric
      reduced_hyperparameters[[band]] <- reduced_selection
      alpha_diagnostic_rows[[length(alpha_diagnostic_rows) + 1L]] <-
        .brm_alpha_diagnostics(
          reduced_selection,
          .brm_selectable_alphas(
            model$reduced_candidates[[band]], model$alpha_scope,
            model$theta_scope, fixed_alpha = model$fixed_alpha,
            fixed_theta = model$reduced_fixed_theta[[band]]
          ),
          paste0("without_", band)
        )
      fit_diagnostic_rows[[length(fit_diagnostic_rows) + 1L]] <-
        .brm_fit_diagnostics(reduced_metric, paste0("without_", band))
    }
    predictive_leave_one_band_out <- list(
      estimand = paste0(
        "predictive leave-one-band-out outer-OOF delta R2: ",
        "R2_full - R2_without_band; independently retuned; not additive"
      ),
      effects = effects,
      maps = effect_maps,
      reduced_performance = reduced_performance,
      reduced_hyperparameters = reduced_hyperparameters,
      reduced_predictions = reduced_predictions,
      diagnostics = reduced_diagnostics,
      candidate_manifests = model$reduced_candidates,
      provenance = list(
        requested_bands = delta_bands,
        conceptual_model_count = 1L + length(delta_bands),
        matched_outer_fold_ids = vapply(
          model$outer_folds, `[[`, character(1), "id"
        ),
        matched_outer_test_rows = setNames(
          lapply(model$outer_folds, `[[`, "test"),
          vapply(model$outer_folds, `[[`, character(1), "id")
        ),
        preprocessing = paste0(
          "centering and scaling are fit independently on each declared ",
          "inner/outer training split"
        ),
        scoring = "pooled outer-OOF R2 on identical response/test rows",
        candidate_projection = "drop band, renormalize retained theta, retune",
        reduced_models_retuned = TRUE,
        clipping = FALSE,
        additive_partition = FALSE
      )
    )
  }
  if (!is.null(weights)) weights <- weights[paste0("voxel_", active_ids)]
  alpha_diagnostics <- do.call(rbind, alpha_diagnostic_rows)
  rownames(alpha_diagnostics) <- NULL
  fit_diagnostics <- do.call(rbind, fit_diagnostic_rows)
  rownames(fit_diagnostics) <- NULL
  selection_diagnostics <- list(
    alpha = alpha_diagnostics,
    fit = fit_diagnostics,
    alpha_grid_origin = model$alpha_grid_origin %||% NA_character_,
    alpha_anchor = model$alpha_anchor %||% NA_real_,
    saturation_share = .brm_alpha_saturation_share,
    r2_floor = .brm_r2_floor
  )
  .brm_warn_alpha_saturation(alpha_diagnostics,
                             selection_diagnostics$alpha_grid_origin)
  .brm_warn_negative_r2(fit_diagnostics)
  out <- list(
    metrics = metrics,
    hyperparameters = selections,
    selection_diagnostics = selection_diagnostics,
    maps = maps,
    predictions = predictions,
    weights = weights,
    diagnostics = diagnostics,
    active_response_ids = active_ids,
    candidate_manifest = model$candidates,
    predictive_leave_one_band_out = predictive_leave_one_band_out,
    provenance = list(
      solver = model$solver,
      solver_plan = runtime_plan,
      outer_fold_ids = vapply(model$outer_folds, `[[`, character(1), "id"),
      seed = model$seed,
      metric = model$metric,
      chunk_manifest = chunk_manifest,
      target_batch_size = batch_size,
      peak_intermediate_dimensions = peak_dimensions,
      max_chunk_result_bytes = max_chunk_result_bytes,
      conceptual_model_count = 1L + length(delta_bands),
      nested_cv_invocations = length(batches) * (1L + length(delta_bands)),
      work_manifest = do.call(rbind, work_rows),
      retention_estimated_bytes = model$retention_estimated_bytes,
      retention_notice = model$retention_notice,
      solver_decomposition_count = if (is.null(cache)) 0L else cache$decomposition_count,
      solver_cache_hits = if (is.null(cache)) 0L else cache$hit_count,
      solver_band_kernel_builds = if (is.null(cache)) 0L else cache$kernel_build_count,
      solver_band_kernel_hits = if (is.null(cache)) 0L else cache$kernel_hit_count,
      solver_cache_evictions = if (is.null(cache)) 0L else cache$eviction_count,
      solver_cache_oversize = if (is.null(cache)) 0L else cache$oversize_count,
      solver_cache_peak_mb = if (is.null(cache)) 0 else cache$peak_bytes / 1024^2,
      solver_cache_limit_mb = model$memory_limit_mb,
      weight_decomposition_count = if (is.null(weight_cache) ||
                                          identical(weight_cache, cache)) {
        0L
      } else weight_cache$decomposition_count,
      weight_cache_hits = if (is.null(weight_cache) ||
                                 identical(weight_cache, cache)) {
        0L
      } else weight_cache$hit_count
    )
  )
  class(out) <- c("banded_ridge_result", "list")
  out
}

#' @export
print.banded_ridge_result <- function(x, ...) {
  cat("banded_ridge_result\n")
  cat("===================\n")
  cat("Responses: ", nrow(x$metrics), "\n", sep = "")
  cat("Solver:    ", x$provenance$solver, "\n", sep = "")
  cat("OOF MSE:   ", sprintf("%.6g", mean(x$metrics$mse)), " (mean)\n", sep = "")
  # results saved before selection_diagnostics existed still print
  diagnostics <- x$selection_diagnostics
  full_fit <- if (is.null(diagnostics$fit)) NULL else {
    diagnostics$fit[diagnostics$fit$model == "full", , drop = FALSE]
  }
  if (!is.null(full_fit) && nrow(full_fit)) {
    cat("OOF R2:    ", sprintf("%.6g", full_fit$median_r2), " (median), ",
        sprintf("%.1f%%", 100 * full_fit$share_r2_positive), " above zero\n",
        sep = "")
  }
  full_alpha <- if (is.null(diagnostics$alpha)) NULL else {
    diagnostics$alpha[diagnostics$alpha$model == "full", , drop = FALSE]
  }
  if (!is.null(full_alpha) && nrow(full_alpha)) {
    if (full_alpha$n_alpha_grid < 2L) {
      cat("Alpha:     fixed at ",
          format(full_alpha$alpha_grid_min, digits = 4L), "\n", sep = "")
    } else {
      cat("Alpha:     modal ", format(full_alpha$modal_alpha, digits = 4L),
          " in [", format(full_alpha$alpha_grid_min, digits = 4L), ", ",
          format(full_alpha$alpha_grid_max, digits = 4L), "]; ",
          sprintf("%.1f%%", 100 * full_alpha$share_at_grid_max),
          " at the ceiling, ",
          sprintf("%.1f%%", 100 * full_alpha$share_at_grid_min),
          " at the floor\n", sep = "")
    }
    if (isTRUE(full_alpha$saturated)) {
      cat("           grid saturated: widen the alpha grid\n")
    }
  }
  if (!is.null(x$predictive_leave_one_band_out)) {
    cat("Dropouts:  ", paste(
      x$predictive_leave_one_band_out$provenance$requested_bands,
      collapse = ", "
    ), " (predictive OOF delta R2)\n", sep = "")
  }
  invisible(x)
}

#' @export
run_searchlight.banded_ridge_model <- function(model_spec, radius, ...) {
  stop(
    "banded_ridge_model does not support overlapping run_searchlight() in v1: ",
    "each response must be estimated exactly once. Use run_banded_ridge(), or ",
    "run_regional() only with a non-overlapping region mask.",
    call. = FALSE
  )
}

#' @export
run_regional.banded_ridge_model <- function(model_spec, region_mask,
                                            backend = c("default", "shard", "auto"),
                                            ...) {
  if (!inherits(region_mask, "NeuroVol") ||
      !identical(dim(region_mask)[1:3], dim(model_spec$dataset$mask)[1:3])) {
    stop("run_regional.banded_ridge_model requires a spatially matching NeuroVol region mask.",
         call. = FALSE)
  }
  active <- model_spec$active_response_ids
  labels <- as.integer(region_mask[active])
  if (anyNA(labels) || any(labels <= 0L)) {
    stop("banded-ridge regional batching requires every active response in exactly one positive region.",
         call. = FALSE)
  }
  partitions <- split(active, labels)
  run_banded_ridge(model_spec, response_partitions = unname(partitions), ...)
}
