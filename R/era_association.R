# Internal helpers for first-order ERA item associations.

#' @keywords internal
#' @noRd
.era_rhs_formula <- function(x, arg) {
  if (is.null(x)) return(NULL)
  if (!inherits(x, "formula") || length(x) != 2L) {
    stop(sprintf("`%s` must be a right-hand-side formula such as `~ vividness`.", arg),
         call. = FALSE)
  }
  x
}

#' @keywords internal
#' @noRd
.era_observation_count <- function(x) {
  dx <- dim(x)
  if (is.null(dx) || length(dx) < 2L) {
    stop("ERA-RSA data must expose an observation dimension.", call. = FALSE)
  }
  as.integer(dx[[length(dx)]])
}

#' Normalize combined or externally split ERA inputs
#'
#' A combined input is split by an explicit phase variable. Pairing is never
#' inferred from row adjacency, so alternating E/R blocks and separated phase
#' blocks follow the same path.
#'
#' @keywords internal
#' @noRd
.era_prepare_inputs <- function(dataset,
                                design,
                                key_var,
                                phase_var = NULL,
                                encoding_level = NULL,
                                retrieval_level = NULL) {
  has_external <- !is.null(dataset$test_data)

  if (has_external) {
    if (is.null(design$test_design)) {
      stop("ERA-RSA external test data require `design$test_design`.", call. = FALSE)
    }
    if (.era_observation_count(dataset$train_data) != nrow(design$train_design) ||
        .era_observation_count(dataset$test_data) != nrow(design$test_design)) {
      stop("ERA-RSA data observations must match the train/test design rows.", call. = FALSE)
    }

    out_design <- design
    if (is.null(out_design$y_test)) {
      out_design$y_test <- factor(parse_variable(key_var, out_design$test_design))
    }
    return(list(
      dataset = dataset,
      design = out_design,
      combined = FALSE,
      encoding_level = encoding_level %||% "encoding",
      retrieval_level = retrieval_level %||% "retrieval"
    ))
  }

  if (is.null(phase_var)) {
    stop(
      "ERA-RSA combined data require `phase_var`, `encoding_level`, and `retrieval_level`.",
      call. = FALSE
    )
  }
  if (!is.null(design$test_design)) {
    stop("ERA-RSA found `design$test_design` but no `dataset$test_data`.", call. = FALSE)
  }
  if (.era_observation_count(dataset$train_data) != nrow(design$train_design)) {
    stop("ERA-RSA combined data observations must match `design$train_design` rows.",
         call. = FALSE)
  }

  phase <- parse_variable(phase_var, design$train_design)
  if (anyNA(phase)) {
    stop("ERA-RSA `phase_var` contains missing values.", call. = FALSE)
  }
  observed <- if (is.factor(phase)) {
    levels(droplevels(phase))
  } else {
    unique(as.character(phase))
  }
  if (is.null(encoding_level) || is.null(retrieval_level)) {
    if (length(observed) != 2L) {
      stop("ERA-RSA must be given explicit encoding/retrieval levels when `phase_var` has other than two levels.",
           call. = FALSE)
    }
    if (is.null(encoding_level)) encoding_level <- observed[[1L]]
    if (is.null(retrieval_level)) retrieval_level <- observed[[2L]]
  }
  if (identical(as.character(encoding_level), as.character(retrieval_level))) {
    stop("ERA-RSA encoding and retrieval levels must differ.", call. = FALSE)
  }

  phase_chr <- as.character(phase)
  enc_idx <- which(phase_chr == as.character(encoding_level))
  ret_idx <- which(phase_chr == as.character(retrieval_level))
  if (!length(enc_idx) || !length(ret_idx)) {
    stop("ERA-RSA encoding/retrieval levels must each select at least one observation.",
         call. = FALSE)
  }

  if (!inherits(dataset, "mvpa_image_dataset")) {
    stop(
      paste0(
        "Automatic ERA-RSA phase splitting currently supports image datasets. ",
        "For other dataset backends, supply pre-split train/test data."
      ),
      call. = FALSE
    )
  }
  enc_data <- neuroim2::sub_vector(dataset$train_data, enc_idx)
  ret_data <- neuroim2::sub_vector(dataset$train_data, ret_idx)
  split_dataset <- mvpa_dataset(enc_data, test_data = ret_data, mask = dataset$mask)

  split_design <- design
  split_design$train_design <- design$train_design[enc_idx, , drop = FALSE]
  split_design$test_design <- design$train_design[ret_idx, , drop = FALSE]
  all_y <- design$cv_labels %||% parse_variable(key_var, design$train_design)
  split_design$cv_labels <- subset_y(all_y, enc_idx)
  split_design$targets <- subset_y(design$targets %||% all_y, enc_idx)
  split_design$y_train <- split_design$cv_labels
  split_design$y_test <- subset_y(all_y, ret_idx)
  if (!is.null(design$block_var) && length(design$block_var) == nrow(design$train_design)) {
    split_design$block_var <- subset_y(design$block_var, enc_idx)
  }
  split_design$split_by <- NULL
  split_design$split_groups <- NULL

  list(
    dataset = split_dataset,
    design = split_design,
    combined = TRUE,
    encoding_level = encoding_level,
    retrieval_level = retrieval_level,
    encoding_indices = enc_idx,
    retrieval_indices = ret_idx
  )
}

#' @keywords internal
#' @noRd
.era_pairing_info <- function(design, key_var, pairing = c("average", "one_to_one")) {
  pairing <- match.arg(pairing)
  key_enc <- as.character(parse_variable(key_var, design$train_design))
  key_ret <- as.character(parse_variable(key_var, design$test_design))
  if (anyNA(key_enc) || anyNA(key_ret) || any(!nzchar(key_enc)) || any(!nzchar(key_ret))) {
    stop("ERA-RSA item keys must be non-missing and non-empty in both phases.", call. = FALSE)
  }

  enc_counts <- table(key_enc)
  ret_counts <- table(key_ret)
  enc_set <- names(enc_counts)
  ret_set <- names(ret_counts)
  common <- sort(intersect(enc_set, ret_set))
  if (length(common) < 2L) {
    stop("ERA-RSA requires at least two item keys represented in both phases.", call. = FALSE)
  }

  if (identical(pairing, "one_to_one")) {
    duplicate_enc <- names(enc_counts)[enc_counts != 1L]
    duplicate_ret <- names(ret_counts)[ret_counts != 1L]
    missing_enc <- setdiff(ret_set, enc_set)
    missing_ret <- setdiff(enc_set, ret_set)
    if (length(duplicate_enc) || length(duplicate_ret) || length(missing_enc) || length(missing_ret)) {
      detail <- c(
        if (length(duplicate_enc)) paste0("encoding count != 1: ", paste(duplicate_enc, collapse = ", ")),
        if (length(duplicate_ret)) paste0("retrieval count != 1: ", paste(duplicate_ret, collapse = ", ")),
        if (length(missing_enc)) paste0("missing from encoding: ", paste(missing_enc, collapse = ", ")),
        if (length(missing_ret)) paste0("missing from retrieval: ", paste(missing_ret, collapse = ", "))
      )
      stop(
        paste0("ERA-RSA one-to-one pairing failed (", paste(detail, collapse = "; "), ")."),
        call. = FALSE
      )
    }
  }

  list(
    pairing = pairing,
    key_enc = key_enc,
    key_ret = key_ret,
    common_keys = common,
    enc_rows = lapply(common, function(key) which(key_enc == key)),
    ret_rows = lapply(common, function(key) which(key_ret == key)),
    enc_counts = enc_counts,
    ret_counts = ret_counts
  )
}

#' Build aligned item prototypes without reparsing design metadata
#'
#' The row maps are invariant across searchlight spheres and are prepared once
#' by \code{.era_pairing_info()}. Strict one-to-one input therefore reduces to
#' direct row indexing; repeated-item input retains the established group-mean
#' estimand.
#'
#' @keywords internal
#' @noRd
.era_item_prototypes <- function(Xenc, Xret, pairing_info) {
  if (nrow(Xenc) != length(pairing_info$key_enc) ||
      nrow(Xret) != length(pairing_info$key_ret)) {
    stop("ERA-RSA ROI observations do not match the prepared item pairing.",
         call. = FALSE)
  }

  keys <- pairing_info$common_keys
  if (identical(pairing_info$pairing, "one_to_one")) {
    enc_idx <- vapply(pairing_info$enc_rows, `[[`, integer(1L), 1L)
    ret_idx <- vapply(pairing_info$ret_rows, `[[`, integer(1L), 1L)
    E <- Xenc[enc_idx, , drop = FALSE]
    R <- Xret[ret_idx, , drop = FALSE]
  } else {
    mean_rows <- function(X, rows) {
      out <- matrix(NA_real_, nrow = length(rows), ncol = ncol(X))
      for (i in seq_along(rows)) {
        idx <- rows[[i]]
        if (length(idx) == 1L) {
          out[i, ] <- X[idx, ]
        } else {
          out[i, ] <- colMeans(X[idx, , drop = FALSE], na.rm = FALSE)
        }
      }
      out
    }
    E <- mean_rows(Xenc, pairing_info$enc_rows)
    R <- mean_rows(Xret, pairing_info$ret_rows)
  }

  rownames(E) <- rownames(R) <- keys
  list(E = E, R = R, keys = keys)
}

#' @keywords internal
#' @noRd
.era_metric_labels <- function(x, arg) {
  out <- sanitize(make.names(x, unique = FALSE))
  if (any(!nzchar(out)) || anyDuplicated(out)) {
    stop(sprintf("`%s` produces empty or colliding metric names after sanitization.", arg),
         call. = FALSE)
  }
  out
}

#' Collapse retrieval metadata to one constant value per item
#'
#' Pattern averaging is compatible with association metadata only when each
#' retrieval-side variable is constant within an item. This avoids silently
#' changing a trial-level estimand by averaging ratings.
#'
#' @keywords internal
#' @noRd
.era_constant_item_column <- function(x, keys, common_keys, label) {
  first <- match(common_keys, keys)
  out <- x[first]
  for (i in seq_along(common_keys)) {
    vals <- x[keys == common_keys[[i]]]
    present <- vals[!is.na(vals)]
    distinct <- unique(as.character(present))
    if (length(distinct) > 1L) {
      stop(sprintf(
        "ERA-RSA association variable `%s` is not constant within retrieval item `%s`.",
        label, common_keys[[i]]
      ), call. = FALSE)
    }
    if (length(present)) out[[i]] <- present[[1L]]
  }
  out
}

#' @keywords internal
#' @noRd
.era_association_item_data <- function(design, key_var, common_keys, formulas) {
  formulas <- Filter(Negate(is.null), formulas)
  if (!length(formulas)) return(NULL)

  needed <- unique(unlist(lapply(formulas, all.vars), use.names = FALSE))
  missing_cols <- setdiff(needed, names(design$test_design))
  if (length(missing_cols)) {
    stop(sprintf(
      "ERA-RSA association variable(s) not found in the retrieval design: %s.",
      paste(missing_cols, collapse = ", ")
    ), call. = FALSE)
  }

  keys <- as.character(parse_variable(key_var, design$test_design))
  out <- data.frame(.era_key = common_keys, stringsAsFactors = FALSE)
  for (nm in needed) {
    out[[nm]] <- .era_constant_item_column(
      design$test_design[[nm]], keys, common_keys, nm
    )
  }
  rownames(out) <- common_keys
  out
}

#' @keywords internal
#' @noRd
.era_prepare_correlates <- function(formula, item_data) {
  formula <- .era_rhs_formula(formula, "era_correlates")
  if (is.null(formula)) return(NULL)
  tt <- stats::terms(formula)
  labels <- attr(tt, "term.labels")
  if (!length(labels)) {
    stop("`era_correlates` must select at least one numeric retrieval variable.",
         call. = FALSE)
  }
  if (!all(labels %in% names(item_data))) {
    stop("`era_correlates` currently accepts direct numeric column names only.",
         call. = FALSE)
  }
  bad <- labels[!vapply(item_data[labels], is.numeric, logical(1L))]
  if (length(bad)) {
    stop(sprintf("`era_correlates` variables must be numeric: %s.", paste(bad, collapse = ", ")),
         call. = FALSE)
  }
  metric_labels <- .era_metric_labels(labels, "era_correlates")
  values <- as.data.frame(item_data[labels], optional = TRUE)
  rownames(values) <- item_data$.era_key
  list(values = values, labels = labels, metric_labels = metric_labels)
}

#' @keywords internal
#' @noRd
.era_prepare_association <- function(formula, effects, item_data) {
  formula <- .era_rhs_formula(formula, "era_association")
  if (is.null(formula)) {
    if (!is.null(effects)) {
      stop("`era_effects` requires `era_association`.", call. = FALSE)
    }
    return(NULL)
  }
  effects <- .era_rhs_formula(effects, "era_effects")
  if (is.null(effects)) {
    stop("`era_association` requires `era_effects` to identify directional focal terms.",
         call. = FALSE)
  }

  tt <- stats::terms(formula)
  if (!isTRUE(attr(tt, "intercept") == 1L)) {
    stop("`era_association` must include an intercept.", call. = FALSE)
  }
  term_labels <- attr(tt, "term.labels")
  if (!length(term_labels)) {
    stop("`era_association` must contain at least one model term.", call. = FALSE)
  }

  mf <- stats::model.frame(
    tt,
    data = item_data,
    na.action = stats::na.pass,
    drop.unused.levels = FALSE
  )
  X <- stats::model.matrix(tt, data = mf)
  if (nrow(X) != nrow(item_data)) {
    stop("ERA-RSA association model matrix lost item rows.", call. = FALSE)
  }
  rownames(X) <- item_data$.era_key

  effect_labels <- attr(stats::terms(effects), "term.labels")
  missing_effects <- setdiff(effect_labels, term_labels)
  if (!length(effect_labels) || length(missing_effects)) {
    stop(sprintf(
      "`era_effects` terms must occur exactly in `era_association`%s.",
      if (length(missing_effects)) paste0(": ", paste(missing_effects, collapse = ", ")) else ""
    ), call. = FALSE)
  }

  assignment <- attr(X, "assign")
  effect_columns <- integer(length(effect_labels))
  for (i in seq_along(effect_labels)) {
    term_id <- match(effect_labels[[i]], term_labels)
    cols <- which(assignment == term_id)
    if (length(cols) != 1L) {
      stop(sprintf(
        paste0(
          "ERA-RSA focal term `%s` has %d model-matrix columns. ",
          "Signed part-r requires a one-degree-of-freedom term; supply an explicit contrast column."
        ),
        effect_labels[[i]], length(cols)
      ), call. = FALSE)
    }
    effect_columns[[i]] <- cols[[1L]]
  }

  list(
    matrix = X,
    effect_labels = effect_labels,
    effect_metric_labels = .era_metric_labels(effect_labels, "era_effects"),
    effect_columns = effect_columns
  )
}

#' @keywords internal
#' @noRd
.era_cross_similarity <- function(E, R, method = c("pearson", "spearman"), min_voxels = 3L) {
  method <- match.arg(method)
  S <- suppressWarnings(stats::cor(
    t(E), t(R), method = method, use = "pairwise.complete.obs"
  ))
  S <- as.matrix(S)
  finite_counts <- (is.finite(E) * 1L) %*% t(is.finite(R) * 1L)
  S[finite_counts < min_voxels] <- NA_real_
  rownames(S) <- rownames(E)
  colnames(S) <- rownames(R)
  S
}

#' Normalize finite item patterns for correlation products
#'
#' @keywords internal
#' @noRd
.era_normalize_similarity_rows <- function(X, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  X <- as.matrix(X)
  if (identical(method, "spearman")) {
    X <- t(apply(X, 1L, rank, ties.method = "average"))
  }
  X <- sweep(X, 1L, rowMeans(X), "-")
  norms <- sqrt(rowSums(X * X))
  valid <- is.finite(norms) & norms > 0
  Z <- matrix(0, nrow = nrow(X), ncol = ncol(X), dimnames = dimnames(X))
  if (any(valid)) {
    Z[valid, ] <- X[valid, , drop = FALSE] / norms[valid]
  }
  list(values = Z, valid = valid)
}

#' Compute per-retrieval matched similarity and item-specific background
#'
#' For finite data this uses normalized pattern sums and is linear in the
#' number of items for the item-level outputs. A full K-by-K matrix is formed
#' only when identification metrics request it. Non-finite data use the exact
#' pairwise-complete matrix implementation so missingness semantics are not
#' approximated.
#'
#' @keywords internal
#' @noRd
.era_item_similarity_scores <- function(E,
                                        R,
                                        method = c("pearson", "spearman"),
                                        min_voxels = 2L,
                                        need_matrix = FALSE,
                                        item_block = NULL) {
  method <- match.arg(method)
  E <- as.matrix(E)
  R <- as.matrix(R)
  if (nrow(E) != nrow(R) || ncol(E) != ncol(R)) {
    stop("ERA-RSA item prototypes must have matching dimensions.", call. = FALSE)
  }
  K <- nrow(E)
  if (K < 2L) {
    stop("ERA-RSA item specificity requires at least two matched items.", call. = FALSE)
  }

  mean_finite <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) mean(x) else NA_real_
  }
  block <- if (is.null(item_block)) NULL else as.character(item_block)
  if (!is.null(block) && length(block) != K) {
    stop("ERA-RSA `item_block` must align with the matched item set.", call. = FALSE)
  }

  all_finite <- all(is.finite(E)) && all(is.finite(R))
  if (!all_finite) {
    S <- .era_cross_similarity(E, R, method = method, min_voxels = min_voxels)
    matched <- diag(S)
    background <- vapply(seq_len(K), function(i) mean_finite(S[-i, i]), numeric(1L))
    same_specificity <- diff_specificity <- rep(NA_real_, K)
    if (!is.null(block)) {
      for (i in seq_len(K)) {
        if (is.na(block[[i]]) || !is.finite(matched[[i]])) next
        same_idx <- which(!is.na(block) & block == block[[i]] & seq_len(K) != i)
        diff_idx <- which(!is.na(block) & block != block[[i]])
        same_bg <- mean_finite(S[same_idx, i])
        diff_bg <- mean_finite(S[diff_idx, i])
        if (is.finite(same_bg)) same_specificity[[i]] <- matched[[i]] - same_bg
        if (is.finite(diff_bg)) diff_specificity[[i]] <- matched[[i]] - diff_bg
      }
    }
    return(list(
      matrix = if (isTRUE(need_matrix)) S else NULL,
      matched = matched,
      background = background,
      specificity = matched - background,
      same_block_specificity = same_specificity,
      diff_block_specificity = diff_specificity
    ))
  }

  enc <- .era_normalize_similarity_rows(E, method = method)
  ret <- .era_normalize_similarity_rows(R, method = method)
  if (ncol(E) < min_voxels) {
    enc$valid[] <- FALSE
    ret$valid[] <- FALSE
  }
  valid_pair <- enc$valid & ret$valid
  matched <- rep(NA_real_, K)
  if (any(valid_pair)) {
    matched[valid_pair] <- rowSums(
      enc$values[valid_pair, , drop = FALSE] *
        ret$values[valid_pair, , drop = FALSE]
    )
  }

  enc_sum <- colSums(enc$values[enc$valid, , drop = FALSE])
  enc_count <- sum(enc$valid)
  background <- rep(NA_real_, K)
  valid_background <- ret$valid & (enc_count - as.integer(enc$valid) > 0L)
  if (any(valid_background)) {
    totals <- as.vector(ret$values[valid_background, , drop = FALSE] %*% enc_sum)
    matched_part <- ifelse(
      enc$valid[valid_background], matched[valid_background], 0
    )
    background[valid_background] <-
      (totals - matched_part) /
      (enc_count - as.integer(enc$valid[valid_background]))
  }

  same_specificity <- diff_specificity <- rep(NA_real_, K)
  if (!is.null(block)) {
    block_levels <- unique(block[!is.na(block)])
    block_id <- match(block, block_levels)
    block_sums <- matrix(0, nrow = length(block_levels), ncol = ncol(E))
    block_counts <- integer(length(block_levels))
    for (g in seq_along(block_levels)) {
      idx <- which(enc$valid & block_id == g)
      block_counts[[g]] <- length(idx)
      if (length(idx)) {
        block_sums[g, ] <- colSums(enc$values[idx, , drop = FALSE])
      }
    }
    for (i in seq_len(K)) {
      g <- block_id[[i]]
      if (is.na(g) || !ret$valid[[i]] || !is.finite(matched[[i]])) next
      same_n <- block_counts[[g]] - as.integer(enc$valid[[i]])
      if (same_n > 0L) {
        same_sum <- block_sums[g, ] - if (enc$valid[[i]]) enc$values[i, ] else 0
        same_specificity[[i]] <- matched[[i]] - sum(ret$values[i, ] * same_sum) / same_n
      }
      diff_n <- enc_count - block_counts[[g]]
      if (diff_n > 0L) {
        diff_sum <- enc_sum - block_sums[g, ]
        diff_specificity[[i]] <- matched[[i]] - sum(ret$values[i, ] * diff_sum) / diff_n
      }
    }
  }

  S <- NULL
  if (isTRUE(need_matrix)) {
    S <- tcrossprod(enc$values, ret$values)
    S[!outer(enc$valid, ret$valid, `&`)] <- NA_real_
    rownames(S) <- rownames(E)
    colnames(S) <- rownames(R)
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

#' @keywords internal
#' @noRd
.era_correlate_metrics <- function(similarity, keys, spec, method, min_complete) {
  if (is.null(spec)) return(numeric())
  idx <- match(keys, rownames(spec$values))
  out <- numeric()
  for (i in seq_along(spec$labels)) {
    value <- spec$values[[i]][idx]
    ok <- is.finite(similarity) & is.finite(value)
    n <- sum(ok)
    estimate <- NA_real_
    if (n >= min_complete && stats::sd(similarity[ok]) > 0 && stats::sd(value[ok]) > 0) {
      estimate <- suppressWarnings(stats::cor(similarity[ok], value[ok], method = method))
    }
    stem <- paste0("era_", spec$metric_labels[[i]])
    out[[paste0(stem, "_cor")]] <- estimate
    out[[paste0(stem, "_n")]] <- n
  }
  unlist(out, use.names = TRUE)
}

#' Compute signed semi-partial correlations for focal association terms
#'
#' For each one-df focal term, residualize its model-matrix column against all
#' remaining columns and correlate that residual with ERA similarity. This is
#' signed part-r; its square equals the term's incremental R-squared, which is
#' intentionally not emitted as a map.
#'
#' @keywords internal
#' @noRd
.era_association_metrics <- function(similarity, keys, spec, min_complete) {
  if (is.null(spec)) return(numeric())
  idx <- match(keys, rownames(spec$matrix))
  X <- spec$matrix[idx, , drop = FALSE]
  keep <- is.finite(similarity) & apply(X, 1L, function(row) all(is.finite(row)))
  n <- sum(keep)

  values <- setNames(rep(NA_real_, length(spec$effect_columns)),
                     paste0("era_assoc_part_r_", spec$effect_metric_labels))
  df_resid <- NA_real_
  if (n >= min_complete) {
    Xk <- X[keep, , drop = FALSE]
    yk <- similarity[keep]
    full_rank <- qr(Xk)$rank
    df_resid <- n - full_rank
    if (df_resid > 0L && stats::sd(yk) > 0) {
      for (i in seq_along(spec$effect_columns)) {
        col <- spec$effect_columns[[i]]
        x <- Xk[, col]
        Z <- Xk[, -col, drop = FALSE]
        x_resid <- stats::lm.fit(x = Z, y = x)$residuals
        residual_norm <- sqrt(sum(x_resid^2))
        x_norm <- sqrt(sum((x - mean(x))^2))
        estimable <- residual_norm > sqrt(.Machine$double.eps) * max(1, x_norm)
        if (all(is.finite(x_resid)) && isTRUE(estimable) && stats::sd(x_resid) > 0) {
          values[[i]] <- suppressWarnings(stats::cor(yk, x_resid))
        }
      }
    }
  }

  c(values, era_assoc_n = n, era_assoc_df_resid = df_resid)
}

#' @keywords internal
#' @noRd
.era_correlate_metric_names <- function(spec) {
  if (is.null(spec)) return(character())
  unlist(lapply(spec$metric_labels, function(label) {
    c(paste0("era_", label, "_cor"), paste0("era_", label, "_n"))
  }), use.names = FALSE)
}

#' @keywords internal
#' @noRd
.era_association_metric_names <- function(spec) {
  if (is.null(spec)) return(character())
  c(
    paste0("era_assoc_part_r_", spec$effect_metric_labels),
    "era_assoc_n",
    "era_assoc_df_resid"
  )
}
