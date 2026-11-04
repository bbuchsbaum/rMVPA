#' Joint nonmatch endpoint decomposition
#'
#' Fit a joint additive row/column model using only nonmatching cells within
#' explicitly specified complete competitor blocks. Templates (encoding or
#' predicted patterns) are rows; retrieval patterns are columns.
#'
#' @param encoding,retrieval Numeric item-by-feature matrices. Feature columns
#'   must have identical identities/order. Rows are aligned by item keys.
#' @param item_key Unique, nonmissing item keys. A vector explicitly declares the
#'   same row order in both inputs; a list with `encoding` and `retrieval` vectors
#'   declares separate orders. If NULL, both inputs must have row names (row and
#'   column names for `similarity_matrix`). No row-order pairing is inferred.
#' @param competitor_block Block labels, either in encoding-key order or named
#'   by item key. Required: use an explicit constant for one common pool. For run
#'   intersections, construct labels from phase-scoped encoding and retrieval runs.
#' @param similarity_matrix Optional square, pre-link score matrix instead of
#'   patterns, with encoding rows and retrieval columns. Scores must be finite
#'   within each eligible block; cross-block cells are unused. Probabilities and
#'   ranks do not inherit the linear similarity decomposition.
#' @param template_source Provenance label: `"observed_encoding"`,
#'   `"held_out_prediction"`, or `"supplied_scores"`. A label is a declaration,
#'   not independent verification of held-out ancestry.
#' @param encoding_run,retrieval_run Optional phase-scoped run labels, aligned
#'   like `competitor_block`. Retained as metadata; they do not change the pool.
#' @param return_matrices Retain within-block similarity matrices? Default FALSE
#'   uses normalized pattern sums without an item-by-item allocation.
#' @param invalid `"status"` returns NA components with a block reason;
#'   `"error"` stops on invalid blocks. Malformed pairing always errors.
#' @param min_voxels Minimum number of common feature columns, at least two.
#' @param tolerance Relative tolerance for constant patterns and effectively
#'   constant score/behavior vectors. Default `sqrt(.Machine$double.eps)`.
#' @param pair_weights Reserved for a generalized solver. Any non-NULL value is
#'   rejected: only complete, equally weighted square blocks are supported.
#' @return An `era_endpoints` list with item-keyed `items`, block sizes, statuses,
#'   identity errors and weights in `blocks`, `metadata`, and optional `matrices`.
#'   Each valid item satisfies `total = shared + encoding + retrieval + excess`.
#' @details
#' Patterns are spatially centered and unit normalized (spatial Pearson
#' correlation). No across-item mean pattern is subtracted. Nonfinite or constant
#' patterns invalidate their entire block; no pairwise-complete fallback or
#' automatic competitor removal is performed. Blocks of fewer than three items
#' are nonidentifiable. The diagonal never enters endpoint fitting.
#'
#' For n items, k=n-1, off-diagonal row/column means r and c, and off-diagonal
#' grand mean mu, set u=r-mu and v=c-mu. The joint least-squares effects are
#' a=(u+v/k)/(1-1/k^2) and b=(v+u/k)/(1-1/k^2). One-sided column specificity
#' equals n/(n-1)*a + excess, not excess alone.
#' @examples
#' set.seed(92)
#' E <- matrix(rnorm(6 * 10), 6, 10)
#' R <- E + matrix(rnorm(6 * 10), 6, 10)
#' z <- era_endpoint_decomposition(E, R, item_key = letters[1:6],
#'                                 competitor_block = rep("pool", 6))
#' summarize_era_endpoints(z)$table
#' @md
#' @export
era_endpoint_decomposition <- function(
    encoding = NULL, retrieval = NULL, item_key = NULL, competitor_block,
    similarity_matrix = NULL,
    template_source = if (is.null(similarity_matrix)) "observed_encoding" else "supplied_scores",
    encoding_run = NULL, retrieval_run = NULL, return_matrices = FALSE,
    invalid = c("status", "error"), min_voxels = 2L,
    tolerance = sqrt(.Machine$double.eps), pair_weights = NULL) {
  invalid <- match.arg(invalid)
  template_source <- match.arg(template_source,
    c("observed_encoding", "held_out_prediction", "supplied_scores"))
  .endpoint_tolerance(tolerance)
  .endpoint_flag(return_matrices, "return_matrices")
  min_voxels <- .era_count_scalar(min_voxels, "min_voxels", 2L)
  if (!is.null(pair_weights)) {
    stop("Weighted competitor designs are unsupported; use complete equally weighted blocks.",
         call. = FALSE)
  }
  scores <- !is.null(similarity_matrix)
  if (scores) {
    if (!is.null(encoding) || !is.null(retrieval)) {
      stop("Supply patterns or similarity_matrix, not both.", call. = FALSE)
    }
    M <- .endpoint_matrix(similarity_matrix, "similarity_matrix")
    if (nrow(M) != ncol(M)) stop("similarity_matrix must be square.", call. = FALSE)
    nr <- nrow(M)
    keys <- .endpoint_pair_keys(item_key, rownames(M), colnames(M), nr, nr)
    M <- M[, match(keys$encoding, keys$retrieval), drop = FALSE]
  } else {
    E <- .endpoint_matrix(encoding, "encoding")
    R <- .endpoint_matrix(retrieval, "retrieval")
    if (ncol(E) != ncol(R) ||
        !identical(colnames(E), colnames(R))) {
      stop("Pattern feature columns must have matching identities and order.", call. = FALSE)
    }
    keys <- .endpoint_pair_keys(item_key, rownames(E), rownames(R), nrow(E), nrow(R))
    R <- R[match(keys$encoding, keys$retrieval), , drop = FALSE]
    nr <- nrow(E)
  }
  key <- keys$encoding
  block <- as.character(.endpoint_align(competitor_block, key, "competitor_block"))
  if (anyNA(block) || any(!nzchar(block))) {
    stop("competitor_block must identify every item explicitly.", call. = FALSE)
  }
  items <- data.frame(item_key = key, competitor_block = block,
                      total = NA_real_, shared = NA_real_, encoding = NA_real_,
                      retrieval = NA_real_, excess = NA_real_, n_items = NA_integer_,
                      status = NA_character_, template_source = template_source,
                      stringsAsFactors = FALSE)
  for (nm in c("encoding_run", "retrieval_run")) {
    value <- get(nm)
    if (!is.null(value)) items[[nm]] <- .endpoint_align(value, key, nm)
  }
  groups <- split(seq_len(nr), factor(block, levels = unique(block)))
  matrices <- if (return_matrices) vector("list", length(groups)) else NULL
  if (return_matrices) names(matrices) <- names(groups)
  blocks <- vector("list", length(groups))
  for (g in seq_along(groups)) {
    idx <- groups[[g]]
    n <- length(idx)
    status <- if (n < 3L) "too_few_items" else "ok"
    B <- NULL
    if (status == "ok" && scores) {
      B <- M[idx, idx, drop = FALSE]
      if (any(!is.finite(B))) status <- "nonfinite_scores"
    } else if (status == "ok") {
      if (ncol(E) < min_voxels) {
        status <- "too_few_voxels"
      } else if (any(!is.finite(E[idx, ])) || any(!is.finite(R[idx, ]))) {
        status <- "nonfinite_patterns"
      } else {
        en <- .endpoint_normalize(E[idx, , drop = FALSE], tolerance)
        re <- .endpoint_normalize(R[idx, , drop = FALSE], tolerance)
        if (!all(en$valid) || !all(re$valid)) status <- "constant_pattern"
      }
    }
    identity_error <- NA_real_
    if (status == "ok") {
      if (scores || return_matrices) {
        if (!scores) B <- tcrossprod(en$values, re$values)
        d <- diag(B)
        off <- B
        diag(off) <- 0
        row_total <- rowSums(off)
        col_total <- colSums(off)
      } else {
        d <- rowSums(en$values * re$values)
        row_total <- as.vector(en$values %*% colSums(re$values)) - d
        col_total <- as.vector(re$values %*% colSums(en$values)) - d
      }
      k <- n - 1L
      mu <- sum(row_total) / (n * k)
      u <- row_total / k - mu
      v <- col_total / k - mu
      a <- (u + v / k) / (1 - 1 / k^2)
      b <- (v + u / k) / (1 - 1 / k^2)
      e <- d - mu - a - b
      items[idx, c("total", "shared", "encoding", "retrieval", "excess")] <-
        cbind(d, mu, a, b, e)
      identity_error <- max(abs(d - (mu + a + b + e)))
      if (return_matrices) {
        dimnames(B) <- list(key[idx], key[idx])
        matrices[[g]] <- B
      }
    } else if (invalid == "error") {
      stop(sprintf("Endpoint block '%s': %s (n=%d); require >=3 complete valid pairs.",
                   names(groups)[g], status, n), call. = FALSE)
    }
    items$n_items[idx] <- n
    items$status[idx] <- status
    blocks[[g]] <- data.frame(competitor_block = names(groups)[g], n_items = n,
      n_nonmatch = n * (n - 1L), status = status, identity_error = identity_error)
  }
  blocks <- do.call(rbind, blocks)
  valid_n <- sum(blocks$n_items[blocks$status == "ok"])
  blocks$raw_weight <- ifelse(blocks$status == "ok", blocks$n_items / valid_n, 0)
  structure(list(items = items, blocks = blocks, matrices = matrices,
    metadata = list(orientation = "encoding/template rows; retrieval columns",
      similarity = if (scores) "supplied_pre_link_scores" else "spatial_pearson",
      fit_pairs = "nonmatch", method = "joint_additive_complete_block",
      template_source = template_source, pairing = "one_to_one",
      weighting = "equal items; block weights proportional to valid block size",
      tolerance = tolerance, min_voxels = min_voxels)), class = "era_endpoints")
}

#' Summarize endpoint performance or behavioral covariance
#'
#' Raw means average all items in valid complete blocks, retaining the shared
#' baseline. Covariance contributions and separate correlations use exactly the
#' same complete rows and nuisance projection for every component and behavior.
#'
#' @param x An `era_endpoints` object.
#' @param summary `"raw_mean"` or `"covariance_contributions"`.
#' @param behavior Numeric behavior vector, named by item key or in `x$items`
#'   order. Required for covariance contributions.
#' @param nuisance Numeric design matrix, optionally with item-key row names.
#'   An intercept is always included. Missing values exclude the same rows from
#'   every component; redundant columns are handled by a common QR projection.
#' @param adjust_blocks Include competitor-block intercepts? Default TRUE.
#'   If FALSE, shared-baseline covariance is retained explicitly.
#' @param cor_method Method for the separately named component correlations:
#'   `"pearson"` or `"spearman"`. Spearman ranks each component and behavior
#'   before the common nuisance projection. Contributions always use unranked
#'   scores and Pearson scaling; there is no additive Spearman decomposition.
#' @param min_complete Minimum common rows for association, at least three.
#' @param tolerance Constant-vector tolerance; defaults to the fit's tolerance.
#' @return A list with `table`, counts, item-keyed exclusion reasons, status and
#'   identity error. Association tables contain signed `contribution` and
#'   separately normalized `correlation` columns, plus correlation statuses.
#'   Inferential fields are NA: group inference requires independent participants.
#' @details
#' Contributions divide every component-behavior residual dot product by the
#' same norm(total)*norm(behavior). Their sum (including shared when needed)
#' equals the residual total Pearson correlation. They are not variance shares
#' or percentages and can cancel or exceed one in magnitude. Endpoint means are
#' structural zeros only for complete block-weighted raw averaging, not for a
#' behavioral subset. Endpoint RMS is descriptive and includes estimation noise.
#' @md
#' @export
summarize_era_endpoints <- function(
    x, summary = c("raw_mean", "covariance_contributions"), behavior = NULL,
    nuisance = NULL, adjust_blocks = TRUE, cor_method = c("pearson", "spearman"),
    min_complete = 3L, tolerance = x$metadata$tolerance) {
  if (!inherits(x, "era_endpoints")) stop("x must be an era_endpoints object.", call. = FALSE)
  summary <- match.arg(summary)
  cor_method <- match.arg(cor_method)
  .endpoint_tolerance(tolerance)
  .endpoint_flag(adjust_blocks, "adjust_blocks")
  min_complete <- .era_count_scalar(min_complete, "min_complete", 3L)
  items <- x$items
  key <- items$item_key
  components <- c("total", "shared", "encoding", "retrieval", "excess")
  values <- as.matrix(items[, components])
  reason <- ifelse(items$status == "ok", "", items$status)
  if (summary == "raw_mean") {
    if (!is.null(behavior) || !is.null(nuisance)) {
      stop("raw_mean uses complete fitted blocks, without behavior/nuisance subsetting.", call. = FALSE)
    }
    ok <- reason == ""
    means <- if (any(ok)) colMeans(values[ok, , drop = FALSE]) else rep(NA_real_, 5L)
    rms <- if (any(ok)) sqrt(colMeans(values[ok, c("encoding", "retrieval"), drop = FALSE]^2)) else c(NA_real_, NA_real_)
    tab <- data.frame(quantity = c(components, "encoding_rms", "retrieval_rms"),
      estimate = unname(c(means, rms)), structural_zero = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE),
      descriptive = c(FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE),
      zero_test = c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE),
      statistic = NA_real_, p_value = NA_real_)
    return(list(table = tab, n = sum(ok), n_total = length(ok), n_excluded = sum(!ok),
      n_blocks = sum(x$blocks$status == "ok"), blocks = x$blocks,
      exclusions = data.frame(item_key = key, included = ok, reason = reason),
      identity_error = if (any(ok)) unname(abs(means[1] - means[2] - means[5])) else NA_real_,
      status = if (any(ok)) "ok" else "no_valid_blocks"))
  }
  behavior <- .endpoint_align(behavior, key, "behavior")
  if (!is.numeric(behavior)) stop("behavior must be numeric.", call. = FALSE)
  Z <- .endpoint_nuisance(nuisance, key)
  reason <- .endpoint_reason(reason, !is.finite(behavior), "nonfinite_behavior")
  reason <- .endpoint_reason(reason, rowSums(!is.finite(Z)) > 0L, "nonfinite_nuisance")
  ok <- reason == ""
  n <- sum(ok)
  tab <- data.frame(quantity = components, contribution = NA_real_, correlation = NA_real_,
    correlation_status = "insufficient_rows", structural_zero = FALSE,
    statistic = NA_real_, p_value = NA_real_)
  status <- "insufficient_rows"
  identity_error <- NA_real_
  nuisance_rank <- NA_integer_
  shared_absorbed <- NA
  if (n >= min_complete) {
    Z <- cbind(intercept = 1, Z[ok, , drop = FALSE])
    if (adjust_blocks) {
      block <- items$competitor_block[ok]
      Z <- cbind(Z, outer(block, unique(block), `==`) * 1)
    }
    qrZ <- qr(Z)
    nuisance_rank <- qrZ$rank
    original <- cbind(values[ok, , drop = FALSE], behavior = behavior[ok])
    residual <- qr.resid(qrZ, original)
    norms <- sqrt(colSums(residual^2))
    stable <- vapply(seq_len(ncol(original)), function(j) {
      .endpoint_nonconstant(residual[, j], tolerance, original[, j])
    }, logical(1))
    shared_absorbed <- unname(norms[2L] <= 100 * .Machine$double.eps *
      max(1, sqrt(sum(original[, 2L]^2))))
    tab$structural_zero[2L] <- shared_absorbed
    status <- if (n - nuisance_rank < 2L) "insufficient_residual_df" else if (!stable[6L]) {
      "constant_behavior"
    } else if (!stable[1L]) "constant_total" else "ok"
    tab$correlation_status <- status
    if (status == "ok") {
      denom <- norms[1L] * norms[6L]
      tab$contribution <- as.vector(crossprod(residual[, 6L], residual[, 1:5, drop = FALSE])) / denom
      # Keep even tiny covariances: thresholding individual contributions can
      # destroy the identity when total variation is small.
      identity_error <- abs(tab$contribution[1L] - sum(tab$contribution[-1L]))
    }
    corr_residual <- residual
    corr_stable <- stable
    if (cor_method == "spearman") {
      # Check original variation first: ranking numerical roundoff must
      # never manufacture an association for a constant excess or endpoint.
      corr_stable <- vapply(seq_len(ncol(original)), function(j)
        .endpoint_nonconstant(original[, j], tolerance), logical(1))
      corr_residual <- qr.resid(qrZ, apply(original, 2L, rank, ties.method = "average"))
    }
    for (j in seq_len(5L)) {
      why <- if (n - nuisance_rank < 2L) "insufficient_residual_df" else if (!corr_stable[6L]) {
        "constant_behavior"
      } else if (!corr_stable[j]) "constant_component" else "ok"
      if (why == "ok" && (!.endpoint_nonconstant(corr_residual[, j], tolerance) ||
                           !.endpoint_nonconstant(corr_residual[, 6L], tolerance))) {
        why <- "constant_after_projection"
      }
      tab$correlation_status[j] <- why
      if (why == "ok") tab$correlation[j] <- sum(corr_residual[, j] * corr_residual[, 6L]) /
        sqrt(sum(corr_residual[, j]^2) * sum(corr_residual[, 6L]^2))
    }
  }
  list(table = tab, n = n, n_total = length(ok), n_excluded = sum(!ok),
    n_blocks = length(unique(items$competitor_block[ok])), nuisance_rank = nuisance_rank,
    residual_df = n - nuisance_rank, shared_absorbed = shared_absorbed,
    identity_error = identity_error, status = status, correlation_method = cor_method,
    exclusions = data.frame(item_key = key, included = ok, reason = reason),
    subset_component_means = if (n) colMeans(values[ok, , drop = FALSE]) else
      setNames(rep(NA_real_, 5L), components))
}

#' Configure opt-in endpoint scoring for ERA or feature RSA
#'
#' This object supplies item identity, competitor pools and optional behavioral
#' summaries to `era_rsa_model(endpoints=)` or `feature_rsa_model(endpoints=)`.
#' It changes scoring only, not feature-model fitting, tuning or folds.
#' @param item_key Unique keys corresponding to ERA common items or feature
#'   design labels. Repeated items/assessments are unsupported.
#' @param competitor_block,encoding_run,retrieval_run See
#'   [era_endpoint_decomposition()]. Competitor pools are never inferred from
#'   cross-validation folds; feature scoring refuses any pool crossing outer folds.
#' @param behavior,nuisance,adjust_blocks,cor_method,min_complete,tolerance See
#'   [summarize_era_endpoints()]. Behavior is optional and numeric.
#' @param return_items Retain item tables and diagnostic summaries in regional
#'   `$fits[[i]]$endpoints`? Default FALSE; searchlights emit scalars only.
#' @return An `era_endpoint_spec` with validated, item-keyed metadata. Scalar
#'   names use `era_endpoint_` or `feature_endpoint_` followed by `mean_total`,
#'   `mean_shared`, `mean_encoding`, `mean_retrieval`, `mean_excess`,
#'   `encoding_rms`, `retrieval_rms`, `n`, `n_excluded`, `n_blocks`, and
#'   `identity_error`. With behavior, `assoc_n`, `assoc_n_excluded`,
#'   `assoc_identity_error`, `contribution_<component>` and `cor_<component>`
#'   are added. These scalars use the standard spatial export contract.
#' @md
#' @export
era_endpoint_spec <- function(item_key, competitor_block, behavior = NULL,
    nuisance = NULL, adjust_blocks = TRUE, cor_method = c("pearson", "spearman"),
    min_complete = 3L, tolerance = sqrt(.Machine$double.eps), return_items = FALSE,
    encoding_run = NULL, retrieval_run = NULL) {
  key <- .endpoint_keys(item_key, "item_key")
  .endpoint_tolerance(tolerance)
  .endpoint_flag(return_items, "return_items")
  .endpoint_flag(adjust_blocks, "adjust_blocks")
  block <- as.character(.endpoint_align(competitor_block, key, "competitor_block"))
  if (anyNA(block) || any(!nzchar(block))) stop("competitor_block must identify every item.", call. = FALSE)
  if (!is.null(behavior)) {
    behavior <- .endpoint_align(behavior, key, "behavior")
    if (!is.numeric(behavior)) stop("behavior must be numeric.", call. = FALSE)
  } else if (!is.null(nuisance)) stop("nuisance requires behavior.", call. = FALSE)
  structure(list(item_key = key, competitor_block = block,
    behavior = behavior, nuisance = .endpoint_nuisance(nuisance, key),
    adjust_blocks = adjust_blocks, cor_method = match.arg(cor_method),
    min_complete = .era_count_scalar(min_complete, "min_complete", 3L),
    tolerance = tolerance, return_items = return_items,
    encoding_run = if (!is.null(encoding_run)) .endpoint_align(encoding_run, key, "encoding_run"),
    retrieval_run = if (!is.null(retrieval_run)) .endpoint_align(retrieval_run, key, "retrieval_run")),
    class = "era_endpoint_spec")
}

.endpoint_keys <- function(x, label) {
  if (is.null(x) || !is.atomic(x) || !is.null(dim(x))) stop(paste(label, "must provide explicit item keys."), call. = FALSE)
  x <- as.character(x)
  if (!length(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x)) {
    stop(paste(label, "must contain unique, nonmissing, nonempty keys."), call. = FALSE)
  }
  x
}

.endpoint_pair_keys <- function(item_key, enc_names, ret_names, ne, nr) {
  if (is.null(item_key)) item_key <- list(encoding = enc_names, retrieval = ret_names)
  if (!is.list(item_key)) item_key <- list(encoding = item_key, retrieval = item_key)
  e <- .endpoint_keys(item_key$encoding, "encoding item keys")
  r <- .endpoint_keys(item_key$retrieval, "retrieval item keys")
  if (length(e) != ne || length(r) != nr || !setequal(e, r)) {
    stop("Require the same one-to-one item keys in both inputs; no recycling or aggregation.", call. = FALSE)
  }
  list(encoding = e, retrieval = r)
}

.endpoint_align <- function(x, key, label) {
  if (is.null(x) || !is.atomic(x) || !is.null(dim(x)) || length(x) != length(key)) {
    stop(paste(label, "must have one value per item key."), call. = FALSE)
  }
  if (!is.null(names(x))) {
    supplied <- .endpoint_keys(names(x), paste(label, "names"))
    if (!setequal(supplied, key)) stop(paste(label, "names must match item keys."), call. = FALSE)
    x <- x[match(key, supplied)]
  }
  unname(x)
}

.endpoint_matrix <- function(x, label) {
  if (!is.matrix(x) || !is.numeric(x) || !nrow(x)) {
    stop(paste(label, "must be a nonempty numeric matrix."), call. = FALSE)
  }
  x
}

.endpoint_nuisance <- function(x, key) {
  if (is.null(x)) return(matrix(numeric(), length(key), 0L))
  x <- .endpoint_matrix(x, "nuisance")
  if (nrow(x) != length(key)) stop("nuisance must have one row per item key.", call. = FALSE)
  if (!is.null(rownames(x))) {
    rn <- .endpoint_keys(rownames(x), "nuisance row names")
    if (!setequal(rn, key)) stop("nuisance row names must match item keys.", call. = FALSE)
    x <- x[match(key, rn), , drop = FALSE]
  }
  x
}

.endpoint_flag <- function(x, label) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) stop(paste(label, "must be TRUE or FALSE."), call. = FALSE)
}

.endpoint_tolerance <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1) {
    stop("tolerance must be a finite number between zero and one.", call. = FALSE)
  }
}

.endpoint_normalize <- function(X, tolerance) {
  # Scale first to avoid overflow; validity is invariant to pattern units.
  scale <- apply(abs(X), 1L, max)
  scaled <- X / ifelse(scale > 0, scale, 1)
  centered <- scaled - rowMeans(scaled)
  norms <- sqrt(rowSums(centered^2))
  valid <- norms > tolerance * sqrt(ncol(X))
  list(values = centered / ifelse(valid, norms, 1), valid = valid)
}

.endpoint_nonconstant <- function(x, tolerance, reference = x) {
  sqrt(mean((x - mean(x))^2)) > tolerance * max(1, max(abs(reference)))
}

.endpoint_reason <- function(reason, bad, label) {
  reason[bad] <- ifelse(nzchar(reason[bad]), paste(reason[bad], label, sep = ";"), label)
  reason
}

.endpoint_validate_spec <- function(spec, keys) {
  if (is.null(spec)) return(invisible(NULL))
  if (!inherits(spec, "era_endpoint_spec")) stop("endpoints must be created by era_endpoint_spec().", call. = FALSE)
  keys <- .endpoint_keys(keys, "model item keys")
  if (!setequal(keys, spec$item_key)) stop("Endpoint item keys must match model item keys exactly.", call. = FALSE)
  invisible(NULL)
}

.endpoint_metric_names <- function(spec, prefix) {
  if (is.null(spec)) return(character())
  nms <- c(paste0("mean_", c("total", "shared", "encoding", "retrieval", "excess")),
    "encoding_rms", "retrieval_rms", "n", "n_excluded", "n_blocks", "identity_error")
  if (!is.null(spec$behavior)) nms <- c(nms, "assoc_n", "assoc_n_excluded", "assoc_identity_error",
    paste0("contribution_", c("total", "shared", "encoding", "retrieval", "excess")),
    paste0("cor_", c("total", "shared", "encoding", "retrieval", "excess")))
  paste0(prefix, nms)
}

.endpoint_score <- function(E, R, keys, spec, prefix, source, fold_id = NULL,
                            min_voxels = 2L) {
  .endpoint_validate_spec(spec, keys)
  idx <- match(keys, spec$item_key)
  block <- spec$competitor_block[idx]
  if (!is.null(fold_id)) {
    if (length(fold_id) != length(keys) || anyNA(fold_id) ||
        any(vapply(split(fold_id, block), function(x) length(unique(x)) != 1L, logical(1)))) {
      stop("Feature endpoint competitor blocks must be nested within outer assessment folds.", call. = FALSE)
    }
  }
  # Input identity is supplied by the model's established item/fold ledger.
  dimnames(E) <- dimnames(R) <- NULL
  z <- era_endpoint_decomposition(E, R, item_key = keys, competitor_block = block,
    template_source = source, tolerance = spec$tolerance, min_voxels = min_voxels,
    encoding_run = spec$encoding_run[idx], retrieval_run = spec$retrieval_run[idx])
  if (!is.null(fold_id)) z$items$fold_id <- fold_id
  raw <- summarize_era_endpoints(z)
  metrics <- c(raw$table$estimate, raw$n, raw$n_excluded, raw$n_blocks, raw$identity_error)
  assoc <- NULL
  if (!is.null(spec$behavior)) {
    nuisance <- spec$nuisance[idx, , drop = FALSE]
    rownames(nuisance) <- NULL
    assoc <- summarize_era_endpoints(z, "covariance_contributions",
      behavior = spec$behavior[idx], nuisance = nuisance,
      adjust_blocks = spec$adjust_blocks, cor_method = spec$cor_method,
      min_complete = spec$min_complete)
    metrics <- c(metrics, assoc$n, assoc$n_excluded, assoc$identity_error,
                 assoc$table$contribution, assoc$table$correlation)
  }
  names(metrics) <- .endpoint_metric_names(spec, prefix)
  list(metrics = metrics, decomposition = z, raw = raw, association = assoc)
}
