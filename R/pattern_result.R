# Restriction is a marginal model, not a crop of the whole-domain weights.
.pattern_restrict_fit <- function(fit, region) {
  keep <- which(fit$feature_index %in% region)
  out <- fit
  out$A <- fit$A[keep, , drop = FALSE]
  out$feature_index <- fit$feature_index[keep]
  out$noise <- if (length(keep)) .noise_restrict(fit$noise, keep) else new_pattern_noise(numeric(0), type = "diag")
  out$precision_A <- .noise_apply_precision(out$noise, out$A)
  out$G <- crossprod(out$A, out$precision_A)
  out$x_transform$mu <- fit$x_transform$mu[keep]
  if (!is.null(fit$x_transform$sd)) out$x_transform$sd <- fit$x_transform$sd[keep]
  out
}

.pattern_regions <- function(regions, p) {
  if (!is.list(regions) || !length(regions) || is.null(names(regions)) ||
      any(!nzchar(names(regions))) || anyDuplicated(names(regions))) {
    stop("regions must be a named list of input-column positions.", call. = FALSE)
  }
  lapply(regions, function(v) {
    if (is.logical(v)) {
      if (length(v) != p || anyNA(v)) stop("Logical regions must have one entry per input column.", call. = FALSE)
      v <- which(v)
    }
    if (!is.numeric(v) || any(!is.finite(v)) || any(v < 1 | v > p) ||
        any(v != as.integer(v)) || anyDuplicated(v) || !length(v)) {
      stop("Each region must contain unique valid input-column positions.", call. = FALSE)
    }
    as.integer(v)
  })
}

#' Evaluate regional access under a whole-brain pattern fit
#'
#' @param result A \code{pattern_global_result} retaining fold fits.
#' @param regions Named list of input-column positions or logical masks.
#'   These are dataset matrix columns, not voxel IDs. Regions may overlap.
#' @param independent_roi Optional named list of fold-resolved
#'   \code{pattern_ledger}s from independent ROI models on exactly the same
#'   rows and folds. Truth, target type, partition, and baseline must agree.
#' @return A table with region, metric, \code{whole_brain}, and
#'   \code{local_restricted} columns (and \code{independent_roi} if supplied).
#'   Regional pooled and fold-resolved ledgers are retained as attributes.
#' @details Each region uses the whole-brain task representation and the exact
#'   marginal covariance \code{Psi_RR}; neither patterns nor noise are refitted.
#'   All preprocessing is columnwise. Repeated assessments are averaged per
#'   observation before scoring, just as for the whole-brain ledger. R-squared
#'   uses fold training means. This measures locally accessible information
#'   under the shared model, not optimal performance of a separately fitted ROI.
#'   Finite-sample regional accuracy need not be below whole-brain accuracy.
#' @export
local_performance <- function(result, regions, independent_roi = NULL) {
  if (!inherits(result, "pattern_global_result") || is.null(result$fold_fits)) {
    stop("local_performance requires run_global(..., return_fits = TRUE).", call. = FALSE)
  }
  regions <- .pattern_regions(regions, result$n_features)
  ref <- result$fold_ledger
  dataset <- result$model_spec$dataset
  X <- if (ref$partition == "external") .pattern_test_matrix(dataset) else get_feature_matrix(dataset)
  if (!is.null(independent_roi) && (!is.list(independent_roi) ||
      !setequal(names(independent_roi), names(regions)))) {
    stop("independent_roi must supply one named ledger per region.", call. = FALSE)
  }
  whole <- .pattern_score_ledger(result$ledger)
  pooled <- folded <- vector("list", length(regions))
  names(pooled) <- names(folded) <- names(regions)
  tables <- lapply(names(regions), function(name) {
    ledger <- ref
    for (k in seq_along(result$fold_fits)) {
      rows <- which(ref$fold == k)
      if (is.null(result$fold_fits[[k]])) next
      fit <- .pattern_restrict_fit(result$fold_fits[[k]], regions[[name]])
      pred <- predict(fit, X[ref$observation[rows], , drop = FALSE],
                      type = if (ref$type == "categorical") "prob" else "decode")
      if (ref$type == "categorical") pred <- .pattern_pad_probs(pred, colnames(ref$prediction))
      ledger$prediction[rows, ] <- pred
    }
    pl <- .pattern_pool_ledger(ledger)
    folded[[name]] <<- ledger; pooled[[name]] <<- pl
    score <- .pattern_score_ledger(pl)
    tab <- data.frame(region = name, metric = names(score), whole_brain = unname(whole[names(score)]),
                      local_restricted = unname(score))
    if (!is.null(independent_roi)) {
      il <- independent_roi[[name]]
      fields <- c("fold", "observation", "truth", "partition", "type", "baseline")
      if (!inherits(il, "pattern_ledger") || !all(vapply(fields, function(f) isTRUE(all.equal(il[[f]], ref[[f]])), logical(1))) ||
          !identical(dim(il$prediction), dim(ref$prediction)) ||
          !identical(colnames(il$prediction), colnames(ref$prediction)) || any(!is.finite(il$prediction))) {
        stop("Independent ROI ledger must match assessment rows, folds, truth, targets, and baseline.", call. = FALSE)
      }
      if (ref$type == "categorical" && (any(il$prediction < 0 | il$prediction > 1) ||
          any(abs(rowSums(il$prediction) - 1) > 1e-8))) {
        stop("Independent ROI probabilities must be nonnegative and sum to one.", call. = FALSE)
      }
      tab$independent_roi <- unname(.pattern_score_ledger(.pattern_pool_ledger(il))[names(score)])
    }
    tab
  })
  out <- do.call(rbind, tables)
  rownames(out) <- NULL
  attr(out, "ledgers") <- pooled
  attr(out, "fold_ledgers") <- folded
  out
}

.pattern_subspace <- function(A, tol = 1e-10) {
  if (!nrow(A)) return(matrix(numeric(0), 0L, 0L))
  s <- svd(A, nv = 0)
  keep <- s$d > tol * max(s$d[1], .Machine$double.eps)
  s$u[, keep, drop = FALSE]
}

#' Compare spatial subspaces across retained folds
#'
#' @param object A global result retaining at least two fits.
#' @return One row per fold pair, with effective ranks, number of common
#'   retained features, principal angles in radians (a list column), and
#'   \code{overlap = sum(cos(angles)^2) / max(rank1, rank2)}. The normalization
#'   penalizes unequal ranks. Empty subspaces have undefined (NA) overlap.
#' @details Comparison uses original feature units on the intersection of
#'   columns retained by both folds. It is invariant to orthogonal component
#'   rotations. It measures spatial subspaces, not matched components or
#'   statistical evidence for a supported rank.
#' @export
component_stability <- function(object) {
  fits <- object$fold_fits
  idx <- if (is.null(fits)) integer(0) else which(!vapply(fits, is.null, logical(1)))
  if (!inherits(object, "pattern_global_result") || length(idx) < 2L) {
    stop("component_stability requires at least two retained fold fits.", call. = FALSE)
  }
  pairs <- utils::combn(idx, 2L)
  rows <- lapply(seq_len(ncol(pairs)), function(k) {
    i <- pairs[1, k]; j <- pairs[2, k]
    a <- fits[[i]]; b <- fits[[j]]
    common <- intersect(a$feature_index, b$feature_index)
    Qa <- .pattern_subspace(.pattern_map_values(a, "forward")[match(common, a$feature_index), , drop = FALSE])
    Qb <- .pattern_subspace(.pattern_map_values(b, "forward")[match(common, b$feature_index), , drop = FALSE])
    ra <- ncol(Qa); rb <- ncol(Qb)
    cs <- if (ra && rb) pmin(1, pmax(0, svd(crossprod(Qa, Qb), nu = 0, nv = 0)$d)) else numeric(0)
    data.frame(fold1 = i, fold2 = j, rank1 = ra, rank2 = rb,
               n_common = length(common), overlap = if (ra && rb) sum(cs^2) / max(ra, rb) else NA_real_,
               angles = I(list(acos(cs))))
  })
  do.call(rbind, rows)
}

#' Save a global pattern result and its invariant maps
#'
#' @inheritParams save_results
#' @return Invisibly, a named list of written paths. The RDS contains the
#'   complete result, including ledgers and any retained fits. Scalar refit maps
#'   are also written using the dataset's image format. Without a refit only
#'   the result RDS is written. Multibasis results retain vectors in RDS only.
#' @export
save_results.pattern_global_result <- function(x, dir,
    level = c("standard", "minimal", "complete"), stack = c("none", "auto", "vec"),
    fname = "pattern.nii.gz", include = NULL, dtype = NULL,
    overwrite = FALSE, quiet = FALSE) {
  level <- match.arg(level); stack <- match.arg(stack)
  .ensure_dir(dir)
  paths <- list(root = dir)
  paths$result <- .unique_path(file.path(dir, "pattern_result.rds"), overwrite)
  saveRDS(x, paths$result)
  if (!is.null(x$refit) && !inherits(x$model_spec$dataset, "mvpa_multibasis_image_dataset")) {
    maps <- lapply(c("signal_sd", "conditional_info"), function(type) model_patterns(x, type))
    names(maps) <- c("signal_sd", "conditional_info")
    pseudo <- structure(list(results = maps, metrics = names(maps)), class = c("searchlight_result", "list"))
    paths$maps <- save_results.searchlight_result(pseudo, dir, level = "minimal", stack = stack,
                  fname = fname, dtype = dtype, overwrite = overwrite, quiet = quiet)
  }
  invisible(paths)
}
