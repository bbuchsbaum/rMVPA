#' Permutation-Based Inference for Searchlight MVPA
#'
#' This module provides permutation testing infrastructure for searchlight MVPA,
#' including covariate-conditioned null distributions, stratified subsampling,
#' and FDR correction.  It wraps existing \pkg{rMVPA} infrastructure without
#' modifying any core files.
#'
#' @name permutation_searchlight
#' @keywords internal
#' @importFrom stats quantile
NULL


# ---------------------------------------------------------------------------
# 1. permutation_control() — Configuration object
# ---------------------------------------------------------------------------

#' Create a Permutation Control Object
#'
#' Specifies all tuning parameters for
#' \code{\link{run_permutation_searchlight}}.
#'
#' @param n_perm Integer >= 1. Number of permutations to run.
#' @param shuffle Character. How to permute labels: \code{"within_block"}
#'   (default) shuffles within each block, \code{"circular_shift"} shifts the
#'   label sequence within each block, \code{"global"} shuffles all labels
#'   ignoring block structure.
#' @param null_method Character. How to build the null distribution:
#'   \code{"adjusted"} conditions on covariate bins (default),
#'   \code{"global"} uses one global null.
#' @param adjust_by Character. Which covariate(s) to condition on:
#'   \code{"nfeatures"} (default), \code{"redundancy"}, or \code{"both"}.
#' @param n_bins Integer >= 2. Number of quantile bins for covariate
#'   stratification.
#' @param subsample Fraction (0, 1] or integer count of searchlight centers
#'   to use for permutation runs.  Only used when
#'   \code{perm_strategy = "iterate"} (see below).
#' @param stratify_subsample Logical. If \code{TRUE} (default), subsample
#'   centers proportionally from nfeatures quantile bins.
#' @param correction Character. Multiple-comparison correction: \code{"fdr"}
#'   (Benjamini-Hochberg, default) or \code{"none"}.
#' @param diagnose Logical. If \code{TRUE} (default), run null diagnostics.
#' @param seed Optional integer random seed.
#' @param perm_strategy Character. Controls how each permutation pass is
#'   executed.  Two strategies are available; neither contains any
#'   engine-specific branching:
#'
#'   \describe{
#'     \item{\code{"iterate"} (default)}{Each permutation runs
#'       \code{\link{mvpa_iterate}} on a \strong{subsampled} set of centers.
#'       This is the universal, safe path: it works with every model type
#'       and every searchlight engine because it goes through the generic
#'       per-ROI iterator.
#'
#'       \strong{When to use}: slow classifiers, large brains, limited
#'       compute.  The \code{subsample} parameter controls how many centers
#'       are evaluated per permutation, giving 5--20\eqn{\times} speedup
#'       over a full-brain pass.
#'
#'       Null pool size: \code{n_perm * n_subsampled_centers}.}
#'
#'     \item{\code{"searchlight"}}{Each permutation runs
#'       \code{\link{run_searchlight}} on the \strong{full brain}, then
#'       extracts metric values at every center.  Because the call goes
#'       through the standard \code{run_searchlight} dispatch, it
#'       automatically benefits from any fast engine the model qualifies
#'       for (e.g.\ SWIFT, dual-LDA) as well as any user-defined
#'       \code{run_searchlight.<class>} method.
#'
#'       \strong{When to use}: models with a fast searchlight engine, or
#'       when you want the richest possible null distribution.  Since the
#'       full brain is computed anyway, \emph{all} centers contribute to the
#'       null (the \code{subsample} parameter is ignored and a note is
#'       logged).
#'
#'       Null pool size: \code{n_perm * all_centers}.}
#'   }
#'
#' @return An S3 object of class \code{"permutation_control"}.
#'
#' @examples
#' # Default: subsampled iterator (safe for any model)
#' pc <- permutation_control(n_perm = 100, shuffle = "within_block",
#'                           subsample = 0.1, seed = 42L)
#'
#' # Full-brain strategy (benefits from fast engines, richer null)
#' pc2 <- permutation_control(n_perm = 20, perm_strategy = "searchlight",
#'                            seed = 42L)
#' print(pc2)
#'
#' @export
permutation_control <- function(
  n_perm              = 5L,
  shuffle             = c("within_block", "circular_shift", "global"),
  null_method         = c("adjusted", "global"),
  adjust_by           = c("nfeatures", "redundancy", "both"),
  n_bins              = 5L,
  subsample           = 0.1,
  stratify_subsample  = TRUE,
  correction          = c("fdr", "none"),
  diagnose            = TRUE,
  seed                = NULL,
  perm_strategy       = c("iterate", "searchlight")
) {
  shuffle       <- match.arg(shuffle)
  null_method   <- match.arg(null_method)
  adjust_by     <- match.arg(adjust_by)
  correction    <- match.arg(correction)
  perm_strategy <- match.arg(perm_strategy)

  assertthat::assert_that(is.numeric(n_perm), length(n_perm) == 1L,
                          n_perm >= 1, msg = "n_perm must be a positive integer >= 1")
  assertthat::assert_that(is.numeric(n_bins), length(n_bins) == 1L,
                          n_bins >= 2, msg = "n_bins must be an integer >= 2")
  assertthat::assert_that(is.numeric(subsample), length(subsample) == 1L,
                          subsample > 0, msg = "subsample must be > 0")
  if (!is.null(seed)) {
    assertthat::assert_that(is.numeric(seed) || is.integer(seed),
                            msg = "seed must be numeric or NULL")
  }

  structure(
    list(
      n_perm             = as.integer(n_perm),
      shuffle            = shuffle,
      null_method        = null_method,
      adjust_by          = adjust_by,
      n_bins             = as.integer(n_bins),
      subsample          = subsample,
      stratify_subsample = isTRUE(stratify_subsample),
      correction         = correction,
      diagnose           = isTRUE(diagnose),
      seed               = seed,
      perm_strategy      = perm_strategy
    ),
    class = "permutation_control"
  )
}

#' @export
print.permutation_control <- function(x, ...) {
  cat("Permutation Control Settings\n")
  cat("  perm_strategy     :", x$perm_strategy, "\n")
  cat("  n_perm            :", x$n_perm, "\n")
  cat("  shuffle           :", x$shuffle, "\n")
  cat("  null_method       :", x$null_method, "\n")
  cat("  adjust_by         :", x$adjust_by, "\n")
  cat("  n_bins            :", x$n_bins, "\n")
  if (identical(x$perm_strategy, "iterate")) {
    cat("  subsample         :", x$subsample, "\n")
    cat("  stratify_subsample:", x$stratify_subsample, "\n")
  } else {
    cat("  subsample         : (ignored - full brain per perm)\n")
  }
  cat("  correction        :", x$correction, "\n")
  cat("  diagnose          :", x$diagnose, "\n")
  cat("  seed              :", if (is.null(x$seed)) "NULL" else x$seed, "\n")
  invisible(x)
}


# ---------------------------------------------------------------------------
# 2. permute_labels() — Label shuffling
# ---------------------------------------------------------------------------

#' Permute Training Labels in an MVPA Design
#'
#' Returns a copy of \code{design} with \code{y_train}, \code{cv_labels},
#' and \code{targets} replaced by permuted versions.  \code{block_var} is
#' never permuted.
#'
#' @param design An \code{mvpa_design} object.
#' @param method Character. One of \code{"within_block"},
#'   \code{"circular_shift"}, or \code{"global"}.
#' @param seed Optional integer seed; RNG state is restored on exit.
#'
#' @return A modified \code{mvpa_design} with permuted labels.
#'
#' @export
permute_labels <- function(design,
                           method = c("within_block", "circular_shift", "global"),
                           seed   = NULL) {
  method <- match.arg(method)
  assertthat::assert_that(inherits(design, "mvpa_design"),
                          msg = "design must be an mvpa_design object")

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    set.seed(seed)
    on.exit({
      if (!is.null(old_seed)) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
  }

  labels <- design$cv_labels
  n      <- length(labels)

  perm_idx <- switch(method,
    global = sample.int(n),
    within_block = {
      block_var <- design$block_var
      if (is.null(block_var)) {
        sample.int(n)
      } else {
        idx <- seq_len(n)
        for (blk in unique(block_var)) {
          pos        <- which(block_var == blk)
          idx[pos]   <- pos[sample.int(length(pos))]
        }
        idx
      }
    },
    circular_shift = {
      block_var <- design$block_var
      if (is.null(block_var)) {
        shift  <- sample.int(max(1L, n - 1L), 1L)
        (seq_len(n) - 1L + shift) %% n + 1L
      } else {
        idx <- seq_len(n)
        for (blk in unique(block_var)) {
          pos <- which(block_var == blk)
          m   <- length(pos)
          if (m > 1L) {
            shift    <- sample.int(m - 1L, 1L)
            idx[pos] <- pos[(seq_len(m) - 1L + shift) %% m + 1L]
          }
        }
        idx
      }
    }
  )

  perm_design            <- design
  perm_design$cv_labels  <- labels[perm_idx]
  perm_design$y_train    <- perm_design$cv_labels
  perm_design$targets    <- if (!is.null(design$targets)) design$targets[perm_idx]
                            else perm_design$cv_labels
  perm_design
}


# ---------------------------------------------------------------------------
# 3. subsample_centers() — Stratified center selection
# ---------------------------------------------------------------------------

#' Subsample Searchlight Centers
#'
#' Selects a representative subset of searchlight centers for permutation
#' runs, optionally stratifying by the number of features per center.
#'
#' @param dataset An \code{mvpa_dataset}.
#' @param searchlight A searchlight iterator (list of integer voxel-index
#'   vectors) as returned by \code{get_searchlight()}.
#' @param n_centers Optional integer count of centers to select.  If
#'   \code{NULL}, derived from \code{fraction}.
#' @param fraction Numeric in (0, 1]. Fraction of all centers to select
#'   when \code{n_centers} is \code{NULL}.
#' @param stratify_by Character. Covariate used for stratification
#'   (\code{"nfeatures"} only for now).
#' @param redundancy_map Optional named numeric vector of per-center
#'   redundancy values (not used unless \code{stratify_by = "redundancy"}).
#' @param seed Optional integer seed.
#'
#' @return A list with:
#'   \describe{
#'     \item{center_ids}{Integer vector of selected center IDs.}
#'     \item{vox_list}{Subsetted searchlight list.}
#'     \item{covariates}{A \code{data.frame} with at least column
#'       \code{nfeatures}.}
#'   }
#'
#' @export
subsample_centers <- function(dataset,
                              searchlight,
                              n_centers      = NULL,
                              fraction       = 0.1,
                              stratify_by    = "nfeatures",
                              redundancy_map = NULL,
                              seed           = NULL) {
  all_ids   <- get_center_ids(dataset)
  nfeatures <- vapply(searchlight, length, integer(1L))
  total     <- length(all_ids)

  if (is.null(n_centers)) {
    n_centers <- max(1L, round(total * fraction))
  }
  n_centers <- as.integer(n_centers)

  if (n_centers >= total) {
    return(list(
      center_ids = all_ids,
      vox_list   = searchlight,
      covariates = data.frame(nfeatures = nfeatures)
    ))
  }

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    set.seed(seed)
    on.exit({
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }, add = TRUE)
  }

  # Stratified sampling by nfeatures quantile bins
  breaks   <- unique(quantile(nfeatures, probs = seq(0, 1, length.out = 6L)))
  bins     <- findInterval(nfeatures, breaks, rightmost.closed = TRUE)
  bins     <- pmax(1L, bins)

  tab        <- table(bins)
  bin_labels <- as.integer(names(tab))
  selected   <- integer(0)

  for (b in bin_labels) {
    in_bin <- which(bins == b)
    prop   <- length(in_bin) / total
    k      <- max(1L, round(n_centers * prop))
    k      <- min(k, length(in_bin))
    selected <- c(selected, in_bin[sample.int(length(in_bin), k)])
  }

  # Adjust to exactly n_centers
  if (length(selected) > n_centers) {
    selected <- selected[seq_len(n_centers)]
  } else if (length(selected) < n_centers) {
    remaining <- setdiff(seq_len(total), selected)
    extra     <- min(n_centers - length(selected), length(remaining))
    if (extra > 0L) {
      selected <- c(selected, remaining[sample.int(length(remaining), extra)])
    }
  }

  selected <- sort(selected)

  list(
    center_ids = all_ids[selected],
    vox_list   = searchlight[selected],
    covariates = data.frame(nfeatures = nfeatures[selected])
  )
}


# ---------------------------------------------------------------------------
# 4. compute_local_redundancy() — Autocorrelation proxy
# ---------------------------------------------------------------------------

#' Compute Local Redundancy for Each Searchlight Center
#'
#' For each searchlight sphere, computes the mean absolute correlation
#' between the center voxel (first column of the extracted data matrix) and
#' all other voxels in the sphere.  This serves as an autocorrelation proxy
#' that can be used as a covariate when building the null distribution.
#'
#' @param dataset An \code{mvpa_dataset}.
#' @param searchlight A searchlight iterator (list of integer voxel-index
#'   vectors).
#'
#' @return Named numeric vector of redundancy values, one per center.
#'
#' @export
compute_local_redundancy <- function(dataset, searchlight) {
  all_ids <- get_center_ids(dataset)
  n       <- length(searchlight)
  out     <- numeric(n)

  for (i in seq_len(n)) {
    tryCatch({
      samp  <- get_samples(dataset, searchlight[i])
      mat   <- as.matrix(neuroim2::values(samp$sample[[1]]))
      if (ncol(mat) < 2L) {
        out[i] <- 0
      } else {
        center_col <- mat[, 1L, drop = TRUE]
        cors <- apply(mat[, -1L, drop = FALSE], 2L, function(v) {
          suppressWarnings(abs(stats::cor(center_col, v, use = "complete.obs")))
        })
        out[i] <- mean(cors, na.rm = TRUE)
      }
    }, error = function(e) {
      out[i] <<- NA_real_
    })
  }

  names(out) <- as.character(all_ids[seq_len(n)])
  out
}


# ---------------------------------------------------------------------------
# 5. build_adjusted_null() — Quantile-binned null construction
# ---------------------------------------------------------------------------

#' Build a Covariate-Adjusted Null Distribution
#'
#' Bins null metric values by a covariate (e.g., number of features per
#' center) so that per-voxel p-values can be computed relative to a
#' matched null.
#'
#' @param null_values Numeric vector of null metric values from permutation
#'   runs.
#' @param covariates A \code{data.frame} with an \code{nfeatures} column
#'   (same length as \code{null_values}).
#' @param n_bins Integer >= 2. Number of quantile bins.
#' @param method Character. \code{"adjusted"} (default) uses covariate
#'   bins; \code{"global"} places all values in one bin.
#'
#' @return An S3 object of class \code{"adjusted_null"}.
#'
#' @keywords internal
build_adjusted_null <- function(null_values,
                                covariates,
                                n_bins = 5L,
                                method = c("adjusted", "global")) {
  method <- match.arg(method)

  if (method == "global" || n_bins <= 1L) {
    bin_assignment <- rep(1L, length(null_values))
    breaks         <- c(-Inf, Inf)
    n_bins_actual  <- 1L
  } else {
    probs  <- seq(0, 1, length.out = n_bins + 1L)
    breaks <- unique(quantile(covariates$nfeatures, probs = probs, na.rm = TRUE))
    bin_assignment <- findInterval(covariates$nfeatures, breaks,
                                   rightmost.closed = TRUE)
    bin_assignment <- pmax(1L, pmin(bin_assignment, length(breaks) - 1L))
    n_bins_actual  <- length(breaks) - 1L
  }

  bin_nulls <- vector("list", n_bins_actual)
  bin_stats <- vector("list", n_bins_actual)

  for (b in seq_len(n_bins_actual)) {
    vals <- null_values[bin_assignment == b]
    bin_nulls[[b]] <- sort(vals)
    bin_stats[[b]] <- list(
      n    = length(vals),
      mean = mean(vals, na.rm = TRUE),
      sd   = stats::sd(vals, na.rm = TRUE)
    )
  }

  structure(
    list(
      bin_nulls      = bin_nulls,
      bin_stats      = bin_stats,
      breaks         = breaks,
      n_bins         = n_bins_actual,
      method         = method,
      bin_assignment = bin_assignment
    ),
    class = "adjusted_null"
  )
}

#' @export
print.adjusted_null <- function(x, ...) {
  cat("Adjusted Null Distribution\n")
  cat("  method :", x$method, "\n")
  cat("  n_bins :", x$n_bins, "\n")
  for (b in seq_len(x$n_bins)) {
    st <- x$bin_stats[[b]]
    cat(sprintf("  bin %d: n=%d, mean=%.4f, sd=%.4f\n",
                b, st$n, st$mean, st$sd))
  }
  invisible(x)
}


# ---------------------------------------------------------------------------
# 6. score_observed() — P-value computation
# ---------------------------------------------------------------------------

#' Compute Permutation P-Values
#'
#' For each observed metric value, assigns it to the matching null bin and
#' computes a conservative one-sided p-value.
#'
#' @param observed_values Numeric vector of observed metric values (one per
#'   searchlight center across the full brain).
#' @param adjusted_null An \code{"adjusted_null"} object from
#'   \code{\link{build_adjusted_null}}.
#' @param covariates_full A \code{data.frame} with an \code{nfeatures}
#'   column for all voxels (same length as \code{observed_values}).
#'
#' @return Numeric vector of p-values in \eqn{[0, 1]}.
#'
#' @keywords internal
score_observed <- function(observed_values, adjusted_null, covariates_full) {
  n   <- length(observed_values)
  p   <- numeric(n)
  adj <- adjusted_null

  if (adj$method == "global" || adj$n_bins == 1L) {
    null_bin <- adj$bin_nulls[[1L]]
    for (i in seq_len(n)) {
      obs_i <- observed_values[i]
      if (is.na(obs_i)) {
        p[i] <- NA_real_
      } else {
        p[i] <- (sum(null_bin >= obs_i, na.rm = TRUE) + 1L) /
                (length(null_bin) + 1L)
      }
    }
  } else {
    bin_assignment <- findInterval(
      covariates_full$nfeatures,
      adj$breaks,
      rightmost.closed = TRUE
    )
    bin_assignment <- pmax(1L, pmin(bin_assignment, adj$n_bins))

    for (i in seq_len(n)) {
      obs_i <- observed_values[i]
      if (is.na(obs_i)) {
        p[i] <- NA_real_
        next
      }
      b        <- bin_assignment[i]
      null_bin <- adj$bin_nulls[[b]]
      if (length(null_bin) == 0L) {
        p[i] <- 1
      } else {
        p[i] <- (sum(null_bin >= obs_i, na.rm = TRUE) + 1L) /
                (length(null_bin) + 1L)
      }
    }
  }

  p
}


# ---------------------------------------------------------------------------
# 7. diagnose_null() — Stationarity diagnostics
# ---------------------------------------------------------------------------

#' Diagnose Null Distribution Stationarity
#'
#' Checks whether the null metric values are systematically associated with
#' covariates, which would indicate that covariate adjustment is needed.
#'
#' @param null_values Numeric vector of null metric values.
#' @param covariates A \code{data.frame} with at least an \code{nfeatures}
#'   column.
#' @param n_perm Integer. Number of permutations used (informational).
#'
#' @return An S3 object of class \code{"null_diagnostics"}.
#'
#' @keywords internal
diagnose_null <- function(null_values, covariates, n_perm) {
  checks <- list()

  # Check 1: correlation with nfeatures
  if (!is.null(covariates$nfeatures)) {
    ct <- tryCatch(
      stats::cor.test(null_values, covariates$nfeatures, method = "spearman",
                      exact = FALSE),
      error = function(e) NULL
    )
    if (!is.null(ct)) {
      flag <- !is.na(ct$p.value) && ct$p.value < 0.01
      checks[["nfeatures_cor"]] <- list(
        covariate = "nfeatures",
        rho       = unname(ct$estimate),
        p_value   = ct$p.value,
        flagged   = flag,
        message   = if (flag) {
          "Null correlates with nfeatures (p < 0.01); covariate adjustment recommended."
        } else {
          "No significant correlation with nfeatures."
        }
      )
    }
  }

  # Check 2: correlation with redundancy (if present)
  if (!is.null(covariates$redundancy)) {
    ct2 <- tryCatch(
      stats::cor.test(null_values, covariates$redundancy, method = "spearman",
                      exact = FALSE),
      error = function(e) NULL
    )
    if (!is.null(ct2)) {
      flag2 <- !is.na(ct2$p.value) && ct2$p.value < 0.01
      checks[["redundancy_cor"]] <- list(
        covariate = "redundancy",
        rho       = unname(ct2$estimate),
        p_value   = ct2$p.value,
        flagged   = flag2,
        message   = if (flag2) {
          "Null correlates with redundancy (p < 0.01); consider adjusting."
        } else {
          "No significant correlation with redundancy."
        }
      )
    }
  }

  structure(
    list(checks = checks, n_perm = n_perm),
    class = "null_diagnostics"
  )
}

#' @export
print.null_diagnostics <- function(x, ...) {
  cat("Null Distribution Diagnostics (n_perm =", x$n_perm, ")\n")
  cat(strrep("-", 50), "\n")
  for (nm in names(x$checks)) {
    chk <- x$checks[[nm]]
    flag_str <- if (isTRUE(chk$flagged)) "[FLAGGED]" else "[OK]"
    cat(sprintf("  %-20s %s\n", chk$covariate, flag_str))
    cat(sprintf("    rho=%.3f  p=%.4f\n", chk$rho, chk$p_value))
    cat(sprintf("    %s\n", chk$message))
  }
  cat(strrep("-", 50), "\n")
  invisible(x)
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Extract performance values from an mvpa_iterate result tibble
#' @keywords internal
#' @noRd
.extract_perf_values <- function(result_tbl, metric = NULL) {
  if (is.null(result_tbl) || nrow(result_tbl) == 0L) {
    return(numeric(0))
  }

  good_rows <- !isTRUE(result_tbl$error)
  if (is.logical(result_tbl$error)) {
    good_rows <- !result_tbl$error
  }

  perf_list <- result_tbl$performance[good_rows]
  ids_good  <- result_tbl$id[good_rows]

  vals <- vapply(perf_list, function(p) {
    if (is.null(p) || length(p) == 0L) return(NA_real_)
    if (!is.null(metric) && metric %in% names(p)) {
      return(as.numeric(p[[metric]]))
    }
    as.numeric(p[[1L]])
  }, numeric(1L))

  names(vals) <- as.character(ids_good)
  vals
}

#' Extract metric values from a searchlight_result at specified center IDs
#'
#' Converts a \code{searchlight_result} (as returned by
#' \code{\link{run_searchlight}}) into a named numeric vector of metric
#' values at the requested center IDs.
#'
#' This is the bridge that lets the \code{"searchlight"} permutation
#' strategy consume results from \emph{any} engine (SWIFT, dual-LDA,
#' legacy, or user-defined \code{run_searchlight.<class>} methods)
#' without engine-specific branching.  It supports:
#'
#' \itemize{
#'   \item A standard \code{searchlight_result} with a \code{$results}
#'     list of maps.
#'   \item Custom list-like outputs exposing a \code{$results} element.
#'   \item A direct numeric vector (preferably named by center IDs).
#' }
#'
#' @param sl_result A \code{searchlight_result} object.
#' @param center_ids Integer vector of center IDs to extract.
#' @param metric Character name of the metric to extract.  If \code{NULL}
#'   or not found in \code{sl_result$results}, the first metric is used.
#'
#' @return A named numeric vector (names = center IDs as character).
#'   Centers where extraction fails are \code{NA}.
#'
#' @keywords internal
#' @noRd
.extract_values_from_searchlight_result <- function(sl_result, center_ids,
                                                    metric = NULL) {
  out <- rep(NA_real_, length(center_ids))
  names(out) <- as.character(center_ids)

  if (is.null(sl_result)) {
    return(out)
  }

  extract_by_id <- function(x) {
    if (is.null(x)) return(NULL)
    tryCatch({
      if (!is.null(names(x))) {
        x[as.character(center_ids)]
      } else if (length(x) == length(center_ids)) {
        x
      } else {
        x[center_ids]
      }
    }, error = function(e) NULL)
  }

  if (is.numeric(sl_result)) {
    raw <- extract_by_id(sl_result)
    if (is.null(raw)) return(out)
    out <- as.numeric(raw)
    names(out) <- as.character(center_ids)
    return(out)
  }

  results_obj <- NULL
  if (is.list(sl_result) && !is.null(sl_result$results)) {
    results_obj <- sl_result$results
  } else if (inherits(sl_result, "searchlight_result")) {
    results_obj <- sl_result$results
  } else if (is.list(sl_result)) {
    # Fallback for custom list-like outputs that directly expose metric maps.
    results_obj <- sl_result
  }

  if (is.null(results_obj) || length(results_obj) == 0L) {
    return(out)
  }

  metric_map <- if (!is.null(metric) && metric %in% names(results_obj)) {
    results_obj[[metric]]
  } else {
    results_obj[[1L]]
  }

  raw <- extract_by_id(metric_map)
  if (is.null(raw)) return(out)

  out <- as.numeric(raw)
  names(out) <- as.character(center_ids)
  out
}

#' Partition variadic arguments for iterate vs searchlight calls
#' @keywords internal
#' @noRd
.partition_permutation_dots <- function(dots) {
  if (length(dots) == 0L) {
    return(list(iterate = list(), ignored_iterate = character(0)))
  }

  dot_names <- names(dots)
  if (is.null(dot_names)) {
    dot_names <- rep("", length(dots))
  }

  named <- nzchar(dot_names)
  iter_formals <- names(formals(mvpa_iterate))

  iterate <- dots[named & dot_names %in% iter_formals]
  ignored <- unique(dot_names[named & !(dot_names %in% iter_formals)])

  list(
    iterate = iterate,
    ignored_iterate = ignored
  )
}

#' Align performance values to a full center ID vector (NAs for missing)
#' @keywords internal
#' @noRd
.align_to_ids <- function(vals, all_ids) {
  if (length(vals) == length(all_ids) && is.null(names(vals))) {
    out <- as.numeric(vals)
    names(out) <- as.character(all_ids)
    return(out)
  }

  out          <- rep(NA_real_, length(all_ids))
  names(out)   <- as.character(all_ids)
  matched      <- match(names(vals), as.character(all_ids))
  valid        <- !is.na(matched)
  out[matched[valid]] <- vals[valid]
  out
}


# ---------------------------------------------------------------------------
# 8. run_permutation_searchlight() — Main orchestrator
# ---------------------------------------------------------------------------

#' Run Permutation Searchlight Inference
#'
#' Computes permutation-based p-values for a searchlight MVPA result by
#' running the analysis on permuted labels, building a covariate-adjusted
#' null distribution, and mapping p-values back to all brain voxels.
#'
#' @section Permutation strategy:
#'
#' The \code{perm_strategy} field in \code{perm_ctrl} determines how each
#' permutation pass is executed.  The two strategies share the same
#' downstream pipeline (null construction, p-value scoring, FDR
#' correction) — they differ only in \emph{how} null metric values are
#' produced.  Neither strategy contains any engine-specific branching.
#'
#' \describe{
#'   \item{\strong{\code{"iterate"}} (default)}{
#'     Calls \code{\link{mvpa_iterate}} on a \strong{subsampled} set of
#'     centers.  This is the generic per-ROI iterator that works with
#'     every model type and every dataset class.  The \code{subsample}
#'     parameter in \code{perm_ctrl} controls how many centers are
#'     evaluated per permutation.
#'
#'     \strong{Null pool size}: \code{n_perm * n_subsampled_centers}.
#'
#'     \strong{Best for}: slow classifiers, large brains, limited compute.
#'     The subsampling gives a 5--20\eqn{\times} speedup over a full-brain
#'     pass.
#'   }
#'
#'   \item{\strong{\code{"searchlight"}}}{
#'     Calls \code{\link{run_searchlight}} on the \strong{full brain} for
#'     each permutation.  Because the call goes through the standard
#'     \code{run_searchlight} S3 dispatch, it automatically benefits from
#'     any fast engine the model qualifies for (e.g.\ SWIFT, dual-LDA)
#'     \emph{and} from any user-defined \code{run_searchlight.<class>}
#'     method.  No engine-specific code exists here — it is purely the
#'     standard \code{run_searchlight} call.
#'
#'     Since the full brain is computed anyway, \strong{all} centers
#'     contribute to the null distribution (the \code{subsample} parameter
#'     is ignored and a note is logged).  This yields a richer null and
#'     therefore better-calibrated p-values.
#'
#'     \strong{Null pool size}: \code{n_perm * all_centers}.
#'
#'     \strong{Best for}: models with a fast searchlight engine, or when
#'     you want the richest possible null distribution.
#'   }
#' }
#'
#' @param model_spec An \code{mvpa_model} or compatible model specification.
#' @param observed Optional pre-computed searchlight result (output of
#'   \code{run_searchlight()}).  If \code{NULL}, it is computed internally.
#'   Can also be a named numeric vector of observed metric values indexed by
#'   center ID.
#' @param radius Numeric. Searchlight radius in mm (default 8).
#' @param method Character. Searchlight method passed to
#'   \code{run_searchlight()} when \code{observed} is \code{NULL} or when
#'   \code{perm_strategy = "searchlight"}.
#' @param perm_ctrl A \code{\link{permutation_control}} object.
#' @param metric Character. Which performance metric to use for inference.
#'   If \code{NULL} (default), the first metric is used.
#' @param ... Additional arguments forwarded to \code{run_searchlight()}
#'   (observed pass and \code{"searchlight"} permutations).  For
#'   \code{"iterate"} permutations, only arguments that are formal
#'   parameters of \code{\link{mvpa_iterate}} are forwarded; other keys
#'   are ignored for that path to avoid argument-mismatch failures
#'   (e.g., \code{engine = "legacy"} is meaningful for
#'   \code{run_searchlight} but not for \code{mvpa_iterate}).
#'
#' @return A \code{permutation_result} S3 object containing:
#'   \describe{
#'     \item{p_map}{Spatial map of raw p-values.}
#'     \item{p_adj_map}{Spatial map of FDR-adjusted p-values (if requested).}
#'     \item{p_values}{Numeric vector of raw p-values (all centers).}
#'     \item{p_adjusted}{Numeric vector of adjusted p-values.}
#'     \item{observed}{The observed searchlight result.}
#'     \item{diagnostics}{A \code{"null_diagnostics"} object (if
#'       \code{diagnose = TRUE}).}
#'     \item{perm_ctrl}{The \code{permutation_control} used.}
#'     \item{metric}{Metric name used for inference.}
#'     \item{perm_strategy}{The strategy that was actually used.}
#'   }
#'
#' @examples
#' \donttest{
#'   ds    <- gen_sample_dataset(c(5, 5, 5), 20, blocks = 2, nlevels = 2)
#'   cval  <- blocked_cross_validation(ds$design$block_var)
#'   mdl   <- load_model("sda_notune")
#'   mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification",
#'                       crossval = cval)
#'
#'   # Strategy 1: subsampled iterator (default, universal)
#'   pc1   <- permutation_control(n_perm = 10, subsample = 0.2, seed = 1L)
#'   res1  <- run_permutation_searchlight(mspec, radius = 3, perm_ctrl = pc1)
#'
#'   # Strategy 2: full-brain via run_searchlight (engine-aware)
#'   pc2   <- permutation_control(n_perm = 5, perm_strategy = "searchlight",
#'                                seed = 1L)
#'   res2  <- run_permutation_searchlight(mspec, radius = 3, perm_ctrl = pc2)
#' }
#'
#' @export
run_permutation_searchlight <- function(
  model_spec,
  observed   = NULL,
  radius     = 8,
  method     = c("standard", "randomized"),
  perm_ctrl  = permutation_control(),
  metric     = NULL,
  ...
) {
  method <- match.arg(method)
  strategy <- perm_ctrl$perm_strategy %||% "iterate"
  dots <- list(...)
  dot_parts <- .partition_permutation_dots(dots)

  assertthat::assert_that(inherits(perm_ctrl, "permutation_control"),
                          msg = "perm_ctrl must be a permutation_control object")
  assertthat::assert_that(
    strategy %in% c("iterate", "searchlight"),
    msg = "perm_ctrl$perm_strategy must be one of: 'iterate', 'searchlight'"
  )

  # Step 1: run observed searchlight if not supplied
  if (is.null(observed)) {
    futile.logger::flog.info("Running observed searchlight (radius = %g mm) ...", radius)
    observed <- do.call(
      run_searchlight,
      c(list(model_spec = model_spec, radius = radius, method = method), dots)
    )
  }

  # Step 2: build searchlight iterator (needed for feature counts)
  futile.logger::flog.info("Building searchlight iterator ...")
  sl <- get_searchlight(model_spec$dataset, "standard", radius)

  all_ids       <- get_center_ids(model_spec$dataset)
  nfeatures_all <- vapply(sl, length, integer(1L))
  covariates_full <- data.frame(nfeatures = nfeatures_all)

  # Step 3: optionally compute redundancy
  redundancy_map <- NULL
  if (perm_ctrl$adjust_by %in% c("redundancy", "both")) {
    futile.logger::flog.info("Computing local redundancy (this may take a while) ...")
    redundancy_map <- compute_local_redundancy(model_spec$dataset, sl)
    covariates_full$redundancy <- redundancy_map[as.character(all_ids)]
  }

  # Step 4: subsample centers (only used by "iterate" strategy)
  #
  # For "searchlight" strategy, subsampling is skipped because the full
  # brain is computed per permutation anyway.  All centers contribute to
  # the null, yielding a richer distribution.
  sub <- NULL
  if (identical(strategy, "iterate")) {
    futile.logger::flog.info("Subsampling searchlight centers ...")
    frac  <- if (perm_ctrl$subsample <= 1) perm_ctrl$subsample else NULL
    n_cnt <- if (perm_ctrl$subsample > 1)  as.integer(perm_ctrl$subsample) else NULL

    sub <- subsample_centers(
      dataset        = model_spec$dataset,
      searchlight    = sl,
      n_centers      = n_cnt,
      fraction       = if (is.null(frac)) 0.1 else frac,
      stratify_by    = "nfeatures",
      redundancy_map = redundancy_map,
      seed           = perm_ctrl$seed
    )
    futile.logger::flog.info("Using %d / %d centers for permutation runs.",
                             length(sub$center_ids), length(all_ids))
  } else {
    futile.logger::flog.info(
      paste0("Strategy = 'searchlight': full brain computed per permutation. ",
             "'subsample' parameter ignored; all %d centers contribute to null."),
      length(all_ids))
  }

  # ------------------------------------------------------------------
  # Step 5: run permutations
  # ------------------------------------------------------------------
  #
  # The two strategies are kept in separate code paths for clarity, but
  # they share the same contract: after the loop, `null_values` and
  # `null_nfeatures` are numeric vectors of equal length, ready for
  # build_adjusted_null().
  #
  # Neither path contains any engine-specific logic.  The "searchlight"
  # path calls run_searchlight(), which handles engine dispatch
  # internally via S3 methods.  The "iterate" path calls mvpa_iterate(),
  # which is the universal per-ROI iterator.
  #
  # Dots policy:
  # - Searchlight path receives all user dots (for run_searchlight dispatch).
  # - Iterate path only receives dots that are formal mvpa_iterate args.
  #   Unsupported dots (e.g., engine = "legacy") are ignored for iterate
  #   to avoid random failures from argument mismatch.
  # ------------------------------------------------------------------

  if (length(dot_parts$ignored_iterate) > 0L) {
    futile.logger::flog.debug(
      "Ignoring iterate-unsupported args in ...: %s",
      paste(dot_parts$ignored_iterate, collapse = ", ")
    )
  }

  null_values    <- numeric(0)
  null_nfeatures <- numeric(0)

  for (i in seq_len(perm_ctrl$n_perm)) {
    perm_seed <- if (!is.null(perm_ctrl$seed)) perm_ctrl$seed + i else NULL

    futile.logger::flog.info("Permutation %d / %d (strategy: %s) ...",
                             i, perm_ctrl$n_perm, strategy)

    perm_design <- permute_labels(
      model_spec$design,
      method = perm_ctrl$shuffle,
      seed   = perm_seed
    )

    perm_spec        <- model_spec
    perm_spec$design <- perm_design

    if (identical(strategy, "iterate")) {
      # ---- "iterate" strategy ----
      # Run the generic per-ROI iterator on subsampled centers only.
      # Works with any model type; no engine dispatch.
      perm_result <- tryCatch(
        do.call(
          mvpa_iterate,
          c(
            list(
              mod_spec = perm_spec,
              vox_list = sub$vox_list,
              ids = sub$center_ids,
              verbose = FALSE
            ),
            dot_parts$iterate
          )
        ),
        error = function(e) {
          futile.logger::flog.warn("Permutation %d failed: %s", i, e$message)
          NULL
        }
      )

      if (is.null(perm_result) || nrow(perm_result) == 0L) next

      perm_vals <- .extract_perf_values(perm_result, metric = metric)
      if (length(perm_vals) == 0L) next

      null_values    <- c(null_values, perm_vals)
      null_nfeatures <- c(null_nfeatures,
                          sub$covariates$nfeatures[
                            match(names(perm_vals), as.character(sub$center_ids))
                          ])

    } else {
      # ---- "searchlight" strategy ----
      # Run the full run_searchlight() pipeline, which dispatches to
      # whatever engine the model qualifies for (SWIFT, dual-LDA,
      # legacy, or any custom run_searchlight.<class> method).
      # Then extract metric values at ALL centers from the spatial map.
      perm_sl <- tryCatch(
        do.call(
          run_searchlight,
          c(list(model_spec = perm_spec, radius = radius, method = method), dots)
        ),
        error = function(e) {
          futile.logger::flog.warn("Permutation %d (searchlight) failed: %s",
                                   i, e$message)
          NULL
        }
      )

      if (is.null(perm_sl)) next

      perm_vals <- .extract_values_from_searchlight_result(
        perm_sl, all_ids, metric = metric
      )

      # Drop NAs (centers that failed or produced no result)
      valid     <- !is.na(perm_vals)
      if (sum(valid) == 0L) next

      null_values    <- c(null_values, perm_vals[valid])
      null_nfeatures <- c(null_nfeatures, nfeatures_all[valid])
    }
  }

  if (length(null_values) == 0L) {
    stop("No valid null values collected. All permutations may have failed.")
  }

  null_covariates <- data.frame(nfeatures = null_nfeatures)

  # Determine metric name used
  if (is.null(metric)) {
    # Try to infer from observed
    if (inherits(observed, "searchlight_result") && length(observed$metrics) > 0L) {
      metric <- observed$metrics[[1L]]
    } else {
      metric <- "metric"
    }
  }

  # Step 6: optional diagnostics on the null
  diag_result <- NULL
  if (perm_ctrl$diagnose) {
    futile.logger::flog.info("Running null diagnostics ...")
    diag_result <- diagnose_null(null_values, null_covariates,
                                 n_perm = perm_ctrl$n_perm)
    print(diag_result)
  }

  # Step 7: build adjusted null
  futile.logger::flog.info("Building adjusted null distribution (%s, %d bins) ...",
                           perm_ctrl$null_method, perm_ctrl$n_bins)
  adj_null <- build_adjusted_null(
    null_values  = null_values,
    covariates   = null_covariates,
    n_bins       = perm_ctrl$n_bins,
    method       = perm_ctrl$null_method
  )

  # Step 8: extract observed values for all centers
  obs_vals <- if (is.numeric(observed)) {
    .align_to_ids(observed, all_ids)
  } else {
    .extract_values_from_searchlight_result(observed, all_ids, metric = metric)
  }

  # Step 9: compute p-values for all centers
  futile.logger::flog.info("Computing p-values for %d centers ...", length(all_ids))
  p_values <- score_observed(obs_vals, adj_null, covariates_full)

  # Step 10: FDR correction
  p_adjusted <- if (perm_ctrl$correction == "fdr") {
    stats::p.adjust(p_values, method = "BH")
  } else {
    p_values
  }

  # Step 11: build spatial maps
  futile.logger::flog.info("Building p-value spatial maps ...")
  p_map <- tryCatch(
    build_output_map(model_spec$dataset, p_values, all_ids),
    error = function(e) {
      futile.logger::flog.warn("build_output_map for p_map failed: %s", e$message)
      NULL
    }
  )

  p_adj_map <- if (perm_ctrl$correction != "none") {
    tryCatch(
      build_output_map(model_spec$dataset, p_adjusted, all_ids),
      error = function(e) {
        futile.logger::flog.warn("build_output_map for p_adj_map failed: %s", e$message)
        NULL
      }
    )
  } else {
    NULL
  }

  n_sig <- sum(p_adjusted < 0.05, na.rm = TRUE)
  futile.logger::flog.info("Done. %d centers significant at FDR < 0.05 (%s).",
                           n_sig, perm_ctrl$correction)

  structure(
    list(
      p_map        = p_map,
      p_adj_map    = p_adj_map,
      p_values     = p_values,
      p_adjusted   = p_adjusted,
      observed     = observed,
      diagnostics  = diag_result,
      adj_null     = adj_null,
      perm_ctrl    = perm_ctrl,
      metric       = metric,
      perm_strategy = strategy,
      all_ids      = all_ids,
      n_perm_used  = perm_ctrl$n_perm,
      n_null_vals  = length(null_values)
    ),
    class = "permutation_result"
  )
}


# ---------------------------------------------------------------------------
# 9. S3 print / summary methods
# ---------------------------------------------------------------------------

#' @export
print.permutation_result <- function(x, ...) {
  cat("Permutation Searchlight Result\n")
  cat(strrep("=", 40), "\n")
  cat("  Metric           :", x$metric, "\n")
  cat("  Strategy         :", x$perm_strategy %||%
        x$perm_ctrl$perm_strategy %||% "iterate", "\n")
  cat("  Permutations used:", x$n_perm_used, "\n")
  cat("  Null values      :", x$n_null_vals, "\n")
  cat("  Total centers    :", length(x$all_ids), "\n")
  cat("  Correction       :", x$perm_ctrl$correction, "\n")
  n_sig_raw <- sum(x$p_values < 0.05, na.rm = TRUE)
  n_sig_adj <- sum(x$p_adjusted < 0.05, na.rm = TRUE)
  cat(sprintf("  Significant (raw p < 0.05)    : %d\n", n_sig_raw))
  cat(sprintf("  Significant (adj p < 0.05)    : %d\n", n_sig_adj))
  invisible(x)
}

#' @export
summary.permutation_result <- function(object, ...) {
  cat("Permutation Searchlight Summary\n")
  cat(strrep("=", 50), "\n")
  cat("  Metric           :", object$metric, "\n")
  cat("  Strategy         :", object$perm_strategy %||%
        object$perm_ctrl$perm_strategy %||% "iterate", "\n")
  cat("  Permutations     :", object$n_perm_used, "\n")
  cat("  Null values      :", object$n_null_vals, "\n")
  cat("  Shuffle method   :", object$perm_ctrl$shuffle, "\n")
  cat("  Null method      :", object$perm_ctrl$null_method, "\n")
  cat("  Bins             :", object$perm_ctrl$n_bins, "\n")
  cat("  Correction       :", object$perm_ctrl$correction, "\n")

  p   <- object$p_values
  adj <- object$p_adjusted
  cat("\n  P-value distribution (raw):\n")
  cat("   ", paste(names(quantile(p, na.rm = TRUE)),
                   round(quantile(p, na.rm = TRUE), 4), sep = "=",
                   collapse = "  "), "\n")
  cat("  P-value distribution (adjusted):\n")
  cat("   ", paste(names(quantile(adj, na.rm = TRUE)),
                   round(quantile(adj, na.rm = TRUE), 4), sep = "=",
                   collapse = "  "), "\n")

  thresholds <- c(0.05, 0.01, 0.001)
  cat("\n  Significant centers:\n")
  for (thr in thresholds) {
    cat(sprintf("    adj p < %.3f : %d\n", thr,
                sum(adj < thr, na.rm = TRUE)))
  }

  if (!is.null(object$diagnostics)) {
    cat("\n")
    print(object$diagnostics)
  }

  invisible(object)
}
