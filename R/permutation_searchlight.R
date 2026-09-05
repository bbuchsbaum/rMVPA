#' Permutation-Based Inference for Searchlight MVPA
#'
#' This module provides permutation testing infrastructure for searchlight MVPA,
#' including covariate-conditioned null distributions, stratified subsampling,
#' and FDR correction.  It wraps existing \pkg{rMVPA} infrastructure; a model
#' type takes part by providing a \code{\link{permute_labels}} method for its
#' design class.
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
#'   ignoring block structure. Circular shifts include zero in each block;
#'   leaving some or all blocks unchanged is a valid null draw.
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
#' @param rsa_null Null hypothesis for \code{rsa_model}: \code{"individual"}
#'   (default) tests marginal associations for correlation models, or the
#'   coefficient of a regression with only one design predictor. Regression
#'   models with multiple predictors (including nuisance predictors) require
#'   \code{"joint"}: no association with \emph{any} design predictor. Raw item
#'   permutations cannot test an individual regression coefficient conditional
#'   on the others. Ignored for other model classes.
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
  perm_strategy       = c("iterate", "searchlight"),
  rsa_null            = c("individual", "joint")
) {
  shuffle       <- match.arg(shuffle)
  null_method   <- match.arg(null_method)
  adjust_by     <- match.arg(adjust_by)
  correction    <- match.arg(correction)
  perm_strategy <- match.arg(perm_strategy)
  rsa_null      <- match.arg(rsa_null)

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
      perm_strategy      = perm_strategy,
      rsa_null           = rsa_null
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
  cat("  rsa_null          :", x$rsa_null %||% "individual", "\n")
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

#' Permute Labels in a Design for Permutation Testing
#'
#' Returns a copy of \code{design} in which the link between the brain data
#' and the quantity being decoded or modelled has been broken by a random
#' relabelling that respects the design's block structure.  This is the
#' null-generating step of \code{\link{run_permutation_searchlight}}.
#'
#' The generic dispatches on the design class:
#' \describe{
#'   \item{\code{mvpa_design}}{Shuffles trial labels: \code{y_train},
#'     \code{cv_labels}, and \code{targets} are replaced by permuted versions.
#'     \code{block_var} is never permuted.}
#'   \item{\code{rsa_design}}{Permutes \emph{item} labels, not RDM entries.
#'     The entries of an RDM are not independent observations: every item
#'     takes part in \code{n - 1} pairs, so shuffling the vectorised RDM would
#'     destroy that dependence and give an anti-conservative null.
#'     Relabelling items keeps it intact.  The permutation is stored on the
#'     design as \code{item_perm}; \code{train_model.rsa_model} applies it to
#'     the rows of the neural pattern matrix before the brain RDM is computed.
#'     That is equivalent to permuting the rows and columns of every model RDM
#'     by the inverse permutation, so the model matrix and any cached kernel
#'     stay untouched.}
#'   \item{\code{pair_rsa_design}}{As \code{rsa_design} in within-domain mode.
#'     For \code{pairs = "between"} the two item sets are permuted
#'     independently, each within its own blocks (\code{block_var_a},
#'     \code{block_var_b}), and the second permutation is stored as
#'     \code{item_perm_b}.}
#' }
#'
#' RSA item permutations break the brain's association with \emph{all}
#' design predictors, including nuisance predictors. For regression with
#' multiple predictors they generate a joint no-association null, not a null
#' for an individual coefficient conditional on the other predictors. See
#' the \code{rsa_null} setting in \code{\link{permutation_control}}.
#'
#' Circular shifts sample every rotation, including zero, independently in
#' each block. Identity draws and draws that move only some blocks belong to
#' the null distribution and are retained.
#'
#' When an RSA design excludes within-block pairs (the default whenever a
#' \code{block_var} is supplied and some pairs fall within a block), items are
#' exchangeable only within a block, and \code{method = "global"} is refused:
#' the null would compare neural patterns from the same run while the observed
#' statistic never does.  If every block holds a single item, no within-block
#' shuffle can move anything and the method errors; when no pairs are excluded
#' in that situation, \code{method = "global"} is the valid choice.
#'
#' @param design An \code{mvpa_design}, \code{rsa_design}, or
#'   \code{pair_rsa_design} object.
#' @param method Character. One of \code{"within_block"} (default),
#'   \code{"circular_shift"}, or \code{"global"}.
#' @param seed Optional integer seed; RNG state is restored on exit.
#'
#' @return A modified copy of \code{design} with the same class.
#'
#' @export
permute_labels <- function(design,
                           method = c("within_block", "circular_shift", "global"),
                           seed   = NULL) {
  UseMethod("permute_labels")
}

#' @rdname permute_labels
#' @export
permute_labels.default <- function(design,
                                   method = c("within_block", "circular_shift", "global"),
                                   seed   = NULL) {
  stop(sprintf(
    "permute_labels() has no method for designs of class <%s>. Supported: mvpa_design, rsa_design, pair_rsa_design.",
    paste(class(design), collapse = "/")
  ), call. = FALSE)
}

#' @rdname permute_labels
#' @export
permute_labels.mvpa_design <- function(design,
                                       method = c("within_block", "circular_shift", "global"),
                                       seed   = NULL) {
  method <- match.arg(method)
  labels <- design$cv_labels

  perm_idx <- .with_perm_seed(seed, function() {
    .permutation_index(length(labels), design$block_var, method)
  })

  perm_design            <- design
  perm_design$cv_labels  <- labels[perm_idx]
  perm_design$y_train    <- perm_design$cv_labels
  perm_design$targets    <- if (!is.null(design$targets)) design$targets[perm_idx]
                            else perm_design$cv_labels
  perm_design
}

#' @rdname permute_labels
#' @export
permute_labels.rsa_design <- function(design,
                                      method = c("within_block", "circular_shift", "global"),
                                      seed   = NULL) {
  method <- match.arg(method)
  .check_rsa_shuffle(design, method)

  n_items <- .rsa_design_n_items(design)
  if (!.blocks_permutable(design$block_var, n_items, method)) {
    .stop_not_permutable(design)
  }
  design$item_perm <- .with_perm_seed(seed, function() {
    .permutation_index(n_items, design$block_var, method)
  })
  design
}

#' @rdname permute_labels
#' @export
permute_labels.pair_rsa_design <- function(design,
                                           method = c("within_block", "circular_shift", "global"),
                                           seed   = NULL) {
  method <- match.arg(method)
  .check_rsa_shuffle(design, method)

  between <- identical(design$pair_kind %||% "within", "between")
  if (between && !is.null(design$include) && is.null(design$block_var_b)) {
    stop(paste0(
      "This between-domain pair_rsa_design excludes within-block pairs but does not ",
      "record `block_var_b`; rebuild it with the current pair_rsa_design() before permuting."
    ), call. = FALSE)
  }

  permutable <- .blocks_permutable(design$block_var, design$n_a, method) ||
    (between && .blocks_permutable(design$block_var_b, design$n_b, method))
  if (!permutable) {
    .stop_not_permutable(design)
  }

  perms <- .with_perm_seed(seed, function() {
    a <- .permutation_index(design$n_a, design$block_var, method)
    b <- if (between) .permutation_index(design$n_b, design$block_var_b, method) else NULL
    list(a = a, b = b)
  })

  design$item_perm <- perms$a
  if (between) design$item_perm_b <- perms$b
  design
}

# Shared helpers -------------------------------------------------------------

#' @keywords internal
#' @noRd
.permutation_index <- function(n, block_var = NULL,
                               method = c("within_block", "circular_shift", "global")) {
  method <- match.arg(method)
  if (n < 1L) return(integer(0))

  switch(method,
    global = sample.int(n),
    within_block = {
      if (is.null(block_var)) {
        sample.int(n)
      } else {
        idx <- seq_len(n)
        for (blk in unique(block_var)) {
          pos      <- which(block_var == blk)
          idx[pos] <- pos[sample.int(length(pos))]
        }
        idx
      }
    },
    circular_shift = {
      if (is.null(block_var)) {
        shift <- sample.int(n, 1L) - 1L
        (seq_len(n) - 1L + shift) %% n + 1L
      } else {
        idx <- seq_len(n)
        for (blk in unique(block_var)) {
          pos <- which(block_var == blk)
          m   <- length(pos)
          if (m > 1L) {
            shift    <- sample.int(m, 1L) - 1L
            idx[pos] <- pos[(seq_len(m) - 1L + shift) %% m + 1L]
          }
        }
        idx
      }
    }
  )
}

#' @keywords internal
#' @noRd
.with_perm_seed <- function(seed, fn) {
  if (is.null(seed)) {
    return(fn())
  }
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
  fn()
}

#' @keywords internal
#' @noRd
.check_rsa_shuffle <- function(design, method) {
  excludes_pairs <- !is.null(design$include) && !all(design$include)
  if (identical(method, "global") && excludes_pairs) {
    stop(paste0(
      "shuffle = 'global' is not valid for an RSA design that excludes within-block pairs: ",
      "relabelling items across blocks would put same-block neural pairs into the null ",
      "that the observed statistic never uses. Use 'within_block' or 'circular_shift'."
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' Can a block-respecting shuffle move anything?
#'
#' \code{FALSE} when fewer than two items exist, or when every block holds a
#' single item so that \code{"within_block"} and \code{"circular_shift"} are
#' forced to the identity.  \code{"global"} ignores blocks and is permutable
#' whenever there are two or more items.
#'
#' @keywords internal
#' @noRd
.blocks_permutable <- function(block_var, n, method) {
  if (n < 2L) return(FALSE)
  if (identical(method, "global") || is.null(block_var)) return(TRUE)
  any(table(block_var) > 1L)
}

#' @keywords internal
#' @noRd
.stop_not_permutable <- function(design) {
  excludes_pairs <- !is.null(design$include) && !all(design$include)
  hint <- if (excludes_pairs) {
    "The design excludes within-block pairs, so a global shuffle is not valid either; an item-permutation null is not available for this design."
  } else {
    "No pairs are excluded, so items are exchangeable across blocks: use method = 'global'."
  }
  stop(paste0(
    "No within-block item permutation is possible: every block holds a single item ",
    "(or there are fewer than two items). ", hint
  ), call. = FALSE)
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
.extract_perf_values <- function(result_tbl, metric = NULL, strict = FALSE) {
  mat <- .extract_perf_matrix(result_tbl, metrics = metric, strict = strict)
  vals <- as.numeric(mat[, 1L])
  names(vals) <- rownames(mat)
  vals
}

#' Extract one value from a per-ROI performance entry
#'
#' Performance entries are named numeric vectors for most model types and
#' one-row matrices with column names for \code{rsa_model}.  Both are read by
#' name.  With \code{strict = FALSE} a metric that is absent falls back to the
#' first entry (the historical behaviour for an inferred metric); with
#' \code{strict = TRUE} it yields \code{NA}, so a user-named metric that does
#' not exist surfaces as an error downstream instead of silently scoring the
#' wrong quantity.
#'
#' @keywords internal
#' @noRd
.perf_entry_value <- function(p, metric = NULL, strict = FALSE) {
  if (is.null(p) || length(p) == 0L) return(NA_real_)
  nm <- if (is.matrix(p)) colnames(p) else names(p)
  if (!is.null(metric)) {
    if (!is.null(nm) && metric %in% nm) {
      v <- if (is.matrix(p)) p[1L, metric] else p[[metric]]
      return(as.numeric(v)[1L])
    }
    if (strict) return(NA_real_)
  }
  as.numeric(p[[1L]])[1L]
}

#' Extract several metrics from an mvpa_iterate result table
#'
#' @return A numeric matrix with one row per successful center (rownames are
#'   center ids) and one column per requested metric.  When \code{metrics} is
#'   \code{NULL} a single positional column named \code{"metric"} is returned.
#' @keywords internal
#' @noRd
.extract_perf_matrix <- function(result_tbl, metrics = NULL, strict = FALSE) {
  col_names <- if (is.null(metrics)) "metric" else metrics
  keys      <- if (is.null(metrics)) list(NULL) else as.list(metrics)
  k         <- length(col_names)

  empty <- matrix(numeric(0), nrow = 0L, ncol = k, dimnames = list(NULL, col_names))
  if (is.null(result_tbl) || nrow(result_tbl) == 0L) {
    return(empty)
  }

  good_rows <- !isTRUE(result_tbl$error)
  if (is.logical(result_tbl$error)) {
    good_rows <- !result_tbl$error
  }

  perf_list <- result_tbl$performance[good_rows]
  ids_good  <- result_tbl$id[good_rows]
  if (length(perf_list) == 0L) {
    return(empty)
  }

  vals <- vapply(perf_list, function(p) {
    vapply(keys, function(m) .perf_entry_value(p, m, strict), numeric(1L))
  }, numeric(k))

  matrix(vals, nrow = length(perf_list), ncol = k, byrow = TRUE,
         dimnames = list(as.character(ids_good), col_names))
}

#' Extract several metrics from a searchlight_result at given center ids
#'
#' @return A numeric matrix, rows = \code{center_ids}, columns = metrics.
#' @keywords internal
#' @noRd
.searchlight_result_matrix <- function(sl_result, center_ids, metrics, strict = FALSE) {
  cols <- lapply(metrics, function(m) {
    .extract_values_from_searchlight_result(sl_result, center_ids, metric = m,
                                            strict = strict)
  })
  mat <- do.call(cbind, cols)
  dimnames(mat) <- list(as.character(center_ids), metrics)
  mat
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
                                                    metric = NULL,
                                                    strict = FALSE) {
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
  } else if (strict && !is.null(metric)) {
    return(out)
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
#' @section RSA null hypothesis:
#' Raw item permutations break associations with every design predictor,
#' including nuisance predictors. With \code{regtype = "lm"} or
#' \code{"rfit"} and multiple design predictors, they do not preserve the
#' effects of the other predictors when testing one coefficient. Such runs
#' are refused by default, even when \code{metric} selects a single output.
#' Use \code{permutation_control(rsa_null = "joint")} only to test the joint
#' no-association null. Each output metric then supplies a statistic for that
#' same joint null; significance does not establish that the corresponding
#' predictor has a unique effect. Conditional coefficient inference is not
#' implemented. This restriction also covers semi-partial and constrained
#' regression statistics.
#'
#' Correlation models test marginal associations, without adjusting for the
#' other RDMs, and single-predictor regressions remain available with the
#' default \code{rsa_null = "individual"}. RSA results record
#' \code{rsa_null}, \code{null_hypothesis}, and \code{null_predictors}.
#' FDR adjustment is performed separately within each metric's spatial map;
#' it does not correct across metrics.
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
#' @param model_spec An \code{mvpa_model}, \code{rsa_model}, or compatible
#'   model specification whose design has a \code{\link{permute_labels}}
#'   method.  For \code{rsa_model} the permutation relabels items, so the null
#'   carries the same item-level dependence as the observed statistic.
#' @param observed Optional pre-computed searchlight result (output of
#'   \code{run_searchlight()}).  If \code{NULL}, it is computed internally.
#'   Can also be a named numeric vector of observed metric values indexed by
#'   center ID.
#' @param radius Numeric. Searchlight radius in mm (default 8).
#' @param method Character. Searchlight method passed to
#'   \code{run_searchlight()} when \code{observed} is \code{NULL} or when
#'   \code{perm_strategy = "searchlight"}.
#' @param perm_ctrl A \code{\link{permutation_control}} object.
#' @param metric Character vector naming the performance metric(s) to test.
#'   If \code{NULL} (default), the first metric is used, except for
#'   \code{rsa_model} specifications, where every model predictor is tested:
#'   a permutation pass returns all predictors at once, so scoring them
#'   uses the same permutation passes. The null hypothesis is controlled by
#'   \code{perm_ctrl$rsa_null}, not by how many metrics are selected. A metric that
#'   is named explicitly but absent from the results is an error rather than
#'   a silent fallback to the first metric.
#' @param ... Additional arguments forwarded to \code{run_searchlight()}
#'   (observed pass and \code{"searchlight"} permutations).  For
#'   \code{"iterate"} permutations, only arguments that are formal
#'   parameters of \code{\link{mvpa_iterate}} are forwarded; other keys
#'   are ignored for that path to avoid argument-mismatch failures
#'   (e.g., \code{engine = "legacy"} is meaningful for
#'   \code{run_searchlight} but not for \code{mvpa_iterate}).
#'
#' @return When one metric is tested, a \code{permutation_result} S3 object.
#'   When several are tested (the default for \code{rsa_model}), a
#'   \code{permutation_result_set}: a named list of \code{permutation_result}
#'   objects, one per metric, all scored against the same permutations.  Each
#'   \code{permutation_result} contains:
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
#'     \item{rsa_null, null_hypothesis, null_predictors}{For RSA, the null
#'       scope, its description, and the predictors covered by that null.
#'       These fields are \code{NULL} for other model classes.}
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
  # Fail before the observed pass, inspecting all fitted predictors rather
  # than only the requested metric or the predictors tagged as models.
  rsa_null <- .resolve_rsa_permutation_null(model_spec, perm_ctrl)

  # Step 1: run observed searchlight if not supplied
  if (is.null(observed)) {
    futile.logger::flog.info("Running observed searchlight (radius = %g mm) ...", radius)
    observed <- do.call(
      run_searchlight,
      c(list(model_spec = model_spec, radius = radius, method = method), dots)
    )
  }

  # Which metrics to score.  Resolved after the observed pass so an inferred
  # metric can be read off the searchlight_result.
  metric_info   <- .resolve_permutation_metrics(model_spec, observed, metric)
  metrics       <- metric_info$metrics
  metric_strict <- metric_info$strict

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

  null_values    <- matrix(numeric(0), nrow = 0L, ncol = length(metrics),
                           dimnames = list(NULL, metrics))
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

      perm_mat <- .extract_perf_matrix(perm_result, metrics, strict = metric_strict)
      keep     <- rowSums(!is.na(perm_mat)) > 0L
      if (!any(keep)) next
      perm_mat <- perm_mat[keep, , drop = FALSE]

      null_values    <- rbind(null_values, perm_mat)
      null_nfeatures <- c(null_nfeatures,
                          sub$covariates$nfeatures[
                            match(rownames(perm_mat), as.character(sub$center_ids))
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

      perm_mat <- .searchlight_result_matrix(perm_sl, all_ids, metrics,
                                             strict = metric_strict)
      keep     <- rowSums(!is.na(perm_mat)) > 0L
      if (!any(keep)) next

      null_values    <- rbind(null_values, perm_mat[keep, , drop = FALSE])
      null_nfeatures <- c(null_nfeatures, nfeatures_all[keep])
    }
  }

  if (nrow(null_values) == 0L) {
    stop("No valid null values collected. All permutations may have failed.")
  }

  # Step 6: observed values for every center and metric
  obs_mat <- if (is.numeric(observed)) {
    matrix(.align_to_ids(observed, all_ids), ncol = 1L,
           dimnames = list(as.character(all_ids), metrics))
  } else {
    .searchlight_result_matrix(observed, all_ids, metrics, strict = metric_strict)
  }

  # Steps 7-11 run once per metric against the shared null pool
  results <- lapply(metrics, function(m) {
    .score_permutation_metric(
      metric          = m,
      null_col        = null_values[, m],
      null_nfeatures  = null_nfeatures,
      obs_vals        = obs_mat[, m],
      all_ids         = all_ids,
      covariates_full = covariates_full,
      perm_ctrl       = perm_ctrl,
      model_spec      = model_spec,
      observed        = observed,
      strategy        = strategy,
      rsa_null        = rsa_null
    )
  })
  names(results) <- metrics

  if (length(results) == 1L) {
    return(results[[1L]])
  }
  structure(results, class = c("permutation_result_set", "list"))
}


#' Validate the hypothesis supported by raw RSA item permutations
#' @keywords internal
#' @noRd
.resolve_rsa_permutation_null <- function(model_spec, perm_ctrl) {
  if (!inherits(model_spec, "rsa_model")) return(NULL)
  mode <- match.arg(perm_ctrl$rsa_null %||% "individual", c("individual", "joint"))
  if (identical(mode, "individual") &&
      model_spec$regtype %in% c("lm", "rfit") &&
      length(model_spec$design$model_mat) > 1L) {
    stop(paste0(
      "Raw RSA item permutations cannot test an individual regression coefficient ",
      "conditional on other design predictors (including nuisance predictors). ",
      "Use permutation_control(rsa_null = 'joint') only for the joint null of ",
      "no association with any design predictor; conditional coefficient inference ",
      "is not implemented. Selecting one `metric` does not change this restriction."
    ), call. = FALSE)
  }
  mode
}


#' Resolve which metrics a permutation run scores
#'
#' An explicit \code{metric} wins.  Otherwise \code{rsa_model} tests every
#' model predictor (its per-ROI output already holds all of them), and any
#' other model type tests the first metric of the observed result.
#'
#' Metrics that come from the user or from the model spec are \emph{strict}:
#' they must be present by name in every result, or the run fails.  Only the
#' legacy path, where the metric is inferred from the observed result or is
#' the positional placeholder \code{"metric"}, keeps the first-entry fallback.
#'
#' @return A list with \code{metrics} (character) and \code{strict} (logical).
#' @keywords internal
#' @noRd
.resolve_permutation_metrics <- function(model_spec, observed, metric) {
  strict <- TRUE
  if (!is.null(metric)) {
    if (!is.character(metric) || length(metric) == 0L || anyNA(metric) ||
        any(!nzchar(metric))) {
      stop("`metric` must be NULL or a character vector of metric names.",
           call. = FALSE)
    }
    metrics <- unique(metric)
  } else if (inherits(model_spec, "rsa_model")) {
    metrics <- model_spec$design$model_predictors
    if (is.null(metrics) || length(metrics) == 0L) {
      metrics <- names(model_spec$design$model_mat)
    }
  } else if (inherits(observed, "searchlight_result") &&
             length(observed$metrics) > 0L) {
    metrics <- observed$metrics[[1L]]
    strict  <- FALSE
  } else {
    metrics <- "metric"
    strict  <- FALSE
  }

  if (is.numeric(observed) && length(metrics) > 1L) {
    stop(paste0(
      "`observed` is a numeric vector, which carries a single metric; ",
      "supply `metric` as one name."
    ), call. = FALSE)
  }

  # A strict metric must exist in the observed result; fail before any
  # permutation is run rather than after all of them come back empty, and
  # never let a missing observed map fall back to the first one.
  if (strict && inherits(observed, "searchlight_result")) {
    available <- names(observed$results)
    if (is.null(available)) available <- as.character(observed$metrics)
    missing <- setdiff(metrics, available)
    if (length(available) > 0L && length(missing) > 0L) {
      stop(sprintf(
        "metric '%s' not found in the observed searchlight result; available: %s.",
        paste(missing, collapse = "', '"), paste(available, collapse = ", ")
      ), call. = FALSE)
    }
  }
  list(metrics = metrics, strict = strict)
}


#' Score one metric against its null pool (steps 6-11 of the pipeline)
#'
#' @keywords internal
#' @noRd
.score_permutation_metric <- function(metric, null_col, null_nfeatures, obs_vals,
                                      all_ids, covariates_full, perm_ctrl,
                                      model_spec, observed, strategy, rsa_null = NULL) {
  ok <- !is.na(null_col)
  if (!any(ok)) {
    stop(sprintf("No valid null values were collected for metric '%s'.", metric),
         call. = FALSE)
  }
  if (all(is.na(obs_vals))) {
    stop(sprintf("No observed values were found for metric '%s'.", metric),
         call. = FALSE)
  }
  null_col        <- null_col[ok]
  null_covariates <- data.frame(nfeatures = null_nfeatures[ok])

  # Step 6: optional diagnostics on the null
  diag_result <- NULL
  if (perm_ctrl$diagnose) {
    futile.logger::flog.info("Running null diagnostics for '%s' ...", metric)
    diag_result <- diagnose_null(null_col, null_covariates,
                                 n_perm = perm_ctrl$n_perm)
    print(diag_result)
  }

  # Step 7: build adjusted null
  futile.logger::flog.info("Building adjusted null distribution for '%s' (%s, %d bins) ...",
                           metric, perm_ctrl$null_method, perm_ctrl$n_bins)
  adj_null <- build_adjusted_null(
    null_values  = null_col,
    covariates   = null_covariates,
    n_bins       = perm_ctrl$n_bins,
    method       = perm_ctrl$null_method
  )

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
  futile.logger::flog.info("Done '%s'. %d centers significant at FDR < 0.05 (%s).",
                           metric, n_sig, perm_ctrl$correction)

  null_hypothesis <- if (identical(rsa_null, "joint")) {
    "Joint no association with any design predictor, including nuisance predictors."
  } else if (identical(rsa_null, "individual")) {
    sprintf("No marginal association with '%s' under item exchangeability.", metric)
  } else {
    NULL
  }
  null_predictors <- if (identical(rsa_null, "joint")) {
    names(model_spec$design$model_mat)
  } else if (identical(rsa_null, "individual")) {
    metric
  } else {
    NULL
  }

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
      rsa_null     = rsa_null,
      null_hypothesis = null_hypothesis,
      null_predictors = null_predictors,
      perm_strategy = strategy,
      all_ids      = all_ids,
      n_perm_used  = perm_ctrl$n_perm,
      n_null_vals  = length(null_col)
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
  if (!is.null(x$null_hypothesis)) cat("  Null hypothesis  :", x$null_hypothesis, "\n")
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
print.permutation_result_set <- function(x, ...) {
  first <- x[[1L]]
  cat("Permutation Searchlight Results (", length(x), " metrics)\n", sep = "")
  cat(strrep("=", 40), "\n")
  if (identical(first$rsa_null, "joint")) {
    cat("  Null hypothesis  :", first$null_hypothesis, "\n")
  } else if (identical(first$rsa_null, "individual")) {
    cat("  Null hypotheses  : individual marginal associations\n")
  }
  cat("  Strategy         :", first$perm_strategy, "\n")
  cat("  Permutations used:", first$n_perm_used, "\n")
  cat("  Null values      :", first$n_null_vals, "\n")
  cat("  Total centers    :", length(first$all_ids), "\n")
  cat("  Correction       :", first$perm_ctrl$correction, "\n\n")
  cat(sprintf("  %-24s %12s %12s\n", "metric", "raw p<0.05", "adj p<0.05"))
  for (m in names(x)) {
    r <- x[[m]]
    cat(sprintf("  %-24s %12d %12d\n", m,
                sum(r$p_values < 0.05, na.rm = TRUE),
                sum(r$p_adjusted < 0.05, na.rm = TRUE)))
  }
  invisible(x)
}

#' @export
summary.permutation_result_set <- function(object, ...) {
  for (m in names(object)) {
    cat("\n[", m, "]\n", sep = "")
    summary(object[[m]], ...)
  }
  invisible(object)
}

#' @export
summary.permutation_result <- function(object, ...) {
  cat("Permutation Searchlight Summary\n")
  cat(strrep("=", 50), "\n")
  cat("  Metric           :", object$metric, "\n")
  if (!is.null(object$null_hypothesis)) cat("  Null hypothesis  :", object$null_hypothesis, "\n")
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
