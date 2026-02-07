#' Temporal/ordinal nuisance RDM (trial-level)
#'
#' Creates a temporal or ordinal nuisance representational dissimilarity matrix (RDM)
#' for use in RSA analyses. This function generates various kernels to model temporal
#' proximity effects in neuroimaging data and can return either similarity- or
#' distance-like outputs.
#'
#' @param index numeric or integer vector representing trial order or time (length N observations)
#' @param block optional vector (length N) of run/block identifiers
#' @param kernel character string specifying the kernel type, one of:
#'   \describe{
#'     \item{adjacent}{Binary kernel for immediate neighbors within specified width}
#'     \item{boxcar}{Binary kernel including diagonal up to specified width}
#'     \item{linear}{Linear distance based on lag}
#'     \item{poly}{Polynomial distance with specified power}
#'     \item{exp}{Exponential similarity with specified lambda}
#'     \item{gauss}{Gaussian similarity with specified sigma}
#'   }
#' @param width integer window for "adjacent"/"boxcar" kernels (default 1)
#' @param power exponent for "poly" kernel (default 1)
#' @param lambda decay constant for "exp" kernel (in units; default chosen automatically if NULL)
#' @param sigma standard deviation for "gauss" kernel (in units; default chosen automatically if NULL)
#' @param within_blocks_only logical; if TRUE, zero nuisance across blocks (default TRUE)
#' @param wrap logical; if TRUE, treat index as circular (default FALSE)
#' @param units one of \code{"auto"}, \code{"sec"}, \code{"TR"}, or \code{"index"}; used to choose sensible defaults for decay parameters
#' @param TR repetition time in seconds (optional, but recommended when \code{units="sec"} or \code{units="TR"})
#' @param metric return type: \code{"distance"} (increasing with lag) or \code{"similarity"} (decreasing with lag). Default \code{"distance"}.
#' @param normalize character string specifying normalization method:
#'   \describe{
#'     \item{rank}{Rank transform (ties averaged)}
#'     \item{z}{Z-score normalization}
#'     \item{none}{No normalization}
#'   }
#' @param as_dist logical; if TRUE return a \code{dist} object, otherwise return matrix (default TRUE)
#' 
#' @return A \code{dist} object or symmetric matrix (N x N) with 0 on the diagonal,
#'   representing temporal/ordinal relationships between observations
#'   
#' @details
#' This function creates temporal nuisance RDMs for modeling carry-over effects,
#' scanner drift, or other temporal confounds in fMRI data. The resulting RDM
#' can be included as a nuisance regressor in RSA models to account for temporal
#' proximity effects while preserving statistical power. By default, it returns
#' a distance-like quantity that increases with temporal separation (\code{metric="distance"}).
#'
#' @examples
#' # Create temporal RDM for 20 trials across 4 runs
#' trial_index <- 1:20
#' run_labels <- rep(1:4, each = 5)
#' 
#' # Exponential decay kernel within runs only
#' temp_rdm <- temporal_rdm(trial_index, block = run_labels, 
#'                          kernel = "exp", lambda = 2,
#'                          within_blocks_only = TRUE)
#' 
#' # Use in RSA design
#' \dontrun{
#' rdes <- rsa_design(~ task_rdm + temporal_rdm(trial_idx, block=run, kernel="adjacent"),
#'                    data = list(task_rdm = my_task_rdm,
#'                               trial_idx = seq_len(n_trials),
#'                               run = run_ids),
#'                    block_var = ~ run,
#'                    keep_intra_run = TRUE)
#' }
#' 
#' @export
#' @importFrom stats dist sd
temporal_rdm <- function(index, 
                        block = NULL,
                        kernel = c("adjacent", "boxcar", "linear", "poly", "exp", "gauss"),
                        width = 1L, 
                        power = 1, 
                        lambda = NULL, 
                        sigma = NULL,
                        within_blocks_only = TRUE, 
                        wrap = FALSE,
                        units = c("auto", "sec", "TR", "index"),
                        TR = NULL,
                        metric = c("similarity", "distance"),
                        normalize = c("rank", "z", "none"),
                        as_dist = TRUE) {
  
  kernel <- match.arg(kernel)
  normalize <- match.arg(normalize)
  units <- match.arg(units)
  metric <- match.arg(metric)
  index <- as.numeric(index)
  n <- length(index)
  
  stopifnot(length(block) %in% c(0, n))
  
  # Compute pairwise lags |delta|
  lag <- abs(outer(index, index, "-"))
  
  if (isTRUE(wrap)) {
    lag <- pmin(lag, n - lag)  # circular distance
  }
  
  # Determine sensible defaults for exp/gauss BEFORE computing kernel
  if (units == "auto") {
    units <- infer_units(index, TR)
  }
  if (kernel %in% c("exp", "gauss")) {
    dparams <- default_decay_params(index, TR, kernel, lambda, sigma, units)
    if (is.null(lambda)) lambda <- dparams$lambda
    if (is.null(sigma))  sigma  <- dparams$sigma
  }

  # Apply kernel once using effective parameters
  K <- switch(kernel,
    adjacent = (lag > 0 & lag <= width) * 1,
    boxcar   = (lag <= width) * 1,
    linear   = lag,
    poly     = lag^power,
    exp      = exp(-lag / lambda),
    gauss    = exp(-(lag^2) / (2 * sigma^2))
  )
  
  diag(K) <- 0
  
  # Constrain to within blocks if requested (mask cross-run as NA until finalization)
  if (!is.null(block) && isTRUE(within_blocks_only)) {
    same_block <- outer(block, block, "==")
    K[!same_block] <- NA_real_
  }

  # Extract lower triangle
  dv <- as.vector(as.dist(K))

  # Convert to distance if requested (do this BEFORE handling NA so masked pairs become 0 later)
  if (metric == "distance") {
    if (kernel %in% c("adjacent", "boxcar", "exp", "gauss")) {
      # Handle edge case where all values are NA (e.g., all pairs masked by within_blocks_only)
      mx <- if (all(is.na(dv))) 1 else max(dv, na.rm = TRUE)
      dv <- if (is.finite(mx) && mx > 0) 1 - (dv / mx) else 1 - dv
    } else {
      # linear/poly already distance-like
    }
  }

  # Apply normalization (ignore masked pairs = NA, normalize only non-NA), then set NA->0
  if (normalize == "rank") {
    nz <- !is.na(dv)
    dv[nz] <- rank(dv[nz], ties.method = "average")
  } else if (normalize == "z") {
    nz <- !is.na(dv)
    m <- mean(dv[nz])
    s <- sd(dv[nz])
    dv[nz] <- if (s > 0) (dv[nz] - m) / s else dv[nz] * 0
  }

  # Finally, masked pairs contribute zero nuisance
  dv[is.na(dv)] <- 0

  # Return in requested format
  if (as_dist) {
    out <- structure(dv,
                    class = "dist",
                    Size = n,
                    call = match.call(),
                    method = paste("temporal", kernel, metric, sep = ":"))
    out
  } else {
    M <- matrix(0, n, n)
    M[lower.tri(M)] <- dv
    M <- M + t(M)
    M
  }
}


#' Temporal nuisance RDM at condition level (for MS-ReVE)
#'
#' Creates a condition-level temporal nuisance RDM for use with MS-ReVE/contrast_rsa_model.
#' This function reduces trial-level temporal relationships to condition-level relationships.
#'
#' @param mvpa_design an mvpa_design object containing Y (condition labels) and block_var
#' @param time_idx numeric/integer vector of temporal indices, length = nrow(mvpa_design$train_design)
#' @param reduce character string specifying reduction method:
#'   \describe{
#'     \item{min}{Minimum lag between any pair of trials from two conditions}
#'     \item{mean}{Average lag between all pairs of trials}
#'     \item{median}{Median lag between all pairs of trials}
#'     \item{nn_min}{Minimum nearest-neighbor distance}
#'   }
#' @param kernel character string specifying kernel type (see \code{\link{temporal_rdm}})
#' @param units one of \code{"auto"}, \code{"sec"}, \code{"TR"}, or \code{"index"}; used to select sensible defaults for decay
#' @param TR repetition time in seconds (optional)
#' @param metric return type: \code{"distance"} or \code{"similarity"} (default \code{"distance"})
#' @param within_blocks_only logical; if TRUE, ignore pairs spanning different blocks (default TRUE)
#' @param ... additional kernel parameters passed to kernel functions (width, power, lambda, sigma)
#' 
#' @return K x K symmetric matrix (0 diagonal) aligned to levels(mvpa_design$Y)
#' 
#' @details
#' This function is designed for MS-ReVE analyses where temporal confounds need to be
#' modeled at the condition level rather than trial level. It computes aggregate temporal
#' relationships between conditions based on the temporal structure of individual trials.
#'
#' @examples
#' \dontrun{
#' # Create temporal nuisance for MS-ReVE
#' temp_K <- temporal_nuisance_for_msreve(
#'   mvpa_design = mvpa_des,
#'   time_idx = seq_len(nrow(mvpa_des$train_design)),
#'   reduce = "min",
#'   kernel = "exp", 
#'   lambda = 3,
#'   within_blocks_only = TRUE
#' )
#' 
#' # Use in msreve_design
#' msreve_des <- msreve_design(
#'   mvpa_design = mvpa_des,
#'   contrast_matrix = C_mat,
#'   nuisance_rdms = list(temp_decay = temp_K)
#' )
#' }
#' 
#' @export
temporal_nuisance_for_msreve <- function(mvpa_design, 
                                         time_idx,
                                         reduce = c("min", "mean", "median", "nn_min"),
                                         kernel = c("adjacent", "boxcar", "linear", 
                                                   "poly", "exp", "gauss"),
                                         units = c("auto", "sec", "TR", "index"),
                                         TR = NULL,
                                         metric = c("distance", "similarity"),
                                         within_blocks_only = TRUE,
                                         ...) {
  
  reduce <- match.arg(reduce)
  kernel <- match.arg(kernel)
  units <- match.arg(units)
  metric <- match.arg(metric)
  
  stopifnot(length(time_idx) == nrow(mvpa_design$train_design))
  
  # Get condition labels and blocks
  cond_vec <- mvpa_design$y_train
  if (is.null(cond_vec) && !is.null(mvpa_design$Y)) cond_vec <- mvpa_design$Y
  if (is.null(cond_vec) && !is.null(mvpa_design$train_design)) {
    # Try common column names "cond" or first factor column
    td <- mvpa_design$train_design
    if ("cond" %in% names(td)) cond_vec <- factor(td$cond)
    if (is.null(cond_vec)) {
      fac_cols <- names(td)[vapply(td, is.factor, logical(1))]
      if (length(fac_cols) >= 1) cond_vec <- factor(td[[fac_cols[1]]])
    }
  }
  if (is.null(cond_vec)) stop("temporal_nuisance_for_msreve: mvpa_design$y_train is NULL and no fallback found")
  cond <- as.character(cond_vec)
  blocks <- mvpa_design$block_var
  lev <- levels(cond_vec)
  K <- length(lev)
  
  # Initialize output matrix
  M <- matrix(0, K, K, dimnames = list(lev, lev))
  
  # Split indices by condition
  split_idx <- split(seq_along(cond), cond)
  
  # Helper function to aggregate lags between two conditions
  agg_lag <- function(i, j) {
    idx_i <- split_idx[[i]]
    idx_j <- split_idx[[j]]
    
    if (within_blocks_only && !is.null(blocks)) {
      # Create mask for same-block pairs
      ok <- outer(blocks[idx_i], blocks[idx_j], "==")
      if (!any(ok, na.rm = TRUE)) return(Inf)
      
      # Compute temporal differences
      tdiff <- abs(outer(time_idx[idx_i], time_idx[idx_j], "-"))
      tdiff <- tdiff[ok]
    } else {
      tdiff <- abs(outer(time_idx[idx_i], time_idx[idx_j], "-"))
    }
    
    if (!length(tdiff)) return(Inf)
    
    # Apply reduction method
    switch(reduce,
      min    = min(tdiff, na.rm = TRUE),
      mean   = mean(tdiff, na.rm = TRUE),
      median = stats::median(tdiff, na.rm = TRUE),
      nn_min = min(apply(tdiff, 1, min, na.rm = TRUE), na.rm = TRUE)
    )
  }
  
  # Compute lag matrix
  Lag <- matrix(0, K, K)
  for (a in seq_len(K - 1)) {
    for (b in (a + 1):K) {
      Lag[a, b] <- Lag[b, a] <- agg_lag(lev[a], lev[b])
    }
  }
  diag(Lag) <- 0
  
  # Get additional parameters with defaults
  dots <- list(...)
  width <- dots$width %||% 1
  power <- dots$power %||% 1
  lambda <- dots$lambda %||% NULL
  sigma <- dots$sigma %||% NULL

  # Sensible defaults if needed
  if (units == "auto") units <- infer_units(time_idx, TR)
  dparams <- default_decay_params(time_idx, TR, kernel, lambda, sigma, units)
  if (is.null(lambda)) lambda <- dparams$lambda
  if (is.null(sigma))  sigma  <- dparams$sigma
  
  # Apply kernel
  Kmat <- switch(kernel,
    adjacent = (Lag > 0 & Lag <= width) * 1,
    boxcar   = (Lag <= width) * 1,
    linear   = Lag,
    poly     = Lag^power,
    exp      = exp(-Lag / lambda),
    gauss    = exp(-(Lag^2) / (2 * sigma^2))
  )
  
  diag(Kmat) <- 0
  Kmat[is.infinite(Lag)] <- NA_real_

  # Convert to distance if requested
  if (metric == "distance") {
    if (kernel %in% c("adjacent", "boxcar", "exp", "gauss")) {
      mx <- max(Kmat[lower.tri(Kmat)], na.rm = TRUE)
      if (is.finite(mx) && mx > 0) {
        Kmat <- 1 - (Kmat / mx)
      } else {
        Kmat <- 1 - Kmat
      }
      diag(Kmat) <- 0
    }
  }
  # Replace NAs (e.g., cross-block or Inf lags) with 0 to indicate no nuisance
  Kmat[is.na(Kmat)] <- 0
  
  Kmat
}


#' Temporal RDM wrapper for formula usage
#'
#' Convenience wrapper for \code{\link{temporal_rdm}} that simplifies usage in RSA formulas.
#'
#' @param index numeric or integer vector representing trial order or time
#' @param block optional vector of run/block identifiers
#' @param ... additional parameters passed to \code{\link{temporal_rdm}}
#' @param as_dist logical; if TRUE return a dist object (default TRUE)
#' 
#' @return A dist object or matrix representing temporal relationships
#' 
#' @details
#' This function provides a shorter name for use in RSA design formulas.
#' It calls \code{temporal_rdm} with the same parameters.
#'
#' @examples
#' \dontrun{
#' # Use directly in RSA formula
#' rdes <- rsa_design(
#'   ~ task_rdm + temporal(trial_index, block=run, kernel="adjacent", width=2),
#'   data = list(task_rdm = task_rdm, trial_index = 1:100, run = run_ids),
#'   block_var = ~ run
#' )
#' }
#' 
#' @export
#' @seealso \code{\link{temporal_rdm}}
temporal <- function(index, block = NULL, ..., as_dist = TRUE) {
  temporal_rdm(index = index, block = block, ..., as_dist = as_dist)
}


# Internal helper for NULL-coalescing operator
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


# Internal helpers ---------------------------------------------------------

#' @keywords internal
infer_units <- function(index, TR) {
  is_intish <- all(abs(index - round(index)) < 1e-6)
  if (!is_intish) return("sec")
  if (!is.null(TR)) return("TR")
  "index"
}

#' @keywords internal
default_decay_params <- function(index, TR, kernel, lambda, sigma, units) {
  out <- list(lambda = lambda, sigma = sigma)
  if (kernel %in% c("exp", "gauss", "adjacent", "boxcar")) {
    idx_sorted <- sort(unique(index))
    step <- if (length(idx_sorted) > 1) stats::median(diff(idx_sorted), na.rm = TRUE) else 1
    base_scale <- step * 2
    if (is.null(lambda)) out$lambda <- ifelse(isTRUE(is.finite(base_scale)), base_scale, 1)
    if (is.null(sigma))  out$sigma  <- ifelse(isTRUE(is.finite(base_scale)), base_scale, 1)
  }
  out
}

#' @keywords internal
to_distance <- function(dv) {
  mx <- max(dv, na.rm = TRUE)
  if (is.finite(mx) && mx > 0) 1 - (dv / mx) else 1 - dv
}


# Convenience wrappers -----------------------------------------------------

#' Temporal RDM from onsets (sugar)
#'
#' Convenience wrapper that accepts event onsets (and optional run labels) and
#' delegates to \code{temporal_rdm}. This is useful in RSA formulas.
#'
#' @param onsets numeric vector of event onsets (typically seconds or TRs)
#' @param run optional vector of run/block identifiers
#' @param ... additional parameters passed to \code{temporal_rdm}
#' @param units units for \code{onsets}: one of \code{"auto"}, \code{"sec"}, \code{"TR"}, or \code{"index"}
#' @param TR repetition time in seconds (optional)
#' @param as_dist logical; if TRUE return a dist object (default TRUE)
#'
#' @return A dist object or matrix representing temporal relationships
#' @export
temporal_from_onsets <- function(onsets, run = NULL, ..., units = c("auto", "sec", "TR", "index"), TR = NULL, as_dist = TRUE) {
  units <- match.arg(units)
  temporal_rdm(index = onsets, block = run, units = units, TR = TR, ..., as_dist = as_dist)
}


#' HRF-overlap temporal confound RDM
#'
#' Builds a temporal confound RDM based on predicted HRF overlap between events.
#' Each event is convolved with a canonical HRF and pairwise similarity is computed
#' as either normalized overlap (cosine similarity) or Pearson correlation. The
#' result can be returned as a distance-like quantity suitable for regression.
#'
#' @param onsets numeric vector of event onsets in seconds
#' @param durations optional numeric vector of event durations in seconds (scalar or length N). Default 0 (impulse).
#' @param run optional vector of run/block identifiers
#' @param TR repetition time in seconds (required)
#' @param hrf character: one of \code{"spm"}, \code{"glover"}, or \code{"gamma"}
#' @param oversampling integer oversampling factor relative to TR (default 16)
#' @param length_s numeric length of HRF kernel in seconds (default 32)
#' @param similarity one of \code{"overlap"} (cosine similarity) or \code{"corr"}
#' @param metric return type: \code{"distance"} or \code{"similarity"} (default \code{"distance"})
#' @param within_blocks_only logical; if TRUE zero-out cross-run entries (default TRUE)
#' @param normalize one of \code{"rank"}, \code{"z"}, or \code{"none"}
#' @param as_dist logical; if TRUE return a \code{dist} object (default TRUE)
#'
#' @return A \code{dist} object or symmetric matrix (N x N)
#' @export
#' @importFrom stats cor convolve sd as.dist
#' @importFrom utils head
temporal_hrf_overlap <- function(onsets,
                                 durations = NULL,
                                 run = NULL,
                                 TR,
                                 hrf = c("spm", "glover", "gamma"),
                                 oversampling = 16L,
                                 length_s = 32,
                                 similarity = c("overlap", "corr"),
                                 metric = c("distance", "similarity"),
                                 within_blocks_only = TRUE,
                                 normalize = c("z", "rank", "none"),
                                 as_dist = TRUE) {

  hrf <- match.arg(hrf)
  similarity <- match.arg(similarity)
  metric <- match.arg(metric)
  normalize <- match.arg(normalize)

  onsets <- as.numeric(onsets)
  n <- length(onsets)
  if (is.null(TR) || !is.finite(TR) || TR <= 0) stop("TR must be provided and > 0 for temporal_hrf_overlap().")
  if (is.null(durations)) durations <- rep(0, n)
  if (length(durations) == 1) durations <- rep(durations, n)
  stopifnot(length(durations) == n)
  stopifnot(length(run) %in% c(0, n))

  # Time grid
  dt <- TR / oversampling
  t0 <- min(onsets)
  on <- onsets - t0
  t_end <- max(on + durations) + length_s
  Tlen <- as.integer(ceiling(t_end / dt)) + 1L
  tt <- seq(0, by = dt, length.out = Tlen)

  # HRF kernel
  h <- switch(hrf,
    spm   = hrf_spm(tt),
    glover= hrf_glover(tt),
    gamma = hrf_gamma(tt)
  )
  if (sum(h) > 0) h <- h / sum(h)

  # Build event timecourses
  R <- matrix(0, Tlen, n)
  for (i in seq_len(n)) {
    onset_i <- on[i]
    dur_i <- durations[i]
    idx_on <- max(1L, as.integer(round(onset_i / dt)) + 1L)
    idx_off <- max(idx_on, as.integer(round((onset_i + dur_i) / dt)) + 1L)
    e <- numeric(Tlen)
    e[idx_on:min(idx_off, Tlen)] <- 1
    conv <- stats::convolve(e, rev(h), type = "open")
    R[, i] <- head(conv, Tlen)
  }

  # Similarity matrix
  if (similarity == "overlap") {
    norms <- sqrt(colSums(R^2))
    norms[norms == 0] <- 1
    Rn <- sweep(R, 2, norms, "/")
    S <- crossprod(Rn)
  } else {
    S <- stats::cor(R)
    S[is.na(S)] <- 0
  }
  diag(S) <- 0

  if (!is.null(run) && isTRUE(within_blocks_only)) {
    same_block <- outer(run, run, "==")
    S[!same_block] <- NA_real_
  }

  if (metric == "distance") {
    M <- 1 - S
    diag(M) <- 0
  } else {
    M <- S
  }

  dv <- as.vector(stats::as.dist(M))
  # Normalize on non-NA then set NA->0 so masked cross-run remain zero
  if (normalize == "rank") {
    nz <- !is.na(dv)
    dv[nz] <- rank(dv[nz], ties.method = "average")
  } else if (normalize == "z") {
    nz <- !is.na(dv)
    m <- mean(dv[nz])
    s <- stats::sd(dv[nz])
    dv[nz] <- if (s > 0) (dv[nz] - m) / s else dv[nz] * 0
  }
  dv[is.na(dv)] <- 0

  if (as_dist) {
    structure(dv,
              class = "dist",
              Size = n,
              call = match.call(),
              method = paste("temporal_hrf", hrf, similarity, metric, sep = ":"))
  } else {
    Out <- matrix(0, n, n)
    Out[lower.tri(Out)] <- dv
    Out + t(Out)
  }
}

# Canonical HRFs
hrf_gamma <- function(t, a = 6, b = 0.9) {
  t <- pmax(t, 0)
  (t^(a - 1) * exp(-t / b)) / ((b^a) * gamma(a))
}

hrf_spm <- function(t) {
  p1 <- hrf_gamma(t, a = 6, b = 0.9)
  p2 <- hrf_gamma(t, a = 12, b = 0.9)
  y <- p1 - 0.35 * p2
  y[y < 0] <- 0
  y
}

hrf_glover <- function(t) {
  a1 <- 6; b1 <- 0.9
  a2 <- 12; b2 <- 0.9
  c  <- 0.35
  y <- hrf_gamma(t, a1, b1) - c * hrf_gamma(t, a2, b2)
  y[y < 0] <- 0
  y
}


#' Build multiple temporal confounds from a spec
#'
#' Convenience function to produce a named list of temporal nuisance RDMs from
#' a compact specification. Items can be kernel-based (\code{temporal_rdm}) or
#' HRF-overlap based (\code{temporal_hrf_overlap}).
#'
#' @param spec a named list. Each element is a list describing one confound.
#'  For kernel-based items, include at least \code{kernel}. For HRF items, set
#'  \code{kind="hrf"} and provide \code{TR}. Common fields: \code{within_blocks_only}, \code{normalize}, \code{metric}.
#' @param onsets numeric vector of event onsets (seconds or TRs depending on \code{units})
#' @param run optional vector of run/block identifiers
#' @param units one of \code{"auto"}, \code{"sec"}, \code{"TR"}, or \code{"index"}
#' @param TR repetition time in seconds (optional for kernel, required for HRF)
#' @param as_dist logical; if TRUE return \code{dist} objects (default TRUE)
#'
#' @return A named list of \code{dist} objects (or matrices if \code{as_dist=FALSE})
#' @export
temporal_confounds <- function(spec, onsets, run = NULL, units = c("auto", "sec", "TR", "index"), TR = NULL, as_dist = TRUE) {
  units <- match.arg(units)
  stopifnot(is.list(spec))
  out <- list()
  if (is.null(names(spec))) names(spec) <- paste0("conf", seq_along(spec))
  for (nm in names(spec)) {
    item <- spec[[nm]]
    kind <- item$kind %||% "kernel"
    if (identical(kind, "hrf")) {
      out[[nm]] <- temporal_hrf_overlap(
        onsets = onsets,
        durations = item$durations %||% NULL,
        run = run,
        TR = item$TR %||% TR,
        hrf = item$hrf %||% "spm",
        oversampling = item$oversampling %||% 16L,
        length_s = item$length_s %||% 32,
        similarity = item$similarity %||% "overlap",
        metric = item$metric %||% "distance",
        within_blocks_only = item$within_blocks_only %||% TRUE,
        normalize = item$normalize %||% "rank",
        as_dist = as_dist
      )
    } else {
      out[[nm]] <- temporal_rdm(
        index = onsets,
        block = run,
        kernel = item$kernel %||% "exp",
        width = item$width %||% 1L,
        power = item$power %||% 1,
        lambda = item$lambda %||% NULL,
        sigma = item$sigma %||% NULL,
        within_blocks_only = item$within_blocks_only %||% TRUE,
        wrap = item$wrap %||% FALSE,
        units = units,
        TR = item$TR %||% TR,
        metric = item$metric %||% "distance",
        normalize = item$normalize %||% "rank",
        as_dist = as_dist
      )
    }
  }
  out
}


#' MS-ReVE: build temporal nuisance RDMs from a spec
#'
#' Convenience helper that produces a named list of condition-level temporal
#' nuisance RDMs for MS-ReVE analyses given an \code{mvpa_design}, a time index
#' (trial onsets or order), and a compact specification. Each entry can be a
#' kernel-based nuisance (delegates to \code{temporal_nuisance_for_msreve}) or
#' an HRF-overlap-based nuisance built trial-wise then reduced to condition level.
#'
#' @param mvpa_design an \code{mvpa_design} object
#' @param time_idx numeric/integer time index per trial (seconds, TRs, or ordinal)
#' @param spec a named list of entries. For kernel entries: include \code{kernel}
#'   and optional \code{width,power,lambda,sigma,units,TR,metric}. For HRF entries,
#'   set \code{kind="hrf"} and include \code{TR} plus optional \code{durations,hrf,oversampling,length_s,similarity,metric}.
#' @param reduce reduction method for condition-level aggregation (\code{"min"},\code{"mean"},\code{"median"},\code{"nn_min"})
#' @param within_blocks_only logical; if TRUE ignore across-run pairs
#'
#' @return a named list of K x K matrices aligned to levels(mvpa_design$Y)
#' @export
msreve_temporal_confounds <- function(mvpa_design,
                                      time_idx,
                                      spec,
                                      reduce = c("min", "mean", "median", "nn_min"),
                                      within_blocks_only = TRUE) {

  stopifnot(inherits(mvpa_design, "mvpa_design"))
  stopifnot(length(time_idx) == nrow(mvpa_design$train_design))
  stopifnot(is.list(spec))
  if (is.null(names(spec))) names(spec) <- paste0("conf", seq_along(spec))

  reduce <- match.arg(reduce)
  cond_vec <- mvpa_design$y_train
  if (is.null(cond_vec) && !is.null(mvpa_design$Y)) cond_vec <- mvpa_design$Y
  if (is.null(cond_vec) && !is.null(mvpa_design$train_design)) {
    td <- mvpa_design$train_design
    if ("cond" %in% names(td)) cond_vec <- factor(td$cond)
    if (is.null(cond_vec)) {
      fac_cols <- names(td)[vapply(td, is.factor, logical(1))]
      if (length(fac_cols) >= 1) cond_vec <- factor(td[[fac_cols[1]]])
    }
  }
  if (is.null(cond_vec)) stop("msreve_temporal_confounds: mvpa_design$y_train is NULL and no fallback found")
  cond <- as.character(cond_vec)
  blocks <- mvpa_design$block_var
  lev <- levels(cond_vec)
  K <- length(lev)

  # helper: reduce a trial-by-trial matrix to condition-level
  reduce_trial_mat <- function(M) {
    Mout <- matrix(0, K, K, dimnames = list(lev, lev))
    split_idx <- split(seq_along(cond), cond)
    for (a in seq_len(K - 1)) {
      for (b in (a + 1):K) {
        ia <- split_idx[[lev[a]]]
        ib <- split_idx[[lev[b]]]
        if (within_blocks_only && !is.null(blocks)) {
          ok <- outer(blocks[ia], blocks[ib], "==")
          if (!any(ok, na.rm = TRUE)) { Mout[a, b] <- Mout[b, a] <- 0; next }
          vals <- M[ia, ib, drop = FALSE][ok]
        } else {
          vals <- M[ia, ib, drop = FALSE]
        }
        if (!length(vals)) { Mout[a, b] <- Mout[b, a] <- 0; next }
        Mout[a, b] <- Mout[b, a] <- switch(reduce,
          min    = min(vals, na.rm = TRUE),
          mean   = mean(vals, na.rm = TRUE),
          median = stats::median(vals, na.rm = TRUE),
          nn_min = { # nearest neighbor min across rows
            if (length(ia) == 0 || length(ib) == 0) 0 else {
              tmp <- M[ia, ib, drop = FALSE]
              min(apply(tmp, 1, min, na.rm = TRUE), na.rm = TRUE)
            }
          }
        )
      }
    }
    diag(Mout) <- 0
    Mout
  }

  out <- list()
  for (nm in names(spec)) {
    item <- spec[[nm]]
    kind <- item$kind %||% "kernel"
    if (identical(kind, "hrf")) {
      # trial-level HRF distance or similarity, then reduce
      Mtrial <- temporal_hrf_overlap(
        onsets = time_idx,
        durations = item$durations %||% NULL,
        run = blocks,
        TR = item$TR,
        hrf = item$hrf %||% "spm",
        oversampling = item$oversampling %||% 16L,
        length_s = item$length_s %||% 32,
        similarity = item$similarity %||% "overlap",
        metric = item$metric %||% "distance",
        within_blocks_only = within_blocks_only,
        normalize = item$normalize %||% "rank",
        as_dist = FALSE
      )
      out[[nm]] <- reduce_trial_mat(Mtrial)
    } else {
      # kernel-based condition-level nuisance
      out[[nm]] <- temporal_nuisance_for_msreve(
        mvpa_design = mvpa_design,
        time_idx = time_idx,
        reduce = item$reduce %||% reduce,
        kernel = item$kernel %||% "exp",
        units = item$units %||% "auto",
        TR = item$TR %||% NULL,
        metric = item$metric %||% "distance",
        within_blocks_only = within_blocks_only,
        width = item$width %||% 1L,
        power = item$power %||% 1,
        lambda = item$lambda %||% NULL,
        sigma = item$sigma %||% NULL
      )
    }
  }
  out
}
