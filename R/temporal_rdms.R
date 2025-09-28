#' Temporal/ordinal nuisance RDM (trial-level)
#'
#' Creates a temporal or ordinal nuisance representational dissimilarity matrix (RDM) 
#' for use in RSA analyses. This function generates various kernels to model temporal
#' proximity effects in neuroimaging data.
#'
#' @param index numeric or integer vector representing trial order or time (length N observations)
#' @param block optional vector (length N) of run/block identifiers
#' @param kernel character string specifying the kernel type, one of:
#'   \itemize{
#'     \item{"adjacent"}{Binary kernel for immediate neighbors within specified width}
#'     \item{"boxcar"}{Binary kernel including diagonal up to specified width}
#'     \item{"linear"}{Linear decay based on lag}
#'     \item{"poly"}{Polynomial decay with specified power}
#'     \item{"exp"}{Exponential decay with specified lambda}
#'     \item{"gauss"}{Gaussian decay with specified sigma}
#'   }
#' @param width integer window for "adjacent"/"boxcar" kernels (default 1)
#' @param power exponent for "poly" kernel (default 1)
#' @param lambda decay constant for "exp" kernel (default 1)
#' @param sigma standard deviation for "gauss" kernel (default 1)
#' @param within_blocks_only logical; if TRUE, zero nuisance across blocks (default TRUE)
#' @param wrap logical; if TRUE, treat index as circular (default FALSE)
#' @param normalize character string specifying normalization method:
#'   \itemize{
#'     \item{"rank"}{Rank transform (ties averaged)}
#'     \item{"z"}{Z-score normalization}
#'     \item{"none"}{No normalization}
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
#' proximity effects while preserving statistical power.
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
                        lambda = 1, 
                        sigma = 1,
                        within_blocks_only = TRUE, 
                        wrap = FALSE,
                        normalize = c("rank", "z", "none"),
                        as_dist = TRUE) {
  
  kernel <- match.arg(kernel)
  normalize <- match.arg(normalize)
  index <- as.numeric(index)
  n <- length(index)
  
  stopifnot(length(block) %in% c(0, n))
  
  # Compute pairwise lags |delta|
  lag <- abs(outer(index, index, "-"))
  
  if (isTRUE(wrap)) {
    lag <- pmin(lag, n - lag)  # circular distance
  }
  
  # Apply kernel
  K <- switch(kernel,
    adjacent = (lag > 0 & lag <= width) * 1,
    boxcar   = (lag <= width) * 1,
    linear   = lag,
    poly     = lag^power,
    exp      = exp(-lag / lambda),
    gauss    = exp(-(lag^2) / (2 * sigma^2))
  )
  
  diag(K) <- 0
  
  # Constrain to within blocks if requested
  if (!is.null(block) && isTRUE(within_blocks_only)) {
    same_block <- outer(block, block, "==")
    K[!same_block] <- 0
  }
  
  # Extract lower triangle for normalization
  if (as_dist) {
    dv <- as.vector(as.dist(K))
  } else {
    dv <- K[lower.tri(K)]
  }
  
  # Apply normalization
  if (normalize == "rank") {
    dv <- rank(dv, ties.method = "average")
  } else if (normalize == "z") {
    m <- mean(dv)
    s <- sd(dv)
    dv <- if (s > 0) (dv - m) / s else dv * 0
  }
  
  # Return in requested format
  if (as_dist) {
    # Create a dist object with the normalized values
    attr_d <- attributes(as.dist(K))
    out <- structure(dv, 
                    class = "dist",
                    Size = n,
                    call = match.call(),
                    method = paste("temporal", kernel))
    out
  } else {
    # Reconstruct full matrix
    K[lower.tri(K)] <- dv
    K[upper.tri(K)] <- t(K)[upper.tri(K)]
    K
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
#'   \itemize{
#'     \item{"min"}{Minimum lag between any pair of trials from two conditions}
#'     \item{"mean"}{Average lag between all pairs of trials}
#'     \item{"median"}{Median lag between all pairs of trials}
#'     \item{"nn_min"}{Minimum nearest-neighbor distance}
#'   }
#' @param kernel character string specifying kernel type (see \code{\link{temporal_rdm}})
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
                                         within_blocks_only = TRUE,
                                         ...) {
  
  reduce <- match.arg(reduce)
  kernel <- match.arg(kernel)
  
  stopifnot(length(time_idx) == nrow(mvpa_design$train_design))
  
  # Get condition labels and blocks
  cond <- as.character(mvpa_design$Y)
  blocks <- mvpa_design$block_var
  lev <- levels(mvpa_design$Y)
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
  lambda <- dots$lambda %||% 1
  sigma <- dots$sigma %||% 1
  
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
  Kmat[is.infinite(Lag)] <- 0
  
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
