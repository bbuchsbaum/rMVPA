# Representational Connectivity (ReNA-RC): repnet_model

#' Representational connectivity model (ReNA-RC)
#'
#' For each ROI or searchlight, computes representational connectivity between
#' the ROI RDM and a seed RDM, optionally controlling for confound RDMs (block,
#' lag, behavior, etc.).
#'
#' @param dataset An \code{mvpa_dataset}.
#' @param design An \code{mvpa_design}.
#' @param repnet_des Output of \code{\link{repnet_design}}.
#' @param distfun Distance function for ROI RDM (e.g. \code{cordist(method = "pearson")}).
#' @param simfun Character: similarity metric between ROI and seed RDMs
#'   (\code{"pearson"} or \code{"spearman"}).
#' @param ... Extra fields stored on the model spec.
#'
#' @return A model spec of class \code{"repnet_model"} compatible with
#'   \code{run_regional()} and \code{run_searchlight()}.
#' @export
repnet_model <- function(dataset,
                         design,
                         repnet_des,
                         distfun = cordist(method = "pearson"),
                         simfun = c("pearson", "spearman"),
                         ...) {

  assertthat::assert_that(inherits(dataset, "mvpa_dataset"))
  assertthat::assert_that(inherits(design,  "mvpa_design"))

  simfun <- match.arg(simfun)

  if (is.character(distfun)) {
    distfun <- create_dist(distfun)
  }
  assertthat::assert_that(
    inherits(distfun, "distfun"),
    msg = "distfun must be a 'distfun' or a shortcut understood by create_dist()."
  )

  create_model_spec(
    "repnet_model",
    dataset       = dataset,
    design        = design,
    key           = repnet_des$key,
    seed_rdm      = repnet_des$seed_rdm,
    confound_rdms = repnet_des$confound_rdms,
    distfun       = distfun,
    simfun        = simfun,
    # We compute per-ROI performance directly in process_roi
    compute_performance = TRUE,
    return_predictions  = FALSE,
    ...
  )
}


#' Per-ROI processing for repnet_model
#'
#' Computes representational connectivity metrics for a single ROI/searchlight.
#'
#' @keywords internal
#' @export
process_roi.repnet_model <- function(mod_spec,
                                     roi,
                                     rnum,
                                     center_global_id = NA,
                                     ...) {

  X   <- as.matrix(neuroim2::values(roi$train_roi))
  ind <- neuroim2::indices(roi$train_roi)

  n_obs <- nrow(X)
  n_vox <- ncol(X)

  if (n_obs < 2L || n_vox < 1L) {
    return(tibble::tibble(
      result          = list(NULL),
      indices         = list(ind),
      performance     = list(NULL),
      id              = rnum,
      error           = TRUE,
      error_message   = "repnet_model: insufficient observations or voxels.",
      warning         = TRUE,
      warning_message = "repnet_model: insufficient observations or voxels."
    ))
  }

  key <- mod_spec$key
  if (length(key) != n_obs) {
    return(tibble::tibble(
      result          = list(NULL),
      indices         = list(ind),
      performance     = list(NULL),
      id              = rnum,
      error           = TRUE,
      error_message   = "repnet_model: key length mismatch with ROI data.",
      warning         = TRUE,
      warning_message = "repnet_model: key length mismatch with ROI data."
    ))
  }

  # Item-level prototypes
  E_full <- group_means(X, margin = 1, group = key)
  items_roi <- rownames(E_full)

  # Align with seed RDM
  seed_mat <- mod_spec$seed_rdm
  items_seed <- rownames(seed_mat)
  common_items <- intersect(items_roi, items_seed)
  K <- length(common_items)
  futile.logger::flog.debug(
    sprintf("repnet_model: items(roi)=%d, items(seed)=%d, common=%d",
            length(items_roi), length(items_seed), K)
  )

  if (K < 3L) {
    return(tibble::tibble(
      result          = list(NULL),
      indices         = list(ind),
      performance     = list(NULL),
      id              = rnum,
      error           = TRUE,
      error_message   = "repnet_model: need at least 3 common items to define RDM connectivity.",
      warning         = TRUE,
      warning_message = "repnet_model: too few common items."
    ))
  }

  common_items <- sort(common_items)
  E <- E_full[common_items, , drop = FALSE]
  seed_sub <- seed_mat[common_items, common_items, drop = FALSE]

  # Drop constant voxels
  sd_E <- apply(E, 2L, stats::sd, na.rm = TRUE)
  keep_vox <- sd_E > 0
  if (!any(keep_vox)) {
    return(tibble::tibble(
      result          = list(NULL),
      indices         = list(ind),
      performance     = list(NULL),
      id              = rnum,
      error           = TRUE,
      error_message   = "repnet_model: all voxels have zero variance in this ROI.",
      warning         = TRUE,
      warning_message = "repnet_model: all voxels have zero variance."
    ))
  }
  E <- E[, keep_vox, drop = FALSE]

  # ROI RDM and seed vector
  D_roi <- pairwise_dist(mod_spec$distfun, E)
  d_roi  <- as.numeric(D_roi[lower.tri(D_roi)])
  d_seed <- as.numeric(seed_sub[lower.tri(seed_sub)])

  # Raw connectivity (simple RSA-style similarity)
  conn_raw <- suppressWarnings(
    stats::cor(d_roi, d_seed, method = mod_spec$simfun, use = "complete.obs")
  )
  conn_na <- is.na(conn_raw)

  # Build regression predictors (seed + confounds)
  mm <- list(seed = d_seed)
  if (!is.null(mod_spec$confound_rdms) && length(mod_spec$confound_rdms) > 0L) {
    for (nm in names(mod_spec$confound_rdms)) {
      M <- mod_spec$confound_rdms[[nm]]
      # Require labeled RDMs that include all common_items to avoid silent misalignment
      if (is.null(rownames(M)) || is.null(colnames(M)) ||
          !all(common_items %in% rownames(M)) || !all(common_items %in% colnames(M))) {
        msg <- sprintf("repnet_model: confound RDM '%s' lacks required item labels for current ROI.", nm)
        return(tibble::tibble(
          result          = list(NULL),
          indices         = list(ind),
          performance     = list(NULL),
          id              = rnum,
          error           = TRUE,
          error_message   = msg,
          warning         = TRUE,
          warning_message = msg
        ))
      }
      M_sub <- M[common_items, common_items, drop = FALSE]
      mm[[nm]] <- as.numeric(M_sub[lower.tri(M_sub)])
    }
  }

  beta_vec <- numeric(0)
  sp_vec   <- numeric(0)

  # Multiple regression + semipartials (if available)
  df <- as.data.frame(mm, stringsAsFactors = FALSE)
  fit <- try(stats::lm(d_roi ~ ., data = df), silent = TRUE)
  if (!inherits(fit, "try-error")) {
    cf <- stats::coef(summary(fit))
    if (nrow(cf) > 1L) {
      keep <- rownames(cf) != "(Intercept)"
      beta <- cf[keep, "Estimate"]
      names(beta) <- rownames(cf)[keep]

      # If confounds are present, compute a partial (residual-residual) effect for seed
      if ("seed" %in% names(mm) && length(mm) > 1L) {
        conf_df <- as.data.frame(mm[names(mm) != "seed"], stringsAsFactors = FALSE)
        resid_seed <- try(stats::lm(mm$seed ~ ., data = conf_df)$residuals, silent = TRUE)
        if (!inherits(resid_seed, "try-error")) {
          pc <- suppressWarnings(stats::cor(resid_seed, d_roi,
                                            method = mod_spec$simfun,
                                            use = "complete.obs"))
          if (is.finite(pc)) {
            beta["seed"] <- pc
          }
        }
      }

      beta_vec <- beta
    }

    # Semi-partial correlations via run_lm_semipartial (if available)
    if (exists("run_lm_semipartial", mode = "function")) {
      tmp_obj <- list(design = list(model_mat = as.data.frame(mm, stringsAsFactors = FALSE)))
      sp <- try(run_lm_semipartial(dvec = d_roi, obj = tmp_obj), silent = TRUE)
      if (!inherits(sp, "try-error") && is.numeric(sp)) {
        names(sp) <- paste0("sp_", names(sp))
        sp_vec <- sp
      }
    }
  }

  perf <- c(
    n_items  = K,
    n_pairs  = length(d_roi),
    conn_raw = conn_raw
  )
  if (length(beta_vec)) {
    names(beta_vec) <- paste0("beta_", sub("^beta_", "", names(beta_vec)))
    perf <- c(perf, beta_vec)
  }
  if (length(sp_vec)) {
    perf <- c(perf, sp_vec)
  }

  # Flag small-K settings and NA connectivity
  warn_flag <- (K < 5L) || conn_na
  warn_msg <- if (conn_na) {
    "repnet_model: conn_raw is NA; check for zero-variance or all-NA distances."
  } else if (K < 5L) {
    "repnet_model: small item count (K < 5); connectivity estimates may be unstable."
  } else "~"

  tibble::tibble(
    result          = list(NULL),
    indices         = list(ind),
    performance     = list(perf),
    id              = rnum,
    error           = FALSE,
    error_message   = "~",
    warning         = warn_flag,
    warning_message = warn_msg
  )
}
