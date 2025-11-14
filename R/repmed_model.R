# Representational Mediation (ReNA-RM): repmed_model

#' Representational mediation model (ReNA-RM)
#'
#' For each ROI/searchlight, tests whether the ROI RDM mediates the relationship
#' between a predictor RDM (X) and an outcome RDM (Y), optionally including
#' confound RDMs.
#'
#' @param dataset An \code{mvpa_dataset}.
#' @param design An \code{mvpa_design}.
#' @param repmed_des Output of \code{\link{repmed_design}}.
#' @param key_var Column or formula giving item identity in the mvpa_design
#'   (e.g. \code{~ ImageID}).
#' @param distfun Distance function to compute the ROI RDM (e.g. \code{cordist(method = "pearson")} ).
#' @param ... Extra fields stored on the model spec.
#'
#' @return A model spec of class \code{"repmed_model"}.
#' @export
repmed_model <- function(dataset,
                         design,
                         repmed_des,
                         key_var,
                         distfun = cordist(method = "pearson"),
                         ...) {

  assertthat::assert_that(inherits(dataset, "mvpa_dataset"))
  assertthat::assert_that(inherits(design,  "mvpa_design"))

  d <- design$train_design
  key <- parse_variable(key_var, d)
  key_fac <- factor(key)

  if (is.character(distfun)) {
    distfun <- create_dist(distfun)
  }
  assertthat::assert_that(inherits(distfun, "distfun"))

  create_model_spec(
    "repmed_model",
    dataset       = dataset,
    design        = design,
    key           = key_fac,
    items         = repmed_des$items,
    X_rdm         = repmed_des$X_rdm,
    Y_rdm         = repmed_des$Y_rdm,
    confound_rdms = repmed_des$confound_rdms,
    distfun       = distfun,
    compute_performance = TRUE,
    return_predictions  = FALSE,
    ...
  )
}


#' Per-ROI processing for repmed_model
#'
#' Computes mediation paths a, b, c' and the indirect effect in RDM space for
#' a single ROI/searchlight.
#'
#' @keywords internal
#' @export
process_roi.repmed_model <- function(mod_spec,
                                     roi,
                                     rnum,
                                     center_global_id = NA,
                                     ...) {

  Xmat <- as.matrix(neuroim2::values(roi$train_roi))
  ind  <- neuroim2::indices(roi$train_roi)
  n_obs <- nrow(Xmat)
  n_vox <- ncol(Xmat)

  if (n_obs < 2L || n_vox < 1L) {
    return(tibble::tibble(
      result          = list(NULL),
      indices         = list(ind),
      performance     = list(NULL),
      id              = rnum,
      error           = TRUE,
      error_message   = "repmed_model: insufficient observations or voxels.",
      warning         = TRUE,
      warning_message = "repmed_model: insufficient observations or voxels."
    ))
  }

  key   <- mod_spec$key
  items <- mod_spec$items

  if (length(key) != n_obs) {
    return(tibble::tibble(
      result          = list(NULL),
      indices         = list(ind),
      performance     = list(NULL),
      id              = rnum,
      error           = TRUE,
      error_message   = "repmed_model: key length mismatch with ROI data.",
      warning         = TRUE,
      warning_message = "repmed_model: key length mismatch with ROI data."
    ))
  }

  # Item-level mediator prototypes
  M_full <- group_means(Xmat, margin = 1, group = key)
  items_roi <- rownames(M_full)

  # Align items across X, M, Y
  items_X <- rownames(mod_spec$X_rdm)
  items_Y <- rownames(mod_spec$Y_rdm)
  common_items <- Reduce(intersect, list(items, items_roi, items_X, items_Y))
  K <- length(common_items)
  futile.logger::flog.debug(
    sprintf("repmed_model: items(design)=%d, items(roi)=%d, items(X)=%d, items(Y)=%d, common=%d",
            length(items), length(items_roi), length(items_X), length(items_Y), K)
  )

  if (K < 3L) {
    return(tibble::tibble(
      result          = list(NULL),
      indices         = list(ind),
      performance     = list(NULL),
      id              = rnum,
      error           = TRUE,
      error_message   = "repmed_model: need at least 3 common items for mediation.",
      warning         = TRUE,
      warning_message = "repmed_model: too few common items."
    ))
  }

  common_items <- sort(common_items)
  M  <- M_full[common_items, , drop = FALSE]
  Xr <- mod_spec$X_rdm[common_items, common_items, drop = FALSE]
  Yr <- mod_spec$Y_rdm[common_items, common_items, drop = FALSE]

  # Drop constant voxels
  sd_M <- apply(M, 2L, stats::sd, na.rm = TRUE)
  keep_vox <- sd_M > 0
  if (!any(keep_vox)) {
    return(tibble::tibble(
      result          = list(NULL),
      indices         = list(ind),
      performance     = list(NULL),
      id              = rnum,
      error           = TRUE,
      error_message   = "repmed_model: all voxels have zero variance.",
      warning         = TRUE,
      warning_message = "repmed_model: all voxels have zero variance."
    ))
  }
  M <- M[, keep_vox, drop = FALSE]

  # RDMs in lower-tri vector form
  D_M <- pairwise_dist(mod_spec$distfun, M)
  x <- as.numeric(Xr[lower.tri(Xr)])
  y <- as.numeric(Yr[lower.tri(Yr)])
  m <- as.numeric(D_M[lower.tri(D_M)])

  # Confounds: each as vector over lower tri
  conf_mat <- NULL
  if (!is.null(mod_spec$confound_rdms) && length(mod_spec$confound_rdms) > 0L) {
    conf_vecs <- lapply(mod_spec$confound_rdms, function(Cr) {
      Csub <- Cr[common_items, common_items, drop = FALSE]
      as.numeric(Csub[lower.tri(Csub)])
    })
    conf_mat <- do.call(cbind, conf_vecs)
    colnames(conf_mat) <- names(mod_spec$confound_rdms)
  }

  # Build data frames for mediation paths
  df_a <- if (is.null(conf_mat)) {
    data.frame(m = m, x = x)
  } else {
    cbind(data.frame(m = m, x = x), as.data.frame(conf_mat))
  }

  df_b <- if (is.null(conf_mat)) {
    data.frame(y = y, m = m, x = x)
  } else {
    cbind(data.frame(y = y, m = m, x = x), as.data.frame(conf_mat))
  }

  # Path c (total effect): y ~ x (+ confounds)
  df_c <- if (is.null(conf_mat)) {
    data.frame(y = y, x = x)
  } else {
    cbind(data.frame(y = y, x = x), as.data.frame(conf_mat))
  }

  # Path a: m ~ x (+ confounds)
  a <- NA_real_
  a_fit <- try(stats::lm(m ~ ., data = df_a), silent = TRUE)
  if (!inherits(a_fit, "try-error")) {
    cf_a <- stats::coef(summary(a_fit))
    if ("x" %in% rownames(cf_a)) {
      a <- cf_a["x", "Estimate"]
    }
  }

  # Path b / c': y ~ m + x (+ confounds)
  b <- NA_real_
  cprime <- NA_real_
  c_tot <- NA_real_
  b_fit <- try(stats::lm(y ~ ., data = df_b), silent = TRUE)
  if (!inherits(b_fit, "try-error")) {
    cf_b <- stats::coef(summary(b_fit))
    if ("m" %in% rownames(cf_b)) {
      b <- cf_b["m", "Estimate"]
    }
    if ("x" %in% rownames(cf_b)) {
      cprime <- cf_b["x", "Estimate"]
    }
  }

  c_fit <- try(stats::lm(y ~ ., data = df_c), silent = TRUE)
  if (!inherits(c_fit, "try-error")) {
    cf_c <- stats::coef(summary(c_fit))
    if ("x" %in% rownames(cf_c)) {
      c_tot <- cf_c["x", "Estimate"]
    }
  }

  indirect <- a * b

  # Naive Sobel test (if enough information)
  sobel_z <- NA_real_
  sobel_p <- NA_real_
  if (!inherits(a_fit, "try-error") &&
      !inherits(b_fit, "try-error") &&
      is.finite(a) && is.finite(b)) {

    cf_a <- stats::coef(summary(a_fit))
    cf_b <- stats::coef(summary(b_fit))

    if ("x" %in% rownames(cf_a) && "m" %in% rownames(cf_b)) {
      sa <- cf_a["x", "Std. Error"]
      sb <- cf_b["m", "Std. Error"]
      se_ind <- sqrt(b^2 * sa^2 + a^2 * sb^2)
      if (is.finite(se_ind) && se_ind > .Machine$double.eps) {
        sobel_z <- indirect / se_ind
        sobel_p <- 2 * stats::pnorm(-abs(sobel_z))
      }
    }
  }

  perf <- c(
    n_items      = K,
    n_pairs      = length(x),
    med_a        = a,
    med_b        = b,
    med_cprime   = cprime,
    med_c        = c_tot,
    med_indirect = indirect,
    med_sobel_z  = sobel_z,
    med_sobel_p  = sobel_p
  )

  # Flag small-K settings as potentially unstable
  warn_flag <- K < 5L
  warn_msg  <- if (warn_flag) "repmed_model: small item count (K < 5); mediation estimates may be unstable." else "~"

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
