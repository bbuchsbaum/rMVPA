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
#' @examples
#' ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 12, blocks = 3, nlevels = 2)
#' items <- as.character(sort(unique(ds$design$train_design$.rownum)))
#' # Predictor (X) and outcome (Y) RDMs over items
#' X <- as.matrix(dist(matrix(rnorm(length(items) * 2), ncol = 2)))
#' Y <- as.matrix(dist(matrix(rnorm(length(items) * 2), ncol = 2)))
#' rownames(X) <- colnames(X) <- rownames(Y) <- colnames(Y) <- items
#' repmed_des <- repmed_design(items = items, X_rdm = X, Y_rdm = Y)
#' model <- repmed_model(ds$dataset, ds$design, repmed_des, key_var = ~ .rownum)
#' res <- run_regional(model, ds$dataset$mask)
#' res$performance_table
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
    .fixed_design = .repmed_fixed_design(key_fac, repmed_des$items, repmed_des$X_rdm,
                                         repmed_des$Y_rdm, repmed_des$confound_rdms),
    ...
  )
}


# Confound RDMs as columns over the lower triangle of the common items.
#' @keywords internal
#' @noRd
.repmed_conf_matrix <- function(confound_rdms, common_items) {
  if (is.null(confound_rdms) || length(confound_rdms) == 0L) return(NULL)
  conf_vecs <- lapply(confound_rdms, function(Cr) {
    Csub <- Cr[common_items, common_items, drop = FALSE]
    as.numeric(Csub[lower.tri(Csub)])
  })
  conf_mat <- do.call(cbind, conf_vecs)
  colnames(conf_mat) <- names(confound_rdms)
  conf_mat
}

# Path a and path c share the design [1, x, confounds] and neither depends on
# the ROI (path c also has no ROI term), so both are fitted once at
# construction. fit_roi() reuses them only when the current inputs are identical.
# NULL when the design cannot be built; fit_roi() then reports as before.
#' @keywords internal
#' @noRd
.repmed_fixed_design <- function(key, items, X_rdm, Y_rdm, confound_rdms) {
  tryCatch({
    if (is.null(key)) return(NULL)
    items_key <- rownames(group_means(matrix(0, nrow = length(key), ncol = 1L),
                                      margin = 1, group = key))
    common_items <- sort(Reduce(intersect, list(items, items_key,
                                                rownames(X_rdm), rownames(Y_rdm))))
    if (length(common_items) < 3L) return(NULL)
    Xr <- X_rdm[common_items, common_items, drop = FALSE]
    Yr <- Y_rdm[common_items, common_items, drop = FALSE]
    x <- as.numeric(Xr[lower.tri(Xr)])
    y <- as.numeric(Yr[lower.tri(Yr)])
    P_xc <- cbind(x = x, .repmed_conf_matrix(confound_rdms, common_items))
    X_xc <- .rep_X(P_xc)
    cache_a <- list(X = X_xc, decomp = if (all(is.finite(X_xc))) qr(X_xc) else NULL)
    fit_c <- .rep_path_fit(y, P_xc, cache_a)
    c_tot <- if (!is.null(fit_c) && "x" %in% names(fit_c$coef)) unname(fit_c$coef["x"]) else NA_real_
    list(a = cache_a, c = list(X = X_xc, y = y, c_tot = c_tot))
  }, error = function(e) NULL)
}


#' fit_roi method for repmed_model
#'
#' Computes mediation paths a, b, c' and the indirect effect in RDM space for
#' a single ROI/searchlight.
#'
#' @rdname fit_roi
#' @method fit_roi repmed_model
#' @export
fit_roi.repmed_model <- function(model, roi_data, context, ...) {

  Xmat  <- roi_data$train_data
  ind   <- roi_data$indices
  n_obs <- nrow(Xmat)
  n_vox <- ncol(Xmat)

  if (n_obs < 2L || n_vox < 1L) {
    return(roi_result(
      metrics       = NULL,
      indices       = ind,
      id            = context$id,
      error         = TRUE,
      error_message = "repmed_model: insufficient observations or voxels."
    ))
  }

  key   <- model$key
  items <- model$items

  if (length(key) != n_obs) {
    return(roi_result(
      metrics       = NULL,
      indices       = ind,
      id            = context$id,
      error         = TRUE,
      error_message = "repmed_model: key length mismatch with ROI data."
    ))
  }

  # Item-level mediator prototypes
  M_full <- group_means(Xmat, margin = 1, group = key)
  items_roi <- rownames(M_full)

  # Align items across X, M, Y
  items_X <- rownames(model$X_rdm)
  items_Y <- rownames(model$Y_rdm)
  common_items <- Reduce(intersect, list(items, items_roi, items_X, items_Y))
  K <- length(common_items)
  .log_debug(
    sprintf("repmed_model: items(design)=%d, items(roi)=%d, items(X)=%d, items(Y)=%d, common=%d",
            length(items), length(items_roi), length(items_X), length(items_Y), K)
  )

  if (K < 3L) {
    return(roi_result(
      metrics       = NULL,
      indices       = ind,
      id            = context$id,
      error         = TRUE,
      error_message = "repmed_model: need at least 3 common items for mediation."
    ))
  }

  common_items <- sort(common_items)
  M  <- M_full[common_items, , drop = FALSE]
  Xr <- model$X_rdm[common_items, common_items, drop = FALSE]
  Yr <- model$Y_rdm[common_items, common_items, drop = FALSE]

  # Drop constant voxels
  sd_M <- apply(M, 2L, stats::sd, na.rm = TRUE)
  keep_vox <- sd_M > 0
  if (!any(keep_vox)) {
    return(roi_result(
      metrics       = NULL,
      indices       = ind,
      id            = context$id,
      error         = TRUE,
      error_message = "repmed_model: all voxels have zero variance."
    ))
  }
  M <- M[, keep_vox, drop = FALSE]

  # RDMs in lower-tri vector form
  D_M <- pairwise_dist(model$distfun, M)
  x <- as.numeric(Xr[lower.tri(Xr)])
  y <- as.numeric(Yr[lower.tri(Yr)])
  m <- as.numeric(D_M[lower.tri(D_M)])

  # Confounds: each as vector over lower tri
  conf_mat <- .repmed_conf_matrix(model$confound_rdms, common_items)

  cache <- model$.fixed_design
  # Path a: m ~ x (+ confounds). The design is fixed, so its QR is cached.
  P_xc <- cbind(x = x, conf_mat)
  path_a <- .rep_path_fit(m, P_xc, cache = cache$a)
  a <- if (!is.null(path_a) && "x" %in% names(path_a$coef)) unname(path_a$coef["x"]) else NA_real_

  # Path b / c': y ~ m + x (+ confounds). The design contains the ROI's m.
  P_mxc <- cbind(m = m, x = x, conf_mat)
  path_b <- .rep_path_fit(y, P_mxc)
  b <- if (!is.null(path_b) && "m" %in% names(path_b$coef)) unname(path_b$coef["m"]) else NA_real_
  cprime <- if (!is.null(path_b) && "x" %in% names(path_b$coef)) unname(path_b$coef["x"]) else NA_real_

  # Path c (total effect): y ~ x (+ confounds). Independent of the ROI, so it
  # is reused from construction when the inputs are unchanged.
  X_xc <- .rep_X(P_xc)
  c_tot <- if (!is.null(cache$c) && identical(cache$c$X, X_xc) && identical(cache$c$y, y)) {
    cache$c$c_tot
  } else {
    fit_c <- .rep_path_fit(y, P_xc, cache = cache$a)
    if (!is.null(fit_c) && "x" %in% names(fit_c$coef)) unname(fit_c$coef["x"]) else NA_real_
  }

  indirect <- a * b

  # Naive Sobel test (if enough information)
  sobel_z <- NA_real_
  sobel_p <- NA_real_
  if (is.finite(a) && is.finite(b) &&
      !is.null(path_a) && !is.null(path_b) &&
      "x" %in% names(path_a$se) && "m" %in% names(path_b$se)) {
    sa <- unname(path_a$se["x"])
    sb <- unname(path_b$se["m"])
    se_ind <- sqrt(b^2 * sa^2 + a^2 * sb^2)
    if (is.finite(se_ind) && se_ind > .Machine$double.eps) {
      sobel_z <- indirect / se_ind
      sobel_p <- 2 * stats::pnorm(-abs(sobel_z))
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

  roi_result(
    metrics = perf,
    indices = ind,
    id      = context$id
  )
}


#' @rdname output_schema
#' @method output_schema repmed_model
#' @export
output_schema.repmed_model <- function(model) {
  nms <- c("n_items", "n_pairs", "med_a", "med_b", "med_cprime",
           "med_c", "med_indirect", "med_sobel_z", "med_sobel_p")
  schema <- as.list(rep("scalar", length(nms)))
  names(schema) <- nms
  schema
}
