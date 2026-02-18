# Representational Mapping (ReNA-Map): repmap_model

#' Representational mapping model (ReNA-Map)
#'
#' For each ROI/searchlight, fits a reduced-rank regression from seed feature
#' vectors (X) to ROI item-level patterns (Y), summarizing mapping rank and fit.
#'
#' @param dataset An \code{mvpa_dataset}.
#' @param design An \code{mvpa_design}.
#' @param repmap_des Output of \code{\link{repmap_design}}.
#' @param key_var Column or formula giving item identity (e.g. \code{~ ImageID}).
#' @param rank \code{"auto"} for \code{rrpack::cv.rrr} rank selection, integer for
#'   fixed rank, or \code{0} for no mapping (zero map; not an identity transform).
#' @param max_rank Maximum rank to search.
#' @param ridge_lambda Optional ridge penalty lambda for \code{rrpack::rrs.fit}.
#' @param ... Extra fields stored on the model spec.
#'
#' @return A model spec of class \code{"repmap_model"}.
#'
#' @details
#' Internally, the item-level seed features (X) and ROI patterns (Y) are column-centered
#' prior to reduced-rank regression. Returned voxelwise R-squared values are in-sample and may be
#' negative when the mapping underperforms the mean model; this can be useful as a diagnostic
#' rather than an error.
#' @examples
#' \dontrun{
#'   # Requires repmap_design with seed features
#'   # model <- repmap_model(dataset, design, repmap_des, key_var=~ImageID)
#' }
#' @export
repmap_model <- function(dataset,
                         design,
                         repmap_des,
                         key_var,
                         rank = "auto",
                         max_rank = 20,
                         ridge_lambda = NULL,
                         ...) {

  assertthat::assert_that(inherits(dataset, "mvpa_dataset"))
  assertthat::assert_that(inherits(design,  "mvpa_design"))

  d <- design$train_design
  key <- parse_variable(key_var, d)
  key_fac <- factor(key)

  create_model_spec(
    "repmap_model",
    dataset       = dataset,
    design        = design,
    key           = key_fac,
    items         = repmap_des$items,
    seed_features = repmap_des$seed_features,
    rank          = rank,
    max_rank      = max_rank,
    ridge_lambda  = ridge_lambda,
    compute_performance = TRUE,
    return_predictions  = FALSE,
    ...
  )
}


#' Internal helper: reduced-rank regression fit for repmap_model
#'
#' @return A list with fitted repmap RRR components.
#' @keywords internal
.repmap_fit_rrr <- function(X, Y,
                            rank = "auto",
                            max_rank = 20,
                            ridge_lambda = NULL) {

  if (!requireNamespace("rrpack", quietly = TRUE)) {
    futile.logger::flog.warn("repmap_model: package 'rrpack' not available; returning zero map (rank=0).")
    return(list(
      C        = matrix(0, ncol(X), ncol(Y)),
      rank     = 0L,
      singvals = numeric(0)
    ))
  }

  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  max_rank <- max(1L, min(max_rank, p, q, n - 1L))

  if (is.numeric(rank) && as.integer(rank) <= 0L) {
    return(list(
      C        = matrix(0, p, q),
      rank     = 0L,
      singvals = numeric(0)
    ))
  }

  fit <- if (!is.null(ridge_lambda)) {
    rrpack::rrs.fit(
      Y      = Y,
      X      = X,
      nrank  = if (identical(rank, "auto")) max_rank else as.integer(rank),
      lambda = ridge_lambda
    )
  } else if (identical(rank, "auto")) {
    rrpack::cv.rrr(
      Y       = Y,
      X       = X,
      nfold   = min(5, max(2, floor(n / 3))),
      maxrank = max_rank
    )
  } else {
    rrpack::rrr.fit(
      Y     = Y,
      X     = X,
      nrank = as.integer(rank)
    )
  }

  C_hat <- tryCatch({
    if (!is.null(fit$coef)) fit$coef else stats::coef(fit)
  }, error = function(e) matrix(0, nrow = p, ncol = q))

  singvals <- tryCatch(as.numeric(fit$Ad), error = function(e) numeric(0))
  r_used   <- tryCatch({
    if (!is.null(fit$rank)) fit$rank else if (length(singvals)) length(singvals) else NA_integer_
  }, error = function(e) NA_integer_)

  list(
    C        = C_hat,
    rank     = r_used,
    singvals = singvals
  )
}


#' @rdname fit_roi
#' @method fit_roi repmap_model
#' @export
fit_roi.repmap_model <- function(model, roi_data, context, ...) {
  Ymat  <- roi_data$train_data
  ind   <- roi_data$indices
  n_obs <- nrow(Ymat)
  n_vox <- ncol(Ymat)

  if (n_obs < 2L || n_vox < 1L) {
    return(roi_result(
      metrics       = NULL,
      indices       = ind,
      id            = context$id,
      error         = TRUE,
      error_message = "repmap_model: insufficient observations or voxels."
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
      error_message = "repmap_model: key length mismatch with ROI data."
    ))
  }

  # Item-level ROI patterns
  Y_full     <- group_means(Ymat, margin = 1, group = key)
  items_roi  <- rownames(Y_full)

  # Align with seed feature items
  items_seed   <- rownames(model$seed_features)
  common_items <- Reduce(intersect, list(items, items_roi, items_seed))
  K <- length(common_items)
  futile.logger::flog.debug(
    sprintf("repmap_model: items(design)=%d, items(roi)=%d, items(seed)=%d, common=%d",
            length(items), length(items_roi), length(items_seed), K)
  )

  if (K < 3L) {
    return(roi_result(
      metrics       = NULL,
      indices       = ind,
      id            = context$id,
      error         = TRUE,
      error_message = "repmap_model: need at least 3 common items to fit mapping."
    ))
  }

  common_items <- sort(common_items)
  X <- model$seed_features[common_items, , drop = FALSE]
  Y <- Y_full[common_items, , drop = FALSE]

  # Drop constant voxels in Y
  sd_Y     <- apply(Y, 2L, stats::sd, na.rm = TRUE)
  keep_vox <- sd_Y > 0
  if (!any(keep_vox)) {
    return(roi_result(
      metrics       = NULL,
      indices       = ind,
      id            = context$id,
      error         = TRUE,
      error_message = "repmap_model: all voxels have zero variance."
    ))
  }
  Y <- Y[, keep_vox, drop = FALSE]

  # Drop constant columns in X
  sd_X      <- apply(X, 2L, stats::sd, na.rm = TRUE)
  keep_feat <- sd_X > 0
  if (!any(keep_feat)) {
    return(roi_result(
      metrics       = NULL,
      indices       = ind,
      id            = context$id,
      error         = TRUE,
      error_message = "repmap_model: all seed features have zero variance."
    ))
  }
  X <- X[, keep_feat, drop = FALSE]

  # Center X, Y (standard RRR assumption)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  Yc <- scale(Y, center = TRUE, scale = FALSE)

  fit <- .repmap_fit_rrr(
    Xc,
    Yc,
    rank         = model$rank,
    max_rank     = model$max_rank,
    ridge_lambda = model$ridge_lambda
  )

  C_hat  <- fit$C
  r_used <- fit$rank
  svals  <- fit$singvals

  # In-sample mapping fit (variance explained) per voxel
  Y_hat  <- Xc %*% C_hat
  resid  <- Yc - Y_hat

  ss_tot <- colSums(Yc^2)
  ss_res <- colSums(resid^2)
  ss_tot[ss_tot < .Machine$double.eps] <- .Machine$double.eps
  r2_vox <- 1 - ss_res / ss_tot

  mapping_r2_mean <- mean(r2_vox, na.rm = TRUE)
  mapping_r2_med  <- stats::median(r2_vox, na.rm = TRUE)
  frob_norm       <- sqrt(sum(C_hat^2))

  perf <- c(
    n_items       = K,
    n_seed_feats  = ncol(X),
    n_vox         = ncol(Y),
    map_rank      = as.numeric(r_used),
    map_sv1       = if (length(svals)) svals[1] else NA_real_,
    map_sv_mean   = if (length(svals)) mean(svals) else NA_real_,
    map_r2_mean   = mapping_r2_mean,
    map_r2_median = mapping_r2_med,
    map_frob_norm = frob_norm
  )

  roi_result(
    metrics = perf,
    indices = ind,
    id      = context$id
  )
}


#' @rdname output_schema
#' @method output_schema repmap_model
#' @export
output_schema.repmap_model <- function(model) {
  nms <- c("n_items", "n_seed_feats", "n_vox", "map_rank", "map_sv1",
           "map_sv_mean", "map_r2_mean", "map_r2_median", "map_frob_norm")
  schema <- as.list(rep("scalar", length(nms)))
  names(schema) <- nms
  schema
}


