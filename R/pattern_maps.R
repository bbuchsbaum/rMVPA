# Interpretation of a pattern fit. Scalar summaries live in feature space;
# component coordinates are a display choice only.

.pattern_base_fit <- function(object) {
  if (inherits(object, "pattern_view")) return(object$fit)
  if (inherits(object, "pattern_global_result")) {
    if (is.null(object$refit)) stop("A refit is required; run_global(..., refit = TRUE).", call. = FALSE)
    return(object$refit)
  }
  if (!inherits(object, "pattern_fit")) stop("Expected a pattern_fit, pattern_view, or pattern_global_result.", call. = FALSE)
  object
}

.pattern_map_values <- function(fit, type) {
  A <- fit$A
  scale <- fit$x_transform$sd %||% rep(1, nrow(A))
  if (type == "forward") return(A * scale)
  if (type == "weights") return((fit$precision_A %*% .sym_pinv(fit$G)) / scale)
  if (type == "signal_sd") return(sqrt(pmax(rowSums((A %*% fit$target_cov) * A), 0)) * scale)
  # K = Phi^(1/2) (I + Phi^(1/2) G Phi^(1/2))^-1 Phi^(1/2)
  # avoids inverting Phi, including when its rank is deficient.
  eg <- eigen(fit$target_cov, symmetric = TRUE)
  B <- sweep(eg$vectors, 2L, sqrt(pmax(eg$values, 0)), "*")
  K <- B %*% solve(diag(ncol(B)) + crossprod(B, fit$G %*% B), t(B))
  h <- fit$precision_A / sqrt(.noise_precision_diag(fit$noise))
  leverage <- rowSums((h %*% K) * h)
  if (any(leverage < -1e-8 | leverage > 1 + 1e-8)) stop("Invalid conditional-information covariance.", call. = FALSE)
  -0.5 * log1p(-pmin(pmax(leverage, 0), 1))
}

#' Extract forward patterns, calibrated weights, or invariant maps
#'
#' @param fit A \code{pattern_fit}, \code{pattern_view}, or global result with a refit.
#' @param type Quantity to extract.
#' @param dataset Optional dataset for \code{build_output_map}. Global results
#'   supply their dataset automatically. Multibasis image maps are refused
#'   because aggregation across channels changes the estimand; extract vectors.
#' @param ... Reserved.
#' @return A vector or matrix aligned to all input columns (screened columns are
#'   \code{NA}), or an image / list of component images when a dataset is supplied.
#' @details Forward patterns and signal standard deviations are in original
#'   feature units. Weights map original, centred measurements to calibrated
#'   scores. Forward patterns and weights depend on the component basis.
#'   Signal SD and conditional information are coordinate invariant.
#'   Conditional information is in nats under the working Gaussian model for
#'   task scores, not empirical information about class labels. A feature with
#'   zero loading can have positive information through noise cancellation.
#' @export
model_patterns <- function(fit, type = c("forward", "weights", "conditional_info", "signal_sd"),
                           dataset = NULL, ...) {
  type <- match.arg(type)
  original <- fit
  fit <- .pattern_base_fit(fit)
  values <- .pattern_map_values(fit, type)
  if (inherits(original, "pattern_view") && type %in% c("forward", "weights")) {
    values <- values %*% original$Q_b
  }
  if (is.matrix(values)) {
    out <- matrix(NA_real_, fit$p_input, ncol(values))
    out[fit$feature_index, ] <- values
  } else {
    out <- rep(NA_real_, fit$p_input)
    out[fit$feature_index] <- values
  }
  if (is.null(dataset) && inherits(original, "pattern_global_result")) dataset <- original$model_spec$dataset
  if (is.null(dataset)) return(out)
  if (inherits(dataset, "mvpa_multibasis_image_dataset")) {
    stop("Multibasis patterns require separate channel maps; extract vectors from the fit.", call. = FALSE)
  }
  ids <- feature_ids_for_dataset(dataset, fit$p_input)
  if (length(ids) != fit$p_input) stop("Dataset feature count does not match fit.", call. = FALSE)
  if (is.matrix(out)) return(lapply(seq_len(ncol(out)), function(j) build_output_map(dataset, out[, j], ids)))
  build_output_map(dataset, out, ids)
}

#' @rdname model_patterns
#' @param object A fit, view, or global result.
#' @param X_train Unused; maps are implied by the fitted model.
#' @export
model_importance.pattern_fit <- function(object, X_train = NULL,
                                        type = c("signal_sd", "conditional_info"), ...) {
  model_patterns(object, type = match.arg(type), ...)
}

#' @rdname model_patterns
#' @export
model_importance.pattern_view <- model_importance.pattern_fit

#' @rdname model_patterns
#' @export
model_importance.pattern_global_result <- model_importance.pattern_fit

.pattern_rotation <- function(L, rotation) {
  r <- ncol(L)
  if (is.character(rotation) && length(rotation) == 1L) {
    rotation <- match.arg(rotation, c("varimax", "none"))
    if (rotation == "none" || r == 1L || !any(L != 0)) return(diag(r))
    return(unclass(stats::varimax(L, normalize = FALSE)$rotmat))
  }
  if (!is.matrix(rotation) || !is.numeric(rotation) ||
      !identical(dim(rotation), c(r, r)) || any(!is.finite(rotation)) ||
      max(abs(crossprod(rotation) - diag(r))) > 1e-10) {
    stop("Rotation must be an orthogonal rank by rank matrix, 'none', or 'varimax'.", call. = FALSE)
  }
  rotation
}

#' Display patterns in orthogonal spatial and target coordinates
#'
#' @param fit A fitted pattern model or global result with a refit.
#' @param spatial,target Orthogonal matrices, \code{"varimax"}, or \code{"none"}.
#' @return A \code{pattern_view} containing \code{L_b = A Q_b},
#'   \eqn{H = Q_b^T Q_t}, and \code{L_t = C Q_t}, in the fitted feature and
#'   whitened target coordinates. Thus \eqn{L_b H L_t^T} reconstructs
#'   \eqn{A C^T}. Back-transforms, a basis identifier, and the Frobenius
#'   reconstruction error (evaluated without a full feature by target matrix)
#'   are retained.
#' @details The underlying fit is unchanged. Predictions (including calibrated
#'   scores) and invariant maps delegate to that fit and are bit-identical.
#'   Display columns are not separately identified scientific components.
#' @export
rotate_patterns <- function(fit, spatial = "varimax", target = "varimax") {
  fit <- .pattern_base_fit(fit)
  Qb <- .pattern_rotation(fit$A, spatial)
  Qt <- .pattern_rotation(fit$C, target)
  Lb <- fit$A %*% Qb; Lt <- fit$C %*% Qt; H <- crossprod(Qb, Qt)
  # C has orthonormal columns, so this is the Frobenius reconstruction
  # error without materializing the potentially huge features x targets map.
  error <- sqrt(sum((fit$A %*% (Qb %*% H %*% t(Qt) - diag(ncol(fit$A))))^2))
  if (error > 1e-8 * max(1, sqrt(sum(fit$A^2)))) stop("Rotation reconstruction failed.", call. = FALSE)
  structure(list(fit = fit, L_b = Lb, H = H, L_t = Lt, Q_b = Qb, Q_t = Qt,
                 basis_id = digest::digest(list(fit$C, fit$y_transform, Qb, Qt)),
                 back_transforms = list(spatial = t(Qb), target = t(Qt),
                                        x = fit$x_transform, y = fit$y_transform),
                 reconstruction_error = error), class = "pattern_view")
}

#' @rdname rotate_patterns
#' @param object A pattern view.
#' @param ... Arguments passed to \code{predict.pattern_fit}.
#' @export
predict.pattern_view <- function(object, ...) predict(object$fit, ...)

#' Compare empirical held-out Haufe patterns with the forward model
#'
#' @param fit A pattern fit or view.
#' @param X_holdout Independent rows, with all original input columns.
#' @param observation_ids Held-out row identifiers in the same namespace as
#'   \code{fit$training_observation_ids}. Required when the fit carries those IDs.
#' @return Empirical and model patterns on retained features in original units,
#'   relative Frobenius discrepancy, effective score rank, and row identifiers.
#' @details Independence is a caller obligation; an ID overlap is rejected.
#'   For a rank-deficient decoder, the comparison uses the identifiable score
#'   subspace \code{A G+ G}. This is a covariance diagnostic, not an inference test.
#' @export
pattern_haufe <- function(fit, X_holdout, observation_ids = NULL) {
  fit <- .pattern_base_fit(fit)
  X <- as.matrix(X_holdout)
  if (!is.numeric(X) || ncol(X) != fit$p_input || nrow(X) < 2L) stop("Supply at least two held-out rows with all input columns.", call. = FALSE)
  if (!is.null(fit$training_observation_ids) && is.null(observation_ids)) stop("Supply held-out observation_ids to check training overlap.", call. = FALSE)
  if (!is.null(observation_ids)) {
    if (length(observation_ids) != nrow(X) || anyNA(observation_ids) || anyDuplicated(observation_ids)) stop("Invalid held-out observation_ids.", call. = FALSE)
    if (any(observation_ids %in% fit$training_observation_ids)) stop("Held-out observation_ids overlap training rows.", call. = FALSE)
  }
  X <- X[, fit$feature_index, drop = FALSE]
  if (any(!is.finite(X))) stop("Held-out retained features must be finite.", call. = FALSE)
  W <- .pattern_map_values(fit, "weights")
  empirical <- haufe_importance(W = W, X = X)$A
  model <- .pattern_map_values(fit, "forward") %*% .sym_pinv(fit$G) %*% fit$G
  denom <- sqrt(sum(model^2))
  gv <- eigen(fit$G, symmetric = TRUE, only.values = TRUE)$values
  list(empirical = empirical, model = model,
       relative_discrepancy = if (denom > 0) sqrt(sum((empirical - model)^2)) / denom else NA_real_,
       score_rank = sum(gv > 1e-10 * max(gv[1], .Machine$double.eps)),
       feature_index = fit$feature_index, observation_ids = observation_ids)
}
