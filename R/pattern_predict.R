# Prediction from a fitted pattern model.
#
# Everything the model needs from an observation x is the r-vector
#   u = A' Psi^{-1} x
# together with the r x r Gram matrix G = A' Psi^{-1} A. Classification,
# decoding, and calibrated scores are all small computations on (u, G).

# Sufficient statistics for new observations: rows of U are u_i'.
.pattern_sufficient <- function(fit, newdata) {
  newdata <- as.matrix(newdata)
  if (ncol(newdata) == fit$p_input) {
    newdata <- newdata[, fit$feature_index, drop = FALSE]
  } else if (ncol(newdata) != length(fit$feature_index)) {
    stop(sprintf("pattern_fit: newdata has %d columns; expected %d (input) or %d (retained).",
                 ncol(newdata), fit$p_input, length(fit$feature_index)), call. = FALSE)
  }
  Xc <- .pattern_x_transform_apply(fit$x_transform, newdata)
  Xc %*% fit$precision_A
}

# Symmetric pseudo-inverse on the numerically supported eigenspace.
.sym_pinv <- function(G, tol = 1e-10) {
  eg <- eigen((G + t(G)) / 2, symmetric = TRUE)
  keep <- eg$values > tol * max(eg$values[1], .Machine$double.eps)
  if (!any(keep)) return(matrix(0, nrow(G), ncol(G)))
  V <- eg$vectors[, keep, drop = FALSE]
  V %*% (t(V) / eg$values[keep])
}

#' Predict from a fitted pattern model
#'
#' @param object A \code{pattern_fit}.
#' @param newdata Numeric matrix of observations (rows) by features (columns):
#'   either all input features or only the retained ones.
#' @param type What to return: \code{"prob"} (class posterior probabilities;
#'   categorical targets), \code{"class"} (the most probable class),
#'   \code{"decode"} (posterior-mean targets on the original scale; for a
#'   categorical fit these are posterior-mean one-hot codes, which are not
#'   probabilities and do not sum to one, so use \code{"prob"} instead),
#'   \code{"scores"} (calibrated component scores \eqn{z = G^{+} u}), or
#'   \code{"encode"} (predicted brain measurements for supplied
#'   \code{targets}; returned over the retained features, with the retained
#'   column positions in attribute \code{"feature_index"}).
#' @param targets Targets for \code{type = "encode"}: a factor or a numeric
#'   vector/matrix on the original scale.
#' @param ... Ignored.
#' @return A matrix (or factor for \code{type = "class"}).
#' @details
#' Classification uses the Gaussian class-conditional model implied by the
#' fit: \eqn{p(c \mid x) \propto \pi_c \exp(m_c' u - m_c' G m_c / 2)} with
#' \eqn{m_c = C' y_w(c)}. Decoding uses the working prior \eqn{y_w \sim N(0, I)}
#' on the whitened targets, giving \eqn{\hat t = \Phi (I + G \Phi)^{-1} u} and
#' \eqn{\hat y_w = C \hat t}. Neither forms \eqn{G^{-1}}; a fit with no retained
#' signal returns the class priors or the target means.
#' @examples
#' ds <- gen_sample_dataset(c(6, 6, 4), 60, nlevels = 3, blocks = 3)
#' spec <- pattern_model(ds$dataset, ds$design, rank = 1, refit = TRUE)
#' fit <- run_global(spec, refit = TRUE)$refit
#' X <- get_feature_matrix(ds$dataset)
#' head(predict(fit, X, type = "prob"))
#' head(predict(fit, X, type = "scores"))
#' @export
predict.pattern_fit <- function(object, newdata, type = c("prob", "class", "decode", "scores", "encode"),
                                targets = NULL, ...) {
  type <- match.arg(type)
  yt <- object$y_transform

  if (type == "encode") {
    if (is.null(targets)) stop("predict.pattern_fit: 'targets' are required for type = 'encode'.", call. = FALSE)
    Yw <- .pattern_targets_apply(yt, targets)
    Xc <- (Yw %*% object$C) %*% t(object$A)
    xt <- object$x_transform
    if (!is.null(xt$sd)) Xc <- sweep(Xc, 2L, xt$sd, "*")
    out <- sweep(Xc, 2L, xt$mu, "+")
    attr(out, "feature_index") <- object$feature_index
    return(out)
  }

  U <- .pattern_sufficient(object, newdata)
  G <- object$G

  if (type == "scores") {
    return(U %*% .sym_pinv(G))
  }

  if (type %in% c("prob", "class")) {
    if (yt$type != "categorical") {
      stop("predict.pattern_fit: type = 'prob'/'class' requires categorical targets.", call. = FALSE)
    }
    M <- yt$class_codes %*% object$C                     # K x r class means in score space
    quad <- rowSums((M %*% G) * M)                       # m_c' G m_c
    S <- U %*% t(M)                                      # n x K
    S <- sweep(S, 2L, 0.5 * quad - log(yt$priors), "-")
    S <- S - apply(S, 1L, max)
    P <- exp(S)
    P <- P / rowSums(P)
    colnames(P) <- yt$levels
    if (type == "prob") return(P)
    return(factor(yt$levels[max.col(P, ties.method = "first")], levels = yt$levels))
  }

  # decode: posterior mean of the whitened targets, back on the original scale
  Phi <- object$target_cov
  r <- ncol(G)
  That <- U %*% t(solve(diag(r) + G %*% Phi)) %*% Phi
  Yw_hat <- That %*% t(object$C)
  .pattern_targets_invert(yt, Yw_hat)
}
