#' Fit a REMAP decoder for independent target observations
#'
#' Learns a residual reduced-rank transformation from matched source and
#' calibration keys. Supply only training observations to this function.
#' A fit without calibration estimates its templates and whitening entirely
#' from the source domain.
#'
#' @md
#' @param source Numeric matrix of source observations by features.
#' @param source_key One item or class key per source row. Repeated keys are
#'   averaged. For precomputed templates, use their row names (the default).
#' @param calibration Optional numeric matrix of labeled target training rows.
#'   Its keys must be a subset of the source keys. At least two common keys
#'   are required. All source keys are candidates for prediction, including
#'   keys absent from calibration.
#' @param calibration_key Keys for calibration rows, defaulting to row names.
#' @param rank Residual rank, or "auto" to select rank using training keys.
#'   Defaults to one with calibration and zero without. A positive rank or
#'   automatic rank requires the optional rrpack package.
#' @param max_rank Upper bound for automatic rank selection.
#' @param ridge_rrr_lambda Optional nonnegative ridge penalty. With "auto",
#'   the ridge fit uses the maximum feasible rank up to max_rank.
#' @param lambda_grid Candidate residual multipliers between zero and one. A single value
#'   fixes the multiplier. Otherwise training prototype retrieval chooses it,
#'   with ties favoring the smaller correction.
#' @param center Subtract the training feature means before applying whitening.
#' @param shrink_whiten Use shrinkage covariance for whitening (requires corpcor
#'   for multivariate samples of three or more prototypes). A source-only
#'   fit with two prototypes uses a spherical covariance with their average
#'   feature variance.
#' @return A serializable `remap_fit` with training center, whitener, correction,
#'   rank, multiplier and frozen class templates. Use [predict.remap_fit()]
#'   on fresh observations. Rank zero with calibration still uses calibration
#'   for whitening. Omit calibration for a strictly source-only baseline.
#' @examples
#' set.seed(12)
#' key <- rep(letters[1:4], each = 3)
#' x <- matrix(rnorm(60), 12, 5)
#' new_x <- matrix(rnorm(15), 3, 5)
#' source_fit <- fit_remap(x, key, shrink_whiten = FALSE)
#' predict(source_fit, new_x)
#' if (requireNamespace("rrpack", quietly = TRUE)) {
#'   calibration <- x + matrix(rnorm(60, sd = 0.2), 12, 5)
#'   fit <- fit_remap(x, key, calibration, key, rank = 1, shrink_whiten = FALSE)
#'   predict(fit, new_x, type = "scores")
#' }
#' @export
fit_remap <- function(source, source_key = rownames(source),
                      calibration = NULL, calibration_key = rownames(calibration),
                      rank = if (is.null(calibration)) 0L else 1L,
                      max_rank = 20L, ridge_rrr_lambda = NULL,
                      lambda_grid = c(0, .25, .5, .75, 1),
                      center = TRUE, shrink_whiten = TRUE) {
  .remap_matrix(source, "source")
  .remap_keys(source_key, nrow(source), "source_key")
  prototypes <- .nx_rowsum_mean(source, source_key)
  if (nrow(prototypes) < 2L) stop("At least two source keys are required.")
  if (is.null(calibration)) {
    if (!is.null(calibration_key) || !isTRUE(rank == 0))
      stop("Without calibration, supply no calibration_key and use rank = 0.")
    target <- prototypes[FALSE, , drop = FALSE]
  } else {
    .remap_matrix(calibration, "calibration")
    .remap_features(calibration, ncol(prototypes), colnames(prototypes))
    .remap_keys(calibration_key, nrow(calibration), "calibration_key")
    target <- .nx_rowsum_mean(calibration, calibration_key)
    if (nrow(target) < 2L || !all(rownames(target) %in% rownames(prototypes)))
      stop("Calibration requires at least two keys, all present in source.")
  }
  .remap_fit_templates(prototypes, target, rank, max_rank, ridge_rrr_lambda,
                       lambda_grid, center, shrink_whiten)
}

.remap_matrix <- function(x, name) {
  if (!is.matrix(x) || !is.numeric(x) || !nrow(x) || !ncol(x) ||
      any(!is.finite(x))) stop(name, " must be a nonempty finite numeric matrix.")
  if (!is.null(colnames(x)) && (anyNA(colnames(x)) || any(!nzchar(colnames(x))) ||
                                anyDuplicated(colnames(x))))
    stop(name, " feature names must be unique and nonempty.")
}

.remap_features <- function(x, count, names) {
  if (ncol(x) != count || !identical(colnames(x), names))
    stop("Features must match training columns in number, name and order.")
}

.remap_keys <- function(x, n, name) {
  if (!is.atomic(x) || !is.null(dim(x)) || length(x) != n || anyNA(x) ||
      any(!nzchar(as.character(x)))) stop(name, " must supply one nonmissing key per row.")
}

.remap_fit_templates <- function(source, target, rank, max_rank, ridge_rrr_lambda,
                                 lambda_grid, center, shrink_whiten) {
  if (!is.numeric(lambda_grid) || !length(lambda_grid) || any(!is.finite(lambda_grid)) ||
      any(lambda_grid < 0 | lambda_grid > 1)) stop("lambda_grid must contain values in [0, 1].")
  if (!rlang::is_bool(center) || !rlang::is_bool(shrink_whiten))
    stop("center and shrink_whiten must be TRUE or FALSE.")
  paired <- if (nrow(target)) source[rownames(target), , drop = FALSE] else source
  jw <- .joint_whiten(paired, target, shrink_whiten, center)
  rf <- .fit_residual_rrr(jw$Xw_fit, jw$Yw_fit, rank, max_rank, ridge_rrr_lambda)
  lambda <- if (nrow(target)) .select_lambda(jw$Xw_fit, jw$Yw_fit, rf$Delta, lambda_grid)$lambda else 0
  raw <- sweep(source, 2L, jw$mu, "-") %*% jw$W
  templates <- raw + lambda * (raw %*% rf$Delta)
  structure(list(templates = templates, mu = jw$mu, W = jw$W,
                 Delta = rf$Delta, rank = rf$rank, lambda = lambda,
                 singular_values = rf$Ad, features = colnames(source),
                 calibration_keys = rownames(target)), class = "remap_fit")
}

#' Score new target observations with a frozen REMAP decoder
#'
#' @md
#' @param object A fitted object from [fit_remap()].
#' @param newdata Finite numeric matrix of observations by training features.
#'   Features must have the training order, and names must match when supplied
#'   at training. Each observation is transformed with the stored training
#'   mean and whitening matrix and can be scored on its own.
#' @param type Return an observation-by-key matrix of Pearson scores, their
#'   row-wise softmax ("prob"), or the highest-scoring key ("class").
#'   A constant spatial pattern has neutral score zero. Ties select the
#'   first template. Softmax values are relative scores.
#' @param ... Reserved for S3 compatibility.
#' @return A numeric matrix for scores or probabilities, or a factor for classes.
#' @export
predict.remap_fit <- function(object, newdata, type = c("scores", "class", "prob"), ...) {
  type <- match.arg(type)
  .remap_matrix(newdata, "newdata")
  .remap_features(newdata, nrow(object$W), object$features)
  whitened <- sweep(newdata, 2L, object$mu, "-") %*% object$W
  scores <- .remap_scores(whitened, object$templates)
  switch(type, scores = scores, prob = .nx_softmax(scores),
         class = factor(colnames(scores)[max.col(scores, ties.method = "first")],
                        levels = colnames(scores)))
}

.remap_scores <- function(observations, templates) {
  scores <- .nx_row_cor(observations, templates)
  scores[is.na(scores)] <- 0
  dimnames(scores) <- list(rownames(observations), rownames(templates))
  scores
}
