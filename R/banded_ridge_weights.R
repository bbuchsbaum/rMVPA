# Full-data weight matrix for a tuned banded-ridge model ---------------------

#' Refit a banded-ridge model on all rows and return its weight matrix
#'
#' Cross-validated accuracies come from [run_banded_ridge()]; this function
#' answers a different question: *what* each response is tuned to. Every
#' response is refit once on all observations, using the hyperparameters it
#' selected most often across the outer folds, and the coefficients are
#' returned as one features-by-responses matrix. That matrix is the input for
#' contrasts over feature columns and for principal-component summaries of
#' tuning.
#'
#' The weights are descriptive. They are estimated from all rows, so they must
#' not be used to score predictions on those rows.
#'
#' @param model A `banded_ridge_model`.
#' @param result The result of `run_banded_ridge(model)`.
#' @param standardized If `TRUE` (default), coefficients are per standard
#'   deviation of each feature, which makes columns with different units
#'   comparable. If `FALSE`, coefficients are in the units of the feature
#'   matrix.
#' @param target_batch_size Number of responses read and fit at once.
#'
#' @return A list of class `banded_ridge_weights` with
#' \describe{
#'   \item{weights}{Numeric matrix, features in rows and responses in columns.}
#'   \item{set}{Feature-band label of each row.}
#'   \item{voxel_index}{Mask index of each column.}
#'   \item{hyperparameters}{One row per response: the selected candidate, its
#'     `alpha` and band weights, and in how many outer folds it was chosen.}
#'   \item{standardized}{The `standardized` argument.}
#' }
#' If the model has nuisance regressors they are regressed out of features and
#' responses (with an intercept) before the refit.
#'
#' @seealso [banded_ridge_model()], [run_banded_ridge()]
#' @export
banded_ridge_weights <- function(model,
                                 result,
                                 standardized = TRUE,
                                 target_batch_size = model$target_batch_size) {
  if (!inherits(model, "banded_ridge_model")) {
    stop("banded_ridge_weights: model must be created by banded_ridge_model().",
         call. = FALSE)
  }
  selections <- result$hyperparameters
  if (!is.data.frame(selections) ||
      !all(c("voxel_index", "candidate_id") %in% names(selections))) {
    stop("banded_ridge_weights: result must come from run_banded_ridge().",
         call. = FALSE)
  }
  if (!is.logical(standardized) || length(standardized) != 1L ||
      is.na(standardized)) {
    stop("standardized must be TRUE or FALSE.", call. = FALSE)
  }
  active_ids <- model$active_response_ids
  if (!setequal(selections$voxel_index, active_ids)) {
    stop("banded_ridge_weights: result and model describe different responses.",
         call. = FALSE)
  }
  candidates <- model$candidates
  candidate_index <- match(selections$candidate_id, candidates$candidate_id)
  if (anyNA(candidate_index)) {
    stop("banded_ridge_weights: result selects candidates the model does not ",
         "contain.", call. = FALSE)
  }
  # Most frequently selected candidate per response; ties resolve to the
  # earliest candidate in canonical order.
  chosen <- vapply(active_ids, function(id) {
    counts <- table(candidate_index[selections$voxel_index == id])
    as.integer(names(counts)[counts == max(counts)][[1L]])
  }, integer(1))
  support <- vapply(seq_along(active_ids), function(ii) {
    sum(candidate_index[selections$voxel_index == active_ids[[ii]]] ==
          chosen[[ii]])
  }, integer(1))

  X <- model$design$X_train$X
  groups <- as.character(model$design$X_train$set)
  n <- nrow(X)
  rows <- seq_len(n)
  nuisance <- model$nuisance
  if (!is.null(nuisance)) {
    X <- .brt_split_data(X, X[, 1L, drop = FALSE], rows, integer(), nuisance)$X
  }
  batch_size <- .brs_batch_size(target_batch_size, length(active_ids),
                                "target_batch_size")
  weights <- matrix(
    NA_real_, ncol(X), length(active_ids),
    dimnames = list(colnames(X), paste0("voxel_", active_ids))
  )
  field <- if (standardized) "standardized_coefficients" else "coefficients"
  for (batch in .brs_batches(length(active_ids), batch_size)) {
    Y <- .brm_response_matrix(model$dataset, active_ids[batch], n)
    if (!is.null(nuisance)) {
      Y <- .brt_split_data(Y, Y, rows, integer(), nuisance)$Y
    }
    for (cc in unique(chosen[batch])) {
      columns <- which(chosen[batch] == cc)
      args <- .brt_candidate_args(candidates, cc)
      fit <- .banded_ridge_fit_direct(
        X, Y[, columns, drop = FALSE], groups = groups,
        alpha = args$alpha, theta = args$theta
      )
      weights[, batch[columns]] <- fit[[field]]
    }
  }
  if (is.null(rownames(weights))) {
    rownames(weights) <- paste0(groups, "_", stats::ave(
      seq_along(groups), groups, FUN = seq_along
    ))
  }
  theta <- candidates$theta[chosen, , drop = FALSE]
  colnames(theta) <- paste0("theta_", candidates$group_names)
  hyperparameters <- data.frame(
    voxel_index = active_ids,
    candidate_id = candidates$candidate_id[chosen],
    alpha = candidates$alpha[chosen],
    theta,
    outer_folds_selected = support,
    stringsAsFactors = FALSE, row.names = NULL
  )
  structure(
    list(weights = weights, set = groups, voxel_index = active_ids,
         hyperparameters = hyperparameters, standardized = standardized),
    class = c("banded_ridge_weights", "list")
  )
}
