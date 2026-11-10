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
#'   \item{contrast_weights}{Contrasts-by-responses matrix: the weights
#'     projected onto the contrasts declared in `feature_sets(contrasts = )`,
#'     or `NULL` when none were declared.}
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
      weights[, batch[columns]] <- .brw_refit(
        model$solver, X, Y[, columns, drop = FALSE], groups, args
      )[[field]]
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
  contrasts <- model$design$X_train$contrasts
  structure(
    list(weights = weights, set = groups, voxel_index = active_ids,
         hyperparameters = hyperparameters, standardized = standardized,
         contrast_weights = if (is.null(contrasts)) NULL else {
           crossprod(contrasts, weights)
         }),
    class = c("banded_ridge_weights", "list")
  )
}

# One full-data fit with the model's own solver, so a feature matrix the model
# handles in dual form is never turned into a features-by-features system.
.brw_refit <- function(solver, X, Y, groups, args) {
  switch(solver,
    svd_primal = .banded_ridge_fit_svd_path(
      X, Y, alphas = args$alpha, theta = args$theta, groups = groups
    )$fits[[1L]],
    dual_kernel = .banded_ridge_fit_dual_path(
      X, Y, alphas = args$alpha, theta = args$theta, groups = groups,
      recover_primal = TRUE
    )$fits[[1L]],
    .banded_ridge_fit_direct(X, Y, groups = groups, alpha = args$alpha,
                             theta = args$theta)
  )
}

#' Permutation null for banded-ridge prediction accuracy
#'
#' Tests the held-out accuracy of every response against a null in which the
#' predicted series is shifted in time against the observed one. Shifting keeps
#' the autocorrelation of both series, which a row shuffle would destroy, so the
#' null is appropriate for continuous time series. No model is refit.
#'
#' The null addresses the full model's accuracy. It does not test a band's
#' unique contribution: shifting the predictions of the full and reduced
#' models together leaves their difference near zero by construction.
#'
#' @param model A `banded_ridge_model` built with `return_predictions = TRUE`.
#' @param result The result of `run_banded_ridge(model)`.
#' @param n_perm Number of permutations.
#' @param metric `"correlation"` (default) or `"r2"`.
#' @param strategy `"circular_shift"` (default) or `"block_shuffle"`.
#' @param block_size Block length for `"block_shuffle"`; defaults to
#'   `floor(sqrt(n))`.
#' @param min_shift Smallest circular shift, in rows. Set it beyond the
#'   autocorrelation length of the data.
#' @param seed Integer seed; the same permutations are applied to every
#'   response.
#'
#' @return A data frame with one row per response: `voxel_index`, the observed
#'   `statistic`, the permutation `p_value` (one-sided, observed value counted
#'   in the null), and `z`, the statistic standardized by the null.
#' @seealso [run_banded_ridge()]
#' @export
banded_ridge_null <- function(model,
                              result,
                              n_perm = 1000L,
                              metric = c("correlation", "r2"),
                              strategy = c("circular_shift", "block_shuffle"),
                              block_size = NULL,
                              min_shift = 1L,
                              seed = 1L) {
  metric <- match.arg(metric)
  strategy <- match.arg(strategy)
  if (!inherits(model, "banded_ridge_model")) {
    stop("banded_ridge_null: model must be created by banded_ridge_model().",
         call. = FALSE)
  }
  predicted <- result$predictions
  if (is.null(predicted)) {
    stop("banded_ridge_null: build the model with return_predictions = TRUE.",
         call. = FALSE)
  }
  if (!is.numeric(n_perm) || length(n_perm) != 1L || !is.finite(n_perm) ||
      n_perm < 1L || n_perm != as.integer(n_perm)) {
    stop("n_perm must be one positive integer.", call. = FALSE)
  }
  active_ids <- model$active_response_ids
  n <- nrow(predicted)
  observed <- result$observed
  if (is.null(observed)) {
    observed <- do.call(cbind, lapply(
      .brs_batches(length(active_ids), model$target_batch_size),
      function(batch) .brm_response_matrix(model$dataset, active_ids[batch], n)
    ))
  }
  statistic <- function(pred) {
    scores <- .banded_ridge_score(observed, pred)
    if (metric == "r2") scores$r2 else scores$correlation
  }
  actual <- statistic(predicted)
  shifts <- .brt_with_seed(seed, lapply(seq_len(n_perm), function(i) {
    .br_perm_indices(n, strategy = strategy, block_size = block_size,
                     min_shift = min_shift)
  }))
  exceed <- numeric(length(actual))
  total <- total_sq <- numeric(length(actual))
  for (index in shifts) {
    shifted <- predicted[index, , drop = FALSE]
    rownames(shifted) <- rownames(predicted)
    value <- statistic(shifted)
    value[is.na(value)] <- -Inf
    exceed <- exceed + (value >= actual)
    finite <- is.finite(value)
    total[finite] <- total[finite] + value[finite]
    total_sq[finite] <- total_sq[finite] + value[finite]^2
  }
  null_mean <- total / n_perm
  null_sd <- sqrt(pmax(total_sq / n_perm - null_mean^2, 0))
  data.frame(
    voxel_index = active_ids,
    statistic = actual,
    p_value = (exceed + 1) / (n_perm + 1),
    z = ifelse(null_sd > 0, (actual - null_mean) / null_sd, NA_real_),
    stringsAsFactors = FALSE
  )
}
