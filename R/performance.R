


#' Calculate the Predicted Class from Probability Matrix
#'
#' This function calculates the predicted class from a matrix of predicted probabilities. The class with the highest probability is selected as the predicted class.
#'
#' @param prob A matrix of predicted probabilities with column names indicating the classes.
#' @return A vector of predicted classes corresponding to the highest probability for each row in the input matrix.
#' @examples
#' prob <- matrix(c(0.2, 0.8,
#'                  0.6, 0.4),
#'                nrow = 2, byrow = TRUE,
#'                dimnames = list(NULL, c("A", "B")))
#' predicted_class(prob)
#' @export
predicted_class <- function(prob) {
  maxid <- max.col(prob, ties.method="random")
  pclass <- colnames(prob)[maxid]
}

#' Calculate Performance Metrics for Regression Result
#'
#' This function calculates performance metrics for a regression result object, including R-squared, Root Mean Squared Error (RMSE), and Spearman correlation.
#'
#' @param x A \code{regression_result} object.
#' @param split_list Optional named list of split index groups for computing metrics on sub-groups in addition to the full result.
#' @param ... extra args (not used).
#' @return A named vector with the calculated performance metrics: R-squared, RMSE, and Spearman correlation.
#' @details
#' The function calculates the following performance metrics for the given regression result object:
#' - R-squared: proportion of variance in the observed data that is predictable from the fitted model.
#' - RMSE: root mean squared error, a measure of the differences between predicted and observed values.
#' - Spearman correlation: a measure of the monotonic relationship between predicted and observed values.
#' @examples
#' \dontrun{
#'   # See performance() generic for examples
#' }
#' @seealso \code{\link{regression_result}}
#' @export
#' @importFrom yardstick rsq_vec rmse_vec
#' @importFrom stats cor
performance.regression_result <- function(x, split_list=NULL,...) {
  metric_fun <- function(res) {
    obs  <- res$observed
    pred <- as.numeric(res$predicted)

    res_rsq <- yardstick::rsq_vec(truth = obs, estimate = pred)
    res_rmse <- yardstick::rmse_vec(truth = obs, estimate = pred)
    res_spearcor <- tryCatch(
      stats::cor(obs, pred, method = "spearman", use = "pairwise.complete.obs"),
      error = function(e) NA_real_
    )

    c(R2 = res_rsq, RMSE = res_rmse, spearcor = res_spearcor)
  }

  apply_metric_with_splits(x, metric_fun, split_groups = split_list)
}


## ----------------------------------------------------------------------
## Split-aware performance helpers
## ----------------------------------------------------------------------

#' @keywords internal
#' @noRd
apply_metric_with_splits <- function(result, metric_fun, split_groups = NULL) {
  stopifnot(is.function(metric_fun))

  base_vals <- metric_fun(result)

  if (is.null(split_groups) || length(split_groups) == 0L) {
    return(base_vals)
  }

  testind <- result$testind

  split_vecs <- lapply(names(split_groups), function(tag) {
    design_idx <- split_groups[[tag]]

    res_idx <- if (!is.null(testind)) {
      which(testind %in% design_idx)
    } else {
      design_idx
    }

    if (!length(res_idx)) {
      vals <- rep(NA_real_, length(base_vals))
      names(vals) <- paste0(names(base_vals), "_", tag)
      return(vals)
    }

    if (inherits(result, "regression_result")) {
      sub_res <- result
      sub_res$observed  <- result$observed[res_idx]
      sub_res$predicted <- result$predicted[res_idx]

      if (!is.null(result$testind)) {
        sub_res$testind <- result$testind[res_idx]
      }

      if (!is.null(result$test_design)) {
        sub_res$test_design <- result$test_design[res_idx, , drop = FALSE]
      }
    } else {
      sub_res <- sub_result(result, res_idx)
    }

    vals <- tryCatch(
      metric_fun(sub_res),
      error = function(e) {
        warning(sprintf("performance metric failed for split '%s': %s", tag, e$message))
        rep(NA_real_, length(base_vals))
      }
    )

    if (is.null(names(vals))) {
      names(vals) <- names(base_vals)
    }
    names(vals) <- paste0(names(vals), "_", tag)
    vals
  })

  c(base_vals, unlist(split_vecs))
}


#' Apply Custom Performance Metric to Prediction Result
#'
#' This function applies a user-supplied performance metric to a prediction
#' result object. The custom function should take a single result object
#' (e.g., a \code{classification_result}) and return a named numeric vector
#' of scalar metrics.
#'
#' @param x The prediction result object.
#' @param custom_fun The function used to compute performance metrics,
#'   i.e., \code{custom_fun(x)}; must return a numeric vector (named if
#'   multiple metrics are returned).
#' @param split_list An optional named list of splitting groups. If provided,
#'   metrics are computed on the full result and separately within each group,
#'   with group-specific metrics suffixed by \code{"_<group>"}.
#' @return A named numeric vector with the calculated custom performance
#'   metric(s).
#' @details
#' When a design includes \code{split_by}, the corresponding split groups
#' (typically \code{mvpa_design$split_groups}) can be passed as
#' \code{split_list}. If the result object contains \code{testind}, indices
#' in \code{split_list} are interpreted in the design space and mapped to
#' rows of the result.
#' @examples
#' cres <- binary_classification_result(
#'   observed  = factor(c("A", "B")),
#'   predicted = factor(c("A", "A")),
#'   probs = matrix(c(0.9, 0.1,
#'                    0.6, 0.4),
#'                  ncol = 2, byrow = TRUE,
#'                  dimnames = list(NULL, c("A", "B")))
#' )
#' acc_fun <- function(x) c(acc = mean(x$observed == x$predicted))
#' custom_performance(cres, acc_fun)
#' @export
custom_performance <- function(x, custom_fun, split_list = NULL) {
  if (!is.function(custom_fun)) {
    stop("custom_fun must be a function taking a single result object.")
  }

  if (is.null(split_list) || length(split_list) == 0L) {
    return(custom_fun(x))
  }

  base_vals <- custom_fun(x)
  testind <- x$testind

  split_vecs <- lapply(names(split_list), function(tag) {
    design_idx <- split_list[[tag]]

    res_idx <- if (!is.null(testind)) {
      which(testind %in% design_idx)
    } else {
      design_idx
    }

    if (!length(res_idx)) {
      vals <- rep(NA_real_, length(base_vals))
      names(vals) <- paste0(names(base_vals), "_", tag)
      return(vals)
    }

    sub_res <- sub_result(x, res_idx)
    vals <- custom_fun(sub_res)

    if (is.null(names(vals))) {
      names(vals) <- names(base_vals)
    }
    names(vals) <- paste0(names(vals), "_", tag)
    vals
  })

  c(base_vals, unlist(split_vecs))
}

#' @rdname merge_results-methods
#' @method merge_results binary_classification_result
#' @export
merge_results.binary_classification_result <- function(obj, ...) {
  rlist <- list(obj, ...)
  probs <- Reduce("+", lapply(rlist, function(result) result$probs))/length(rlist)
  
  mc <- max.col(probs)
  predicted <- levels(obj$observed)[mc]
  binary_classification_result(observed = obj$observed,
                               predicted = predicted,
                               probs = probs,
                               testind = obj$testind,
                               test_design = obj$test_design,
                               predictor = obj$predictor)
}

#' @rdname merge_results-methods
#' @method merge_results regression_result
#' @export
merge_results.regression_result <- function(obj, ...) {
  rlist <- list(obj, ...)
  pred <- Reduce("+", lapply(rlist, function(result) result$predicted))/length(rlist)
  regression_result(observed = obj$observed,
                    predicted = pred,
                    testind = obj$testind,
                    test_design = obj$test_design,
                    predictor = obj$predictor)
}



#' @export
prob_observed.binary_classification_result <- function(x) {
  x$probs[cbind(seq(1,nrow(x$probs)),as.integer(x$observed))]
}

#' @export
prob_observed.multiway_classification_result <- function(x) {
  x$probs[cbind(seq(1,nrow(x$probs)),as.integer(x$observed))]
}

#' @export
prob_observed.regression_result <- function(x) {
  # Regression results don't have probabilities, return NULL
  # This allows the searchlight combiner to skip prob_observed for regression
  NULL
}

#' @rdname merge_results-methods
#' @method merge_results multiway_classification_result
#' @export
merge_results.multiway_classification_result <- function(obj, ...) {
  
  rlist <- list(obj, ...)
  #ds <- sapply(rlist, function(x) nrow(x$probs))
  
  probs <- Reduce("+", lapply(rlist, function(result) result$probs))/length(rlist)
  mc <- max.col(probs)
  predicted <- levels(obj$observed)[mc]
  
  multiway_classification_result(observed = obj$observed,
                                 predicted = predicted,
                                 probs = probs,
                                 testind = obj$testind,
                                 test_design = obj$test_design,
                                 predictor = obj$predictor)
}

#' @export
performance.binary_classification_result <- function(x, split_list=NULL,...) {
  stopifnot(length(x$observed) == length(x$predicted))
  metric_fun <- function(res) {
    binary_perf(res$observed, res$predicted, res$probs)
  }
  apply_metric_with_splits(x, metric_fun, split_list)
}


#' @export
performance.multiway_classification_result <- function(x, split_list=NULL, class_metrics=FALSE,...) {
  stopifnot(length(x$observed) == length(x$predicted))

  metric_fun <- function(res) {
    multiclass_perf(res$observed, res$predicted, res$probs, class_metrics)
  }
  apply_metric_with_splits(x, metric_fun, split_list)
}

#' @keywords internal
#' @noRd
combinedACC <- function(Pred, Obs) {
  levs <- levels(as.factor(Obs))
  maxind <- apply(Pred, 1, which.max)
  pclass <- levs[maxind]
  sum(pclass == Obs)/length(pclass)
  
}


#' @keywords internal
#' @importFrom yardstick accuracy_vec roc_auc_vec
binary_perf <- function(observed, predicted, probs) {
  # Ensure observed is a factor with levels in a consistent order
  lvls <- levels(observed)
  if (length(lvls) != 2) stop("binary_perf expects 2 levels in observed.")

  # Ensure predicted is a factor with the same levels as observed
  predicted_factor <- factor(predicted, levels = lvls)

  # Positive class is the second level; prefer named columns if available
  if (!is.null(colnames(probs)) && all(lvls %in% colnames(probs))) {
    prob_positive_class <- probs[, lvls[2]]
  } else if (ncol(probs) >= 2) {
    prob_positive_class <- probs[, 2]
  } else {
    prob_positive_class <- probs[, 1]
  }

  res_acc <- yardstick::accuracy_vec(truth = observed, estimate = predicted_factor)

  tab <- table(observed)
  if (length(tab) < 2L || any(tab == 0L)) {
    # Yardstick would warn about missing control/event; return NA without calling it.
    res_auc <- NA_real_
  } else {
    res_auc <- suppressWarnings(
      yardstick::roc_auc_vec(
        truth       = observed,
        estimate    = prob_positive_class,
        event_level = "second"
      )
    )
  }

  auc_centered <- if (is.na(res_auc)) NA_real_ else (2 * res_auc - 1)  # chance-centered, range [-1,1]
  c(Accuracy = res_acc, AUC = auc_centered)
}

#' @keywords internal
#' @importFrom yardstick accuracy_vec roc_auc_vec
multiclass_perf <- function(observed, predicted, probs, class_metrics=FALSE) {
  lvls <- levels(observed)
  predicted_factor <- factor(predicted, levels = lvls)

  acc <- yardstick::accuracy_vec(truth = observed, estimate = predicted_factor)
  
  # Calculate per-class AUC using one-vs-rest approach (matching original logic)
  aucres <- sapply(seq_along(lvls), function(i) {
    lev <- lvls[i]
    pos <- observed == lev
    pclass <- probs[,i]
    pother <- rowMeans(probs[,-i, drop=FALSE])
    # Original uses pclass - pother as the score
    score <- pclass - pother
    binary_truth <- factor(ifelse(pos, "positive", "negative"),
                           levels = c("negative", "positive"))

    # Skip AUC when one of the levels is absent to avoid yardstick warnings.
    tab <- table(binary_truth)
    if (length(tab) < 2L || any(tab == 0L)) {
      return(NA_real_)
    }

    suppressWarnings(
      tryCatch(
        yardstick::roc_auc_vec(
          truth       = binary_truth,
          estimate    = score,
          event_level = "second"
        ),
        error = function(e) NA_real_
      )
    )
  })
  
  names(aucres) <- paste0("AUC_", colnames(probs))

  auc_centered <- 2 * aucres - 1
  mean_auc_centered <- mean(auc_centered, na.rm = TRUE)
  
  metrics <- c(Accuracy = acc, AUC = mean_auc_centered)
  
  if (class_metrics) {
    names(auc_centered) <- paste0("AUC_", colnames(probs))
    c(metrics, auc_centered)
  } else {
    metrics
  }
}
  
