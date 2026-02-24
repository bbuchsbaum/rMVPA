#' Construct a Standardized ROI Result
#'
#' Creates a uniform return type for per-ROI analysis. Every
#' \code{\link{fit_roi}} method should return an \code{roi_result}.
#'
#' @param metrics Named numeric vector of performance metrics.
#' @param indices Integer vector of voxel/vertex indices for this ROI.
#' @param id Scalar ROI identifier (center voxel ID or region number).
#' @param result Optional detailed result object (e.g.,
#'   \code{classification_result}).
#' @param error Logical; \code{TRUE} if this ROI failed.
#' @param error_message Character error description, or \code{"~"} if no error.
#' @param warning Logical; \code{TRUE} if a warning was raised.
#' @param warning_message Character warning description, or \code{"~"} if none.
#' @return A list of class \code{"roi_result"}.
#'
#' @examples
#' # Successful ROI result
#' res <- roi_result(
#'   metrics = c(accuracy = 0.85, AUC = 0.9),
#'   indices = 1:10,
#'   id = 42
#' )
#' res$metrics
#'
#' # Error ROI result
#' err <- roi_result(
#'   metrics = NULL,
#'   indices = 1:10,
#'   id = 42,
#'   error = TRUE,
#'   error_message = "Too few features"
#' )
#' err$error
#'
#' @export
roi_result <- function(metrics, indices, id,
                       result = NULL,
                       error = FALSE,
                       error_message = "~",
                       warning = FALSE,
                       warning_message = "~") {
  structure(
    list(
      metrics = metrics,
      indices = indices,
      id = id,
      result = result,
      error = error,
      error_message = error_message,
      warning = warning,
      warning_message = warning_message
    ),
    class = "roi_result"
  )
}

#' @export
#' @method print roi_result
print.roi_result <- function(x, ...) {
  if (x$error) {
    cat(sprintf("roi_result [ERROR] id=%s: %s\n", x$id, x$error_message))
  } else {
    cat(sprintf("roi_result id=%s, %d indices, %d metrics\n",
                x$id, length(x$indices), length(x$metrics)))
    if (length(x$metrics) > 0) {
      cat("  metrics:", paste(names(x$metrics), round(x$metrics, 4),
                              sep = "=", collapse = ", "), "\n")
    }
  }
  invisible(x)
}

#' Convert an roi_result to the tibble format expected by combiners
#'
#' Translates the structured \code{roi_result} into the single-row tibble
#' format that existing combiner functions (\code{combine_standard},
#' \code{combine_rsa_standard}, etc.) expect from \code{process_roi}.
#'
#' @param res An \code{roi_result} object.
#' @return A single-row tibble with columns \code{result}, \code{indices},
#'   \code{performance}, \code{id}, \code{error}, \code{error_message}.
#' @keywords internal
#' @noRd
roi_result_to_tibble <- function(res) {
  stopifnot(inherits(res, "roi_result"))
  if (res$error) {
    tibble::tibble(
      result = list(res$result),
      indices = list(res$indices),
      performance = list(NULL),
      id = res$id,
      error = TRUE,
      error_message = res$error_message,
      warning = res$warning %||% TRUE,
      warning_message = res$warning_message %||% res$error_message
    )
  } else {
    tibble::tibble(
      result = list(res$result),
      indices = list(res$indices),
      performance = list(res$metrics),
      id = res$id,
      error = FALSE,
      error_message = "~",
      warning = res$warning %||% FALSE,
      warning_message = res$warning_message %||% "~"
    )
  }
}


#' Split an iteration results tibble into good and bad rows
#'
#' @param results A data frame (typically from \code{mvpa_iterate}) with an
#'   \code{error} column.
#' @return A list with \code{good} and \code{bad} data frames.
#' @keywords internal
#' @noRd
split_results <- function(results) {
  if (!is.data.frame(results) || !"error" %in% names(results)) {
    empty <- if (is.data.frame(results)) results[0, , drop = FALSE] else tibble::tibble()
    return(list(good = results, bad = empty))
  }
  list(
    good = results[!results$error, , drop = FALSE],
    bad  = results[results$error, , drop = FALSE]
  )
}


#' Log a summary of errors from an iteration results tibble
#'
#' Replaces the many copy-pasted error-reporting blocks throughout the
#' codebase with a single shared helper.
#'
#' @param results A data frame with \code{error} and optionally
#'   \code{error_message} columns.
#' @param context A short label used in log messages (e.g. \code{"run_regional"}).
#' @return The number of errors (invisibly).
#' @keywords internal
#' @noRd
summarize_errors <- function(results, context = "mvpa") {
  if (!is.data.frame(results) || !"error" %in% names(results)) {
    return(invisible(0L))
  }
  n_total   <- nrow(results)
  n_errors  <- sum(results$error, na.rm = TRUE)
  n_success <- n_total - n_errors

  futile.logger::flog.info("%s: %d ROIs processed (success=%d, errors=%d)",
                           context, n_total, n_success, n_errors)

  if (n_errors > 0 && "error_message" %in% names(results)) {
    futile.logger::flog.warn("%s: %d of %d ROIs failed (%.1f%%)",
                             context, n_errors, n_total,
                             100 * n_errors / max(n_total, 1L))
    bad <- results[results$error, , drop = FALSE]
    err_tab <- sort(table(bad$error_message), decreasing = TRUE)
    for (i in seq_len(min(length(err_tab), 5L))) {
      futile.logger::flog.warn("  - [%d ROIs] %s", err_tab[i], names(err_tab)[i])
    }
  }
  invisible(n_errors)
}
