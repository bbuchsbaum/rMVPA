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
    .roi_error_rows(
      id = res$id %||% NA_integer_,
      error_message = res$error_message,
      indices = list(res$indices),
      warning_message = if (is.null(res$warning_message) ||
                            identical(res$warning_message, "~")) {
        res$error_message
      } else {
        res$warning_message
      },
      result = list(res$result)
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

#' Per-ROI error rows
#'
#' The one constructor for failed-ROI rows, used by both iteration backends,
#' the cross-validation paths and the size filter, so every failure carries the
#' same columns: result, indices, performance, id, error, error_message,
#' error_trace, warning, warning_message. Vectorised over `id`.
#'
#' `error_message` must not contain ROI ids or tracebacks, so identical
#' failures group in [summarize_errors()]; a traceback goes in `error_trace`.
#' @param id ROI ids (one row each).
#' @param error_message Error message (recycled).
#' @param indices `NULL` or a list of index vectors, one per row (recycled).
#' @param warning_message Warning message (recycled).
#' @param error_trace Optional traceback text (recycled).
#' @param result `NULL` or a list of per-row results (recycled).
#' @keywords internal
#' @noRd
.roi_error_rows <- function(id, error_message, indices = NULL,
                            warning_message = error_message,
                            error_trace = NA_character_, result = NULL) {
  n <- length(id)
  as_list_col <- function(x) {
    if (is.null(x)) return(rep(list(NULL), n))
    stopifnot(is.list(x), length(x) %in% c(1L, n))
    rep_len(x, n)
  }
  tibble::tibble(
    result = as_list_col(result),
    indices = as_list_col(indices),
    performance = rep(list(NULL), n),
    id = id,
    error = rep(TRUE, n),
    error_message = rep_len(as.character(error_message), n),
    error_trace = rep_len(as.character(error_trace), n),
    warning = rep(TRUE, n),
    warning_message = rep_len(as.character(warning_message), n)
  )
}

#' Truncated traceback text for an error row, or NA when unavailable
#' @keywords internal
#' @noRd
.roi_error_trace <- function(max_chars = 500L) {
  tryCatch({
    raw <- paste(utils::capture.output(rlang::trace_back()), collapse = "\n")
    if (nchar(raw) > max_chars) {
      paste0(substr(raw, 1L, max_chars), "\n... [truncated]")
    } else {
      raw
    }
  }, error = function(...) NA_character_)
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
  failed <- results$error %in% TRUE
  list(
    good = results[!failed, , drop = FALSE],
    bad  = .with_error_trace(results[failed, , drop = FALSE])
  )
}

# Bad-results tables always carry error_trace (after error_message), also
# when empty or built only from rows of producers without a trace.
#' @keywords internal
#' @noRd
.with_error_trace <- function(bad) {
  if ("error_trace" %in% names(bad) || !"error_message" %in% names(bad)) {
    return(bad)
  }
  bad$error_trace <- rep(NA_character_, nrow(bad))
  pos <- match("error_message", names(bad))
  bad[, append(setdiff(names(bad), "error_trace"), "error_trace", after = pos),
      drop = FALSE]
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
    # `%in% TRUE` treats NA flags as not failed. error_message is kept free of
    # ROI ids (the id is a column), so identical failures group together.
    bad_msg <- as.character(results$error_message[which(results$error %in% TRUE)])
    bad_msg[is.na(bad_msg)] <- "<no message>"
    err_tab <- sort(table(bad_msg), decreasing = TRUE)
    for (i in seq_len(min(length(err_tab), 5L))) {
      futile.logger::flog.warn("  - [%d ROIs] %s", err_tab[i], names(err_tab)[i])
    }
  }
  invisible(n_errors)
}
