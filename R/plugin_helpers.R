#' Build Mock ROI Data for Plugin Tests
#'
#' Constructs a minimal ROI payload compatible with \code{\link{fit_roi}}.
#' This is useful for unit-testing plugin methods without running
#' \code{run_searchlight()} or \code{run_regional()}.
#'
#' @param train_data Optional matrix of training data (observations x features).
#' @param test_data Optional matrix of test data (observations x features).
#' @param indices Optional integer feature indices. Defaults to
#'   \code{seq_len(ncol(train_data))}.
#' @param n_train Number of training observations to simulate when
#'   \code{train_data} is \code{NULL}.
#' @param n_features Number of features to simulate when \code{train_data}
#'   is \code{NULL}.
#' @param n_test Number of test observations to simulate when \code{test_data}
#'   is \code{NULL} and \code{n_test > 0}.
#' @param seed Optional RNG seed used only when simulating data.
#' @param train_roi Optional ROI object to include in the returned list.
#' @param test_roi Optional ROI object to include in the returned list.
#'
#' @return A list with fields \code{train_data}, \code{test_data},
#'   \code{indices}, \code{train_roi}, and \code{test_roi}.
#' @export
mock_roi_data <- function(train_data = NULL,
                          test_data = NULL,
                          indices = NULL,
                          n_train = 20L,
                          n_features = 10L,
                          n_test = 0L,
                          seed = NULL,
                          train_roi = NULL,
                          test_roi = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (is.null(train_data)) {
    n_train <- as.integer(n_train)
    n_features <- as.integer(n_features)
    if (length(n_train) != 1L || is.na(n_train) || n_train < 1L) {
      stop("mock_roi_data: `n_train` must be a positive integer.", call. = FALSE)
    }
    if (length(n_features) != 1L || is.na(n_features) || n_features < 1L) {
      stop("mock_roi_data: `n_features` must be a positive integer.", call. = FALSE)
    }
    train_data <- matrix(stats::rnorm(n_train * n_features), nrow = n_train, ncol = n_features)
  } else {
    train_data <- as.matrix(train_data)
  }

  if (is.null(test_data)) {
    n_test <- as.integer(n_test)
    if (length(n_test) != 1L || is.na(n_test) || n_test < 0L) {
      stop("mock_roi_data: `n_test` must be a non-negative integer.", call. = FALSE)
    }
    if (n_test > 0L) {
      test_data <- matrix(stats::rnorm(n_test * ncol(train_data)), nrow = n_test, ncol = ncol(train_data))
    }
  } else {
    test_data <- as.matrix(test_data)
  }

  if (!is.null(test_data) && ncol(test_data) != ncol(train_data)) {
    stop(
      sprintf(
        "mock_roi_data: `test_data` must have %d columns to match `train_data`.",
        ncol(train_data)
      ),
      call. = FALSE
    )
  }

  if (is.null(indices)) {
    indices <- seq_len(ncol(train_data))
  }
  indices <- as.integer(indices)
  if (length(indices) != ncol(train_data)) {
    stop(
      sprintf(
        "mock_roi_data: length(indices) must equal number of features (%d).",
        ncol(train_data)
      ),
      call. = FALSE
    )
  }

  list(
    train_data = train_data,
    test_data = test_data,
    indices = indices,
    train_roi = train_roi,
    test_roi = test_roi
  )
}

#' Build a Mock \code{fit_roi} Context
#'
#' Creates the standard \code{context} list passed to \code{\link{fit_roi}}
#' methods.
#'
#' @param design Optional design object. If \code{NULL}, a small synthetic
#'   \code{\link{mvpa_design}} is created.
#' @param cv_spec Optional cross-validation specification.
#' @param id ROI identifier to include in the context.
#' @param center_global_id Optional global center index.
#'
#' @return A list with fields \code{design}, \code{cv_spec}, \code{id}, and
#'   \code{center_global_id}.
#' @export
mock_context <- function(design = NULL,
                         cv_spec = NULL,
                         id = 1L,
                         center_global_id = NA_integer_) {
  if (is.null(design)) {
    dtab <- data.frame(
      condition = factor(rep(c("A", "B"), each = 5)),
      block = rep(1:2, each = 5)
    )
    design <- mvpa_design(dtab, y_train = ~condition, block_var = ~block)
  }

  list(
    design = design,
    cv_spec = cv_spec,
    id = as.integer(id),
    center_global_id = as.integer(center_global_id)
  )
}

#' Validate a Plugin Model Contract
#'
#' Executes one \code{\link{fit_roi}} call and checks that the model returns
#' a valid \code{\link{roi_result}} with metrics that match
#' \code{\link{output_schema}} (when present).
#'
#' @param model_spec A model specification object.
#' @param roi_data ROI payload passed to \code{fit_roi}. Defaults to
#'   \code{\link{mock_roi_data}()}.
#' @param context Context list passed to \code{fit_roi}. If \code{NULL}, a
#'   context is constructed from \code{model_spec}.
#' @param check_schema Logical; if \code{TRUE}, enforce metric name/width
#'   agreement with \code{output_schema(model_spec)} when schema is defined.
#'
#' @return An object of class \code{plugin_validation_result}.
#' @export
validate_plugin_model <- function(model_spec,
                                  roi_data = mock_roi_data(),
                                  context = NULL,
                                  check_schema = TRUE) {
  if (!inherits(model_spec, "model_spec")) {
    stop("validate_plugin_model: `model_spec` must inherit from 'model_spec'.", call. = FALSE)
  }

  if (is.null(context)) {
    context <- mock_context(
      design = model_spec$design,
      cv_spec = model_spec$crossval,
      id = 1L,
      center_global_id = NA_integer_
    )
  }

  if (is.null(context$design)) {
    context$design <- model_spec$design
  }
  if (is.null(context$cv_spec)) {
    context$cv_spec <- model_spec$crossval
  }
  if (is.null(context$id) || length(context$id) != 1L) {
    context$id <- 1L
  }
  if (is.null(context$center_global_id)) {
    context$center_global_id <- NA_integer_
  }

  res <- fit_roi(model_spec, roi_data, context)
  if (!inherits(res, "roi_result")) {
    stop(
      sprintf("validate_plugin_model: fit_roi.%s must return an roi_result.",
              class(model_spec)[1]),
      call. = FALSE
    )
  }

  if (isTRUE(res$error)) {
    stop(
      sprintf("validate_plugin_model: fit_roi returned an error: %s", res$error_message),
      call. = FALSE
    )
  }

  metrics <- unlist(res$metrics, recursive = TRUE, use.names = TRUE)
  if (!is.numeric(metrics)) {
    stop("validate_plugin_model: `roi_result$metrics` must be numeric.", call. = FALSE)
  }

  if (is.null(names(metrics)) || any(!nzchar(names(metrics)))) {
    stop(
      "validate_plugin_model: `roi_result$metrics` must be a named numeric vector.",
      call. = FALSE
    )
  }

  schema <- output_schema(model_spec)
  if (isTRUE(check_schema) && !is.null(schema)) {
    expected_names <- schema_metric_names(schema)
    if (length(metrics) != length(expected_names)) {
      stop(
        sprintf(
          "validate_plugin_model: schema/metric width mismatch (schema=%d, metrics=%d).",
          length(expected_names), length(metrics)
        ),
        call. = FALSE
      )
    }
    if (!identical(names(metrics), expected_names)) {
      stop(
        sprintf(
          paste(
            "validate_plugin_model: schema/metric name mismatch.",
            "schema: %s.",
            "metrics: %s."
          ),
          paste(expected_names, collapse = ", "),
          paste(names(metrics), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  structure(
    list(
      valid = TRUE,
      model_class = class(model_spec)[1],
      metric_names = names(metrics),
      metric_count = length(metrics),
      schema = schema,
      result = res
    ),
    class = "plugin_validation_result"
  )
}

#' Evaluate One ROI with rMVPA Cross-Validation Helpers
#'
#' Public wrapper around rMVPA's internal ROI cross-validation helpers.
#' This allows plugin authors to reuse standard fold execution without
#' calling non-exported functions.
#'
#' @param model_spec A model specification object.
#' @param roi_data ROI payload as used by \code{\link{fit_roi}}.
#' @param context Context list (defaults to \code{\link{mock_context}()} based on
#'   \code{model_spec}).
#' @param mode One of \code{"auto"}, \code{"internal"}, or \code{"external"}.
#'   \code{"auto"} chooses \code{"external"} when \code{has_test_set(model_spec)}
#'   is true, otherwise \code{"internal"}.
#' @param return One of \code{"roi_result"} (default) or \code{"row"}.
#'   \code{"row"} returns the raw tibble row produced by the CV helper.
#' @param ... Additional arguments passed to the underlying CV helper.
#'
#' @return A \code{\link{roi_result}} by default, or a tibble row when
#'   \code{return = "row"}.
#' @export
cv_evaluate_roi <- function(model_spec,
                            roi_data,
                            context = NULL,
                            mode = c("auto", "internal", "external"),
                            return = c("roi_result", "row"),
                            ...) {
  if (!inherits(model_spec, "model_spec")) {
    stop("cv_evaluate_roi: `model_spec` must inherit from 'model_spec'.", call. = FALSE)
  }
  if (!is.list(roi_data) || is.null(roi_data$train_data)) {
    stop("cv_evaluate_roi: `roi_data` must be a list containing `train_data`.", call. = FALSE)
  }

  mode <- match.arg(mode)
  return <- match.arg(return)

  if (is.null(context)) {
    context <- mock_context(
      design = model_spec$design,
      cv_spec = model_spec$crossval,
      id = 1L,
      center_global_id = NA_integer_
    )
  }

  train_data <- as.matrix(roi_data$train_data)
  test_data <- if (!is.null(roi_data$test_data)) as.matrix(roi_data$test_data) else NULL

  indices <- roi_data$indices
  if (is.null(indices)) {
    if (!is.null(roi_data$train_roi)) {
      indices <- neuroim2::indices(roi_data$train_roi)
    } else {
      indices <- seq_len(ncol(train_data))
    }
  }
  indices <- as.integer(indices)

  if (mode == "auto") {
    mode <- if (has_test_set(model_spec)) "external" else "internal"
  }

  if (mode == "external" && is.null(test_data) && is.null(roi_data$test_roi)) {
    stop(
      "cv_evaluate_roi: external mode requires `roi_data$test_data` or `roi_data$test_roi`.",
      call. = FALSE
    )
  }
  if (mode == "external" && !isTRUE(has_test_set(model_spec))) {
    stop(
      "cv_evaluate_roi: external mode requires a model spec with external test labels (has_test_set(model_spec) == TRUE).",
      call. = FALSE
    )
  }

  roi_obj <- list(train_roi = roi_data$train_roi, test_roi = roi_data$test_roi)
  roi_id <- as.integer(context$id %||% 1L)
  center_global_id <- as.integer(context$center_global_id %||% NA_integer_)

  row <- switch(
    mode,
    internal = internal_crossval(
      model_spec,
      roi = roi_obj,
      id = roi_id,
      center_global_id = center_global_id,
      x_all = train_data,
      ind = indices,
      ...
    ),
    external = external_crossval(
      model_spec,
      roi = roi_obj,
      id = roi_id,
      center_global_id = center_global_id,
      xtrain_mat = train_data,
      xtest_mat = test_data,
      ind = indices,
      ...
    )
  )

  if (identical(return, "row")) {
    return(row)
  }

  if (!is.data.frame(row) || nrow(row) == 0L) {
    stop("cv_evaluate_roi: CV helper did not return a result row.", call. = FALSE)
  }

  row1 <- row[1, , drop = FALSE]
  is_error <- isTRUE(row1$error[[1]])
  err_msg <- if ("error_message" %in% names(row1)) row1$error_message[[1]] else "~"
  warn_flag <- if ("warning" %in% names(row1)) isTRUE(row1$warning[[1]]) else FALSE
  warn_msg <- if ("warning_message" %in% names(row1)) row1$warning_message[[1]] else "~"

  if (is_error) {
    return(roi_result(
      metrics = NULL,
      indices = indices,
      id = roi_id,
      result = if ("result" %in% names(row1)) row1$result[[1]] else NULL,
      error = TRUE,
      error_message = err_msg %||% "CV evaluation failed",
      warning = warn_flag,
      warning_message = warn_msg %||% err_msg
    ))
  }

  roi_result(
    metrics = if ("performance" %in% names(row1)) row1$performance[[1]] else NULL,
    indices = indices,
    id = roi_id,
    result = if ("result" %in% names(row1)) row1$result[[1]] else NULL,
    warning = warn_flag,
    warning_message = warn_msg %||% "~"
  )
}

#' Validate a Model Specification for Plugin Readiness
#'
#' Performs structural checks on a model spec and optionally executes a
#' one-ROI dry run.
#'
#' @param model_spec A model specification object.
#' @param require_schema Logical; if \code{TRUE}, fail when
#'   \code{output_schema(model_spec)} is \code{NULL}.
#' @param dry_run Logical; if \code{TRUE}, run \code{\link{validate_plugin_model}}
#'   on mock ROI/context inputs.
#' @param roi_data ROI payload used for \code{dry_run}.
#' @param context Context used for \code{dry_run}.
#'
#' @return An object of class \code{model_spec_validation_result}.
#' @export
validate_model_spec <- function(model_spec,
                                require_schema = FALSE,
                                dry_run = TRUE,
                                roi_data = mock_roi_data(),
                                context = NULL) {
  checks <- list()
  add_check <- function(name, status, message) {
    checks[[length(checks) + 1L]] <<- list(name = name, status = status, message = message)
  }

  if (!inherits(model_spec, "model_spec")) {
    add_check("class_chain", "fail", "Object does not inherit from 'model_spec'.")
    return(structure(
      list(valid = FALSE, model_class = class(model_spec)[1] %||% "<unknown>", checks = checks),
      class = "model_spec_validation_result"
    ))
  }

  cls <- class(model_spec)
  primary_class <- cls[1]
  add_check("class_chain", "pass", sprintf("Class chain: %s", paste(cls, collapse = " -> ")))

  if (identical(primary_class, "model_spec")) {
    add_check("primary_class", "warn", "Primary class is 'model_spec'; plugin class should be first.")
  } else {
    add_check("primary_class", "pass", sprintf("Primary class is '%s'.", primary_class))
  }

  fit_method_class <- NULL
  for (cls_i in cls) {
    if (!is.null(utils::getS3method("fit_roi", cls_i, optional = TRUE))) {
      fit_method_class <- cls_i
      break
    }
  }
  if (is.null(fit_method_class)) {
    add_check("fit_roi_method", "fail", sprintf("No fit_roi method found in class chain for '%s'.", primary_class))
  } else {
    add_check("fit_roi_method", "pass", sprintf("Found fit_roi.%s.", fit_method_class))
  }

  schema <- tryCatch(output_schema(model_spec), error = function(e) e)
  if (inherits(schema, "error")) {
    add_check("output_schema", "fail", sprintf("output_schema failed: %s", conditionMessage(schema)))
  } else if (is.null(schema)) {
    if (isTRUE(require_schema)) {
      add_check("output_schema", "fail", "output_schema returned NULL but require_schema=TRUE.")
    } else {
      add_check("output_schema", "warn", "output_schema returned NULL (legacy combiner path).")
    }
  } else {
    schema_ok <- tryCatch({
      schema_metric_names(schema)
      TRUE
    }, error = function(e) {
      add_check("output_schema", "fail", sprintf("Invalid output_schema: %s", conditionMessage(e)))
      FALSE
    })
    if (isTRUE(schema_ok)) {
      add_check("output_schema", "pass", sprintf("output_schema declares %d metric column(s).",
                                                 length(schema_metric_names(schema))))
    }
  }

  if (isTRUE(dry_run) && !is.null(fit_method_class)) {
    dry <- tryCatch({
      validate_plugin_model(
        model_spec = model_spec,
        roi_data = roi_data,
        context = context,
        check_schema = !is.null(schema)
      )
    }, error = function(e) e)

    if (inherits(dry, "error")) {
      add_check("dry_run", "fail", sprintf("Dry run failed: %s", conditionMessage(dry)))
    } else {
      add_check("dry_run", "pass", "Dry run fit_roi validation succeeded.")
    }
  }

  any_fail <- any(vapply(checks, function(ch) identical(ch$status, "fail"), logical(1)))
  structure(
    list(valid = !any_fail, model_class = primary_class, checks = checks),
    class = "model_spec_validation_result"
  )
}

#' @keywords internal
#' @noRd
.plugin_preflight <- function(model_spec, context = "analysis") {
  if (!inherits(model_spec, "model_spec")) {
    return(invisible(NULL))
  }

  vres <- validate_model_spec(
    model_spec = model_spec,
    require_schema = FALSE,
    dry_run = FALSE
  )
  if (isTRUE(vres$valid)) {
    return(invisible(vres))
  }

  fail_checks <- Filter(function(ch) identical(ch$status, "fail"), vres$checks)
  detail <- if (length(fail_checks) == 0L) {
    "unknown validation failure"
  } else {
    paste(
      vapply(
        fail_checks,
        function(ch) sprintf("- %s: %s", ch$name, ch$message),
        character(1)
      ),
      collapse = "\n"
    )
  }

  stop(
    sprintf(
      "%s: plugin contract preflight failed.\n%s",
      context, detail
    ),
    call. = FALSE
  )
}

#' @export
#' @method print model_spec_validation_result
print.model_spec_validation_result <- function(x, ...) {
  n_pass <- sum(vapply(x$checks, function(ch) identical(ch$status, "pass"), logical(1)))
  n_warn <- sum(vapply(x$checks, function(ch) identical(ch$status, "warn"), logical(1)))
  n_fail <- sum(vapply(x$checks, function(ch) identical(ch$status, "fail"), logical(1)))

  cat(sprintf(
    "<model_spec_validation_result> class=%s valid=%s (pass=%d warn=%d fail=%d)\n",
    x$model_class, if (isTRUE(x$valid)) "TRUE" else "FALSE", n_pass, n_warn, n_fail
  ))
  invisible(x)
}

#' @export
#' @method print plugin_validation_result
print.plugin_validation_result <- function(x, ...) {
  cat(sprintf(
    "<plugin_validation_result> class=%s, metrics=%d\n",
    x$model_class,
    x$metric_count
  ))
  if (length(x$metric_names) > 0L) {
    cat("  metric_names:", paste(x$metric_names, collapse = ", "), "\n")
  }
  invisible(x)
}
