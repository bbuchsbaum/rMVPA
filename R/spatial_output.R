#' Metric Schema Constructors
#'
#' Helpers for defining typed metric entries in \code{\link{output_schema}}
#' methods for plugin models.
#'
#' \code{schema_scalar()} declares a single metric column.
#' \code{schema_vector(n)} declares \code{n} metric columns expanded as
#' \code{metric.1}, \code{metric.2}, ..., \code{metric.n}.
#'
#' These constructors are optional. Legacy string specs
#' (\code{"scalar"}, \code{"vector[N]"}) remain supported.
#'
#' @param n Positive integer length for vector-valued metrics.
#' @return An object of class \code{rmvpa_metric_spec}.
#' @examples
#' schema_scalar()
#' schema_vector(3)
#' @name schema_metric_spec
NULL

#' @rdname schema_metric_spec
#' @export
schema_scalar <- function() {
  structure(list(type = "scalar", size = 1L), class = "rmvpa_metric_spec")
}

#' @rdname schema_metric_spec
#' @export
schema_vector <- function(n) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("schema_vector: `n` must be a positive integer.", call. = FALSE)
  }
  structure(list(type = "vector", size = n), class = "rmvpa_metric_spec")
}

#' @keywords internal
#' @noRd
.as_metric_schema_spec <- function(spec, metric) {
  if (inherits(spec, "rmvpa_metric_spec")) {
    type <- spec$type
    size <- spec$size
    if (!is.character(type) || length(type) != 1L || !type %in% c("scalar", "vector")) {
      stop(sprintf(
        "Invalid typed schema for metric '%s': `type` must be 'scalar' or 'vector'.",
        metric
      ), call. = FALSE)
    }
    if (!is.numeric(size) || length(size) != 1L || is.na(size)) {
      stop(sprintf(
        "Invalid typed schema for metric '%s': `size` must be a positive integer.",
        metric
      ), call. = FALSE)
    }
    size <- as.integer(size)
    if (type == "scalar" && size != 1L) {
      stop(sprintf(
        "Invalid typed schema for metric '%s': scalar metrics must have size 1.",
        metric
      ), call. = FALSE)
    }
    if (type == "vector" && size < 1L) {
      stop(sprintf(
        "Invalid typed schema for metric '%s': vector size must be >= 1.",
        metric
      ), call. = FALSE)
    }
    return(list(type = type, size = size))
  }

  if (is.character(spec) && length(spec) == 1L && !is.na(spec) && identical(spec, "scalar")) {
    return(list(type = "scalar", size = 1L))
  }

  if (is.character(spec) && length(spec) == 1L && !is.na(spec) &&
      grepl("^vector\\[[1-9][0-9]*\\]$", spec)) {
    n <- as.integer(sub("vector\\[([1-9][0-9]*)\\]", "\\1", spec))
    return(list(type = "vector", size = n))
  }

  stop(sprintf(
    paste(
      "Unknown schema type '%s' for metric '%s'.",
      "Expected 'scalar', 'vector[N]', schema_scalar(), or schema_vector(N)."
    ),
    as.character(spec), metric
  ), call. = FALSE)
}

#' @keywords internal
#' @noRd
.normalize_output_schema <- function(schema) {
  if (is.character(schema) && !is.null(names(schema))) {
    schema <- as.list(schema)
  }
  if (!is.list(schema)) {
    stop("output_schema must return a list.", call. = FALSE)
  }
  nms <- names(schema)
  if (is.null(nms) || any(!nzchar(nms))) {
    stop("output_schema must be a named list with non-empty metric names.", call. = FALSE)
  }
  if (anyDuplicated(nms)) {
    dups <- unique(nms[duplicated(nms)])
    stop(sprintf("output_schema has duplicate metric names: %s",
                 paste(dups, collapse = ", ")), call. = FALSE)
  }

  out <- lapply(seq_along(schema), function(i) {
    .as_metric_schema_spec(schema[[i]], metric = nms[[i]])
  })
  names(out) <- nms
  out
}

#' @keywords internal
schema_metric_names <- function(schema) {
  schema_specs <- .normalize_output_schema(schema)

  out <- character()
  for (metric in names(schema_specs)) {
    spec <- schema_specs[[metric]]
    if (identical(spec$type, "scalar")) {
      out <- c(out, metric)
    } else {
      out <- c(out, paste0(metric, ".", seq_len(spec$size)))
    }
  }

  out
}

#' @keywords internal
coerce_performance_vector <- function(x, row_idx) {
  if (is.null(x)) {
    stop(sprintf("performance row %d is NULL.", row_idx), call. = FALSE)
  }

  if (is.matrix(x)) {
    if (nrow(x) != 1L) {
      stop(sprintf("performance row %d must have exactly 1 row (got %d).",
                   row_idx, nrow(x)), call. = FALSE)
    }
    vec <- as.numeric(x[1, , drop = TRUE])
    nms <- colnames(x)
    if (!is.null(nms) && length(nms) == length(vec)) {
      names(vec) <- nms
    }
    return(vec)
  }

  raw <- unlist(x, recursive = TRUE, use.names = TRUE)
  if (!is.atomic(raw)) {
    stop(sprintf("performance row %d could not be flattened to an atomic vector.", row_idx),
         call. = FALSE)
  }
  if (!is.numeric(raw)) {
    raw_chr <- as.character(raw)
    vec <- suppressWarnings(as.numeric(raw_chr))
    if (any(is.na(vec) & !is.na(raw_chr))) {
      stop(sprintf("performance row %d contains non-numeric values.", row_idx),
           call. = FALSE)
    }
    names(vec) <- names(raw)
    return(vec)
  }

  vec <- as.numeric(raw)
  names(vec) <- names(raw)
  vec
}

#' @keywords internal
coerce_performance_matrix <- function(perf_list) {
  if (!is.list(perf_list) || length(perf_list) == 0L) {
    stop("No performance metrics available to combine.", call. = FALSE)
  }

  perf_vecs <- lapply(seq_along(perf_list), function(i) {
    vec <- coerce_performance_vector(perf_list[[i]], i)
    if (is.null(names(vec)) || length(names(vec)) != length(vec)) {
      names(vec) <- names(unlist(perf_list[[i]], recursive = TRUE, use.names = TRUE))
    }
    vec
  })

  widths <- vapply(perf_vecs, length, integer(1))
  if (length(unique(widths)) != 1L) {
    stop(sprintf("Inconsistent metric lengths across ROIs: %s",
                 paste(widths, collapse = ", ")), call. = FALSE)
  }

  has_complete_names <- vapply(
    perf_vecs,
    function(v) !is.null(names(v)) && length(names(v)) == length(v) &&
      all(nzchar(names(v))),
    logical(1)
  )

  perf_names <- NULL
  if (all(has_complete_names)) {
    perf_names <- names(perf_vecs[[1]])
    mismatch_idx <- which(!vapply(perf_vecs, function(v) identical(names(v), perf_names), logical(1)))
    if (length(mismatch_idx) > 0L) {
      stop(sprintf("Inconsistent metric names across ROIs (first mismatch at row %d).",
                   mismatch_idx[[1]]), call. = FALSE)
    }
  }

  perf_mat <- do.call(rbind, lapply(perf_vecs, unname))
  if (!is.null(perf_names)) {
    colnames(perf_mat) <- perf_names
  }
  perf_mat
}

#' Build spatial output maps from an output schema
#'
#' Uses a named schema list to select columns from a performance matrix
#' and create one spatial map per metric.
#'
#' @param schema Named list; names are metric labels, values are
#'   \code{"scalar"} or \code{"vector[N]"}.
#' @param perf_mat Numeric matrix with rows = ROIs and columns = metrics.
#' @param dataset The dataset object (used by \code{\link{build_output_map}}).
#' @param ids Integer vector of ROI center IDs corresponding to rows of
#'   \code{perf_mat}.
#' @return A named list of spatial map objects.
#' @keywords internal
build_maps_from_schema <- function(schema, perf_mat, dataset, ids) {
  perf_mat <- as.matrix(perf_mat)
  schema_specs <- .normalize_output_schema(schema)
  expected_names <- schema_metric_names(schema)
  if (ncol(perf_mat) != length(expected_names)) {
    stop(
      sprintf(
        "Schema/performance width mismatch: schema declares %d metric columns but ROI results produced %d.",
        length(expected_names), ncol(perf_mat)
      ),
      call. = FALSE
    )
  }

  maps <- list()
  col_idx <- 1L

  for (nm in names(schema_specs)) {
    spec <- schema_specs[[nm]]
    if (identical(spec$type, "scalar")) {
      maps[[nm]] <- build_output_map(dataset, perf_mat[, col_idx], ids)
      col_idx <- col_idx + 1L
    } else if (identical(spec$type, "vector")) {
      for (k in seq_len(spec$size)) {
        sub_name <- paste0(nm, ".", k)
        maps[[sub_name]] <- build_output_map(dataset, perf_mat[, col_idx], ids)
        col_idx <- col_idx + 1L
      }
    } else {
      stop(sprintf("Unknown normalized schema type '%s' for metric '%s'", spec$type, nm))
    }
  }

  maps
}


#' Schema-driven combiner for searchlight results
#'
#' Drop-in replacement for \code{\link{combine_standard}} that uses the
#' model's \code{\link{output_schema}} to build output maps. Falls back to
#' \code{combine_standard} if the schema is \code{NULL}.
#'
#' @inheritParams combine_standard
#' @return A \code{searchlight_result} object.
#' @importFrom futile.logger flog.error
#' @keywords internal
combine_schema_standard <- function(model_spec, good_results, bad_results) {
  schema <- output_schema(model_spec)

  if (is.null(schema)) {
    return(combine_standard(model_spec, good_results, bad_results))
  }

  if (nrow(good_results) == 0) {
    futile.logger::flog.error(
      "No successful results to combine (schema combiner). %d ROIs failed.",
      nrow(bad_results)
    )
    stop("No valid results for standard searchlight: all ROIs failed to process")
  }

  ind <- unlist(good_results$id)
  expected_names <- schema_metric_names(schema)
  perf_mat <- coerce_performance_matrix(good_results$performance)

  if (ncol(perf_mat) != length(expected_names)) {
    stop(
      sprintf(
        "Schema/performance width mismatch: schema declares %d metric columns but ROI results produced %d.",
        length(expected_names), ncol(perf_mat)
      ),
      call. = FALSE
    )
  }

  perf_names <- colnames(perf_mat)
  if (!is.null(perf_names) && !identical(perf_names, expected_names)) {
    stop(
      sprintf(
        paste(
          "Schema/performance metric name mismatch.",
          "Schema expects: %s.",
          "ROI metrics are: %s."
        ),
        paste(expected_names, collapse = ", "),
        paste(perf_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  output_maps <- build_maps_from_schema(schema, perf_mat, model_spec$dataset, ind)

  metric_names <- names(output_maps)

  structure(
    list(
      results = output_maps,
      n_voxels = estimate_mask_size(model_spec$dataset),
      active_voxels = count_active_voxels(model_spec$dataset, ind, perf_mat),
      metrics = metric_names
    ),
    class = c("searchlight_result", "list")
  )
}
