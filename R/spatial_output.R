#' Build spatial output maps from an output schema
#'
#' Uses a named schema list to select columns from a performance matrix
#' and create one spatial map per metric.  Each schema entry is either
#' \code{"scalar"} (one column) or \code{"vector[N]"} (N consecutive columns).
#'
#' @param schema Named list; names are metric labels, values are
#'   \code{"scalar"} or \code{"vector[N]"}.
#' @param perf_mat Numeric matrix with rows = ROIs and columns = metrics.
#' @param dataset The dataset object (used by \code{\link{build_output_map}}).
#' @param ids Integer vector of ROI center IDs corresponding to rows of
#'   \code{perf_mat}.
#' @return A named list of spatial map objects (one per schema entry for
#'   scalar, N per entry for vector types).
#' @keywords internal
build_maps_from_schema <- function(schema, perf_mat, dataset, ids) {
  perf_mat <- as.matrix(perf_mat)
  maps <- list()
  col_idx <- 1L

  for (nm in names(schema)) {
    spec <- schema[[nm]]

    if (spec == "scalar") {
      maps[[nm]] <- build_output_map(dataset, perf_mat[, col_idx], ids)
      col_idx <- col_idx + 1L
    } else if (grepl("^vector\\[\\d+\\]$", spec)) {
      n <- as.integer(sub("vector\\[(\\d+)\\]", "\\1", spec))
      for (k in seq_len(n)) {
        sub_name <- paste0(nm, ".", k)
        maps[[sub_name]] <- build_output_map(dataset, perf_mat[, col_idx], ids)
        col_idx <- col_idx + 1L
      }
    } else {
      stop(sprintf("Unknown schema type '%s' for metric '%s'", spec, nm))
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
  perf_mat <- do.call(rbind, good_results$performance)

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
