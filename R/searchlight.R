#' Wrap output results
#'
#' This function wraps the output results of the performance matrix into a list
#' of spatial objects (NeuroVol or NeuroSurface) for each column 
#' in the performance matrix, and structures it as a searchlight_result.
#'
#' @keywords internal
#' @param perf_mat A performance matrix (voxels/vertices x metrics) containing classifier results.
#' @param dataset A dataset object containing the dataset information (including mask and type).
#' @param ids An integer vector of voxel/vertex indices corresponding to the rows of `perf_mat`.
#'   These are typically global indices into the mask space for volumetric data, or vertex numbers for surface data.
#' @return A `searchlight_result` object containing
#'   \itemize{
#'     \item `results`: Named spatial maps for each metric.
#'     \item `n_voxels`: Total number of voxels/vertices defined by the mask.
#'     \item `active_voxels`: Number of voxels/vertices with results.
#'     \item `metrics`: Character vector of metric names.
#'   }
wrap_out <- function(perf_mat, dataset, ids = NULL) {
  validate_wrap_inputs(perf_mat, ids)

  if (is_perf_empty(perf_mat)) {
    return(empty_searchlight_result(dataset))
  }

  metric_names <- colnames(perf_mat)
  if (is.null(metric_names)) {
    metric_names <- paste0("Metric", seq_len(ncol(perf_mat)))
  }

  is_surface_dataset <- inherits(dataset, "mvpa_surface_dataset")
  output_maps <- vector("list", length(metric_names))

  for (idx in seq_along(metric_names)) {
    metric_vector <- perf_mat[, idx]
    output_maps[[idx]] <- if (is_surface_dataset) {
      build_surface_map(dataset, metric_vector, ids)
    } else {
      build_volume_map(dataset, metric_vector, ids)
    }
  }
  names(output_maps) <- metric_names

  structure(
    list(
      results = output_maps,
      n_voxels = estimate_mask_size(dataset),
      active_voxels = count_active_voxels(dataset, ids, perf_mat),
      metrics = metric_names
    ),
    class = c("searchlight_result", "list")
  )
}

validate_wrap_inputs <- function(perf_mat, ids) {
  if (!is.null(ids) && !is.null(perf_mat) && nrow(perf_mat) > 0L && nrow(perf_mat) != length(ids)) {
    stop("Number of rows in `perf_mat` must match the length of `ids` when ids are provided.")
  }
}

is_perf_empty <- function(perf_mat) {
  is.null(perf_mat) || ncol(perf_mat) == 0L
}

empty_searchlight_result <- function(dataset) {
  structure(
    list(
      results = list(),
      n_voxels = estimate_mask_size(dataset),
      active_voxels = 0L,
      metrics = character(0)
    ),
    class = c("searchlight_result", "list")
  )
}

build_surface_map <- function(dataset, metric_vector, ids) {
  if (!requireNamespace("neurosurf", quietly = TRUE)) {
    stop("Package 'neurosurf' is required to handle surface datasets.")
  }
  surface_ids <- resolve_surface_ids(dataset, ids, length(metric_vector))
  neurosurf::NeuroSurface(
    geometry = neurosurf::geometry(dataset$train_data),
    indices = surface_ids,
    data = metric_vector
  )
}

resolve_surface_ids <- function(dataset, ids, vec_len) {
  if (!is.null(ids)) {
    return(ids)
  }
  if (!requireNamespace("neurosurf", quietly = TRUE)) {
    stop("Package 'neurosurf' is required to handle surface datasets.")
  }
  mask <- dataset$mask
  candidate <- NULL
  if (inherits(mask, "NeuroSurface")) {
    geom <- neurosurf::geometry(mask)
    candidate <- neurosurf::nodes(geom)
  } else if (inherits(mask, "SurfaceGeometry")) {
    candidate <- neurosurf::nodes(mask)
  } else if (is.numeric(mask) || is.logical(mask)) {
    candidate <- which(mask != 0)
  } else if (!is.null(dataset$train_data)) {
    geom <- neurosurf::geometry(dataset$train_data)
    candidate <- neurosurf::nodes(geom)
  }
  if (is.null(candidate)) {
    stop("Cannot determine vertex indices for surface data when 'ids' is NULL.")
  }
  if (vec_len != length(candidate)) {
    stop(sprintf("Length of metric_vector (%s) does not match number of vertices (%s) for dense surface map.",
                 vec_len, length(candidate)))
  }
  candidate
}

build_volume_map <- function(dataset, metric_vector, ids) {
  space_obj <- resolve_volume_space(dataset)
  if (is.null(ids)) {
    expected_len <- spatial_dim_product(space_obj)
    if (length(metric_vector) != expected_len) {
      stop(sprintf("Length of metric_vector (%s) does not match total voxels in mask space (%s) for dense volumetric map when ids is NULL.",
                   length(metric_vector), expected_len))
    }
    vol_dims <- spatial_dim_shape(space_obj)
    vol_array <- array(metric_vector, dim = vol_dims)
    return(neuroim2::NeuroVol(vol_array, space_obj))
  }
  neuroim2::NeuroVol(
    data = metric_vector,
    space = space_obj,
    indices = ids
  )
}

resolve_volume_space <- function(dataset) {
  mask <- dataset$mask
  if (inherits(mask, c("NeuroVol", "NeuroVec"))) {
    return(spatial_only_space(neuroim2::space(mask)))
  }
  if (!is.null(mask) && (is.numeric(mask) || is.logical(mask))) {
    if (!is.null(dataset$train_data) && inherits(dataset$train_data, "NeuroVec")) {
      return(spatial_only_space(neuroim2::space(dataset$train_data)))
    }
    stop("Cannot determine volumetric space from numeric/logical mask without train_data.")
  }
  if (!is.null(dataset$train_data) && inherits(dataset$train_data, "NeuroVec")) {
    return(spatial_only_space(neuroim2::space(dataset$train_data)))
  }
  stop("Cannot determine volumetric space from dataset mask or train_data.")
}

spatial_dim_product <- function(space_obj) {
  dims <- spatial_dim_shape(space_obj)
  prod(dims)
}

spatial_dim_shape <- function(space_obj) {
  dims <- dim(space_obj)
  if (length(dims) > 3L) {
    dims <- dims[seq_len(3L)]
  }
  dims
}

spatial_only_space <- function(space_obj) {
  dims <- spatial_dim_shape(space_obj)
  spacing_vals <- neuroim2::spacing(space_obj)[seq_len(length(dims))]
  origin_vals <- neuroim2::origin(space_obj)[seq_len(length(dims))]
  neuroim2::NeuroSpace(dims, spacing = spacing_vals, origin = origin_vals)
}

resolve_volume_mask <- function(mask, spatial_length = NULL) {
  if (inherits(mask, c("NeuroVol", "NeuroVec"))) {
    vals <- neuroim2::values(mask)
    if (is.matrix(vals)) {
      vals <- vals[, 1, drop = TRUE]
    }
    return(as.logical(as.numeric(vals)))
  }
  if (is.numeric(mask) || is.logical(mask)) {
    return(as.logical(mask))
  }
  if (!is.null(spatial_length)) {
    return(rep(TRUE, spatial_length))
  }
  logical(0)
}

estimate_mask_size <- function(dataset) {
  mask <- dataset$mask
  if (inherits(mask, c("NeuroVol", "NeuroVec"))) {
    return(spatial_dim_product(neuroim2::space(mask)))
  }
  if (inherits(mask, "NeuroSurface")) {
    if (!requireNamespace("neurosurf", quietly = TRUE)) {
      stop("Package 'neurosurf' is required to handle surface datasets.")
    }
    geom <- neurosurf::geometry(mask)
    return(length(neurosurf::nodes(geom)))
  }
  if (is.numeric(mask) || is.logical(mask)) {
    return(length(mask))
  }
  if (!is.null(dataset$train_data) && inherits(dataset$train_data, "NeuroVec")) {
    return(spatial_dim_product(neuroim2::space(dataset$train_data)))
  }
  0L
}

count_active_voxels <- function(dataset, ids, perf_mat) {
  if (!is.null(ids)) {
    return(length(ids))
  }
  mask <- dataset$mask
  if (inherits(mask, c("NeuroVol", "NeuroVec"))) {
    mask_values <- neuroim2::values(mask)
    return(sum(mask_values != 0))
  }
  if (inherits(mask, "NeuroSurface")) {
    if (!requireNamespace("neurosurf", quietly = TRUE)) {
      stop("Package 'neurosurf' is required to handle surface datasets.")
    }
    geom <- neurosurf::geometry(mask)
    return(length(neurosurf::nodes(geom)))
  }
  if (is.numeric(mask) || is.logical(mask)) {
    return(sum(mask != 0))
  }
  if (!is.null(perf_mat)) {
    return(nrow(perf_mat))
  }
  0L
}

#' @export
#' @method print searchlight_result
print.searchlight_result <- function(x, ...) {
  # Ensure crayon is available
  if (!requireNamespace("crayon", quietly = TRUE)) {
    # stop("Package 'crayon' is required for pretty printing. Please install it.")
    # Fallback to basic print if crayon is not there, to avoid breaking R CMD check
    has_crayon <- FALSE
  } else {
    has_crayon <- TRUE
  }
  
  # Define color scheme
  header_style <- if (has_crayon) crayon::bold$cyan else function(txt) txt
  section_style <- if (has_crayon) crayon::yellow else function(txt) txt
  info_style <- if (has_crayon) crayon::white else function(txt) txt
  number_style <- if (has_crayon) crayon::green else function(txt) txt
  metric_style <- if (has_crayon) crayon::magenta else function(txt) txt
  type_style <- if (has_crayon) crayon::blue else function(txt) txt # For object types
  
  # Print header
  cat("\n", header_style("Searchlight Analysis Results"), "\n\n")
  
  # Basic information
  cat(section_style("- Coverage"), "\n")
  cat(info_style("  - Voxels/Vertices in Mask: "), number_style(format(x$n_voxels, big.mark=",")), "\n")
  cat(info_style("  - Voxels/Vertices with Results: "), number_style(format(x$active_voxels, big.mark=",")), "\n")
  
  # Performance metrics (now direct spatial objects)
  cat(section_style("- Output Maps (Metrics)"), "\n")
  if (length(x$metrics) > 0 && !is.null(x$results)) {
    for (metric_name in x$metrics) {
      if (metric_name %in% names(x$results)) {
        metric_map <- x$results[[metric_name]]
        map_type <- class(metric_map)[1]
        cat(info_style("  - "), metric_style(metric_name), 
            info_style(" (Type: "), type_style(map_type), info_style(")"), "\n")
        
        # Optionally, add simple summary stats if feasible and desired
        # This requires knowing how to get data from metric_map (NeuroVol/NeuroSurfaceVector)
        # For example (conceptual, needs specific accessors for neuroim2/neurosurf objects):
        # data_values <- try(as.vector(metric_map), silent = TRUE) # or specific accessor
        # if (!inherits(data_values, "try-error") && is.numeric(data_values)) {
        #   data_values_no_na <- data_values[!is.na(data_values) & data_values != 0]
        #   if (length(data_values_no_na) > 0) {
        #     cat(info_style("    - Mean (non-zero): "), number_style(sprintf("%.4f", mean(data_values_no_na))), "\n")
        #     cat(info_style("    - SD   (non-zero): "), number_style(sprintf("%.4f", sd(data_values_no_na))), "\n")
        #     cat(info_style("    - Range(non-zero): "), number_style(sprintf("[%.4f, %.4f]", min(data_values_no_na), max(data_values_no_na))), "\n")
        #   } else {
        #     cat(info_style("    - (No non-zero data for summary stats)"), "\n")
        #   }
        # } else {
        #    cat(info_style("    - (Summary stats not available for this map type)"), "\n")
        # }
      } else {
        cat(info_style("  - "), metric_style(metric_name), info_style(" (Map data not found in results)"), "\n")
      }
    }
  } else {
    cat(info_style("  - No output maps found or metrics list is empty."),"\n")
  }
  
  if (!is.null(x$pobserved)) {
    cat(section_style("\n- Observed Probabilities Map"), "\n")
    cat(info_style("  - Type: "), type_style(class(x$pobserved)[1]), "\n")
  }
  
  cat("\n")
}

#' Combine standard classifier results
#'
#' This function combines the standard classifier results from a good results data frame
#' by binding the performance rows together and optionally computes the observed probabilities.
#'
#' @keywords internal
#' @param model_spec A list containing the model specification
#' @param good_results A data frame containing the successful classifier results
#' @param bad_results A data frame containing the unsuccessful classifier results
#' @return A list containing the combined performance matrix and other information
combine_standard <- function(model_spec, good_results, bad_results) {
  # Add improved error handling with diagnostics
  if (nrow(good_results) == 0) {
    futile.logger::flog.error("No successful results to combine. Examining errors from bad results:")
    
    # Analyze bad results to provide better diagnostics
    if (nrow(bad_results) > 0) {
      # Group error messages and count occurrences
      error_summary <- table(bad_results$error_message)
      futile.logger::flog.error("Error summary from %d failed ROIs:", nrow(bad_results))
      
      for (i in seq_along(error_summary)) {
        error_msg <- names(error_summary)[i]
        count <- error_summary[i]
        futile.logger::flog.error("  - %s: %d occurrences", error_msg, count)
      }
      
      # Show a sample of the first few errors
      sample_size <- min(5, nrow(bad_results))
      futile.logger::flog.error("Sample of first %d errors:", sample_size)
      for (i in 1:sample_size) {
        futile.logger::flog.error("  ROI %s: %s", 
                                 bad_results$id[i], 
                                 bad_results$error_message[i])
      }
    } else {
      futile.logger::flog.error("No error information available.")
    }
    
    stop("No valid results for standard searchlight: all ROIs failed to process")
  }
  
  result <- NULL


  
  # Proceed with combining results
  tryCatch({
    ind <- unlist(good_results$id)
    perf_mat <- good_results %>% dplyr::select(performance) %>% (function(x) do.call(rbind, x[[1]]))
    
    ret <- wrap_out(perf_mat, model_spec$dataset, ind)
    
    # Optionally construct a 4D (space × trial) map of
    # prob_observed for classification results, when available.
    has_results <- any(unlist(purrr::map(good_results$result, function(x) !is.null(x))))
    if (has_results) {
      pob_list <- good_results %>%
        dplyr::select(result) %>%
        dplyr::pull(result) %>%
        purrr::map(~ prob_observed(.))

      if (any(!vapply(pob_list, is.null, logical(1)))) {
        if (inherits(model_spec$dataset, "mvpa_surface_dataset")) {
          # TODO: implement surface-based prob_observed maps if needed.
          futile.logger::flog.debug("pobserved maps for surface datasets are not yet implemented; skipping.")
        } else {
          # Volumetric case: build a SparseNeuroVec where rows are mask
          # voxels and columns are trials (probability of true class).
          space_obj <- resolve_volume_space(model_spec$dataset)
          mask_vec  <- resolve_volume_mask(model_spec$dataset$mask, spatial_dim_product(space_obj))
          mask_idx  <- which(mask_vec)
          n_mask    <- length(mask_idx)

          # Map global voxel index -> row in prob matrix
          row_map <- integer(length(mask_vec))
          row_map[mask_idx] <- seq_len(n_mask)

          # Trials correspond to y_test (or y_train if no explicit test set)
          n_trials <- length(y_test(model_spec$design))
          prob_mat <- matrix(NA_real_, nrow = n_mask, ncol = n_trials)

          for (i in seq_len(nrow(good_results))) {
            p_i <- pob_list[[i]]
            if (is.null(p_i) || length(p_i) == 0L) next

            res_i <- good_results$result[[i]]
            testind <- res_i$testind
            if (is.null(testind)) {
              # Fallback: assume full coverage in order
              testind <- seq_along(p_i)
            }
            if (length(testind) != length(p_i)) {
              futile.logger::flog.warn(
                "combine_standard: length mismatch between testind (%s) and prob_observed (%s); skipping ROI %s.",
                length(testind), length(p_i), as.character(good_results$id[i])
              )
              next
            }

            # Keep only indices within the available trial range
            keep <- testind >= 1L & testind <= n_trials
            if (!any(keep)) next
            if (!all(keep)) {
              p_i     <- p_i[keep]
              testind <- testind[keep]
            }

            trial_vec <- rep(NA_real_, n_trials)
            trial_vec[testind] <- p_i

            center_id <- unlist(good_results$id[i])
            row_idx   <- if (center_id >= 1L && center_id <= length(row_map)) row_map[center_id] else 0L
            if (row_idx == 0L) next

            prob_mat[row_idx, ] <- trial_vec
          }

          if (any(is.finite(prob_mat))) {
            time_dim <- n_trials
            # Extend space with a time dimension if needed
            if (length(dim(space_obj)) == 3L) {
              space_prob <- neuroim2::add_dim(space_obj, time_dim)
            } else if (dim(space_obj)[length(dim(space_obj))] != time_dim) {
              space_prob <- neuroim2::add_dim(spatial_only_space(space_obj), time_dim)
            } else {
              space_prob <- space_obj
            }

            created_map <- neuroim2::SparseNeuroVec(
              data  = prob_mat,
              space = space_prob,
              mask  = mask_vec
            )
            ret$pobserved <- created_map
          } else {
            futile.logger::flog.debug("combine_standard: no finite prob_observed values; skipping pobserved map.")
          }
        }
      }
    }
    
    return(ret)
  }, error = function(e) {
    futile.logger::flog.error("Error combining results: %s", e$message)
    futile.logger::flog.debug("Error details: %s", e$call)
    stop(paste("Failed to combine searchlight results:", e$message))
  })
}

#' Combine RSA standard classifier results
#'
#' This function combines the RSA standard classifier results from a good results data frame
#' by binding the performance rows together.
#'
#' @keywords internal
#' @param model_spec A list containing the model specification.
#' @param good_results A data frame containing the successful classifier results.
#' @param bad_results A data frame containing the unsuccessful classifier results.
#' @return A list containing the combined performance matrix along with other information from the dataset.
combine_rsa_standard <- function(model_spec, good_results, bad_results) {
  # Enhanced error handling with detailed diagnostics
  if (nrow(good_results) == 0) {
    futile.logger::flog.error("No successful results for RSA searchlight. Examining errors from bad results:")
    
    # Analyze bad results to provide better diagnostics
    if (nrow(bad_results) > 0) {
      # Group error messages and count occurrences
      error_summary <- table(bad_results$error_message)
      futile.logger::flog.error("Error summary from %d failed ROIs:", nrow(bad_results))
      
      for (i in seq_along(error_summary)) {
        error_msg <- names(error_summary)[i]
        count <- error_summary[i]
        futile.logger::flog.error("  - %s: %d occurrences", error_msg, count)
      }
      
      # Show a sample of the first few errors
      sample_size <- min(5, nrow(bad_results))
      futile.logger::flog.error("Sample of first %d errors:", sample_size)
      for (i in 1:sample_size) {
        futile.logger::flog.error("  ROI %s: %s", 
                                 bad_results$id[i], 
                                 bad_results$error_message[i])
      }
      
      # Check if there are any common issues
      if (any(grepl("unable to find an inherited method for function 'values'", bad_results$error_message))) {
        futile.logger::flog.error("  - Many errors involve NULL ROIs. This may be caused by insufficient voxels in searchlight regions.")
      }
      if (any(grepl("insufficient data dimensions", bad_results$error_message))) {
        futile.logger::flog.error("  - Many errors involve insufficient data dimensions. Your feature data may be too high-dimensional for small searchlight regions.")
      }
    } else {
      futile.logger::flog.error("No error information available.")
    }
    
    stop("No valid results for RSA searchlight: all ROIs failed to process")
  }
  
  # Check if performance data is available
  if (length(good_results$performance) == 0 || all(sapply(good_results$performance, is.null))) {
    futile.logger::flog.error("Performance metrics missing in all results")
    stop("No valid performance metrics for RSA searchlight")
  }
  
  ind <- unlist(good_results$id)
  
  # Extract the performance matrix using safe error handling
  tryCatch({
    perf_mat <- good_results %>% dplyr::select(performance) %>% (function(x) do.call(rbind, x[[1]]))
    if (is.null(colnames(perf_mat)) || any(colnames(perf_mat) == "")) {
      expected_names <- names(model_spec$design$model_mat)
      if (length(expected_names) == ncol(perf_mat)) {
        colnames(perf_mat) <- expected_names
      }
    }
    ret <- wrap_out(perf_mat, model_spec$dataset, ind)
    return(ret)
  }, error = function(e) {
    futile.logger::flog.error("Error combining RSA results: %s", e$message)
    
    # Try to provide more specific diagnostic information
    if (grepl("requires numeric/complex matrix/vector arguments", e$message)) {
      futile.logger::flog.error("This error often occurs when performance metrics are incompatible across ROIs. Check your performance calculation.")
    } else if (grepl("subscript out of bounds", e$message)) {
      futile.logger::flog.error("This may be caused by inconsistent dimensions in your results. Check that all ROIs return the same metrics.")
    }
    
    stop(paste("Failed to combine RSA results:", e$message))
  })
}

#' Combine Vector RSA standard classifier results
#'
#' This function combines the Vector RSA standard classifier results from a good results data frame
#' by binding the performance rows together.
#'
#' @keywords internal
#' @param model_spec A list containing the model specification.
#' @param good_results A data frame containing the successful classifier results.
#' @param bad_results A data frame containing the unsuccessful classifier results.
#' @return A list containing the combined performance matrix along with other information from the dataset.
combine_vector_rsa_standard <- function(model_spec, good_results, bad_results) {
  ind <- unlist(good_results$id)
  perf_mat <- good_results %>% dplyr::select(performance) %>% (function(x) do.call(rbind, x[[1]]))
  score_mat <- data.frame(sim=rowMeans(perf_mat))
  ret <- wrap_out(score_mat, model_spec$dataset, ind)
  ret
}

#' Combine randomized classifier results
#'
#' This function combines the randomized classifier results from a good results data frame
#' and normalizes the performance matrix by the number of instances for each voxel index.
#'
#' @keywords internal
#' @param model_spec A list containing the model specification.
#' @param good_results A data frame containing the successful classifier results.
#' @param bad_results A data frame containing the unsuccessful classifier results.
#' @return A list containing the combined and normalized performance matrix along with other information from the dataset.
combine_randomized <- function(model_spec, good_results, bad_results=NULL) {
  # Check if we have results
  if (nrow(good_results) == 0 || length(good_results$performance) == 0) {
    futile.logger::flog.error("No valid results for randomized searchlight")
    stop("No valid results for randomized searchlight")
  }
  
  all_ind   <- sort(unlist(good_results$indices))
  ind_count <- table(all_ind)
  ind_set   <- unique(all_ind)
  
  # Prototype performance entry to determine metric count and names
  perf_proto <- good_results$performance[[1]]
  if (is.matrix(perf_proto) || is.data.frame(perf_proto)) {
    ncols        <- ncol(perf_proto)
    metric_names <- colnames(perf_proto)
  } else {
    ncols        <- length(perf_proto)
    metric_names <- names(perf_proto)
  }
  
  # Create sparse matrix to hold results
  perf_mat <- Matrix::sparseMatrix(
    i    = rep(ind_set, ncols),
    j    = rep(seq_len(ncols), each = length(ind_set)),
    x    = rep(0, length(ind_set) * ncols),
    dims = c(length(model_spec$dataset$mask), ncols)
  )
  
  # Process each result
  for (i in seq_len(nrow(good_results))) {
    ind_i   <- good_results$indices[[i]]
    perf_i  <- good_results$performance[[i]]
    if (is.null(ind_i) || is.null(perf_i)) next

    # Coerce performance entry to a numeric vector of length ncols
    perf_vec <- if (is.matrix(perf_i) || is.data.frame(perf_i)) {
      as.numeric(perf_i[1, , drop = FALSE])
    } else {
      as.numeric(perf_i)
    }
    if (length(perf_vec) != ncols) {
      futile.logger::flog.warn("combine_randomized: performance length mismatch in row %d", i)
      next
    }

    tryCatch({
      m <- kronecker(matrix(perf_vec, 1, ncols), rep(1, length(ind_i)))
      perf_mat[ind_i, ] <- perf_mat[ind_i, ] + m
    }, error = function(e) {
      futile.logger::flog.warn("Error processing result %d: %s", i, e$message)
    })
  }

  # Normalize by the count of overlapping searchlights
  perf_mat[ind_set, ] <- sweep(perf_mat[ind_set, , drop = FALSE],
                               1, as.integer(ind_count), FUN = "/")
  
  # Set column names from the performance metrics (fallback to generic names if needed)
  if (is.null(metric_names)) {
    metric_names <- paste0("Metric", seq_len(ncols))
  }
  colnames(perf_mat) <- metric_names
  
  # Wrap and return results
  ret <- wrap_out(perf_mat, model_spec$dataset)
  ret
}

#' Pool classifier results
#'
#' This function pools classifier results collected over a set of overlapping indices.
#'
#' @keywords internal
#' @param ... A variable list of data frames containing classifier results to be pooled.
#' @return A list of merged classifier results.
#' @noRd
pool_results <- function(...) {
  reslist <- list(...)
  is_df <- sapply(reslist, function(res) inherits(res, "data.frame"))
  assertthat::assert_that(all(is_df), msg="pool_results: all arguments must be of type 'data.frame'")
  good_results <- do.call(rbind, reslist)

  ## map every result to the set of indices in that set
  indmap <- do.call(rbind, lapply(1:nrow(good_results), function(i) {
    ind <- good_results$indices[[i]]
    cbind(i, ind)
  }))
  
  
  respsets <- split(indmap[,1], indmap[,2])
  
  merged_results <- purrr::map(respsets, do_merge_results, good_results=good_results)

  return(merged_results)
}



#' Merge searchlight results
#'
#' This function merges searchlight results, combining the first result with the rest of the results.
#'
#' @keywords internal
#' @param r1 A list of indices representing the searchlight results to be merged.
#' @param good_results A data frame containing the valid searchlight results.
#' @return A combined searchlight result object.
do_merge_results <- function(r1, good_results) {
  if (length(r1) > 1) {
    first <- r1[1]
    rest <- r1[2:length(r1)]
    z1 <- good_results$result[[first]]
    z2 <- good_results$result[rest]
    ff <- purrr::partial(merge_results, x=z1)
    do.call(ff, z2)
  } else {
    good_results$result[[r1[1]]]
  }
}

#' Combine randomized searchlight results by pooling
#'
#' This function combines randomized searchlight results by pooling the good results.
#'
#' @keywords internal
#' @param model_spec An object specifying the model used in the searchlight analysis.
#' @param good_results A data frame containing the valid searchlight results.
#' @param bad_results A data frame containing the invalid searchlight results.
#' @return An object containing the combined searchlight results.
pool_randomized <- function(model_spec, good_results, bad_results) {
  if (nrow(good_results) == 0) {
    stop("searchlight: no searchlight samples produced valid results")
  }
  
  
  merged_results <- pool_results(good_results)
  
  # Get prob_observed values - these may be NULL for regression
  pob_list <- merged_results %>% purrr::map(~ prob_observed(.))
  
  # Only process pobserved if we have non-NULL values (i.e., classification results)
  pobserved <- NULL
  if (any(!sapply(pob_list, is.null))) {
    # Filter out NULL values (from regression results)
    pob_list <- pob_list[!sapply(pob_list, is.null)]
    
    if (length(pob_list) > 0) {
      pobserved <- as.data.frame(do.call(cbind, pob_list), stringsAsFactors = FALSE)
      if (is.null(colnames(pobserved)) || any(colnames(pobserved) == "")) {
        colnames(pobserved) <- paste0("class_", seq_len(ncol(pobserved)))
      }
    }
  }
  
  ind_set <- sort(unique(unlist(good_results$indices)))

  all_ids <- which(model_spec$dataset$mask > 0)
  ## if we did not get a result for all voxel ids returned results...
  mask <- if (length(ind_set) != length(all_ids)) {
    mask <- model_spec$dataset$mask
    keep <- all_ids %in% ind_set
    mask[all_ids[!keep]] <- 0
    mask
  } else {
    model_spec$dataset$mask
  }
  
  
  # Only create SparseNeuroVec for pobserved if we have classification results
  if (!is.null(pobserved)) {
    pobserved <- SparseNeuroVec(as.matrix(pobserved), neuroim2::space(mask), mask=as.logical(mask))
  }
  
  #perf_list <- furrr::future_map(merged_results, function(res) compute_performance(model_spec, res))
  perf_list <- purrr::map(merged_results, function(res) compute_performance(model_spec, res))
  
  ncols <- length(perf_list[[1]])
  pmat <- do.call(rbind, perf_list)
  
  perf_mat <- Matrix::sparseMatrix(i=rep(ind_set, ncols), j=rep(1:ncols, each=length(ind_set)), 
                                   x=as.vector(pmat), dims=c(length(model_spec$dataset$mask), ncols))
  
  
  colnames(perf_mat) <- names(perf_list[[1]])
  ret <- wrap_out(perf_mat, model_spec$dataset, ids=NULL) 
  
  # Only add pobserved if it's not NULL (i.e., for classification results)
  if (!is.null(pobserved)) {
    ret$pobserved <- pobserved
  }
  
  ret
}

#' Perform randomized searchlight analysis
#'
#' This function performs randomized searchlight analysis using a specified model, radius, and number of iterations.
#' It can be customized with different MVPA functions, combiners, and permutation options.
#'
#' @keywords internal
#' @param model_spec An object specifying the model to be used in the searchlight analysis.
#' @param radius The radius of the searchlight sphere.
#' @param niter The number of iterations for randomized searchlight.
#' @param mvpa_fun The MVPA function to be used in the searchlight analysis (default is \code{mvpa_iterate}).
#' @param combiner The function to be used to combine results (default is \code{pool_randomized}).
#' @param ... Additional arguments to be passed to the MVPA function.
#'
#' @importFrom futile.logger flog.error flog.info
#' @importFrom dplyr filter bind_rows
#' @importFrom furrr future_map
do_randomized <- function(model_spec, radius, niter, 
                         mvpa_fun=mvpa_iterate, 
                         combiner=pool_randomized, 
                         ...) {
  error=NULL
  total_models <- 0
  total_errors <- 0
  
  futile.logger::flog.info("Starting randomized searchlight analysis:")
  futile.logger::flog.info("- Radius: %s", crayon::blue(radius))
  futile.logger::flog.info("- Iterations: %s", crayon::blue(niter))
  
  ret <- purrr::map(seq(1,niter), function(i) {
    futile.logger::flog.info("\nIteration %s/%s", crayon::blue(i), crayon::blue(niter))
    slight <- get_searchlight(model_spec$dataset, "randomized", radius)
    
    ## hacky
  
    cind <- if (is.integer(slight[[1]])) {
      ## SurfaceSearchlight...
      purrr::map_int(slight, ~ attr(., "center.index"))
    } else {
      purrr::map_int(slight, ~ .@parent_index)
    }
    
    # Pass analysis_type to the mvpa function
    result <- mvpa_fun(model_spec, slight, cind, analysis_type="searchlight", ...)
    
    # Count successful and failed models
    n_success <- sum(!result$error, na.rm=TRUE)
    n_errors <- sum(result$error, na.rm=TRUE)
    total_models <<- total_models + n_success
    total_errors <<- total_errors + n_errors
    
    if (n_errors > 0) {
      futile.logger::flog.debug("- %s ROIs failed in this iteration", n_errors)
    }
    
    result
  })
  
  results <- dplyr::bind_rows(ret)
  good_results <- results %>% dplyr::filter(error == FALSE)
  bad_results <- results %>% dplyr::filter(error == TRUE)
  
  # Final summary with improved formatting
  futile.logger::flog.info("\nSearchlight analysis complete")
  futile.logger::flog.info("- Total Models Fit: %s", crayon::green(total_models))
  if (total_errors > 0) {
    futile.logger::flog.info("- Failed ROIs: %s (%s%%)", 
                            crayon::yellow(total_errors),
                            crayon::yellow(sprintf("%.1f", total_errors/(total_models + total_errors)*100)))
  } else {
    futile.logger::flog.info("- All ROIs processed successfully!")
  }
  
  if (nrow(good_results) == 0) {
    futile.logger::flog.error("No valid results for randomized searchlight")
    stop("No valid results produced")
  }
  
  combiner(model_spec, good_results)
}



#' Perform standard searchlight analysis
#'
#' This function performs standard searchlight analysis using a specified model and radius.
#' It can be customized with different MVPA functions, combiners, and permutation options.
#'
#' @keywords internal
#' @importFrom utils head
#' @param model_spec An object specifying the model to be used in the searchlight analysis.
#' @param radius The radius of the searchlight sphere.
#' @param mvpa_fun The MVPA function to be used in the searchlight analysis (default is \code{mvpa_iterate}).
#' @param combiner The function to be used to combine results (default is \code{combine_standard}).
#' @param ... Additional arguments to be passed to the MVPA function.
do_standard <- function(model_spec, radius, mvpa_fun=mvpa_iterate, combiner=combine_standard, ...) {
  error=NULL
  flog.info("creating standard searchlight")
  t_sl_create <- proc.time()[3]
  slight <- get_searchlight(model_spec$dataset, "standard", radius)
  flog.debug("get_searchlight (standard) took %.3f sec", proc.time()[3] - t_sl_create)
  
  t_iterate <- proc.time()[3]
  cind <- which(model_spec$dataset$mask > 0)
  flog.info("running standard searchlight iterator")
  ret <- mvpa_fun(model_spec, slight, cind, analysis_type="searchlight", ...)
  flog.debug("mvpa_iterate (standard searchlight) took %.3f sec",
             proc.time()[3] - t_iterate)
  good_results <- ret %>% dplyr::filter(!error)
  bad_results <- ret %>% dplyr::filter(error == TRUE)
  
  if (nrow(bad_results) > 0) {
    flog.info(bad_results$error_message)
  }
  
  if (nrow(good_results) == 0) {
    failed <- nrow(bad_results)
    flog.error(
      "no valid results for standard searchlight: %s ROIs failed.",
      failed
    )
    if (failed > 0) {
      unique_errors <- unique(bad_results$error_message)
      msg_sample <- paste(head(unique_errors, 3), collapse = " | ")
      flog.debug("sample error messages: %s", msg_sample)
    }
  }
  
  t_combine <- proc.time()[3]
  out <- combiner(model_spec, good_results, bad_results)
  flog.debug("Combiner '%s' took %.3f sec",
             deparse(substitute(combiner)),
             proc.time()[3] - t_combine)
  out
}



#' A "base" function for searchlight analysis
#'
#' This function implements the generic logic for running a searchlight:
#' \enumerate{
#'   \item Checks \code{radius} and \code{method}.
#'   \item For "standard" searchlight, calls \code{do_standard(...)}.
#'   \item For "randomized", calls \code{do_randomized(...)} with \code{niter} times.
#'   \item Handles the \code{combiner} function or string ("pool", "average").
#' }
#'
#' It does not assume any specific model type, but expects that \code{model_spec}
#' is compatible with \code{do_standard(...)} or \code{do_randomized(...)} in your code.
#'
#' @param model_spec A model specification object (e.g., \code{mvpa_model}, \code{vector_rsa_model}, etc.).
#' @param radius Numeric searchlight radius (1 to 100).
#' @param method Character: "standard" or "randomized".
#' @param niter Number of iterations if \code{method="randomized"}.
#' @param combiner Either a function that combines partial results or a string
#'        ("pool", "average") that selects a built-in combiner.
#' @param ... Additional arguments passed on to \code{do_standard} or \code{do_randomized}.
#'
#' @return The result object from \code{do_standard} or \code{do_randomized} (often a \code{searchlight_result} or similar).
#'
#' @export
run_searchlight_base <- function(model_spec,
                                 radius = 8,
                                 method = c("randomized", "standard"),
                                 niter = 4,
                                 combiner = "average",
                                 ...) {

  
  # 1) Check radius
  if (radius < 1 || radius > 100) {
    stop(paste("radius", radius, "outside allowable range (1-100)"))
  }
  
  # 2) Match method
  method <- match.arg(method)
  
  # If method is randomized, check niter
  if (method == "randomized") {
    assert_that(niter >= 1, msg = "Number of iterations for randomized searchlight must be >= 1")
  }
  
  # 3) Decide combiner if it's "pool" or "average"
  #    (In your code, you might have do_standard/do_randomized handle this logic directly - this is just an example.)
  chosen_combiner <- combiner
  if (!is.function(combiner)) {
    # Default mapping for string-based combiner argument
    if (method == "standard") {
      if (combiner == "average") { # Default for run_searchlight.default standard method
        chosen_combiner <- combine_standard
        attr(chosen_combiner, "name") <- "combine_standard"
      } else if (combiner == "standard") { # Allow explicit string name
        chosen_combiner <- combine_standard
        attr(chosen_combiner, "name") <- "combine_standard"
      } else if (combiner == "rsa_standard") { # if other models need specific string refs
        chosen_combiner <- combine_rsa_standard
        attr(chosen_combiner, "name") <- "combine_rsa_standard"
      } else if (combiner == "vector_rsa_standard") {
        chosen_combiner <- combine_vector_rsa_standard
        attr(chosen_combiner, "name") <- "combine_vector_rsa_standard"
      } else if (combiner == "msreve_standard") {
        chosen_combiner <- combine_msreve_standard
        attr(chosen_combiner, "name") <- "combine_msreve_standard"
      } else {
        stop(paste0("Unknown string combiner '", combiner, "' for method 'standard'."))
      }
    } else if (method == "randomized") {
      if (combiner == "pool") {
        chosen_combiner <- pool_randomized
        attr(chosen_combiner, "name") <- "pool_randomized"
      } else if (combiner == "average") {
        chosen_combiner <- combine_randomized
        attr(chosen_combiner, "name") <- "combine_randomized"
      } else {
        stop(paste0("Unknown string combiner '", combiner, "' for method 'randomized'."))
      }
    } else {
      stop(paste0("Unknown method '", method, "' for resolving string combiner."))
    }
  } else {
    # combiner was supplied as a function; record its name if possible
    attr(chosen_combiner, "name") <- deparse(substitute(combiner))
  }
  
  # Ensure chosen_combiner is actually a function now
  if (!is.function(chosen_combiner)){
      stop(paste0("Internal error: Combiner resolution failed. 'chosen_combiner' is not a function for combiner string: ", combiner, " and method: ", method))
  }

  # print(paste("combiner is", str(chosen_combiner))) # Original debug line
  if (getOption("rMVPA.debug", FALSE)) {
      combiner_name <- attr(chosen_combiner, "name")
      if (is.null(combiner_name)) {
          combiner_name <- deparse(substitute(combiner))
      }
      message(paste("Using combiner:", combiner_name, "for method:", method))
  }
  
  # 4) Dispatch to do_standard or do_randomized
  res <- if (method == "standard") {
    flog.info("Running standard searchlight with radius = %s", radius)
    do_standard(model_spec, radius, combiner = chosen_combiner, ...)
  } else {  # method == "randomized"
    flog.info("Running randomized searchlight with radius = %s and niter = %s", radius, niter)
    do_randomized(model_spec, radius, niter = niter, combiner = chosen_combiner, ...)
  }
  
  res
}

#' Default method for run_searchlight
#'
#' By default, if an object's class does not implement a specific 
#' \code{run_searchlight.<class>} method, this fallback will call
#' \code{run_searchlight_base}.
#'
#' @param model_spec The generic model specification object.
#' @inheritParams run_searchlight_base
#'
#' @export
run_searchlight.default <- function(model_spec, radius = 8, method = c("randomized","standard"),
                                    niter = 4, combiner = "average", ...) {
  run_searchlight_base(
    model_spec    = model_spec,
    radius        = radius,
    method        = method,
    niter         = niter,
    combiner      = combiner,
    ...
  )
}

#' run_searchlight method for vector_rsa_model
#'
#' This sets a custom \code{mvpa_fun} (e.g., \code{vector_rsa_iterate}) or 
#' different combiners for standard vs. randomized, etc.
#'
#' @param model_spec A \code{vector_rsa_model} object.
#' @inheritParams run_searchlight_base
#' @export
run_searchlight.vector_rsa <- function(model_spec,
                                       radius = 8,
                                       method = c("randomized","standard"),
                                       niter = 4,
                                       ...) {
  method <- match.arg(method)
  
  if (method == "standard") {
    flog.info("Running standard vector RSA searchlight (radius = %s)", radius)
    do_standard(model_spec, radius, mvpa_fun = vector_rsa_iterate, combiner = combine_vector_rsa_standard, ...)
  } else {
    flog.info("Running randomized vector RSA searchlight (radius = %s, niter = %s)", radius, niter)
    do_randomized(model_spec, radius, niter = niter, mvpa_fun = vector_rsa_iterate, combiner = combine_randomized, ...)
  }
}



#' Combine MS-ReVE (Contrast RSA) Searchlight Results
#'
#' This function gathers the Q-dimensional performance vectors from each successful
#' searchlight center and combines them into Q separate output maps.
#'
#' @param model_spec The \code{contrast_rsa_model} specification.
#' @param good_results A tibble containing successful results from \code{train_model.contrast_rsa_model}.
#'   Each row corresponds to a searchlight center. Expected columns include \code{id}
#'   (center voxel global index) and \code{performance} (a named numeric vector of length Q).
#' @param bad_results A tibble containing information about failed searchlights (for error reporting).
#'
#' @return A \code{searchlight_result} object containing:
#'   \item{results}{A named list of \code{SparseNeuroVec} or \code{NeuroSurfaceVector} objects,
#'     one for each contrast (Q maps in total).}
#'   \item{...}{Other standard searchlight metadata.}
#' @keywords internal
#' @importFrom neuroim2 SparseNeuroVec space add_dim
#' @importFrom neurosurf NeuroSurfaceVector geometry
#' @importFrom purrr map 
#' @importFrom dplyr bind_cols select pull
#' @export
combine_msreve_standard <- function(model_spec, good_results, bad_results) {
  output_metric <- model_spec$output_metric
  dataset <- model_spec$dataset
  output_maps <- list()
  final_metrics <- NULL
  
  # Handle cases where perf_list might contain NULLs if some searchlights failed
  # but good_results still had rows (e.g. error occurred after performance was NULL)
  # This shouldn't happen if mvpa_iterate filters errors correctly, but as a safeguard:
  valid_perf_list <- Filter(Negate(is.null), good_results$performance)
  if (length(valid_perf_list) == 0 && length(good_results$performance) > 0) {
      stop("All performance results in perf_list are NULL, cannot combine.")
  } else if (length(valid_perf_list) == 0 && length(good_results$performance) == 0) {
      # This case is already handled by the nrow(good_results) == 0 check earlier
      # but kept for logical completeness if perf_list was somehow empty independently.
      stop("Performance results are empty, cannot combine.")
  }
  # Use the filtered list for further processing
  # We also need to filter center_ids to align with valid_perf_list if some were NULL
  # However, good_results should already contain only successful iterations where performance is non-NULL.
  # Let's assume good_results$performance (aliased as perf_list) only has non-NULL valid results.

  first_perf <- valid_perf_list[[1]] # Assumes perf_list is not empty (checked by nrow(good_results))
  
  # Check if the output is a single-value metric (like recon_score or composite)
  # vs. a Q-length vector metric (like beta_delta, beta_only, delta_only, beta_delta_norm)
  if (output_metric %in% c("recon_score", "composite")) {
      # --- Handle single-value metrics --- 
      expected_metric_name <- output_metric
      if (length(first_perf) != 1 || names(first_perf)[1] != expected_metric_name) {
          stop(paste0("Expected single ", expected_metric_name, " metric but found: ", 
                      paste(names(first_perf), collapse=", ")))
      }
      
      # Combine the single values into one vector
      perf_vec <- sapply(valid_perf_list, function(x) x[[expected_metric_name]])
      if (length(perf_vec) != length(good_results$id)) {
         stop(paste0("Mismatch between number of performance scores and center IDs for ", expected_metric_name, "."))
      }
      
      # Create a single output map
      if (inherits(dataset, "mvpa_surface_dataset")) {
          output_maps[[expected_metric_name]] <- build_surface_map(dataset, perf_vec, good_results$id)
      } else {
          output_maps[[expected_metric_name]] <- build_volume_map(dataset, perf_vec, good_results$id)
      }
      final_metrics <- expected_metric_name
      
  } else {
      # --- Handle multi-value (Q-length) metrics --- 
      Q <- length(first_perf)
      contrast_names <- names(first_perf)
      if (is.null(contrast_names) || Q == 0) {
          # This case (Q=0) should ideally be caught earlier, e.g. in train_model if no contrasts
          # or if contrast_names are NULL for a Q-length vector.
          warning("Performance vectors are unnamed or have zero length. Map names will be generic or map creation might fail.")
          if (Q > 0) contrast_names <- paste0("contrast_", 1:Q) else contrast_names <- character(0)
      }
      
      # Ensure all performance vectors in the list have the same length Q
      all_lengths_match <- all(sapply(valid_perf_list, length) == Q)
      if (!all_lengths_match) {
          stop("Inconsistent number of performance metrics (contrasts) across searchlights.")
      }
      
      if (Q > 0) {
          # Bind the list of vectors into a matrix (N_voxels x Q)
          perf_mat_from_list <- do.call(rbind, valid_perf_list)
          if (!is.matrix(perf_mat_from_list)) {
              # Handle case where only one voxel succeeded, or perf_list had single vectors
              perf_mat_from_list <- matrix(valid_perf_list[[1]], nrow=length(valid_perf_list), byrow = TRUE)
          }
          # Ensure correct dimensions if only one center_id
          if (length(good_results$id) == 1 && nrow(perf_mat_from_list) == 1 && ncol(perf_mat_from_list) != Q) {
             # If rbind on single list element transposes, fix it
             if (length(valid_perf_list[[1]]) == Q) { # check if first element was Q-length
                perf_mat_from_list <- matrix(valid_perf_list[[1]], nrow=1, ncol=Q)
             }
          }
          
          if (nrow(perf_mat_from_list) != length(good_results$id) || (!is.null(ncol(perf_mat_from_list)) && ncol(perf_mat_from_list) != Q) ){
             stop(paste0("Dimension mismatch when creating performance matrix from list. Expected N_voxels x Q (", 
                        length(good_results$id), "x", Q, ") but got ", nrow(perf_mat_from_list), "x", ncol(perf_mat_from_list))) 
          }
          colnames(perf_mat_from_list) <- contrast_names
          
          # Create Q output maps
          for (q_idx in 1:Q) {
            q_contrast_name <- contrast_names[q_idx]
            q_perf_vec <- perf_mat_from_list[, q_idx]
            
            if (inherits(dataset, "mvpa_surface_dataset")) {
              output_maps[[q_contrast_name]] <- build_surface_map(dataset, q_perf_vec, good_results$id)
            } else {
              output_maps[[q_contrast_name]] <- build_volume_map(dataset, q_perf_vec, good_results$id)
            }
          }
      } # End if Q > 0
      final_metrics <- contrast_names
  }
  

  # --- Return Structure ---
  structure(
    list(
      results = output_maps, # List of Q maps
      n_voxels = estimate_mask_size(dataset), # Total voxels/vertices in mask space
      active_voxels = length(good_results$id), # Number of voxels with results
      metrics = final_metrics # Names of the contrasts/maps
    ),
    class = c("msreve_searchlight_result", "searchlight_result", "list") # Reuse existing class
  )
}
