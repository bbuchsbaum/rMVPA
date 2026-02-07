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
  # Normalize to a standard matrix type for downstream operations
  if (!is.null(perf_mat)) {
    perf_mat <- as.matrix(perf_mat)
  }

  validate_wrap_inputs(perf_mat, ids)

  if (is_perf_empty(perf_mat)) {
    return(empty_searchlight_result(dataset))
  }

  metric_names <- colnames(perf_mat)
  if (is.null(metric_names)) {
    metric_names <- paste0("Metric", seq_len(ncol(perf_mat)))
  }

  output_maps <- vector("list", length(metric_names))

  for (idx in seq_along(metric_names)) {
    metric_vector <- perf_mat[, idx]
    output_maps[[idx]] <- build_output_map(dataset, metric_vector, ids)
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
      pob_list <- if ("prob_observed" %in% names(good_results)) {
        good_results$prob_observed
      } else {
        good_results %>%
          dplyr::select(result) %>%
          dplyr::pull(result) %>%
          purrr::map(~ prob_observed(.))
      }

      if (any(!vapply(pob_list, is.null, logical(1)))) {
        if (!inherits(model_spec$dataset, "mvpa_image_dataset") ||
            searchlight_scope(model_spec$dataset) != "searchlight") {
          futile.logger::flog.debug("pobserved maps not supported for this dataset type/scope; skipping.")
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
combine_randomized <- function(model_spec, good_results, bad_results=NULL, ...) {
  futile.logger::flog.debug("combine_randomized: Starting with %d ROI results", nrow(good_results))

  # Check if we have results
  if (nrow(good_results) == 0 || length(good_results$performance) == 0) {
    futile.logger::flog.error("No valid results for randomized searchlight")
    stop("No valid results for randomized searchlight")
  }

  # Find first non-NULL performance entry to determine metric count and names
  futile.logger::flog.debug("combine_randomized: Scanning for metric prototype and names")
  perf_proto <- NULL
  metric_names <- NULL
  for (idx in seq_len(nrow(good_results))) {
    perf_candidate <- good_results$performance[[idx]]
    if (!is.null(perf_candidate)) {
      perf_proto <- perf_candidate
      if (is.matrix(perf_proto) || is.data.frame(perf_proto)) {
        ncols <- ncol(perf_proto)
        metric_names <- colnames(perf_proto)
      } else {
        ncols <- length(perf_proto)
        metric_names <- names(perf_proto)
      }
      # If we found names, we're done; otherwise keep looking
      if (!is.null(metric_names)) break
    }
  }

  if (is.null(perf_proto)) {
    futile.logger::flog.error("No valid performance entries found")
    stop("No valid performance entries found")
  }

  futile.logger::flog.debug("combine_randomized: Found %d metrics: %s",
                           ncols,
                           paste(if (is.null(metric_names)) paste0("Metric", seq_len(ncols)) else metric_names, collapse=", "))

  # Accumulate triplets for sparse matrix construction
  futile.logger::flog.debug("combine_randomized: Building triplet lists from ROI results")
  I_list <- list()
  J_list <- list()
  X_list <- list()
  k <- 1L
  skipped <- 0L

  # Process each result
  for (i in seq_len(nrow(good_results))) {
    ind_i   <- good_results$indices[[i]]
    perf_i  <- good_results$performance[[i]]
    if (is.null(ind_i) || is.null(perf_i)) {
      skipped <- skipped + 1L
      next
    }

    # Coerce performance entry to a numeric vector of length ncols
    perf_vec <- if (is.matrix(perf_i) || is.data.frame(perf_i)) {
      as.numeric(perf_i[1, , drop = FALSE])
    } else {
      as.numeric(perf_i)
    }
    if (length(perf_vec) != ncols) {
      futile.logger::flog.warn("combine_randomized: performance length mismatch in row %d", i)
      skipped <- skipped + 1L
      next
    }

    tryCatch({
      len <- length(ind_i)
      # Build triplets: each voxel gets all metric values
      I_list[[k]] <- rep(ind_i, each = ncols)
      J_list[[k]] <- rep.int(seq_len(ncols), times = len)
      X_list[[k]] <- rep(perf_vec, times = len)
      k <- k + 1L
    }, error = function(e) {
      futile.logger::flog.warn("Error processing result %d: %s", i, e$message)
      skipped <- skipped + 1L
    })
  }

  futile.logger::flog.debug("combine_randomized: Processed %d valid ROIs, skipped %d", k - 1L, skipped)

  # Flatten triplet lists
  futile.logger::flog.debug("combine_randomized: Flattening triplet lists")
  I <- unlist(I_list, use.names = FALSE)
  J <- unlist(J_list, use.names = FALSE)
  X <- unlist(X_list, use.names = FALSE)

  futile.logger::flog.debug("combine_randomized: Total triplets: %d (%.1f MB)",
                           length(I),
                           (length(I) * 8 * 3) / 1024^2)

  # Build sparse matrices for sums and counts in one pass
  # sparseMatrix automatically sums duplicate (i,j) entries
  futile.logger::flog.debug("combine_randomized: Constructing sparse values matrix (%d x %d)",
                           length(model_spec$dataset$mask), ncols)
  vals_mat <- Matrix::sparseMatrix(
    i = I,
    j = J,
    x = X,
    dims = c(length(model_spec$dataset$mask), ncols)
  )

  futile.logger::flog.debug("combine_randomized: Constructing sparse counts matrix")
  counts_mat <- Matrix::sparseMatrix(
    i = I,
    j = J,
    x = 1,
    dims = c(length(model_spec$dataset$mask), ncols)
  )

  # Normalize by counts (avoid division by zero)
  futile.logger::flog.debug("combine_randomized: Normalizing by counts")
  perf_mat <- vals_mat / pmax(counts_mat, 1)

  futile.logger::flog.debug("combine_randomized: Result sparsity: %.2f%% (non-zero: %d)",
                           100 * Matrix::nnzero(perf_mat) / prod(dim(perf_mat)),
                           Matrix::nnzero(perf_mat))

  # Set column names from the performance metrics (fallback to generic names if needed)
  if (is.null(metric_names)) {
    metric_names <- paste0("Metric", seq_len(ncols))
  }
  colnames(perf_mat) <- metric_names

  # Wrap and return results
  futile.logger::flog.debug("combine_randomized: Wrapping output")
  ret <- wrap_out(perf_mat, model_spec$dataset)
  futile.logger::flog.debug("combine_randomized: Complete")
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
pool_randomized <- function(model_spec,
                            good_results,
                            bad_results = NULL,
                            chunk_size = NULL,
                            return_pobserved = TRUE,
                            ...) {
  if (nrow(good_results) == 0) {
    stop("searchlight: no searchlight samples produced valid results")
  }
  # Heuristic: aim for ~20 chunks, bounded to avoid extremes.
  # This keeps temporary merged objects small without paying excessive loop overhead.
  nvox <- length(unique(unlist(good_results$indices)))
  chunk_size <- if (is.null(chunk_size) || is.na(chunk_size) || chunk_size <= 0) {
    as.integer(max(1000L, min(50000L, ceiling(nvox / 20))))
  } else {
    max(1L, as.integer(chunk_size))
  }

  # Build voxel -> ROI map (keeps deterministic ordering)
  indmap <- do.call(rbind, lapply(seq_len(nrow(good_results)), function(i) {
    cbind(i, good_results$indices[[i]])
  }))
  respsets <- split(indmap[, 1], indmap[, 2])
  ind_set <- as.integer(names(respsets))

  # Prototype for metrics and optional pobserved
  proto_merge <- do_merge_results(respsets[[1]], good_results)
  proto_perf  <- compute_performance(model_spec, proto_merge)
  perf_names  <- names(proto_perf)
  ncols       <- length(proto_perf)
  if (is.null(perf_names)) {
    perf_names <- paste0("Metric", seq_len(ncols))
  }

  proto_pobs <- if (return_pobserved) prob_observed(proto_merge) else NULL
  keep_pobs  <- return_pobserved && !is.null(proto_pobs)
  n_trials   <- if (keep_pobs) length(proto_pobs) else 0L

  nvox <- length(ind_set)
  trip_len <- nvox * ncols
  I_perf <- integer(trip_len)
  J_perf <- integer(trip_len)
  X_perf <- numeric(trip_len)

  if (keep_pobs) {
    pobs_mat <- matrix(NA_real_, nrow = n_trials, ncol = nvox)
  }

  fill_voxel_block <- function(block_inds, block_offset) {
    for (i in seq_along(block_inds)) {
      voxel_id <- block_inds[i]
      rset     <- respsets[[as.character(voxel_id)]]
      merged   <- do_merge_results(rset, good_results)

      # performance
      perf_vec <- compute_performance(model_spec, merged)
      perf_num <- as.numeric(perf_vec)
      if (length(perf_num) != ncols) {
        futile.logger::flog.warn(
          "pool_randomized: performance length mismatch for voxel %s (expected %s, got %s)",
          voxel_id, ncols, length(perf_num)
        )
        perf_num <- rep(NA_real_, ncols)
      }
      idx_range <- ((block_offset + i - 1L) * ncols + 1L):((block_offset + i) * ncols)
      I_perf[idx_range] <<- voxel_id
      J_perf[idx_range] <<- seq_len(ncols)
      X_perf[idx_range] <<- perf_num

      # pobserved
      if (keep_pobs) {
        pobs_vec <- prob_observed(merged)
        if (is.null(pobs_vec)) {
          pobs_vec <- rep(NA_real_, n_trials)
        } else if (length(pobs_vec) != n_trials) {
          futile.logger::flog.warn(
            "pool_randomized: prob_observed length mismatch for voxel %s (expected %s, got %s)",
            voxel_id, n_trials, length(pobs_vec)
          )
          pobs_vec <- rep(NA_real_, n_trials)
        }
        pobs_mat[, block_offset + i] <<- pobs_vec
      }

      rm(merged)
    }
    invisible()
  }

  # Stream through voxels in chunks
  voxel_indices <- seq_len(nvox)
  block_starts  <- seq(1L, nvox, by = chunk_size)
  for (bs in block_starts) {
    be <- min(bs + chunk_size - 1L, nvox)
    fill_voxel_block(ind_set[bs:be], bs - 1L)
    # Light GC hint; avoid forced collection
    if ((be - bs + 1L) * ncols > 1e6) gc(FALSE)
  }

  perf_mat <- Matrix::sparseMatrix(
    i    = I_perf,
    j    = J_perf,
    x    = X_perf,
    dims = c(length(model_spec$dataset$mask), ncols)
  )
  colnames(perf_mat) <- perf_names

  all_ids <- which(model_spec$dataset$mask > 0)
  mask <- if (length(ind_set) != length(all_ids)) {
    mask <- model_spec$dataset$mask
    keep <- all_ids %in% ind_set
    mask[all_ids[!keep]] <- 0
    mask
  } else {
    model_spec$dataset$mask
  }

  ret <- wrap_out(perf_mat, model_spec$dataset, ids = NULL)

  if (keep_pobs) {
    if (is.null(colnames(pobs_mat)) || any(colnames(pobs_mat) == "")) {
      colnames(pobs_mat) <- paste0("class_", seq_len(ncol(pobs_mat)))
    }
    ret$pobserved <- SparseNeuroVec(as.matrix(pobs_mat), neuroim2::space(mask), mask = as.logical(mask))
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
                         ...,
                         chunk_size = NULL,
                         return_pobserved = TRUE,
                         drop_probs = FALSE,
                         fail_fast = FALSE) {
  error=NULL
  total_models <- 0
  total_errors <- 0
  
  if (drop_probs && identical(combiner, pool_randomized)) {
    stop("drop_probs = TRUE is incompatible with combiner=pool_randomized (probs are required for merging). Choose combiner='average' or FALSE.")
  }
  
  futile.logger::flog.info("Starting randomized searchlight analysis:")
  futile.logger::flog.info("- Radius: %s", crayon::blue(radius))
  futile.logger::flog.info("- Iterations: %s", crayon::blue(niter))
  
  ret <- purrr::map(seq(1,niter), function(i) {
    futile.logger::flog.info("\nIteration %s/%s", crayon::blue(i), crayon::blue(niter))
    futile.logger::flog.debug("do_randomized iter %d: Generating searchlight with radius %d", i, radius)
    slight <- get_searchlight(model_spec$dataset, "randomized", radius)

    ## hacky
    futile.logger::flog.debug("do_randomized iter %d: Got %d ROIs, extracting center indices", i, length(slight))

    cind <- if (is.integer(slight[[1]])) {
      ## SurfaceSearchlight...
      purrr::map_int(slight, ~ attr(., "center.index"))
    } else {
      purrr::map_int(slight, ~ .@parent_index)
    }

    # Pass analysis_type to the mvpa function
    futile.logger::flog.debug("do_randomized iter %d: Calling mvpa_fun with %d ROIs", i, length(slight))

    result <- tryCatch({
      mvpa_fun(model_spec, slight, cind, analysis_type="searchlight",
               drop_probs = drop_probs, fail_fast = fail_fast, ...)
    }, error = function(e) {
      futile.logger::flog.error("do_randomized iter %d: mvpa_fun threw error: %s", i, e$message)
      # Return a tibble with all errors to maintain structure
      tibble::tibble(
        result = rep(list(NULL), length(slight)),
        indices = rep(list(NULL), length(slight)),
        performance = rep(list(NULL), length(slight)),
        id = seq_along(slight),
        error = TRUE,
        error_message = sprintf("Iteration failed: %s", e$message)
      )
    })

    # Defensive check: ensure result is a valid data frame with error column
    if (is.null(result) || !is.data.frame(result) || !"error" %in% names(result)) {
      has_err_col <- if (!is.null(result)) "error" %in% names(result) else FALSE
      futile.logger::flog.error("do_randomized iter %d: mvpa_fun returned invalid result (type: %s, has_error_col: %s)",
                               i, class(result)[1], has_err_col)
      futile.logger::flog.error("This usually indicates an interrupted process or unexpected failure in mvpa_iterate")
      # Instead of aborting, wrap into an all-error row so the run can finish and report
      raw_desc <- tryCatch(
        paste(utils::capture.output(utils::str(result, vec.len = 10, max.level = 1)), collapse = "\n"),
        error = function(...) NA_character_
      )
      msg <- sprintf("Invalid mvpa_fun result (iter %d): type=%s, has_error_col=%s\nraw:%s",
                     i, class(result)[1], has_err_col, raw_desc)
      result <- tibble::tibble(
        result = list(NULL),
        indices = list(NULL),
        performance = list(NULL),
        id = seq_along(cind),
        error = TRUE,
        error_message = msg
      )
    }

    futile.logger::flog.debug("do_randomized iter %d: mvpa_fun returned %d results", i, nrow(result))

    # Count successful and failed models
    n_success <- sum(!result$error, na.rm=TRUE)
    n_errors <- sum(result$error, na.rm=TRUE)
    total_models <<- total_models + n_success
    total_errors <<- total_errors + n_errors

    if (n_errors > 0) {
      futile.logger::flog.debug("- %s ROIs failed in this iteration", n_errors)
    }

    futile.logger::flog.debug("do_randomized iter %d: Complete (success=%d, errors=%d)", i, n_success, n_errors)
    result
  })
  
  futile.logger::flog.debug("do_randomized: Combining %d iteration results", length(ret))
  results <- dplyr::bind_rows(ret)
  futile.logger::flog.debug("do_randomized: Combined into %d total results", nrow(results))
  good_results <- results %>% dplyr::filter(error == FALSE)
  bad_results <- results %>% dplyr::filter(error == TRUE)
  futile.logger::flog.debug("do_randomized: Split into %d good, %d bad results", nrow(good_results), nrow(bad_results))

  # Final summary with improved formatting
  futile.logger::flog.info("\nSearchlight analysis complete")
  futile.logger::flog.info("- Total Models Fit: %s", crayon::green(total_models))
  if (total_errors > 0) {
    futile.logger::flog.info("- Failed ROIs: %s (%s%%)", 
                            crayon::yellow(total_errors),
                            crayon::yellow(sprintf("%.1f", total_errors/(total_models + total_errors)*100)))
    # Surface the most common failure reasons
    top_errs <- sort(table(bad_results$error_message), decreasing = TRUE)
    if (length(top_errs) > 0) {
      top_errs <- head(top_errs, 3L)
      futile.logger::flog.info("- Top failure causes:")
      for (msg in names(top_errs)) {
        futile.logger::flog.info("    - %s (n=%s)", msg, top_errs[[msg]])
      }
    }
  } else {
    futile.logger::flog.info("- All ROIs processed successfully!")
  }
  
  if (nrow(good_results) == 0) {
    futile.logger::flog.warn(
      "do_randomized: all randomized ROIs failed across %d iteration(s); returning NA maps sized to mask.",
      niter
    )
    mask_ids <- which(model_spec$dataset$mask > 0)
    if (length(mask_ids) == 0) {
      stop("do_randomized: all ROIs failed and mask has zero active voxels.")
    }
    metric_name <- "NA_metric"
    na_vec <- rep(NA_real_, length(mask_ids))
    na_map <- if (inherits(model_spec$dataset, "mvpa_surface_dataset")) {
      build_surface_map(model_spec$dataset, na_vec, mask_ids)
    } else {
      build_volume_map(model_spec$dataset, na_vec, mask_ids)
    }
    return(structure(
      list(
        results = setNames(list(na_map), metric_name),
        n_voxels = estimate_mask_size(model_spec$dataset),
        active_voxels = length(mask_ids),
        metrics = metric_name
      ),
      class = c("searchlight_result", "list"),
      bad_results = bad_results
    ))
  }

  futile.logger::flog.debug("do_randomized: Calling combiner function with %d good results", nrow(good_results))
  result <- combiner(model_spec, good_results, bad_results,
                     chunk_size = chunk_size,
                     return_pobserved = return_pobserved)
  futile.logger::flog.debug("do_randomized: Combiner complete, returning results")
  attr(result, "bad_results") <- bad_results
  result
}

#' Perform resampled searchlight analysis
#'
#' Similar to randomized searchlight but uses neuroim2::resampled_searchlight,
#' which draws a fixed total number of centers (`niter`), optionally with a
#' vector of radii. Results are aggregated per voxel using the randomized combiner.
#'
#' @keywords internal
#' @param model_spec An object specifying the model to be used in the searchlight analysis.
#' @param radius The radius (or vector of radii) of the searchlight sphere.
#' @param niter Total number of sampled searchlights.
#' @param mvpa_fun The MVPA function to be used in the searchlight analysis (default is \code{mvpa_iterate}).
#' @param combiner The function to be used to combine results (default is \code{combine_randomized}).
#' @param ... Additional arguments to be passed to the MVPA function.
#' @param drop_probs Logical; drop per-ROI probability matrices after computing metrics (default \code{FALSE}).
#' @param return_pobserved Logical; placeholder for API symmetry with randomized searchlight. Currently ignored because \code{combine_randomized} does not aggregate prob-observed.
do_resampled <- function(model_spec, radius, niter,
                        mvpa_fun = mvpa_iterate,
                        combiner = combine_randomized,
                        ...,
                        drop_probs = FALSE,
                        return_pobserved = FALSE,
                        fail_fast = FALSE) {
  futile.logger::flog.info("Starting resampled searchlight analysis:")
  futile.logger::flog.info("- Radius: %s", crayon::blue(paste(radius, collapse = ", ")))
  futile.logger::flog.info("- Samples: %s", crayon::blue(niter))
  
  if (drop_probs && identical(combiner, pool_randomized)) {
    stop("drop_probs = TRUE is incompatible with combiner=pool_randomized (probs are required for merging). Choose combiner='average' or set drop_probs=FALSE.")
  }

  slight <- get_searchlight(model_spec$dataset, type = "resampled", radius = radius, iter = niter)
  futile.logger::flog.debug("do_resampled: Got %d ROIs, extracting center indices", length(slight))

  cind <- if (is.integer(slight[[1]])) {
    purrr::map_int(slight, ~ attr(., "center.index"))
  } else {
    purrr::map_int(slight, ~ .@parent_index)
  }

  result <- tryCatch({
    mvpa_fun(model_spec, slight, cind, analysis_type = "searchlight",
             drop_probs = drop_probs, fail_fast = fail_fast, ...)
  }, error = function(e) {
    futile.logger::flog.error("do_resampled: mvpa_fun threw error: %s", e$message)
    tibble::tibble(
      result = rep(list(NULL), length(slight)),
      indices = rep(list(NULL), length(slight)),
      performance = rep(list(NULL), length(slight)),
      id = seq_along(slight),
      error = TRUE,
      error_message = sprintf("Iteration failed: %s", e$message)
    )
  })

  if (is.null(result) || !is.data.frame(result) || !"error" %in% names(result)) {
    has_err_col <- if (!is.null(result)) "error" %in% names(result) else FALSE
    futile.logger::flog.error("do_resampled: mvpa_fun returned invalid result (type: %s, has_error_col: %s)",
                             class(result)[1], has_err_col)
    # Wrap into an error-only tibble to keep pipeline alive and surface the issue
    raw_desc <- tryCatch(
      paste(utils::capture.output(utils::str(result, vec.len = 10, max.level = 1)), collapse = "\n"),
      error = function(...) NA_character_
    )
    msg <- sprintf("Invalid mvpa_fun result (resampled): type=%s, has_error_col=%s\nraw:%s",
                   class(result)[1], has_err_col, raw_desc)
    result <- tibble::tibble(
      result = rep(list(NULL), length(cind)),
      indices = rep(list(NULL), length(cind)),
      performance = rep(list(NULL), length(cind)),
      id = seq_along(cind),
      error = TRUE,
      error_message = msg
    )
  }

  n_success <- sum(!result$error, na.rm = TRUE)
  n_errors  <- sum(result$error, na.rm = TRUE)
  futile.logger::flog.debug("do_resampled: mvpa_fun returned %d results (success=%d, errors=%d)",
                            nrow(result), n_success, n_errors)

  good_results <- result %>% dplyr::filter(error == FALSE)
  bad_results  <- result %>% dplyr::filter(error == TRUE)

  futile.logger::flog.info("\nSearchlight analysis complete")
  futile.logger::flog.info("- Total Models Fit: %s", crayon::green(n_success))
  if (n_errors > 0) {
    futile.logger::flog.info("- Failed ROIs: %s (%s%%)",
                             crayon::yellow(n_errors),
                             crayon::yellow(sprintf("%.1f", n_errors/(n_success + n_errors)*100)))
  } else {
    futile.logger::flog.info("- All ROIs processed successfully!")
  }

  if (nrow(good_results) == 0) {
    futile.logger::flog.error("No valid results for resampled searchlight")
    stop("No valid results produced")
  }

  futile.logger::flog.debug("do_resampled: Calling combiner function with %d good results", nrow(good_results))
  res <- combiner(model_spec, good_results, bad_results)
  futile.logger::flog.debug("do_resampled: Combiner complete, returning results")
  attr(res, "bad_results") <- bad_results
  res
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
do_standard <- function(model_spec, radius, mvpa_fun=mvpa_iterate, combiner=combine_standard, ..., k = NULL, drop_probs = FALSE, fail_fast = FALSE) {
  error=NULL
  flog.info("creating standard searchlight")
  t_sl_create <- proc.time()[3]
  slight <- get_searchlight(model_spec$dataset, "standard", radius, k = k)
  flog.debug("get_searchlight (standard) took %.3f sec", proc.time()[3] - t_sl_create)

  t_iterate <- proc.time()[3]
  cind <- get_center_ids(model_spec$dataset)
  atype <- searchlight_scope(model_spec$dataset)
  flog.info("running standard searchlight iterator")
  ret <- mvpa_fun(model_spec, slight, cind, analysis_type=atype,
                  drop_probs = drop_probs, fail_fast = fail_fast, ...)
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
  attr(out, "bad_results") <- bad_results
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
                                 method = c("standard", "randomized", "resampled"),
                                 niter = 4,
                                 combiner = "average",
                                 drop_probs = FALSE,
                                 fail_fast = FALSE,
                                 k = NULL,
                                 ...) {

  
  # 1) Check radius (allow vectors for resampled)
  if (length(radius) == 1L) {
    if (radius < 1 || radius > 100) {
      stop(paste("radius", radius, "outside allowable range (1-100)"))
    }
  } else {
    if (any(radius < 1) || any(radius > 100)) {
      stop("All radii must be within allowable range (1-100)")
    }
  }
  
  # 2) Match method
  method <- match.arg(method)
  
  # If method is randomized, check niter
  if (method %in% c("randomized", "resampled")) {
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
    } else if (method %in% c("randomized", "resampled")) {
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
    do_standard(model_spec, radius, combiner = chosen_combiner, k = k,
                drop_probs = drop_probs, fail_fast = fail_fast, ...)
  } else if (method == "randomized") {
    flog.info("Running randomized searchlight with radius = %s and niter = %s", radius, niter)
    do_randomized(model_spec, radius, niter = niter, combiner = chosen_combiner,
                  drop_probs = drop_probs, fail_fast = fail_fast, ...)
  } else { # resampled
    flog.info("Running resampled searchlight with radius = %s and samples = %s", radius, niter)
    do_resampled(model_spec, radius, niter = niter, combiner = chosen_combiner,
                 drop_probs = drop_probs, fail_fast = fail_fast, ...)
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
run_searchlight.default <- function(model_spec, radius = 8, method = c("standard","randomized","resampled"),
                                    niter = 4, combiner = "average", drop_probs = FALSE,
                                    fail_fast = FALSE, k = NULL, ...) {
  run_searchlight_base(
    model_spec    = model_spec,
    radius        = radius,
    method        = method,
    niter         = niter,
    combiner      = combiner,
    drop_probs    = drop_probs,
    fail_fast     = fail_fast,
    k             = k,
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
                                       method = c("randomized","standard","resampled"),
                                       niter = 4,
                                       drop_probs = FALSE,
                                       fail_fast = FALSE,
                                       ...) {
  method <- match.arg(method)
  
  if (method == "standard") {
    flog.info("Running standard vector RSA searchlight (radius = %s)", radius)
    do_standard(model_spec, radius, mvpa_fun = vector_rsa_iterate, combiner = combine_vector_rsa_standard,
                drop_probs = drop_probs, fail_fast = fail_fast, ...)
  } else if (method == "randomized") {
    flog.info("Running randomized vector RSA searchlight (radius = %s, niter = %s)", radius, niter)
    do_randomized(model_spec, radius, niter = niter, mvpa_fun = vector_rsa_iterate,
                  combiner = combine_randomized, drop_probs = drop_probs, fail_fast = fail_fast, ...)
  } else {
    flog.info("Running resampled vector RSA searchlight (radius = %s, samples = %s)", radius, niter)
    do_resampled(model_spec, radius, niter = niter, mvpa_fun = vector_rsa_iterate,
                 combiner = combine_randomized, drop_probs = drop_probs, fail_fast = fail_fast, ...)
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
