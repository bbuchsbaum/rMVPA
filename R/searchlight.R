#' Wrap output results
#'
#' This function wraps the output results of the performance matrix into a list
#' of SparseNeuroVec objects for each column in the performance matrix.
#'
#' @keywords internal
#' @param perf_mat A performance matrix containing classifier results.
#' @param dataset A dataset object containing the dataset information.
#' @param ids An optional vector of voxel IDs.
#' @return A named list of SparseNeuroVec objects representing the wrapped output results.
wrap_out <- function(perf_mat, dataset, ids=NULL) {
  out <- lapply(1:ncol(perf_mat), function(i) create_searchlight_performance(dataset, perf_mat[,i], ids))
  names(out) <- colnames(perf_mat)
  
  # Add class and metadata
  structure(
    list(
      results = out,
      n_voxels = length(dataset$mask),
      active_voxels = sum(dataset$mask > 0),
      metrics = colnames(perf_mat)
    ),
    class = c("searchlight_result", "list")
  )
}

#' @export
#' @method print searchlight_result
print.searchlight_result <- function(x, ...) {
  # Ensure crayon is available
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty printing. Please install it.")
  }
  
  # Define color scheme
  header_style <- crayon::bold$cyan
  section_style <- crayon::yellow
  info_style <- crayon::white
  number_style <- crayon::green
  metric_style <- crayon::magenta
  
  # Print header
  cat("\n", header_style("█▀▀ Searchlight Analysis Results ▀▀█"), "\n\n")
  
  # Basic information
  cat(section_style("├─ Coverage"), "\n")
  cat(info_style("│  ├─ Total Voxels: "), number_style(format(x$n_voxels, big.mark=",")), "\n")
  cat(info_style("│  └─ Active Voxels: "), number_style(format(x$active_voxels, big.mark=",")), "\n")
  
  # Performance metrics
  cat(section_style("└─ Performance Metrics"), "\n")
  for (metric in x$metrics) {
    results <- x$results[[metric]]
    if (inherits(results, "searchlight_performance")) {
      cat(info_style("   ├─ "), metric_style(metric), "\n")
      cat(info_style("   │  ├─ Mean: "), number_style(sprintf("%.4f", results$summary_stats$mean)), "\n")
      cat(info_style("   │  ├─ SD: "), number_style(sprintf("%.4f", results$summary_stats$sd)), "\n")
      cat(info_style("   │  ├─ Min: "), number_style(sprintf("%.4f", results$summary_stats$min)), "\n")
      cat(info_style("   │  └─ Max: "), number_style(sprintf("%.4f", results$summary_stats$max)), "\n")
    }
  }
  
  if (!is.null(x$pobserved)) {
    cat(section_style("\n└─ Observed Probabilities"), "\n")
    # Add probability summary if needed
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
    
    has_results <- any(unlist(purrr::map(good_results$result, function(x) !is.null(x))))
    ret <- wrap_out(perf_mat, model_spec$dataset, ind)
    
    if (has_results) {
      pobserved <- good_results %>% 
        dplyr::select(result) %>% 
        pull(result) %>% 
        purrr::map(~ prob_observed(.)) %>% 
        bind_cols()
      
      # Create appropriate type of output based on dataset type
      if (inherits(model_spec$dataset, "mvpa_surface_dataset")) {
        # For surface data
        pobserved <- neurosurf::NeuroSurfaceVector(
          geometry = geometry(model_spec$dataset$train_data),
          indices = seq_len(nrow(pobserved)),  # Or appropriate indices
          mat = as.matrix(pobserved)
        )
      } else {
        # For volume data
        pobserved <- SparseNeuroVec(
          as.matrix(pobserved), 
          space(model_spec$dataset$mask), 
          mask=as.logical(model_spec$dataset$mask)
        )
      }
      
      ret$pobserved <- pobserved
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
  
  all_ind <- sort(unlist(good_results$indices))
  ind_count <- table(all_ind)
  ind_set <- unique(all_ind)
  
  # Get the number of performance metrics
  ncols <- length(good_results$performance[[1]])
  
  # Create sparse matrix to hold results
  perf_mat <- Matrix::sparseMatrix(i=rep(ind_set, ncols), j=rep(1:ncols, each=length(ind_set)), 
                                  x=rep(0, length(ind_set)*ncols), 
                                  dims=c(length(model_spec$dataset$mask), ncols))
  
  # Process each result
  for (i in 1:nrow(good_results)) {
    ind <- good_results$indices[[i]]
    if (!is.null(ind) && !is.null(good_results$performance[[i]])) {
      # Use tryCatch to handle errors gracefully
      tryCatch({
        m <- kronecker(matrix(good_results$performance[[i]], 1, ncols), rep(1, length(ind)))
        perf_mat[ind,] <- perf_mat[ind,] + m
      }, error = function(e) {
        futile.logger::flog.warn("Error processing result %d: %s", i, e$message)
      })
    }
  }

  # Normalize by the count of overlapping searchlights
  perf_mat[ind_set,] <- sweep(perf_mat[ind_set,,drop=FALSE], 1, as.integer(ind_count), FUN="/")
  
  # Set column names from the performance metrics
  colnames(perf_mat) <- names(good_results$performance[[1]])
  
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
  check <- sapply(reslist, function(res) inherits(res, "data.frame")) 
  assertthat::assert_that(all(check), msg="pool_results: all arguments must be of type 'data.frame'")
  good_results <- do.call(rbind, reslist)
 
  ## the sorted vector of all voxel indices
  all_ind <- sort(unlist(good_results$indices))
  ## how many instances of each voxel?
  ind_count <- table(all_ind)
  ind_set <- unique(all_ind)
  
  ## map every result to the set of indices in that set
  indmap <- do.call(rbind, lapply(1:nrow(good_results), function(i) {
    ind <- good_results$indices[[i]]
    cbind(i, ind)
  }))
  
  
  respsets <- split(indmap[,1], indmap[,2])
  
  merged_results <- purrr::map(respsets, do_merge_results, good_results=good_results)
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
  pobserved <- merged_results %>% purrr::map( ~ prob_observed(.)) %>% bind_cols()
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
  
  
  pobserved <- SparseNeuroVec(as.matrix(pobserved), neuroim2::space(mask), mask=as.logical(mask))
  
  #perf_list <- furrr::future_map(merged_results, function(res) compute_performance(model_spec, res))
  perf_list <- purrr::map(merged_results, function(res) compute_performance(model_spec, res))
  
  ncols <- length(perf_list[[1]])
  pmat <- do.call(rbind, perf_list)
  
  perf_mat <- Matrix::sparseMatrix(i=rep(ind_set, ncols), j=rep(1:ncols, each=length(ind_set)), 
                                   x=as.vector(pmat), dims=c(length(model_spec$dataset$mask), ncols))
  
  
  colnames(perf_mat) <- names(perf_list[[1]])
  ret <- wrap_out(perf_mat, model_spec$dataset, ids=NULL) 
  ret$pobserved <- pobserved
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
  
  futile.logger::flog.info("🔄 Starting randomized searchlight analysis:")
  futile.logger::flog.info("├─ Radius: %s", crayon::blue(radius))
  futile.logger::flog.info("└─ Iterations: %s", crayon::blue(niter))
  
  ret <- purrr::map(seq(1,niter), function(i) {
    futile.logger::flog.info("\n📊 Iteration %s/%s", crayon::blue(i), crayon::blue(niter))
    slight <- get_searchlight(model_spec$dataset, "randomized", radius)
    
    ## hacky
  
    cind <- if (is.integer(slight[[1]])) {
      ## SurfaceSearchlight...
      purrr::map_int(slight, ~ attr(., "center.index"))
    } else {
      purrr::map_int(slight, ~ .@parent_index)
    }
    
    result <- mvpa_fun(model_spec, slight, cind, ...)
    
    # Count successful and failed models
    n_success <- sum(!result$error, na.rm=TRUE)
    n_errors <- sum(result$error, na.rm=TRUE)
    total_models <<- total_models + n_success
    total_errors <<- total_errors + n_errors
    
    if (n_errors > 0) {
      futile.logger::flog.debug("└─ %s ROIs failed in this iteration", n_errors)
    }
    
    result
  })
  
  results <- dplyr::bind_rows(ret)
  good_results <- results %>% dplyr::filter(error == FALSE)
  bad_results <- results %>% dplyr::filter(error == TRUE)
  
  # Final summary with improved formatting
  futile.logger::flog.info("\n✨ Searchlight Analysis Complete")
  futile.logger::flog.info("├─ Total Models Fit: %s", crayon::green(total_models))
  if (total_errors > 0) {
    futile.logger::flog.info("└─ Failed ROIs: %s (%s%%)", 
                            crayon::yellow(total_errors),
                            crayon::yellow(sprintf("%.1f", total_errors/(total_models + total_errors)*100)))
  } else {
    futile.logger::flog.info("└─ All ROIs processed successfully!")
  }
  
  if (nrow(good_results) == 0) {
    futile.logger::flog.error("❌ No valid results for randomized searchlight")
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
#' @param model_spec An object specifying the model to be used in the searchlight analysis.
#' @param radius The radius of the searchlight sphere.
#' @param mvpa_fun The MVPA function to be used in the searchlight analysis (default is \code{mvpa_iterate}).
#' @param combiner The function to be used to combine results (default is \code{combine_standard}).
#' @param ... Additional arguments to be passed to the MVPA function.
do_standard <- function(model_spec, radius, mvpa_fun=mvpa_iterate, combiner=combine_standard, ...) {
  error=NULL
  flog.info("creating standard searchlight")
  slight <- get_searchlight(model_spec$dataset, "standard", radius)
  
   
  cind <- which(model_spec$dataset$mask > 0)
  flog.info("running standard searchlight iterator")
  ret <- mvpa_fun(model_spec, slight, cind, ...)
  good_results <- ret %>% dplyr::filter(!error)
  bad_results <- ret %>% dplyr::filter(error == TRUE)
  
  if (nrow(bad_results) > 0) {
    flog.info(bad_results$error_message)
  }
  
  if (nrow(good_results) == 0) {
    ## TODO print out some debug information
    flog.error("no valid results for standard searchlight, exiting.")
  }
  
  combiner(model_spec, good_results, bad_results)
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
  #    (In your code, you might have do_standard/do_randomized handle this logic directly—this is just an example.)
  chosen_combiner <- combiner
  if (!is.function(combiner)) {
    if (combiner == "pool") {
      chosen_combiner <- pool_randomized        # or combine_standard, depending on your code
    } else if (combiner == "average") {
      chosen_combiner <- combine_randomized     # or combine_standard, depending on your code
    } else {
      stop("'combiner' must be 'average', 'pool', or a custom function.")
    }
  }
  
  print(paste("combiner is", str(chosen_combiner)))
  
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


#' Create a searchlight performance object
#'
#' @keywords internal
#' @param dataset The dataset object
#' @param perf_vec Performance vector for a single metric
#' @param ids Optional vector of voxel IDs
#' @return A searchlight_performance object
create_searchlight_performance <- function(dataset, perf_vec, ids=NULL) {
  # First use the S3 wrap_output method to create the NeuroVol
  ret <- wrap_output(dataset, perf_vec, ids)
  
  # Get non-zero and non-NA values for statistics
  vals <- perf_vec[perf_vec != 0 & !is.na(perf_vec)]
  
  # Then wrap it in our searchlight_performance structure
  structure(
    list(
      data = ret,
      metric_name = names(perf_vec)[1],
      n_nonzero = sum(perf_vec != 0, na.rm = TRUE),
      summary_stats = list(
        mean = if(length(vals) > 0) mean(vals, na.rm = TRUE) else NA,
        sd = if(length(vals) > 0) sd(vals, na.rm = TRUE) else NA,
        min = if(length(vals) > 0) min(vals, na.rm = TRUE) else NA,
        max = if(length(vals) > 0) max(vals, na.rm = TRUE) else NA
      ),
      indices = ids  # Store the indices for reference
    ),
    class = c("searchlight_performance", "list")
  )
}

#' @export
#' @method print searchlight_performance
print.searchlight_performance <- function(x, ...) {
  # Ensure crayon is available
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty printing. Please install it.")
  }
  
  # Define color scheme
  header_style <- crayon::bold$cyan
  section_style <- crayon::yellow
  info_style <- crayon::white
  number_style <- crayon::green
  metric_style <- crayon::magenta
  
  # Print header
  cat("\n", header_style("█▀▀ Searchlight Performance: "), 
      metric_style(x$metric_name), header_style(" ▀▀█"), "\n\n")
  
  # Data information
  cat(section_style("├─ Data Summary"), "\n")
  cat(info_style("│  ├─ Non-zero Values: "), 
      if(is.null(x$n_nonzero)) crayon::red("NULL") else number_style(format(x$n_nonzero, big.mark=",")), 
      "\n")
  
  # Statistics
  cat(section_style("└─ Statistics"), "\n")
  
  # Helper function to format stats with better NULL handling
  format_stat <- function(val) {
    if (is.null(val) || (length(val) == 0)) {
      crayon::red("No data")
    } else if (is.na(val)) {
      crayon::red("No valid data")
    } else {
      number_style(sprintf("%.4f", val))
    }
  }
  
  # Safely extract stats with NULL checking
  stats <- x$summary_stats
  if (is.null(stats)) {
    stats <- list(mean=NULL, sd=NULL, min=NULL, max=NULL)
  }
  
  cat(info_style("   ├─ Mean: "), format_stat(stats$mean), "\n")
  cat(info_style("   ├─ SD: "), format_stat(stats$sd), "\n")
  cat(info_style("   ├─ Min: "), format_stat(stats$min), "\n")
  cat(info_style("   └─ Max: "), format_stat(stats$max), "\n\n")
}




