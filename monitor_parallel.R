#!/usr/bin/env Rscript

# Real-time monitoring of parallel execution
library(rMVPA)
library(future)

cat("========================================\n")
cat("Real-time Parallel Process Monitoring\n")
cat("========================================\n\n")

# Create a function that logs when it's called and by which process
logged_process_roi <- function(obj, roi, rnum, center_global_id = NA) {
  pid <- Sys.getpid()
  
  # Write to a temp file to track process activity
  log_file <- file.path(tempdir(), paste0("worker_", pid, ".log"))
  cat(format(Sys.time(), "%H:%M:%OS3"), "- ROI", rnum, "\n", 
      file = log_file, append = TRUE)
  
  # Call the original function
  rMVPA:::process_roi(obj, roi, rnum, center_global_id)
}

# Generate dataset
dataset <- gen_sample_dataset(
  D = c(10, 10, 10),
  nobs = 50,
  response_type = "categorical",
  data_mode = "image",
  blocks = 2,
  nlevels = 2
)

cval <- blocked_cross_validation(dataset$design$block_var)
model <- load_model("sda_notune")
mspec <- mvpa_model(
  model = model,
  dataset = dataset$dataset,
  design = dataset$design,
  model_type = "classification",
  crossval = cval
)

# Function to analyze worker patterns
analyze_workers <- function(batch_size) {
  cat("\n========================================\n")
  cat("Testing with batch_size =", batch_size, "\n")
  cat("========================================\n")
  
  # Clean up old logs
  temp_files <- list.files(tempdir(), pattern = "worker_.*\\.log", full.names = TRUE)
  if (length(temp_files) > 0) file.remove(temp_files)
  
  # Set up parallel processing
  plan(multicore, workers = 4)
  
  # Count searchlights
  n_searchlights <- sum(mspec$dataset$mask > 0)
  actual_batch <- ifelse(is.null(batch_size), as.integer(0.1 * n_searchlights), batch_size)
  n_batches <- ceiling(n_searchlights / actual_batch)
  
  cat("\nSetup:\n")
  cat("  Total searchlights:", n_searchlights, "\n")
  cat("  Batch size:", actual_batch, "\n")
  cat("  Number of batches:", n_batches, "\n")
  cat("  Workers requested:", 4, "\n")
  
  # Check the theoretical parallelism
  items_per_batch <- ceiling(n_searchlights / n_batches)
  theoretical_chunks <- ceiling(items_per_batch / 4)  # furrr default chunking
  
  cat("\nTheoretical analysis:\n")
  cat("  Items per batch:", items_per_batch, "\n")
  cat("  Items per worker (if evenly divided):", items_per_batch / 4, "\n")
  
  if (items_per_batch < 4) {
    cat("  ⚠️  WARNING: Not enough items per batch for 4 workers!\n")
    cat("              Only", items_per_batch, "items but 4 workers requested\n")
  }
  
  # Run searchlight with custom processor
  cat("\nRunning searchlight...\n")
  
  start_time <- Sys.time()
  
  # We need to trace the actual execution
  # Let's patch mvpa_iterate temporarily
  orig_mvpa_iterate <- rMVPA:::mvpa_iterate
  
  traced_mvpa_iterate <- function(mod_spec, vox_list, ids, ...) {
    cat("  mvpa_iterate called with", length(ids), "searchlights\n")
    
    # Get batch_size from ...
    args <- list(...)
    bs <- args$batch_size
    if (is.null(bs)) bs <- as.integer(0.1 * length(ids))
    
    cat("  Using batch_size:", bs, "\n")
    cat("  Will create", ceiling(length(ids) / bs), "batches\n")
    
    # Call original
    orig_mvpa_iterate(mod_spec, vox_list, ids, ...)
  }
  
  # Temporarily replace
  assignInNamespace("mvpa_iterate", traced_mvpa_iterate, "rMVPA")
  
  tryCatch({
    if (is.null(batch_size)) {
      results <- run_searchlight(mspec, radius = 3, method = "standard")
    } else {
      results <- run_searchlight(mspec, radius = 3, method = "standard", 
                                batch_size = batch_size)
    }
  }, finally = {
    # Restore original
    assignInNamespace("mvpa_iterate", orig_mvpa_iterate, "rMVPA")
  })
  
  end_time <- Sys.time()
  
  cat("\n✓ Completed in", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n")
  
  # Analyze worker logs
  log_files <- list.files(tempdir(), pattern = "worker_.*\\.log", full.names = TRUE)
  cat("\nWorker analysis:\n")
  cat("  Unique worker processes spawned:", length(log_files), "\n")
  
  if (length(log_files) > 0) {
    for (lf in log_files) {
      pid <- gsub(".*worker_(\\d+)\\.log", "\\1", basename(lf))
      n_rois <- length(readLines(lf))
      cat("  Worker", pid, "processed", n_rois, "ROIs\n")
    }
  }
}

# Test scenarios
analyze_workers(batch_size = NULL)    # Default 10%
analyze_workers(batch_size = 50)      # Small batches
analyze_workers(batch_size = 300)     # Medium batches
analyze_workers(batch_size = 1000)    # Large batches

cat("\n========================================\n")
cat("Summary:\n")
cat("========================================\n")
cat("Worker visibility in 'top' depends on:\n")
cat("1. Batch size - smaller batches = more frequent forking\n")
cat("2. Work per item - heavier computation = longer-lived processes\n")
cat("3. Total items - more items = longer overall runtime\n")
cat("\nWith default 10% batching:\n")
cat("- Large datasets: batches are big, workers live longer\n")
cat("- Small datasets: batches are small, workers are ephemeral\n")
cat("\nWith batch_size=300:\n")
cat("- Forces more frequent batching\n")
cat("- Workers fork/join more often\n")
cat("- More likely to catch them in 'top'\n")