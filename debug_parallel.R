#!/usr/bin/env Rscript

# Debug script to understand parallel processing visibility
library(rMVPA)
library(future)
library(furrr)

# Function to show current future plan details
show_future_info <- function() {
  cat("\n=== Future Plan Information ===\n")
  p <- plan()
  cat("Current plan:", class(p)[1], "\n")
  
  # Get worker info
  if (inherits(p, "multicore")) {
    cat("Workers:", nbrOfWorkers(), "\n")
    cat("Available cores:", availableCores(), "\n")
  }
  
  # Check if we're actually running in parallel
  cat("Sequential:", inherits(plan(), "sequential"), "\n")
  cat("Parallel:", !inherits(plan(), "sequential"), "\n")
  
  # Show future options
  cat("\nFuture options:\n")
  cat("- future.globals.maxSize:", getOption("future.globals.maxSize"), "\n")
  cat("- future.rng.onMisuse:", getOption("future.rng.onMisuse"), "\n")
  
  # Check furrr settings
  cat("\nFurrr chunk settings:\n")
  cat("- future.chunk.size:", getOption("future.chunk.size"), "\n")
  cat("- future.scheduling:", getOption("future.scheduling"), "\n")
}

# Function to monitor parallel execution
test_parallel_visibility <- function(n_items, batch_size, workers = 4) {
  cat("\n=== Testing with", n_items, "items, batch_size =", batch_size, "===\n")
  
  # Set up parallel plan
  plan(multicore, workers = workers)
  show_future_info()
  
  # Calculate batches
  n_batches <- ceiling(n_items / batch_size)
  items_per_batch <- ceiling(n_items / n_batches)
  
  cat("\nBatch calculation:\n")
  cat("- Total items:", n_items, "\n")
  cat("- Batch size:", batch_size, "\n")
  cat("- Number of batches:", n_batches, "\n")
  cat("- Items per batch:", items_per_batch, "\n")
  
  # Simulate what mvpa_iterate does
  cat("\nSimulating batch processing...\n")
  
  for (i in 1:n_batches) {
    start_idx <- (i-1) * items_per_batch + 1
    end_idx <- min(i * items_per_batch, n_items)
    batch_items <- end_idx - start_idx + 1
    
    cat("\nBatch", i, ":", batch_items, "items\n")
    
    # This simulates what happens in run_future.default
    start_time <- Sys.time()
    
    # Create a frame similar to what mvpa_iterate creates
    frame <- data.frame(
      id = start_idx:end_idx,
      item = paste0("item_", start_idx:end_idx)
    )
    
    # Run parallel processing like furrr::future_pmap does
    results <- future_map(1:nrow(frame), function(idx) {
      # Add a small delay to make workers visible
      Sys.sleep(0.1)
      
      # Return process ID to see which worker handled this
      list(
        item = frame$item[idx],
        pid = Sys.getpid(),
        worker = Sys.getpid()
      )
    }, .options = furrr_options(
      seed = TRUE,
      chunk_size = NULL  # Let furrr decide chunking
    ))
    
    end_time <- Sys.time()
    
    # Count unique workers
    pids <- unique(sapply(results, function(x) x$pid))
    cat("  Unique process IDs seen:", length(pids), "\n")
    cat("  Time taken:", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n")
    
    # If only one PID, workers might not be spawning
    if (length(pids) == 1) {
      cat("  ⚠️  Only one process detected - might be running sequentially!\n")
    } else {
      cat("  ✓", length(pids), "parallel processes detected\n")
    }
  }
}

# Test different scenarios
cat("========================================\n")
cat("Testing Parallel Processing Visibility\n")
cat("========================================\n")

# Scenario 1: Large batch (might not spawn workers)
test_parallel_visibility(n_items = 1000, batch_size = 900, workers = 4)

# Scenario 2: Small batch (should spawn workers)
test_parallel_visibility(n_items = 1000, batch_size = 100, workers = 4)

# Scenario 3: Very small batch
test_parallel_visibility(n_items = 1000, batch_size = 20, workers = 4)

cat("\n========================================\n")
cat("Debugging Tips:\n")
cat("========================================\n")
cat("1. If you don't see multiple R processes in 'top', check:\n")
cat("   - Batch size vs number of items per batch\n")
cat("   - Whether futures are actually being created\n")
cat("   - If the work is too small to justify spawning workers\n")
cat("\n2. Use 'htop' instead of 'top' for better visibility\n")
cat("\n3. Check with: ps aux | grep R | grep -v grep\n")
cat("\n4. Monitor with: watch -n 0.5 'ps aux | grep R | wc -l'\n")