#!/usr/bin/env Rscript

# Debug script specifically for searchlight parallelism
library(rMVPA)
library(future)

cat("========================================\n")
cat("Searchlight Parallelism Debugging\n")
cat("========================================\n\n")

# Function to monitor system processes
monitor_r_processes <- function() {
  cmd <- "ps aux | grep '[R] ' | grep -v grep | wc -l"
  as.numeric(system(cmd, intern = TRUE))
}

# Generate realistic dataset
dataset <- gen_sample_dataset(
  D = c(20, 20, 20),  # Realistic brain volume size
  nobs = 100,
  response_type = "categorical",
  data_mode = "image",
  blocks = 3,
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

# Function to test with different configurations
test_config <- function(workers, batch_size = NULL) {
  cat("\n----------------------------------------\n")
  cat("Configuration:\n")
  cat("  Workers:", workers, "\n")
  cat("  Batch size:", ifelse(is.null(batch_size), "default (10%)", batch_size), "\n")
  
  # Set up parallel plan
  if (workers > 1) {
    plan(multicore, workers = workers)
  } else {
    plan(sequential)
  }
  
  # Count searchlights
  mask <- dataset$dataset$mask
  n_searchlights <- sum(mask > 0)
  default_batch <- as.integer(0.1 * n_searchlights)
  actual_batch <- ifelse(is.null(batch_size), default_batch, batch_size)
  n_batches <- ceiling(n_searchlights / actual_batch)
  
  cat("  Total searchlights:", n_searchlights, "\n")
  cat("  Actual batch size:", actual_batch, "\n")
  cat("  Number of batches:", n_batches, "\n")
  cat("  Searchlights per batch:", ceiling(n_searchlights / n_batches), "\n")
  
  # Check if batches are small enough to not spawn workers
  if (actual_batch < workers * 2) {
    cat("  ⚠️  WARNING: Batch size (", actual_batch, ") may be too small\n")
    cat("              for", workers, "workers (need at least", workers * 2, "items)\n")
  }
  
  # Monitor processes during execution
  cat("\nMonitoring R processes during execution:\n")
  initial_procs <- monitor_r_processes()
  cat("  Initial R processes:", initial_procs, "\n")
  
  # Start searchlight in background
  cat("  Starting searchlight...\n")
  
  # We'll capture output to reduce noise
  capture.output({
    if (is.null(batch_size)) {
      results <- run_searchlight(
        mspec,
        radius = 4,
        method = "standard"
      )
    } else {
      results <- run_searchlight(
        mspec,
        radius = 4,
        method = "standard",
        batch_size = batch_size
      )
    }
  })
  
  cat("  ✓ Searchlight completed\n")
}

# Test different scenarios
cat("\nScenario 1: Default batch size with 4 workers\n")
test_config(workers = 4, batch_size = NULL)

cat("\nScenario 2: Small batch size (300) with 4 workers\n")
test_config(workers = 4, batch_size = 300)

cat("\nScenario 3: Very small batch size (50) with 4 workers\n")
test_config(workers = 4, batch_size = 50)

cat("\nScenario 4: Large batch size (2000) with 4 workers\n")
test_config(workers = 4, batch_size = 2000)

cat("\nScenario 5: Sequential (no parallelism)\n")
test_config(workers = 1, batch_size = 300)

cat("\n========================================\n")
cat("Key Insights:\n")
cat("========================================\n")
cat("1. Worker processes are forked on-demand for each batch\n")
cat("2. They terminate after processing their batch\n")
cat("3. Small batches = more frequent forking (visible in top)\n")
cat("4. Large batches = less frequent forking (harder to see)\n")
cat("5. If batch_size < 2*workers, parallelism may not engage\n")
cat("\nTo monitor in real-time, run in another terminal:\n")
cat("  watch -n 0.1 'ps aux | grep \"[R] \" | wc -l'\n")
cat("Or use htop and filter for R processes\n")