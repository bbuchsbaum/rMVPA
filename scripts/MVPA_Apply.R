#!/usr/bin/env Rscript
###############################################################################
# MVPA Apply Script
#
# Apply saved model fits from a prior MVPA_Regional.R run to new data.
#
# How it works:
# - Loads saved config and pre-trained model fits from a prior run
# - Optionally overrides test_data/test_design with --new_data/--new_design
# - Applies the saved fits to predict on the new dataset (NO re-training!)
# - Results are written in the same structure (maps + tables)
#
# Requirements:
# - Prior run must have used --save_predictors flag to save fits
# - Fits are loaded from model_dir/fits/roi_*.rds
#
# Enhancements:
# - Filter specific ROIs with --roinum
# - Subset test data with --test_subset (R expression)
# - Merge design variables into predictions with --coalesce_design_vars
###############################################################################

suppressPackageStartupMessages({
  library(neuroim2)
  library(rMVPA)
  library(optparse)
  library(futile.logger)
  library(io)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# ------------------------- CLI options ---------------------------------------
option_list <- list(
  make_option(c("-f","--model_dir"), type = "character",
              help = "Path to folder containing config.yaml and fits/ subdirectory"),
  make_option(c("-n","--new_data"), type = "character",
              help = "New 4D NIfTI (or surface) dataset to predict on"),
  make_option(c("-d","--new_design"), type = "character",
              help = "Design table for the new dataset (optional, for evaluation)"),
  make_option(c("-o","--output"), type = "character",
              help = "Output folder for predictions (defaults to model_dir + '/apply')"),
  make_option(c("--roinum"), type = "character",
              help = "Comma-separated list of ROI IDs to apply (default: all ROIs)"),
  make_option(c("--test_subset"), type = "character",
              help = "Quoted R expression to subset test_design (e.g., 'run == 1')"),
  make_option(c("--coalesce_design_vars"), action = "store_true", default = FALSE,
              help = "Merge design variables into prediction table")
)

oparser <- OptionParser(usage = "MVPA_Apply.R [options]\n\nApply saved model fits to new data without re-training.",
                        option_list = option_list)
opt <- parse_args(oparser, positional_arguments = TRUE)
args <- opt$options

if (is.null(args$model_dir)) {
  stop("--model_dir is required (must contain config.yaml and fits/)")
}

cfgfile <- file.path(args$model_dir, "config.yaml")
if (!file.exists(cfgfile)) stop("Could not find config.yaml in ", args$model_dir)

fits_dir <- file.path(args$model_dir, "fits")
if (!dir.exists(fits_dir)) {
  stop("No fits/ directory found in ", args$model_dir,
       "\nMake sure the original run used --save_predictors=TRUE")
}

base_cfg <- io::qread(cfgfile)
config <- as.environment(base_cfg)

# Resolve output
outdir <- args$output %||% file.path(args$model_dir, "apply")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Override test data/design if provided
if (!is.null(args$new_data))   config$test_data   <- args$new_data
if (!is.null(args$new_design)) config$test_design <- args$new_design

# Safety: ensure we have test data to predict on
if (is.null(config$test_data)) {
  stop("No test data specified. Use --new_data or ensure config has test_data.")
}

# Normalize flag compatibility
if (!is.null(config$normalize)) config$normalize_samples <- isTRUE(config$normalize)

###############################################################################
# Load Saved Fits
###############################################################################

flog.info("Loading saved model fits from: %s", fits_dir)

fit_files <- list.files(fits_dir, pattern = "^roi_.*\\.rds$", full.names = TRUE)
if (length(fit_files) == 0) {
  stop("No fit files found in ", fits_dir)
}

# Extract ROI numbers from filenames and sort
roi_nums <- as.integer(sub("^roi_(\\d+)\\.rds$", "\\1", basename(fit_files)))
fit_order <- order(roi_nums)
fit_files <- fit_files[fit_order]
roi_nums <- roi_nums[fit_order]

# Filter ROIs if --roinum specified
if (!is.null(args$roinum)) {
  requested_rois <- as.integer(strsplit(args$roinum, ",")[[1]])
  keep_idx <- which(roi_nums %in% requested_rois)

  if (length(keep_idx) == 0) {
    stop("None of the requested ROIs (", args$roinum, ") found in saved fits")
  }

  fit_files <- fit_files[keep_idx]
  roi_nums <- roi_nums[keep_idx]
  flog.info("Filtered to %d ROIs: %s", length(roi_nums), paste(roi_nums, collapse = ", "))
}

# Load all fits
fits <- lapply(fit_files, readRDS)
names(fits) <- sprintf("roi_%03d", roi_nums)

flog.info("Loaded %d ROI fits", length(fits))

###############################################################################
# Load Test Data
###############################################################################

flog.info("Loading test data: %s", config$test_data)

# Load mask
mask_vol <- if (identical(config$data_mode, "image")) {
  as(rMVPA:::load_mask(config), "LogicalNeuroVol")
} else {
  NULL
}

# Load test dataset (override train_data temporarily to load test)
test_dataset <- if (identical(config$data_mode, "image")) {
  test_cfg <- config
  test_cfg$train_data <- config$test_data
  rMVPA:::initialize_image_data(test_cfg, mask_vol)
} else if (identical(config$data_mode, "surface")) {
  test_cfg <- config
  test_cfg$train_data <- config$test_data
  rMVPA:::initialize_surface_data(test_cfg)[[1]]
} else {
  stop("Unsupported data_mode in config: ", config$data_mode)
}

flog.info("Test data loaded: %d observations", nobs(test_dataset))

# Load test design for labeling/evaluation
test_design <- if (!is.null(config$test_design)) {
  test_cfg <- config
  test_cfg$train_design <- config$test_design

  # Handle test_subset if provided
  if (!is.null(args$test_subset)) {
    test_cfg$test_subset <- args$test_subset
  }

  design_obj <- rMVPA:::initialize_design(test_cfg)

  # Apply subset if present
  if (!is.null(args$test_subset)) {
    subset_expr <- parse(text = args$test_subset)
    subset_idx <- eval(subset_expr, envir = design_obj$design_table)

    if (!is.logical(subset_idx)) {
      stop("--test_subset must evaluate to a logical vector")
    }

    flog.info("Applying test_subset: '%s' (keeping %d/%d observations)",
              args$test_subset, sum(subset_idx), length(subset_idx))

    # Subset design
    design_obj$design_table <- design_obj$design_table[subset_idx, , drop = FALSE]
    design_obj$y_train <- design_obj$y_train[subset_idx]
    if (!is.null(design_obj$block_var)) {
      design_obj$block_var <- design_obj$block_var[subset_idx]
    }

    # Subset test data
    if (inherits(test_dataset$train_data, "NeuroVec")) {
      test_dataset$train_data <- neuroim2::sub_vec(test_dataset$train_data, subset_idx)
    } else {
      test_dataset$train_data <- test_dataset$train_data[subset_idx, , drop = FALSE]
    }
  }

  design_obj
} else {
  NULL
}

###############################################################################
# Detect Special Models (REMAP/naive_xdec)
###############################################################################

# Check if any fit requires special process_roi handling
needs_process_roi <- FALSE
if (length(fits) > 0) {
  first_fit <- fits[[1]]
  # REMAP and naive_xdec models have custom process_roi methods
  if (inherits(first_fit, c("remap_rrr_model", "naive_xdec_model", "custom_model"))) {
    needs_process_roi <- TRUE
    flog.warn("Detected special model type requiring process_roi: %s", class(first_fit)[1])
    flog.warn("Current implementation uses generic predict(). Results may differ from run_regional().")
  }
}

###############################################################################
# Apply Fits to New Data
###############################################################################

flog.info("Applying saved fits to new data...")

# Extract region mask
region_mask <- if (identical(config$data_mode, "image")) {
  rMVPA:::load_mask(config)
} else {
  test_dataset$mask
}

# Get unique ROI IDs from mask
roi_ids_mask <- sort(unique(as.vector(region_mask)))
roi_ids_mask <- roi_ids_mask[roi_ids_mask != 0]

# Use the ROI IDs from saved fits (more reliable)
roi_ids <- roi_nums

# Prepare results storage
predictions_by_roi <- list()
performance_by_roi <- list()

for (i in seq_along(fits)) {
  roi_id <- roi_ids[i]
  fit <- fits[[i]]

  flog.info("Processing ROI %d (%d/%d)", roi_id, i, length(fits))

  # Extract ROI data
  roi_mask <- region_mask == roi_id

  if (sum(roi_mask) == 0) {
    warning("ROI ", roi_id, " not found in mask, skipping")
    next
  }

  if (inherits(test_dataset$train_data, "NeuroVec")) {
    roi_data <- neuroim2::series(test_dataset$train_data, which(roi_mask))
  } else {
    # Matrix-based
    roi_indices <- which(roi_mask)
    roi_data <- test_dataset$train_data[, roi_indices, drop = FALSE]
  }

  # Apply normalization if needed
  if (isTRUE(config$normalize_samples)) {
    roi_data <- scale(roi_data)
  }

  # Predict using saved fit
  preds <- tryCatch({
    predict(fit, roi_data)
  }, error = function(e) {
    warning("Failed to predict for ROI ", roi_id, ": ", e$message)
    NULL
  })

  if (!is.null(preds)) {
    predictions_by_roi[[i]] <- preds

    # If we have test_design, compute performance
    if (!is.null(test_design)) {
      obs <- test_design$y_train
      # Compute performance (simplified, matches common metrics)
      if (is.factor(obs) || is.character(obs)) {
        # Classification
        pred_class <- if (inherits(preds, "classification_result")) {
          predicted_class(preds)
        } else if (is.list(preds) && !is.null(preds$class)) {
          preds$class
        } else {
          as.character(preds)
        }
        acc <- mean(pred_class == obs, na.rm = TRUE)
        performance_by_roi[[i]] <- data.frame(
          roinum = roi_id,
          accuracy = acc,
          n_obs = length(obs),
          stringsAsFactors = FALSE
        )
      } else {
        # Regression
        pred_vals <- if (is.numeric(preds)) preds else as.numeric(preds)
        rmse <- sqrt(mean((pred_vals - obs)^2, na.rm = TRUE))
        cor_val <- cor(pred_vals, obs, use = "complete.obs")
        performance_by_roi[[i]] <- data.frame(
          roinum = roi_id,
          rmse = rmse,
          correlation = cor_val,
          n_obs = length(obs),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

###############################################################################
# Write Results
###############################################################################

flog.info("Writing results to: %s", outdir)

# Combine predictions if available
if (length(predictions_by_roi) > 0) {
  # Write prediction table
  pred_table <- do.call(rbind, lapply(seq_along(predictions_by_roi), function(i) {
    preds <- predictions_by_roi[[i]]
    roi_id <- roi_ids[i]

    # Extract predictions into tabular form
    if (inherits(preds, "classification_result")) {
      base_df <- data.frame(
        .rownum = seq_along(preds$observed),
        roinum = roi_id,
        observed = preds$observed,
        predicted = predicted_class(preds),
        correct = predicted_class(preds) == preds$observed,
        stringsAsFactors = FALSE
      )

      # Add pobserved if probs available
      if (!is.null(preds$probs)) {
        base_df$pobserved <- sapply(seq_along(preds$observed), function(j) {
          preds$probs[j, as.character(preds$observed[j])]
        })
        # Add prob_* columns
        base_df <- cbind(base_df, preds$probs)
      }

      base_df
    } else if (is.list(preds) && !is.null(preds$class)) {
      base_df <- data.frame(
        .rownum = seq_len(length(preds$class)),
        roinum = roi_id,
        predicted = preds$class,
        stringsAsFactors = FALSE
      )
      # Add probabilities if available
      if (!is.null(preds$prob)) {
        base_df <- cbind(base_df, preds$prob)
      }
      base_df
    } else {
      data.frame(
        .rownum = seq_along(preds),
        roinum = roi_id,
        predicted = as.numeric(preds),
        stringsAsFactors = FALSE
      )
    }
  }))

  # Coalesce design vars if requested
  if (isTRUE(args$coalesce_design_vars) && !is.null(test_design)) {
    flog.info("Merging design variables into prediction table")

    design_table <- test_design$design_table
    design_table$.rownum <- seq_len(nrow(design_table))

    # Left join design vars into predictions
    pred_table <- merge(pred_table, design_table, by = ".rownum", all.x = TRUE, sort = FALSE)
  }

  utils::write.table(pred_table, file = file.path(outdir, "prediction_table.txt"),
                     row.names = FALSE, quote = FALSE, sep = "\t")
  flog.info("Wrote prediction table: %s", file.path(outdir, "prediction_table.txt"))
}

# Write performance table if available
if (length(performance_by_roi) > 0) {
  perf_table <- do.call(rbind, performance_by_roi)
  utils::write.table(perf_table, file = file.path(outdir, "performance_table.txt"),
                     row.names = FALSE, quote = FALSE, sep = "\t")
  flog.info("Wrote performance table: %s", file.path(outdir, "performance_table.txt"))
}

# Save applied config snapshot
config_applied <- as.list(config)
config_applied$apply_roinum <- args$roinum
config_applied$apply_test_subset <- args$test_subset
config_applied$apply_coalesce_design_vars <- args$coalesce_design_vars
io::qwrite(config_applied, file.path(outdir, "config_applied.yaml"))

cat("\n========================================\n")
cat("MVPA Apply complete!\n")
cat("Results written to: ", outdir, "\n")
cat("Predictions made using saved fits (no re-training)\n")
if (!is.null(args$roinum)) {
  cat("Filtered to ROIs: ", args$roinum, "\n")
}
if (!is.null(args$test_subset)) {
  cat("Applied test_subset: ", args$test_subset, "\n")
}
if (isTRUE(args$coalesce_design_vars)) {
  cat("Merged design variables into predictions\n")
}
cat("========================================\n")

