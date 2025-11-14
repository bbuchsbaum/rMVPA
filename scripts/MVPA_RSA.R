#!/usr/bin/env Rscript
###############################################################################
# MVPA RSA Script
#
# Command-line entry point to run RSA analyses in regional or searchlight modes.
# Supports:
#   - Standard RSA (rsa_model + rsa_design)
#   - Feature RSA (feature_rsa_model + feature_rsa_design)
#   - Vector RSA (vector_rsa_model + vector_rsa_design)
#   - Contrast RSA / MS-ReVE (contrast_rsa_model + msreve_design)
#
# Modeled after MVPA_Regional.R and MVPA_Searchlight.R.
###############################################################################

suppressPackageStartupMessages({
  library(neuroim2)
  library(rMVPA)
  library(optparse)
  library(futile.logger)
  library(io)
})

# ------------------------- CLI options ---------------------------------------
option_list <- list(
  make_option(c("-M","--mode"), type = "character", default = "regional",
              help = "Analysis mode: 'regional' or 'searchlight' [default: %default]"),
  make_option(c("-T","--rsa_type"), type = "character", default = "standard",
              help = "RSA type: 'standard', 'feature', 'vector', or 'contrast' [default: %default]"),
  # Common data/design
  make_option(c("-t","--train_design"), type = "character", help = "Training design table"),
  make_option("--train_data", type = "character", help = "4D NIfTI for training data"),
  make_option(c("-a","--mask"), type = "character", help = "Mask image (regional: labeled ROI map; searchlight: binary mask)"),
  make_option(c("-l","--label_column"), type = "character", help = "Label column in design"),
  make_option(c("-b","--block_column"), type = "character", help = "Block/run column in design"),
  make_option(c("--data_mode"), type = "character", default = "image", help = "'image' or 'surface' [default: %default]"),
  # Searchlight knobs
  make_option(c("-r","--radius"), type = "numeric", default = 6, help = "Searchlight radius (mm) [default: %default]"),
  make_option(c("-s","--sl_type"), type = "character", default = "standard", help = "Searchlight type: 'standard' or 'randomized' [default: %default]"),
  make_option(c("-i","--niter"), type = "integer", default = 8, help = "Randomized searchlight iterations [default: %default]"),
  # Standard RSA inputs (minimal)
  make_option(c("--rsa_formula"), type = "character", default = "~ dismat", help = "RSA formula over RDM names [default: %default]"),
  make_option(c("--rdm_file"), type = "character", help = "Path to a CSV/TSV RDM (items x items) for standard or vector RSA"),
  # Feature RSA inputs
  make_option(c("--feature_matrix"), type = "character", help = "Path to CSV/TSV feature matrix (items x features)"),
  make_option(c("--similarity_matrix"), type = "character", help = "Path to CSV/TSV similarity matrix (items x items)"),
  make_option(c("--feature_rsa_method"), type = "character", default = "pls", help = "Feature RSA method: pls|pca|glmnet [default: %default]"),
  make_option(c("--max_comps"), type = "integer", default = 10, help = "Max components for feature RSA [default: %default]"),
  # Vector RSA inputs
  make_option(c("--rsa_simfun"), type = "character", default = "pearson", help = "RSA similarity fun (pearson|spearman) [default: %default]"),
  # Contrast RSA inputs
  make_option(c("--contrast_matrix"), type = "character", help = "Path to CSV/TSV contrast matrix (K x Q) for MS-ReVE"),
  make_option(c("--estimation_method"), type = "character", default = "average", help = "MS-ReVE estimation: average|L2_norm|crossnobis [default: %default]"),
  make_option(c("--regression_type"), type = "character", default = "lm", help = "RSA regression type: pearson|spearman|lm|rfit|ridge_hkb [default: %default]"),
  # Misc
  make_option(c("-o","--output"), type = "character", help = "Output folder"),
  make_option(c("-c","--config"), type = "character", help = "YAML/R config file (values override CLI defaults)"),
  make_option(c("-p","--ncores"), type = "integer", default = 1, help = "Cores (searchlight randomized may use futures) [default: %default]")
)

oparser <- OptionParser(usage = "MVPA_RSA.R [options]", option_list = option_list)
opt <- parse_args(oparser, positional_arguments = TRUE)
args <- opt$options

# ------------------------- Config merge --------------------------------------
config <- rMVPA:::initialize_configuration(args)
config <- rMVPA:::initialize_standard_parameters(config, args, "rsa")
config$ncores   <- args$ncores %||% config$ncores
config$mode     <- args$mode     %||% config$mode     %||% "regional"
config$rsa_type <- tolower(args$rsa_type %||% config$rsa_type %||% "standard")
config$data_mode<- args$data_mode%||% config$data_mode%||% "image"

# Core files
config$train_design <- args$train_design %||% config$train_design
config$train_data   <- args$train_data   %||% config$train_data
config$mask         <- args$mask         %||% config$mask
config$label_column <- args$label_column %||% config$label_column
config$block_column <- args$block_column %||% config$block_column
config$output       <- args$output       %||% config$output %||% paste0("rsa_", config$label_column %||% "analysis")

# RSA-specific
config$rsa_formula  <- args$rsa_formula  %||% config$rsa_formula %||% "~ dismat"
config$rdm_file     <- args$rdm_file     %||% config$rdm_file
config$feature_matrix <- args$feature_matrix %||% config$feature_matrix
config$similarity_matrix <- args$similarity_matrix %||% config$similarity_matrix
config$feature_rsa_method <- args$feature_rsa_method %||% config$feature_rsa_method %||% "pls"
config$max_comps    <- args$max_comps    %||% config$max_comps %||% 10L
config$rsa_simfun   <- args$rsa_simfun   %||% config$rsa_simfun %||% "pearson"
config$contrast_matrix <- args$contrast_matrix %||% config$contrast_matrix
config$estimation_method <- args$estimation_method %||% config$estimation_method %||% "average"
config$regression_type  <- args$regression_type  %||% config$regression_type  %||% "lm"

# Searchlight knobs
config$radius  <- args$radius  %||% config$radius  %||% 6
config$sl_type <- args$sl_type %||% config$sl_type %||% "standard"
config$niter   <- args$niter   %||% config$niter   %||% 8

flog.info("Configuration:")
print(as.list(config))

# ------------------------- Load design/data ----------------------------------
design_mvpa <- rMVPA:::initialize_design(config)

dataset <- if (config$data_mode == "image") {
  mask_vol <- as(rMVPA:::load_mask(config), "LogicalNeuroVol")
  rMVPA:::initialize_image_data(config, mask_vol)
} else if (config$data_mode == "surface") {
  rMVPA:::initialize_surface_data(config)[[1]]
} else {
  stop("Unknown data_mode: ", config$data_mode)
}

# ------------------------- Helper readers ------------------------------------
read_matrix_file <- function(path) {
  if (is.null(path)) return(NULL)
  if (!file.exists(path)) stop("File not found: ", path)
  mat <- as.matrix(utils::read.table(path, header = TRUE, check.names = FALSE))
  rownames(mat) <- rownames(mat) %||% as.character(seq_len(nrow(mat)))
  mat
}

# ------------------------- Build model spec ----------------------------------
rsa_type <- tolower(config$rsa_type)
model_spec <- switch(rsa_type,
  # Standard RSA
  standard = {
    mats <- list()
    dm <- read_matrix_file(config$rdm_file)
    if (is.null(dm)) stop("For standard RSA, provide --rdm_file (CSV/TSV RDM)")
    mats$dismat <- if (inherits(dm, "dist")) dm else stats::as.dist(dm)
    rdes <- rsa_design(stats::as.formula(config$rsa_formula), mats,
                       block_var = design_mvpa$block_var,
                       split_by = NULL, keep_intra_run = FALSE)
    rsa_model(dataset = dataset, design = rdes,
              distmethod = "spearman", regtype = config$regression_type)
  },
  # Feature RSA
  feature = {
    labels <- as.character(y_train(design_mvpa))
    Fm <- read_matrix_file(config$feature_matrix)
    Sm <- read_matrix_file(config$similarity_matrix)
    if (is.null(Fm) && is.null(Sm)) stop("Provide either --feature_matrix or --similarity_matrix for feature RSA")
    if (!is.null(Fm) && nrow(Fm) != length(labels)) stop("feature_matrix rows must equal number of items/labels")
    if (!is.null(Sm) && nrow(Sm) != ncol(Sm)) stop("similarity_matrix must be square (items x items)")
    fdes <- feature_rsa_design(S = Sm, F = Fm, labels = labels, k = 0, max_comps = config$max_comps,
                               block_var = design_mvpa$block_var)
    # Construct crossval if blocks exist
    crossval <- if (!is.null(design_mvpa$block_var)) {
      blocked_cross_validation(design_mvpa$block_var)
    } else {
      NULL
    }
    feature_rsa_model(dataset = dataset, design = fdes, method = config$feature_rsa_method,
                      crossval = crossval, max_comps = config$max_comps)
  },
  # Vector RSA
  vector = {
    labels <- as.character(y_train(design_mvpa))
    Dm <- read_matrix_file(config$rdm_file)
    if (is.null(Dm)) stop("--rdm_file required for vector RSA")
    vdes <- vector_rsa_design(D = Dm, labels = labels, block_var = design_mvpa$block_var)
    vector_rsa_model(dataset = dataset, design = vdes,
                     distfun = cordist(),
                     rsa_simfun = config$rsa_simfun)
  },
  # Contrast RSA / MS-ReVE
  contrast = {
    C <- read_matrix_file(config$contrast_matrix)
    if (is.null(C)) stop("--contrast_matrix required for contrast RSA (MS-ReVE)")
    msd <- msreve_design(design_mvpa, C)
    contrast_rsa_model(dataset = dataset, design = msd,
                       estimation_method = config$estimation_method,
                       regression_type  = config$regression_type)
  },
  stop("Unsupported rsa_type: ", rsa_type)
)

# ------------------------- Run analysis --------------------------------------
output_dir <- rMVPA:::make_output_dir(config$output)

if (tolower(config$mode) == "regional") {
  # Region mask: labeled volume for image mode; surface mask for surface mode
  rmask <- if (config$data_mode == "image") rMVPA:::load_mask(config) else dataset$mask
  reg_res <- run_regional(model_spec, rmask, return_fits = TRUE)
  # Write maps/tables
  rMVPA::save_results(reg_res$vol_results, dir = output_dir, level = "minimal")
  utils::write.table(reg_res$performance_table, file = file.path(output_dir, "performance_table.txt"),
                     row.names = FALSE, quote = FALSE, sep = "\t")
  if (!is.null(reg_res$prediction_table)) {
    utils::write.table(reg_res$prediction_table, file = file.path(output_dir, "prediction_table.txt"),
                       row.names = FALSE, quote = FALSE, sep = "\t")
  }
} else if (tolower(config$mode) == "searchlight") {
  sl_res <- run_searchlight(model_spec, radius = config$radius,
                            method = config$sl_type, niter = config$niter)
  rMVPA::save_results(sl_res, dir = output_dir, level = "minimal")
} else {
  stop("Unknown mode: ", config$mode)
}

# Save resolved config
cfg_out <- file.path(output_dir, "config.yaml")
io::qwrite(as.list(config), cfg_out)

cat("\nRSA analysis complete. Results saved to:\n", output_dir, "\n")

