#' Install rMVPA Command-Line Wrappers
#'
#' Copies the packaged command-line wrappers into a directory on your
#' \code{PATH}. The installed commands are \code{rmvpa-searchlight} and
#' \code{rmvpa-regional}.
#'
#' @param dest_dir Destination directory for the wrappers.
#' @param overwrite Logical; overwrite existing wrapper files if \code{TRUE}.
#' @param commands Which wrappers to install. Any subset of
#'   \code{c("searchlight", "regional")}.
#'
#' @return Invisibly, a named character vector of installed wrapper paths.
#' @export
install_cli <- function(dest_dir = "~/.local/bin",
                        overwrite = FALSE,
                        commands = c("searchlight", "regional")) {
  command_map <- .cli_command_map()
  commands <- match.arg(commands, names(command_map), several.ok = TRUE)
  commands <- unique(commands)

  exec_dir <- system.file("exec", package = "rMVPA")
  if (!nzchar(exec_dir) && requireNamespace("pkgload", quietly = TRUE)) {
    exec_dir <- tryCatch(pkgload::inst_path("exec"), error = function(e) "")
  }
  if (!nzchar(exec_dir)) {
    stop(
      "install_cli() requires an installed rMVPA package so the packaged executables can be located.",
      call. = FALSE
    )
  }

  dest_dir <- path.expand(dest_dir)
  if (!dir.exists(dest_dir) &&
      !dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create destination directory: ", dest_dir, call. = FALSE)
  }

  installed <- setNames(character(length(commands)), commands)
  for (cmd in commands) {
    src <- file.path(exec_dir, unname(command_map[[cmd]]))
    dest <- file.path(dest_dir, basename(src))

    if (!file.exists(src)) {
      stop("Packaged CLI wrapper not found: ", src, call. = FALSE)
    }
    if (file.exists(dest) && !isTRUE(overwrite)) {
      stop(
        "Refusing to overwrite existing file: ", dest,
        ". Re-run with overwrite = TRUE.",
        call. = FALSE
      )
    }
    if (!file.copy(src, dest, overwrite = overwrite, copy.mode = FALSE)) {
      stop("Failed to copy CLI wrapper to: ", dest, call. = FALSE)
    }
    Sys.chmod(dest, mode = "755")
    installed[[cmd]] <- dest
  }

  path_entries <- strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]
  if (!dest_dir %in% path_entries) {
    message(
      "Installed wrappers to ", dest_dir,
      ". Add this directory to PATH if you want to invoke them directly."
    )
  }

  invisible(installed)
}

#' @keywords internal
#' @noRd
.cli_command_map <- function() {
  c(
    searchlight = "rmvpa-searchlight",
    regional = "rmvpa-regional"
  )
}

#' @keywords internal
#' @noRd
.cli_or <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

#' @keywords internal
#' @noRd
.cli_mode_spec <- function(mode) {
  mode <- match.arg(mode, c("searchlight", "regional"))
  if (identical(mode, "searchlight")) {
    return(list(
      mode = "searchlight",
      command = .cli_command_map()[["searchlight"]],
      legacy_script = "scripts/MVPA_Searchlight.R",
      title = "Run voxelwise searchlight MVPA from the rMVPA workflow API.",
      default_output = "searchlight_results",
      mask_help = paste(
        "Binary brain mask for image analyses.",
        "For surface analyses, pass a comma-separated mask list when needed."
      ),
      backend_help = "Execution backend: default, shard, or auto. Defaults to auto."
    ))
  }

  list(
    mode = "regional",
    command = .cli_command_map()[["regional"]],
    legacy_script = "scripts/MVPA_Regional.R",
    title = "Run region-based MVPA from the rMVPA workflow API.",
    default_output = "regional_results",
    mask_help = paste(
      "Labelled ROI mask or parcellation.",
      "For surface analyses, pass a comma-separated mask list with non-zero region IDs."
    ),
    backend_help = "Execution backend: default, shard, or auto. Defaults to default."
  )
}

#' @keywords internal
#' @noRd
.cli_description <- function(mode) {
  spec <- .cli_mode_spec(mode)
  paste(
    spec$title,
    "CLI flags override config-file values.",
    "Kebab-case flags are preferred; underscore aliases are accepted for compatibility."
  )
}

#' @keywords internal
#' @noRd
.cli_epilogue <- function(mode) {
  spec <- .cli_mode_spec(mode)

  if (identical(mode, "searchlight")) {
    examples <- c(
      paste0(
        spec$command,
        " --train-design design.tsv --train-data betas.nii.gz",
        " --mask brain_mask.nii.gz --label-column condition --radius 6",
        " --model sda_notune --workers 4"
      ),
      paste0(spec$command, " --config analysis.yaml --dry-run"),
      paste0(spec$legacy_script, " --help")
    )
    notes <- c(
      "Use --workers to configure a temporary future plan for this run.",
      "Searchlight type defaults to standard; radius defaults to 8 mm if omitted.",
      "Existing non-empty output directories are rejected unless you pass --overwrite or --skip-if-exists."
    )
  } else {
    examples <- c(
      paste0(
        spec$command,
        " --train-design design.tsv --train-data betas.nii.gz",
        " --mask roi_labels.nii.gz --label-column condition",
        " --model sda_notune --save-predictors"
      ),
      paste0(spec$command, " --config regional.yaml --pool-predictions mean"),
      paste0(spec$legacy_script, " --help")
    )
    notes <- c(
      "Regional analysis requires a labelled ROI mask; non-zero values define regions.",
      "Use --save-predictors if you want fitted ROI models stored with the output bundle.",
      "Use --pool-predictions mean or stack to compute pooled ROI-level predictions."
    )
  }

  paste(
    c(
      "Examples:",
      paste0("  ", examples),
      "",
      "Notes:",
      paste0("  - ", notes)
    ),
    collapse = "\n"
  )
}

#' @keywords internal
#' @noRd
.cli_common_options <- function(mode) {
  spec <- .cli_mode_spec(mode)

  list(
    optparse::make_option(
      c("-c", "--config"),
      type = "character",
      metavar = "FILE",
      default = NULL,
      help = "YAML or R config file. CLI flags override file values."
    ),
    optparse::make_option(
      c("-t", "--train-design"),
      type = "character",
      metavar = "FILE[,FILE2...]",
      default = NULL,
      help = "Training design table(s). Multiple files can be passed as a comma-separated list."
    ),
    optparse::make_option(
      c("--test-design"),
      type = "character",
      metavar = "FILE[,FILE2...]",
      default = NULL,
      help = "Optional held-out test design table(s)."
    ),
    optparse::make_option(
      c("--train-data"),
      type = "character",
      metavar = "FILE[,FILE2...]",
      default = NULL,
      help = "Training data. For surface analyses, pass one file per section as a comma-separated list."
    ),
    optparse::make_option(
      c("--test-data"),
      type = "character",
      metavar = "FILE[,FILE2...]",
      default = NULL,
      help = "Optional held-out test data."
    ),
    optparse::make_option(
      c("-a", "--mask"),
      type = "character",
      metavar = "FILE[,FILE2...]",
      default = NULL,
      help = spec$mask_help
    ),
    optparse::make_option(
      c("-m", "--model"),
      type = "character",
      metavar = "NAME",
      default = NULL,
      help = "Model registry entry. If omitted, the workflow falls back to corsim."
    ),
    optparse::make_option(
      c("--model-type"),
      type = "character",
      metavar = "TYPE",
      default = NULL,
      help = "Problem type: classification or regression. Use regression for numeric targets."
    ),
    optparse::make_option(
      c("-l", "--label-column"),
      type = "character",
      metavar = "NAME",
      default = NULL,
      help = "Target column in the training design. Defaults to labels if omitted."
    ),
    optparse::make_option(
      c("--test-label-column"),
      type = "character",
      metavar = "NAME",
      default = NULL,
      help = "Target column in the test design. Defaults to the training label column."
    ),
    optparse::make_option(
      c("-b", "--block-column"),
      type = "character",
      metavar = "NAME",
      default = NULL,
      help = "Blocking column used to build cross-validation folds."
    ),
    optparse::make_option(
      c("--split-by"),
      type = "character",
      metavar = "EXPR",
      default = NULL,
      help = "Optional grouping variable or expression for split-wise metrics."
    ),
    optparse::make_option(
      c("--train-subset"),
      type = "character",
      metavar = "EXPR",
      default = NULL,
      help = "R expression used to subset the training design before loading data."
    ),
    optparse::make_option(
      c("--test-subset"),
      type = "character",
      metavar = "EXPR",
      default = NULL,
      help = "R expression used to subset the test design before loading data."
    ),
    optparse::make_option(
      c("-g", "--tune-grid"),
      type = "character",
      metavar = "EXPR",
      default = NULL,
      help = "Tuning-grid expression, e.g. 'lambda=c(0.1,0.5,1.0)'."
    ),
    optparse::make_option(
      c("--data-mode"),
      type = "character",
      metavar = "MODE",
      default = NULL,
      help = "Input mode: image or surface. Defaults to image."
    ),
    optparse::make_option(
      c("--link-by"),
      type = "character",
      metavar = "NAME",
      default = NULL,
      help = "Optional pairing key used by cross-decoding models such as remap_rrr."
    ),
    optparse::make_option(
      c("-n", "--normalize"),
      action = "store_true",
      default = NULL,
      help = "Center and scale each sample before fitting."
    ),
    optparse::make_option(
      c("--no-normalize"),
      action = "store_false",
      dest = "normalize",
      default = NULL,
      help = "Disable sample normalization even if the config file enables it."
    ),
    optparse::make_option(
      c("--class-metrics"),
      action = "store_true",
      default = NULL,
      help = "Write per-class metrics for multiclass problems."
    ),
    optparse::make_option(
      c("--no-class-metrics"),
      action = "store_false",
      dest = "class_metrics",
      default = NULL,
      help = "Do not write per-class metrics."
    ),
    optparse::make_option(
      c("--backend"),
      type = "character",
      metavar = "BACKEND",
      default = NULL,
      help = spec$backend_help
    ),
    optparse::make_option(
      c("-p", "--workers"),
      type = "integer",
      metavar = "N",
      default = NULL,
      help = "Number of workers for the temporary future plan. Defaults to 1."
    ),
    optparse::make_option(
      c("--future-plan"),
      type = "character",
      metavar = "PLAN",
      default = "auto",
      help = "Future strategy: auto, sequential, multisession, or multicore. [default: %default]"
    ),
    optparse::make_option(
      c("--preflight"),
      type = "character",
      metavar = "POLICY",
      default = "warn",
      help = "Preflight policy: warn, error, or off. [default: %default]"
    ),
    optparse::make_option(
      c("--seed"),
      type = "integer",
      metavar = "INT",
      default = NULL,
      help = "Set the RNG seed before building the analysis."
    ),
    optparse::make_option(
      c("-o", "--output"),
      type = "character",
      metavar = "DIR",
      default = NULL,
      help = sprintf("Output directory. Defaults to %s.", spec$default_output)
    ),
    optparse::make_option(
      c("--save-level"),
      type = "character",
      metavar = "LEVEL",
      default = "complete",
      help = "Result bundle level passed to save_results(): minimal, standard, or complete. [default: %default]"
    ),
    optparse::make_option(
      c("--overwrite"),
      action = "store_true",
      default = FALSE,
      help = "Allow writing into an existing non-empty output directory."
    ),
    optparse::make_option(
      c("--skip-if-exists"),
      action = "store_true",
      default = FALSE,
      help = "Exit successfully without running when the output directory already exists and is non-empty."
    ),
    optparse::make_option(
      c("--dry-run"),
      action = "store_true",
      default = FALSE,
      help = "Build and validate the analysis, print the resolved configuration, then exit without running."
    ),
    optparse::make_option(
      c("--print-config"),
      action = "store_true",
      default = FALSE,
      help = "Print the resolved configuration before running."
    ),
    optparse::make_option(
      c("--list-models"),
      action = "store_true",
      default = FALSE,
      help = "List available built-in models and exit."
    ),
    optparse::make_option(
      c("--example-config"),
      action = "store_true",
      default = FALSE,
      help = "Print an example YAML configuration for this command and exit."
    ),
    optparse::make_option(
      c("--version"),
      action = "store_true",
      default = FALSE,
      help = "Print the installed rMVPA version and exit."
    ),
    optparse::make_option(
      c("--verbose"),
      action = "store_true",
      default = FALSE,
      help = "Emit progress messages from the analysis runners."
    )
  )
}

#' @keywords internal
#' @noRd
.cli_mode_options <- function(mode) {
  if (identical(mode, "searchlight")) {
    return(list(
      optparse::make_option(
        c("-r", "--radius"),
        type = "double",
        metavar = "MM",
        default = NULL,
        help = "Searchlight radius in millimetres. Defaults to 8."
      ),
      optparse::make_option(
        c("-s", "--searchlight-type"),
        type = "character",
        metavar = "TYPE",
        default = NULL,
        help = "Searchlight type: standard, randomized, or resampled. Defaults to standard."
      ),
      optparse::make_option(
        c("-i", "--iterations"),
        type = "integer",
        metavar = "N",
        default = NULL,
        help = "Iteration count for randomized or resampled searchlights. Defaults to 4."
      ),
      optparse::make_option(
        c("--engine"),
        type = "character",
        metavar = "NAME",
        default = NULL,
        help = "Searchlight engine: auto, legacy, swift, or dual_lda_fast."
      ),
      optparse::make_option(
        c("--batch-size"),
        type = "integer",
        metavar = "N",
        default = NULL,
        help = "Number of searchlights per batch."
      )
    ))
  }

  list(
    optparse::make_option(
      c("--save-predictors"),
      action = "store_true",
      default = NULL,
      help = "Store fitted ROI predictors alongside the result bundle."
    ),
    optparse::make_option(
      c("--no-save-predictors"),
      action = "store_false",
      dest = "save_predictors",
      default = NULL,
      help = "Do not store fitted ROI predictors."
    ),
    optparse::make_option(
      c("--pool-predictions"),
      type = "character",
      metavar = "MODE",
      default = NULL,
      help = "Pool regional predictions using none, mean, or stack."
    ),
    optparse::make_option(
      c("--stack-folds"),
      type = "integer",
      metavar = "N",
      default = NULL,
      help = "Number of folds to use when --pool-predictions=stack."
    ),
    optparse::make_option(
      c("--stack-seed"),
      type = "integer",
      metavar = "INT",
      default = NULL,
      help = "Random seed used when auto-generating stacking folds."
    ),
    optparse::make_option(
      c("--stack-lambda"),
      type = "double",
      metavar = "NUM",
      default = NULL,
      help = "Ridge penalty used when --pool-predictions=stack. Defaults to 1e-3."
    ),
    optparse::make_option(
      c("--coalesce-design-vars"),
      action = "store_true",
      default = FALSE,
      help = "Merge design variables into the regional prediction table."
    )
  )
}

#' @keywords internal
#' @noRd
.cli_parser <- function(mode) {
  spec <- .cli_mode_spec(mode)
  optparse::OptionParser(
    usage = sprintf("%s [options]", spec$command),
    description = .cli_description(mode),
    epilogue = .cli_epilogue(mode),
    option_list = c(.cli_common_options(mode), .cli_mode_options(mode))
  )
}

#' @keywords internal
#' @noRd
.cli_split_arg <- function(x) {
  if (is.null(x) || !nzchar(x)) {
    return(NULL)
  }
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0L) NULL else parts
}

#' @keywords internal
#' @noRd
.cli_normalize_argv <- function(argv) {
  replacements <- c(
    "--train_design" = "--train-design",
    "--test_design" = "--test-design",
    "--train_data" = "--train-data",
    "--test_data" = "--test-data",
    "--label_column" = "--label-column",
    "--test_label_column" = "--test-label-column",
    "--block_column" = "--block-column",
    "--split_by" = "--split-by",
    "--train_subset" = "--train-subset",
    "--test_subset" = "--test-subset",
    "--tune_grid" = "--tune-grid",
    "--data_mode" = "--data-mode",
    "--link_by" = "--link-by",
    "--class_metrics" = "--class-metrics",
    "--save_predictors" = "--save-predictors",
    "--skip_if_folder_exists" = "--skip-if-exists",
    "--type" = "--searchlight-type",
    "--niter" = "--iterations",
    "--ncores" = "--workers",
    "--pthreads" = "--workers"
  )

  out <- argv
  for (i in seq_along(out)) {
    arg <- out[[i]]
    for (old in names(replacements)) {
      if (identical(arg, old)) {
        arg <- replacements[[old]]
        break
      }
      prefix <- paste0(old, "=")
      if (startsWith(arg, prefix)) {
        arg <- paste0(replacements[[old]], "=", substring(arg, nchar(prefix) + 1L))
        break
      }
    }
    out[[i]] <- arg
  }
  out
}

#' @keywords internal
#' @noRd
.cli_standardize_options <- function(options) {
  names(options) <- gsub("-", "_", names(options), fixed = TRUE)
  options
}

#' @keywords internal
#' @noRd
.cli_config_args <- function(mode, options) {
  args <- list(
    config = options$config,
    train_design = .cli_split_arg(options$train_design),
    test_design = .cli_split_arg(options$test_design),
    train_data = .cli_split_arg(options$train_data),
    test_data = .cli_split_arg(options$test_data),
    mask = .cli_split_arg(options$mask),
    model = options$model,
    model_type = options$model_type,
    label_column = options$label_column,
    test_label_column = options$test_label_column,
    block_column = options$block_column,
    split_by = options$split_by,
    train_subset = options$train_subset,
    test_subset = options$test_subset,
    tune_grid = options$tune_grid,
    data_mode = options$data_mode,
    link_by = options$link_by,
    normalize = options$normalize,
    class_metrics = options$class_metrics,
    backend = options$backend,
    output = options$output
  )

  if (!is.null(options$workers)) {
    args$pthreads <- as.integer(options$workers)
  }

  if (identical(mode, "searchlight")) {
    args$radius <- options$radius
    args$type <- options$searchlight_type
    args$niter <- options$iterations
    args$engine <- options$engine
    args$batch_size <- options$batch_size
  } else {
    if (!is.null(options$save_predictors)) {
      args$save_predictors <- isTRUE(options$save_predictors)
      args$return_fits <- isTRUE(options$save_predictors)
    }
  }

  args[!vapply(args, is.null, logical(1))]
}

#' @keywords internal
#' @noRd
.cli_finalize_config <- function(cfg, mode) {
  cfg <- as.list(cfg)
  class(cfg) <- c("rmvpa_config", "list")

  if (is.null(cfg$model)) {
    cfg$model <- "corsim"
  }
  if (is.null(cfg$model_type)) {
    cfg$model_type <- "classification"
  }
  if (is.null(cfg$data_mode)) {
    cfg$data_mode <- "image"
  }
  if (is.null(cfg$label_column)) {
    cfg$label_column <- "labels"
  }
  if (is.null(cfg$class_metrics)) {
    cfg$class_metrics <- TRUE
  }
  if (is.null(cfg$output)) {
    cfg$output <- .cli_mode_spec(mode)$default_output
  }
  if (!is.null(cfg$ncores) && is.null(cfg$pthreads)) {
    cfg$pthreads <- cfg$ncores
  }
  if (!is.null(cfg$searchlight_type) && is.null(cfg$type)) {
    cfg$type <- cfg$searchlight_type
  }

  if (identical(mode, "searchlight")) {
    if (is.null(cfg$radius)) {
      cfg$radius <- 8
    }
    if (is.null(cfg$type)) {
      cfg$type <- "standard"
    }
    if (is.null(cfg$niter)) {
      cfg$niter <- 4L
    }
  } else if (is.null(cfg$return_fits) && !is.null(cfg$save_predictors)) {
    cfg$return_fits <- isTRUE(cfg$save_predictors)
  }

  if (!is.null(cfg$test_design) && is.null(cfg$test_label_column)) {
    cfg$test_label_column <- cfg$label_column
  }

  cfg
}

#' @keywords internal
#' @noRd
.cli_list_models <- function() {
  registry <- get("MVPAModels", envir = asNamespace("rMVPA"))
  models <- ls(envir = registry, all.names = TRUE)
  models <- sort(models)

  info <- data.frame(
    model = models,
    type = vapply(
      models,
      function(name) {
        entry <- registry[[name]]
        if (is.null(entry$type)) "" else as.character(entry$type)
      },
      character(1)
    ),
    packages = vapply(
      models,
      function(name) {
        entry <- registry[[name]]
        libs <- entry$library
        if (is.null(libs) || length(libs) == 0L) "" else paste(libs, collapse = ",")
      },
      character(1)
    ),
    stringsAsFactors = FALSE
  )

  utils::write.table(info, row.names = FALSE, quote = FALSE, sep = "\t")
  invisible(info)
}

#' @keywords internal
#' @noRd
.cli_example_config <- function(mode) {
  if (identical(mode, "searchlight")) {
    return(paste(
      c(
        "train_design: train_design.tsv",
        "train_data: train_betas.nii.gz",
        "mask: brain_mask.nii.gz",
        "label_column: condition",
        "block_column: run",
        "model: sda_notune",
        "model_type: classification",
        "data_mode: image",
        "radius: 6",
        "type: standard",
        "niter: 4",
        "backend: auto",
        "output: searchlight_results",
        "class_metrics: true"
      ),
      collapse = "\n"
    ))
  }

  paste(
    c(
      "train_design: train_design.tsv",
      "train_data: train_betas.nii.gz",
      "mask: roi_labels.nii.gz",
      "label_column: condition",
      "block_column: run",
      "model: sda_notune",
      "model_type: classification",
      "data_mode: image",
      "backend: default",
      "output: regional_results",
      "return_fits: true",
      "pool_predictions: mean"
    ),
    collapse = "\n"
  )
}

#' @keywords internal
#' @noRd
.cli_file_fields <- function(mode) {
  base <- c("train_design", "test_design", "train_data", "test_data", "mask", "config_file")
  if (identical(mode, "searchlight")) base else base
}

#' @keywords internal
#' @noRd
.cli_validate_config <- function(cfg, mode) {
  data_mode <- match.arg(as.character(cfg$data_mode)[1], c("image", "surface"))
  cfg$data_mode <- data_mode

  if (!is.null(cfg$model_type)) {
    cfg$model_type <- match.arg(as.character(cfg$model_type)[1], c("classification", "regression"))
  }
  if (!is.null(cfg$backend)) {
    cfg$backend <- match.arg(as.character(cfg$backend)[1], c("default", "shard", "auto"))
  }

  if (identical(mode, "searchlight")) {
    cfg$type <- match.arg(as.character(cfg$type)[1], c("standard", "randomized", "resampled"))
    if (!is.numeric(cfg$radius) || length(cfg$radius) != 1L || is.na(cfg$radius) || cfg$radius <= 0) {
      stop("Searchlight radius must be a single positive number.", call. = FALSE)
    }
    if (!is.null(cfg$niter) &&
        (!is.numeric(cfg$niter) || length(cfg$niter) != 1L || is.na(cfg$niter) || cfg$niter < 1)) {
      stop("Searchlight iterations must be a positive integer.", call. = FALSE)
    }
    if (!is.null(cfg$engine)) {
      cfg$engine <- match.arg(as.character(cfg$engine)[1], c("auto", "legacy", "swift", "dual_lda_fast"))
    }
  } else {
    if (!is.null(cfg$pool_predictions)) {
      cfg$pool_predictions <- match.arg(
        as.character(cfg$pool_predictions)[1],
        c("none", "mean", "stack")
      )
    }
    if (!is.null(cfg$stack_folds) &&
        (!is.numeric(cfg$stack_folds) || length(cfg$stack_folds) != 1L || is.na(cfg$stack_folds) || cfg$stack_folds < 2)) {
      stop("stack_folds must be a single integer >= 2.", call. = FALSE)
    }
    if (!is.null(cfg$stack_lambda) &&
        (!is.numeric(cfg$stack_lambda) || length(cfg$stack_lambda) != 1L || is.na(cfg$stack_lambda) || cfg$stack_lambda <= 0)) {
      stop("stack_lambda must be a single positive number.", call. = FALSE)
    }
  }

  if (is.null(cfg$train_design)) {
    stop("Missing required input: train_design.", call. = FALSE)
  }
  if (is.null(cfg$train_data)) {
    stop("Missing required input: train_data.", call. = FALSE)
  }
  if (is.null(cfg$mask)) {
    stop("Missing required input: mask.", call. = FALSE)
  }

  has_test_design <- !is.null(cfg$test_design)
  has_test_data <- !is.null(cfg$test_data)
  if (xor(has_test_design, has_test_data)) {
    stop("test_design and test_data must be supplied together.", call. = FALSE)
  }

  for (field in .cli_file_fields(mode)) {
    value <- cfg[[field]]
    if (is.null(value) || !is.character(value)) {
      next
    }
    missing <- value[!file.exists(value)]
    if (length(missing) > 0L) {
      stop(
        sprintf("Path(s) supplied via %s do not exist: %s", field, paste(missing, collapse = ", ")),
        call. = FALSE
      )
    }
  }

  if (!is.character(cfg$output) || length(cfg$output) != 1L || !nzchar(cfg$output)) {
    stop("Output directory must be a non-empty character string.", call. = FALSE)
  }

  cfg
}

#' @keywords internal
#' @noRd
.cli_resolve_workers <- function(options, cfg) {
  workers <- .cli_or(options$workers, .cli_or(cfg$pthreads, .cli_or(cfg$ncores, 1L)))
  workers <- as.integer(workers)[1]
  if (is.na(workers) || workers < 1L) {
    stop("workers must be a positive integer.", call. = FALSE)
  }
  workers
}

#' @keywords internal
#' @noRd
.cli_prepare_output <- function(output_dir, overwrite = FALSE, skip_if_exists = FALSE) {
  if (isTRUE(overwrite) && isTRUE(skip_if_exists)) {
    stop("Use either overwrite or skip_if_exists, not both.", call. = FALSE)
  }

  if (!dir.exists(output_dir)) {
    return("run")
  }

  existing <- list.files(output_dir, all.files = TRUE, no.. = TRUE)
  if (length(existing) == 0L) {
    return("run")
  }

  if (isTRUE(skip_if_exists)) {
    return("skip")
  }

  if (!isTRUE(overwrite)) {
    stop(
      "Output directory already exists and is not empty: ", output_dir,
      ". Re-run with --overwrite or --skip-if-exists.",
      call. = FALSE
    )
  }

  "run"
}

#' @keywords internal
#' @noRd
.cli_apply_future_plan <- function(workers, future_plan = c("auto", "sequential", "multisession", "multicore")) {
  future_plan <- match.arg(future_plan)

  available <- tryCatch(future::availableCores(), error = function(e) NA_integer_)
  if (is.numeric(available) && length(available) == 1L && !is.na(available) && available > 0L) {
    actual_workers <- min(as.integer(workers), as.integer(available))
    if (actual_workers < workers) {
      warning(
        sprintf("Requested %d workers but only %d are available; using %d.", workers, available, actual_workers),
        call. = FALSE
      )
    }
  } else {
    actual_workers <- as.integer(workers)
  }

  strategy <- if (actual_workers <= 1L || identical(future_plan, "sequential")) {
    "sequential"
  } else if (identical(future_plan, "auto")) {
    if (isTRUE(future::supportsMulticore())) "multicore" else "multisession"
  } else if (identical(future_plan, "multicore") && !isTRUE(future::supportsMulticore())) {
    warning("multicore is not supported here; falling back to multisession.", call. = FALSE)
    "multisession"
  } else {
    future_plan
  }

  old_plan <- future::plan()
  if (identical(strategy, "sequential")) {
    future::plan(future::sequential)
  } else if (identical(strategy, "multicore")) {
    future::plan(future::multicore, workers = actual_workers)
  } else {
    future::plan(future::multisession, workers = actual_workers)
  }

  list(
    workers = actual_workers,
    strategy = strategy,
    old_plan = old_plan
  )
}

#' @keywords internal
#' @noRd
.cli_print_config <- function(cfg) {
  cat("Resolved configuration:\n")
  utils::str(as.list(cfg), max.level = 2, give.attr = FALSE)
  invisible(cfg)
}

#' @keywords internal
#' @noRd
.cli_print_dry_run <- function(analysis, cfg, mode, workers, future_plan) {
  .cli_print_config(cfg)
  cat(
    sprintf(
      "\nDry run summary:\n  mode: %s\n  output: %s\n  workers: %d\n  future plan: %s\n",
      mode,
      cfg$output,
      workers,
      future_plan
    )
  )
  print(analysis)

  if (identical(mode, "searchlight") && length(analysis$entries) > 0L) {
    cat("\nEngine eligibility:\n")
    print(searchlight_engines(
      model_spec = analysis$entries[[1]]$model_spec,
      method = cfg$type
    ))
  }

  invisible(analysis)
}

#' @keywords internal
#' @noRd
.cli_runner_args <- function(mode, options, cfg) {
  args <- list(verbose = isTRUE(options$verbose))

  if (identical(mode, "regional")) {
    pool_predictions <- .cli_or(options$pool_predictions, cfg$pool_predictions)
    if (!is.null(pool_predictions)) {
      args$pool_predictions <- pool_predictions
    }

    stack_folds <- .cli_or(options$stack_folds, cfg$stack_folds)
    if (!is.null(stack_folds)) {
      args$stack_folds <- as.integer(stack_folds)
    }

    stack_seed <- .cli_or(options$stack_seed, cfg$stack_seed)
    if (!is.null(stack_seed)) {
      args$stack_seed <- as.integer(stack_seed)
    }

    stack_lambda <- .cli_or(options$stack_lambda, cfg$stack_lambda)
    if (!is.null(stack_lambda)) {
      args$stack_lambda <- as.numeric(stack_lambda)
    }

    if (isTRUE(options$coalesce_design_vars)) {
      args$coalesce_design_vars <- TRUE
    }
  }

  args
}

#' @keywords internal
#' @noRd
run_cli_command <- function(mode = c("searchlight", "regional"),
                            argv = commandArgs(trailingOnly = TRUE)) {
  mode <- match.arg(mode)

  if (!requireNamespace("optparse", quietly = TRUE)) {
    stop(
      "This command requires the 'optparse' package. Install it with install.packages('optparse').",
      call. = FALSE
    )
  }

  parser <- .cli_parser(mode)
  argv <- .cli_normalize_argv(argv)
  if (length(argv) == 0L) {
    optparse::print_help(parser)
    return(invisible(NULL))
  }

  parsed <- optparse::parse_args(parser, args = argv, positional_arguments = TRUE)
  if (length(parsed$args) > 0L) {
    stop("Unexpected positional arguments: ", paste(parsed$args, collapse = " "), call. = FALSE)
  }
  options <- .cli_standardize_options(parsed$options)

  log_levels <- list(
    default = futile.logger::flog.threshold(),
    ROOT = futile.logger::flog.threshold(name = "ROOT"),
    rMVPA = futile.logger::flog.threshold(name = "rMVPA")
  )
  on.exit({
    futile.logger::flog.threshold(log_levels$default)
    futile.logger::flog.threshold(log_levels$ROOT, name = "ROOT")
    futile.logger::flog.threshold(log_levels$rMVPA, name = "rMVPA")
  }, add = TRUE)

  target_log_level <- if (isTRUE(options$verbose)) futile.logger::INFO else futile.logger::WARN
  futile.logger::flog.threshold(target_log_level)
  futile.logger::flog.threshold(target_log_level, name = "ROOT")
  futile.logger::flog.threshold(target_log_level, name = "rMVPA")

  if (isTRUE(options$version)) {
    cat(as.character(utils::packageVersion("rMVPA")), "\n")
    return(invisible(NULL))
  }
  if (isTRUE(options$list_models)) {
    .cli_list_models()
    return(invisible(NULL))
  }
  if (isTRUE(options$example_config)) {
    cat(.cli_example_config(mode), "\n")
    return(invisible(NULL))
  }

  if (!is.null(options$seed)) {
    set.seed(as.integer(options$seed)[1])
  }

  cfg <- do.call(mvpa_config, c(list(mode = mode), .cli_config_args(mode, options)))
  cfg <- .cli_finalize_config(cfg, mode)
  cfg <- .cli_validate_config(cfg, mode)

  workers <- .cli_resolve_workers(options, cfg)
  preflight <- match.arg(as.character(options$preflight)[1], c("warn", "error", "off"))
  save_level <- match.arg(as.character(options$save_level)[1], c("minimal", "standard", "complete"))
  future_plan <- match.arg(as.character(options$future_plan)[1], c("auto", "sequential", "multisession", "multicore"))

  if (!isTRUE(options$dry_run)) {
    output_action <- .cli_prepare_output(
      output_dir = cfg$output,
      overwrite = isTRUE(options$overwrite),
      skip_if_exists = isTRUE(options$skip_if_exists) || isTRUE(cfg$skip_if_folder_exists)
    )
    if (identical(output_action, "skip")) {
      message("Output directory already exists; skipping analysis: ", cfg$output)
      return(invisible(NULL))
    }
  }

  plan_state <- .cli_apply_future_plan(workers = workers, future_plan = future_plan)
  on.exit(future::plan(plan_state$old_plan), add = TRUE)

  analysis <- build_analysis(cfg)
  if (isTRUE(options$dry_run)) {
    return(invisible(.cli_print_dry_run(
      analysis = analysis,
      cfg = cfg,
      mode = mode,
      workers = plan_state$workers,
      future_plan = plan_state$strategy
    )))
  }

  if (isTRUE(options$print_config)) {
    .cli_print_config(cfg)
    cat(
      sprintf(
        "\nRun settings:\n  workers: %d\n  future plan: %s\n\n",
        plan_state$workers,
        plan_state$strategy
      )
    )
  }

  runner_args <- .cli_runner_args(mode, options, cfg)
  result <- do.call(
    run_analysis,
    c(
      list(
        x = analysis,
        preflight = preflight
      ),
      runner_args
    )
  )

  save_results(
    result,
    dir = cfg$output,
    level = save_level,
    overwrite = isTRUE(options$overwrite),
    quiet = FALSE
  )

  message(
    sprintf(
      "%s analysis complete. Results written to: %s",
      tools::toTitleCase(mode),
      cfg$output
    )
  )

  invisible(result)
}
