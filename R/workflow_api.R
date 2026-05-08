#' Build a Public Analysis Configuration
#'
#' Creates a lightweight configuration object for the supported high-level
#' workflow API. Config values can be supplied directly as named arguments or
#' loaded from a YAML or R config file via \code{config}.
#'
#' @param mode Analysis mode: \code{"searchlight"}, \code{"regional"}, or
#'   \code{"global"}.
#' @param config Optional path to a YAML or R config file.
#' @param ... Additional named configuration values.
#'
#' @return An object of class \code{c("rmvpa_config", "list")}.
#' @export
mvpa_config <- function(mode = c("searchlight", "regional", "global"),
                        config = NULL,
                        ...) {
  mode <- match.arg(mode)
  file_cfg <- .read_public_config_file(config)
  dot_cfg <- list(...)
  merged <- utils::modifyList(file_cfg, dot_cfg)
  merged$mode <- mode
  merged$config_file <- config
  structure(merged, class = c("rmvpa_config", "list"))
}

#' Build an Analysis Object
#'
#' Converts an \code{rmvpa_config} into a concrete analysis object with
#' initialized design, datasets, and model specification(s).
#'
#' @param config An \code{rmvpa_config} returned by \code{\link{mvpa_config}},
#'   or a named list containing at least \code{mode}.
#'
#' @return An object of class \code{c("rmvpa_analysis", "list")}.
#' @export
build_analysis <- function(config) {
  config <- .coerce_public_config(config)

  if (!is.null(config$model_spec)) {
    entry_name <- .rmvpa_or(config$name, "")
    return(structure(
      list(
        mode = config$mode,
        config = as.list(config),
        design = config$model_spec$design,
        entries = list(list(
          name = entry_name,
          model_spec = config$model_spec,
          region_mask = config$region_mask
        ))
      ),
      class = c("rmvpa_analysis", "list")
    ))
  }

  cfg_env <- .public_config_env(config)
  cfg_args <- as.list(config)
  analysis_type <- switch(
    config$mode,
    searchlight = "searchlight",
    regional = "mvpa_regional",
    global = "mvpa_global"
  )

  cfg_env <- initialize_standard_parameters(cfg_env, cfg_args, analysis_type)
  cfg_env$tune_grid <- initialize_tune_grid(cfg_args, cfg_env)
  cfg_env$design <- initialize_design(cfg_env)

  design <- cfg_env$design
  feature_selector <- initialize_feature_selection(cfg_env)
  crossval <- initialize_crossval(cfg_env, design)

  built <- .build_analysis_entries(cfg_env, mode = config$mode)
  entries <- lapply(seq_along(built$datasets), function(i) {
    dset <- built$datasets[[i]]
    list(
      name = names(built$datasets)[[i]],
      model_spec = load_mvpa_model(cfg_env, dset, design, crossval, feature_selector),
      region_mask = built$region_masks[[i]]
    )
  })

  structure(
    list(
      mode = config$mode,
      config = as.list(cfg_env),
      design = design,
      entries = entries
    ),
    class = c("rmvpa_analysis", "list")
  )
}

#' Run a Built Analysis
#'
#' Runs a high-level workflow object created by \code{\link{build_analysis}}.
#' The input can also be an \code{rmvpa_config}, in which case it is built
#' first.
#'
#' @param x An \code{rmvpa_analysis} or \code{rmvpa_config}.
#' @param preflight One of \code{"warn"} (default), \code{"error"}, or
#'   \code{"off"}.
#' @param ... Additional arguments forwarded to the underlying runner.
#'
#' @return An object of class \code{c("rmvpa_analysis_run", "list")}.
#' @export
run_analysis <- function(x,
                         preflight = c("warn", "error", "off"),
                         ...) {
  preflight <- match.arg(preflight)
  analysis <- if (inherits(x, "rmvpa_analysis")) x else build_analysis(x)

  results <- vector("list", length(analysis$entries))
  for (i in seq_along(analysis$entries)) {
    entry <- analysis$entries[[i]]
    result <- .run_analysis_entry(
      analysis = analysis,
      entry = entry,
      preflight = preflight,
      ...
    )
    results[[i]] <- .attach_analysis_context(
      result = result,
      analysis = analysis,
      entry = entry,
      preflight = preflight
    )
  }

  names(results) <- vapply(
    seq_along(analysis$entries),
    function(i) .analysis_output_label(analysis$entries[[i]]$name, i),
    character(1)
  )

  structure(
    list(
      analysis = analysis,
      results = results,
      preflight = preflight
    ),
    class = c("rmvpa_analysis_run", "list")
  )
}

#' Save a High-Level Analysis Run
#'
#' Writes all results from an \code{rmvpa_analysis_run} and a root manifest
#' describing the workflow configuration and per-entry outputs.
#'
#' @inheritParams save_results
#' @param x An \code{rmvpa_analysis_run}.
#' @return Invisibly, a nested list of file paths.
#' @export
save_results.rmvpa_analysis_run <- function(x, dir = NULL,
                                            level = c("standard", "minimal", "complete"),
                                            stack = c("none", "auto", "vec"),
                                            fname = "analysis.nii.gz",
                                            include = NULL,
                                            dtype = NULL,
                                            overwrite = FALSE,
                                            quiet = FALSE) {
  level <- match.arg(level)
  stack <- match.arg(stack)
  dir <- .rmvpa_or(dir, .rmvpa_or(x$analysis$config$output, "rmvpa_output"))

  .ensure_dir(dir)
  paths <- list(root = dir, results = list())

  aux_dir <- file.path(dir, "aux")
  .ensure_dir(aux_dir)
  config_rds <- .unique_path(file.path(aux_dir, "analysis_config.rds"), overwrite)
  saveRDS(x$analysis$config, config_rds)
  paths$aux <- c(paths$aux %||% character(), config_rds)

  for (i in seq_along(x$results)) {
    label <- names(x$results)[[i]]
    subdir <- file.path(dir, .slugify(label))
    paths$results[[label]] <- save_results(
      x$results[[i]],
      dir = subdir,
      level = level,
      stack = stack,
      fname = fname,
      include = include,
      dtype = dtype,
      overwrite = overwrite,
      quiet = quiet
    )
  }

  paths$manifest <- .try_write_manifest(
    build = function() list(
      created = as.character(Sys.time()),
      mode = x$analysis$mode,
      preflight_policy = x$preflight,
      config_summary = .summarize_public_config(x$analysis$config),
      entries = lapply(seq_along(x$analysis$entries), function(i) {
        entry <- x$analysis$entries[[i]]
        list(
          name = .analysis_output_label(entry$name, i),
          model_class = class(entry$model_spec)[1],
          dataset_class = class(entry$model_spec$dataset)[1],
          design_class = class(entry$model_spec$design)[1]
        )
      }),
      git_sha = .git_head_sha(),
      session_info = .safe_session_info(),
      files = paths
    ),
    dir = dir, quiet = quiet, paths = paths
  )
  invisible(paths)
}

#' @export
#' @method print rmvpa_analysis
print.rmvpa_analysis <- function(x, ...) {
  cat(
    sprintf(
      "<rmvpa_analysis> mode=%s entries=%d\n",
      x$mode,
      length(x$entries)
    )
  )
  invisible(x)
}

#' @export
#' @method print rmvpa_analysis_run
print.rmvpa_analysis_run <- function(x, ...) {
  cat(
    sprintf(
      "<rmvpa_analysis_run> mode=%s results=%d preflight=%s\n",
      x$analysis$mode,
      length(x$results),
      x$preflight
    )
  )
  invisible(x)
}

#' @keywords internal
#' @noRd
.rmvpa_or <- function(x, y) {
  if (is.null(x)) y else x
}

#' @keywords internal
#' @noRd
.coerce_public_config <- function(config) {
  if (inherits(config, "rmvpa_config")) {
    return(config)
  }

  if (!is.list(config)) {
    stop("build_analysis: `config` must be an rmvpa_config or a named list.", call. = FALSE)
  }

  mode <- config$mode
  if (is.null(mode) || !is.character(mode) || length(mode) != 1L) {
    stop("build_analysis: list config must include a scalar `mode`.", call. = FALSE)
  }

  dots <- config
  dots$mode <- NULL
  do.call(mvpa_config, c(list(mode = mode), dots))
}

#' @keywords internal
#' @noRd
.read_public_config_file <- function(path) {
  if (is.null(path)) {
    return(list())
  }

  if (!file.exists(path)) {
    stop("mvpa_config: config file does not exist: ", path, call. = FALSE)
  }

  if (grepl("\\.(yaml|yml)$", path, ignore.case = TRUE)) {
    if (requireNamespace("yaml", quietly = TRUE)) {
      cfg <- yaml::read_yaml(path)
    } else if (requireNamespace("io", quietly = TRUE)) {
      cfg <- io::qread(path)
    } else {
      stop("mvpa_config: YAML config support requires package 'yaml' or 'io'.", call. = FALSE)
    }
    return(if (is.null(cfg)) list() else as.list(cfg))
  }

  if (grepl("\\.[rR]$", path)) {
    env <- new.env(parent = baseenv())
    sys.source(path, envir = env)
    return(as.list(env, all.names = TRUE))
  }

  stop("mvpa_config: unsupported config file type: ", path, call. = FALSE)
}

#' @keywords internal
#' @noRd
.public_config_env <- function(config) {
  env <- new.env(parent = baseenv())
  cfg <- as.list(config)
  for (nm in names(cfg)) {
    env[[nm]] <- cfg[[nm]]
  }
  env
}

#' @keywords internal
#' @noRd
.build_analysis_entries <- function(config_env, mode) {
  datasets <- NULL
  region_masks <- NULL

  if (identical(config_env$data_mode, "image")) {
    mask_obj <- load_mask(config_env)
    dataset_mask <- if (identical(mode, "regional")) mask_obj else methods::as(mask_obj, "LogicalNeuroVol")
    datasets <- list(initialize_image_data(config_env, dataset_mask))
    names(datasets) <- ""
    region_masks <- list(if (identical(mode, "regional")) mask_obj else NULL)
  } else {
    datasets <- initialize_surface_data(config_env)
    if (is.null(names(datasets))) {
      names(datasets) <- rep.int("", length(datasets))
    }
    region_masks <- lapply(datasets, function(dset) {
      if (identical(mode, "regional")) dset$mask else NULL
    })
  }

  list(datasets = datasets, region_masks = region_masks)
}

#' @keywords internal
#' @noRd
.run_analysis_entry <- function(analysis, entry, preflight, ...) {
  cfg <- analysis$config
  dots <- list(...)

  if (identical(analysis$mode, "searchlight")) {
    call_args <- c(
      list(
        model_spec = entry$model_spec,
        radius = .rmvpa_or(cfg$radius, 8),
        method = .rmvpa_or(cfg$method, .rmvpa_or(cfg$type, "standard")),
        niter = .rmvpa_or(cfg$niter, 4),
        preflight = preflight
      ),
      dots
    )
    if (!is.null(cfg$backend) && is.null(call_args$backend)) {
      call_args$backend <- cfg$backend
    }
    if (!is.null(cfg$engine) && is.null(call_args$engine)) {
      call_args$engine <- cfg$engine
    }
    if (!is.null(cfg$batch_size) && is.null(call_args$batch_size)) {
      call_args$batch_size <- cfg$batch_size
    }
    return(do.call(run_searchlight, call_args))
  }

  if (identical(analysis$mode, "regional")) {
    region_mask <- .rmvpa_or(entry$region_mask, entry$model_spec$dataset$mask)
    call_args <- c(
      list(
        model_spec = entry$model_spec,
        region_mask = region_mask,
        return_fits = isTRUE(.rmvpa_or(cfg$return_fits, cfg$save_predictors)),
        preflight = preflight
      ),
      dots
    )
    if (!is.null(cfg$backend) && is.null(call_args$backend)) {
      call_args$backend <- cfg$backend
    }
    return(do.call(run_regional, call_args))
  }

  call_args <- c(
    list(
      model_spec = entry$model_spec,
      return_fits = isTRUE(cfg$return_fits),
      preflight = preflight
    ),
    dots
  )
  do.call(run_global, call_args)
}

#' @keywords internal
#' @noRd
.analysis_output_label <- function(name, index) {
  if (is.character(name) && length(name) == 1L && nzchar(name)) {
    return(name)
  }
  if (index == 1L) {
    return("analysis")
  }
  sprintf("analysis_%02d", index)
}

#' @keywords internal
#' @noRd
.attach_analysis_context <- function(result, analysis, entry, preflight) {
  attr(result, "analysis_context") <- list(
    mode = analysis$mode,
    entry_name = entry$name,
    preflight_policy = preflight,
    model_class = class(entry$model_spec)[1],
    dataset_class = class(entry$model_spec$dataset)[1],
    design_class = class(entry$model_spec$design)[1],
    config = .summarize_public_config(analysis$config)
  )
  result
}

#' @keywords internal
#' @noRd
.summarize_public_config <- function(config) {
  out <- list()
  cfg <- as.list(config)
  skip <- c(
    "dataset", "design", "model_spec", "full_train_design", "full_test_design",
    "train_design", "test_design", "cross_validation"
  )
  keep <- setdiff(names(cfg), skip)

  for (nm in keep) {
    val <- cfg[[nm]]
    if (is.null(val)) {
      next
    }
    if (is.atomic(val) && length(val) <= 20L) {
      out[[nm]] <- val
      next
    }
    if (is.data.frame(val)) {
      out[[nm]] <- list(
        class = class(val)[1],
        nrow = nrow(val),
        ncol = ncol(val),
        names = colnames(val)
      )
      next
    }
    if (inherits(val, "formula")) {
      out[[nm]] <- Reduce(paste, deparse(val))
      next
    }
    out[[nm]] <- list(class = class(val)[1])
  }

  out
}
