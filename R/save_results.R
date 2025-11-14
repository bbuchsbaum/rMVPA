#' Save searchlight results to disk (S3 generic)
#'
#' @param x     A result object (e.g., `searchlight_result`, list of NeuroVols, etc.)
#' @param dir   Directory to write into (created if needed)
#' @param level One of "minimal", "standard", "complete"
#'              - minimal: maps only
#'              - standard: maps + manifest (default)
#'              - complete: maps + manifest + tables + aux (if present)
#' @param stack One of c("none","auto","vec"). 
#'              - "none": Write individual NIfTI files for each metric (default)
#'              - "auto": Stack into 4D if all maps are compatible, otherwise individual files
#'              - "vec": Force stacking into one 4D NIfTI file
#' @param fname Base filename when writing a 4D file (default "searchlight.nii.gz")
#'              Only used when stack="vec" or stack="auto" with compatible volumes
#' @param include Character vector of extras to include; subset of
#'                c("manifest","tables","aux"). `level` sets sensible defaults.
#' @param dtype Optional data_type for neuroim2 write_* (e.g., "FLOAT","DOUBLE")
#' @param overwrite Logical; if FALSE and a target exists, stop.
#' @param quiet Logical; suppress messages.
#' @return (invisible) a list describing what was written.
#' @export
save_results <- function(x, dir,
                         level = c("standard","minimal","complete"),
                         stack = c("none","auto","vec"),
                         fname = "searchlight.nii.gz",
                         include = NULL,
                         dtype = NULL,
                         overwrite = FALSE,
                         quiet = FALSE) {
  UseMethod("save_results")
}

# ---------- helpers (internal) ----------

.slugify <- function(s) {
  s <- gsub("[^A-Za-z0-9._-]+", "_", s)
  s <- gsub("_+", "_", s)
  s <- sub("^_+", "", s)
  s <- sub("_+$", "", s)
  s
}

.ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

.unique_path <- function(path, overwrite) {
  if (overwrite || !file.exists(path)) return(path)
  root <- tools::file_path_sans_ext(path)
  ext  <- sub("^[^.]*", "", basename(path)) # includes dot(s)
  i <- 1L
  repeat {
    cand <- file.path(dirname(path), paste0(basename(root), "_", i, ext))
    if (!file.exists(cand)) return(cand)
    i <- i + 1L
  }
}

.write_manifest <- function(manifest, dir, quiet) {
  f_yaml <- file.path(dir, "manifest.yaml")
  f_json <- file.path(dir, "manifest.json")
  ok <- FALSE
  if (requireNamespace("yaml", quietly = TRUE)) {
    yaml::write_yaml(manifest, f_yaml); ok <- TRUE
    if (!quiet) message("Wrote manifest: ", f_yaml)
    return(f_yaml)
  }
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(manifest, f_json, pretty = TRUE, auto_unbox = TRUE); ok <- TRUE
    if (!quiet) message("Wrote manifest: ", f_json)
    return(f_json)
  }
  # fallback
  f_rds <- file.path(dir, "manifest.rds")
  saveRDS(manifest, f_rds)
  if (!quiet) message("Wrote manifest: ", f_rds)
  f_rds
}

# heuristics from level
.compute_includes <- function(level, include) {
  level <- match.arg(level, c("minimal", "standard", "complete"))
  base <- switch(level,
                 minimal  = character(0),
                 standard = "manifest",
                 complete = c("manifest","tables","aux"))
  if (is.null(include)) base else unique(c(base, include))
}

# ---------- default method ----------
# Includes handling for lists of NeuroVol/NeuroVec so older pipelines keep working.
#' @export
save_results.default <- function(x, dir,
                                 level = c("standard","minimal","complete"),
                                 stack = c("none","auto","vec"),
                                 fname = "searchlight.nii.gz",
                                 include = NULL,
                                 dtype = NULL,
                                 overwrite = FALSE,
                                 quiet = FALSE) {
  # If it's a list of named neuroim objects, treat like searchlight_result$results
  if (is.list(x) && length(x) > 0L &&
      all(vapply(x, function(o) inherits(o, c("NeuroVol","NeuroVec")) ||
                           inherits(o, c("NeuroSurface","NeuroSurfaceVector")), logical(1)))) {
    pseudo <- structure(list(results = x,
                             n_voxels = NA_integer_,
                             active_voxels = NA_integer_,
                             metrics = names(x)),
                        class = c("searchlight_result","list"))
    return(save_results.searchlight_result(pseudo, dir, level, stack, fname,
                                           include, dtype, overwrite, quiet))
  }
  # Fallback: save .rds
  .ensure_dir(dir)
  out <- file.path(dir, "result.rds")
  out <- .unique_path(out, overwrite)
  saveRDS(x, out)
  if (!quiet) message("Saved R object: ", out)
  invisible(list(files = list(rds = out)))
}

# ---------- searchlight_result method ----------
#' @export
save_results.searchlight_result <- function(x, dir,
                                            level = c("standard","minimal","complete"),
                                            stack = c("none","auto","vec"),
                                            fname = "searchlight.nii.gz",
                                            include = NULL,
                                            dtype = NULL,
                                            overwrite = FALSE,
                                            quiet = FALSE) {

  level   <- match.arg(level)
  stack   <- match.arg(stack)
  include <- .compute_includes(level, include)

  .ensure_dir(dir)
  paths <- list()
  paths$root <- dir

  # subdirs used as needed
  maps_dir     <- file.path(dir, "maps")
  surfaces_dir <- file.path(dir, "surfaces")
  tables_dir   <- file.path(dir, "tables")
  aux_dir      <- file.path(dir, "aux")

  vols   <- list()
  surfs  <- list()
  names_all <- names(x$results) %||% rep.int(paste0("metric", seq_along(x$results)), length(x$results))

  # Partition volumetric vs surface
  for (nm in names_all) {
    obj <- x$results[[nm]]
    if (inherits(obj, c("NeuroVol","NeuroVec"))) {
      vols[[nm]] <- obj
    } else if (inherits(obj, c("NeuroSurface","NeuroSurfaceVector"))) {
      surfs[[nm]] <- obj
    } else {
      # tolerate unknown by saving as RDS
      .ensure_dir(aux_dir)
      f <- .unique_path(file.path(aux_dir, paste0(.slugify(nm), ".rds")), overwrite)
      saveRDS(obj, f)
      paths$aux <- c(paths$aux %||% character(), f)
      if (!quiet) message("Saved unknown map as RDS: ", f)
    }
  }

  # ---- write volumetric maps ----
  if (length(vols)) {
    .ensure_dir(maps_dir)

    write_vol_one <- function(v, file) {
      file <- .unique_path(file, overwrite)
      if (inherits(v, "NeuroVol")) {
        if (!is.null(dtype)) neuroim2::write_vol(v, file, data_type = dtype) else neuroim2::write_vol(v, file)
      } else if (inherits(v, "NeuroVec")) {
        if (!is.null(dtype)) neuroim2::write_vec(v, file, data_type = dtype) else neuroim2::write_vec(v, file)
      } else stop("Not a NeuroVol/NeuroVec")
      file
    }

    can_stack <- length(vols) > 1 &&
      all(vapply(vols, function(v) inherits(v, "NeuroVol"), logical(1))) &&
      {
        ref <- neuroim2::space(vols[[1]])
        all(vapply(vols, function(v) identical(neuroim2::space(v), ref), logical(1)))
      }

    do_stack <- switch(stack,
                       vec  = TRUE,
                       none = FALSE,
                       auto = isTRUE(can_stack))

    if (do_stack) {
      # stack into 4D via concat and write_vec
      v4 <- Reduce(neuroim2::concat, vols)
      out4 <- file.path(maps_dir, fname)
      out4 <- sub("\\.nii(\\.gz)?$", ".nii.gz", out4)  # prefer gz
      out4 <- write_vol_one(v4, out4)
      paths$maps <- c(paths$maps %||% character(), out4)
      paths$map_series_order <- names(vols)
      if (!quiet) message("Wrote 4D NIfTI: ", out4, " (series = ", paste(names(vols), collapse = ", "), ")")
    } else {
      # write one file per metric (3D each)
      for (nm in names(vols)) {
        out3 <- file.path(maps_dir, paste0(.slugify(nm), ".nii.gz"))
        out3 <- write_vol_one(vols[[nm]], out3)
        paths$maps <- c(paths$maps %||% character(), out3)
        if (!quiet) message("Wrote map: ", out3)
      }
      names(paths$maps) <- names(vols)
    }
  }

  # ---- write surface maps ----
  if (length(surfs)) {
    if (!requireNamespace("neurosurf", quietly = TRUE))
      stop("Saving surfaces requires the 'neurosurf' package.")
    .ensure_dir(surfaces_dir)
    for (nm in names(surfs)) {
      base <- file.path(surfaces_dir, .slugify(nm))
      # neurosurf::write_surf_data decides concrete extension(s)
      neurosurf::write_surf_data(surfs[[nm]], base)
      # We cannot predict exact extension(s) here; record base
      paths$surfaces <- c(paths$surfaces %||% character(), base)
      if (!quiet) message("Wrote surface data (base): ", base)
    }
    names(paths$surfaces) <- names(surfs)
  }

  # ---- tables (summaries) ----
  if ("tables" %in% include && length(vols)) {
    .ensure_dir(tables_dir)
    summ <- lapply(names(vols), function(nm) {
      v <- vols[[nm]]
      vals <- try(neuroim2::values(v), silent = TRUE)
      vals <- if (inherits(vals, "try-error")) NA_real_ else as.numeric(vals)
      data.frame(metric = nm,
                 n = sum(is.finite(vals)),
                 mean = mean(vals, na.rm = TRUE),
                 sd = stats::sd(vals, na.rm = TRUE),
                 median = stats::median(vals, na.rm = TRUE),
                 min = suppressWarnings(min(vals, na.rm = TRUE)),
                 max = suppressWarnings(max(vals, na.rm = TRUE)),
                 stringsAsFactors = FALSE)
    })
    tbl <- do.call(rbind, summ)
    fcsv <- file.path(tables_dir, "metric_summary.csv")
    fcsv <- .unique_path(fcsv, overwrite)
    utils::write.csv(tbl, fcsv, row.names = FALSE)
    paths$tables <- c(paths$tables %||% character(), fcsv)
    if (!quiet) message("Wrote summary table: ", fcsv)
  }

  # ---- aux (heavy) ----
  if ("aux" %in% include) {
    heavy <- intersect(c("predictions","folds","cv","design","config","aux"), names(x))
    if (length(heavy)) {
      .ensure_dir(aux_dir)
      for (nm in heavy) {
        f <- file.path(aux_dir, paste0(.slugify(nm), ".rds"))
        f <- .unique_path(f, overwrite)
        saveRDS(x[[nm]], f)
        paths$aux <- c(paths$aux %||% character(), f)
        if (!quiet) message("Wrote aux: ", f)
      }
      names(paths$aux) <- heavy
    }
  }

  # ---- manifest ----
  if ("manifest" %in% include) {
    man <- list(
      created = as.character(Sys.time()),
      rMVPA_version   = tryCatch(as.character(utils::packageVersion("rMVPA")), error = function(e) NA_character_),
      neuroim2_version = tryCatch(as.character(utils::packageVersion("neuroim2")), error = function(e) NA_character_),
      neurosurf_version = if (requireNamespace("neurosurf", quietly = TRUE))
        as.character(utils::packageVersion("neurosurf")) else NA_character_,
      class = class(x),
      metrics = names_all,
      files = paths
    )
    mfile <- .write_manifest(man, dir, quiet)
    paths$manifest <- mfile
  }

  invisible(paths)
}

# ---------- regional_mvpa_result method ----------
#' @export
save_results.regional_mvpa_result <- function(x, dir,
                                               level = c("standard","minimal","complete"),
                                               stack = c("none","auto","vec"),
                                               fname = "regional.nii.gz",
                                               include = NULL,
                                               dtype = NULL,
                                               overwrite = FALSE,
                                               quiet = FALSE) {

  level   <- match.arg(level)
  include <- .compute_includes(level, include)

  .ensure_dir(dir)
  paths <- list()
  paths$root <- dir

  # Use searchlight_result method for vol_results (maps)
  if (!is.null(x$vol_results)) {
    vol_paths <- save_results.searchlight_result(
      x$vol_results, dir, level = "minimal", stack, fname,
      include = NULL, dtype, overwrite, quiet
    )
    paths <- c(paths, vol_paths)
  }

  # Save performance_table
  if (!is.null(x$performance_table)) {
    perf_file <- file.path(dir, "performance_table.txt")
    perf_file <- .unique_path(perf_file, overwrite)
    utils::write.table(x$performance_table, perf_file,
                      row.names = FALSE, quote = FALSE, sep = "\t")
    paths$performance_table <- perf_file
    if (!quiet) message("Wrote performance table: ", perf_file)
  }

  # Save prediction_table
  if (!is.null(x$prediction_table)) {
    pred_file <- file.path(dir, "prediction_table.txt")
    pred_file <- .unique_path(pred_file, overwrite)
    utils::write.table(x$prediction_table, pred_file,
                      row.names = FALSE, quote = FALSE, sep = "\t")
    paths$prediction_table <- pred_file
    if (!quiet) message("Wrote prediction table: ", pred_file)
  }

  # Save fits if present and not minimal level
  if (!is.null(x$fits) && length(x$fits) > 0 && level != "minimal") {
    fits_dir <- file.path(dir, "fits")
    .ensure_dir(fits_dir)

    fit_files <- character(length(x$fits))
    for (i in seq_along(x$fits)) {
      # Get ROI identifier from performance_table if available
      roi_id <- if (!is.null(x$performance_table$roinum)) {
        x$performance_table$roinum[i]
      } else {
        i
      }

      fit_file <- file.path(fits_dir, sprintf("roi_%03d.rds", roi_id))
      fit_file <- .unique_path(fit_file, overwrite)
      saveRDS(x$fits[[i]], fit_file)
      fit_files[i] <- fit_file
    }

    paths$fits_dir <- fits_dir
    paths$fit_files <- fit_files
    if (!quiet) message("Wrote ", length(fit_files), " ROI fits to: ", fits_dir)
  }

  # Manifest
  if ("manifest" %in% include) {
    man <- list(
      created = as.character(Sys.time()),
      rMVPA_version = tryCatch(as.character(utils::packageVersion("rMVPA")),
                               error = function(e) NA_character_),
      class = class(x),
      n_rois = if (!is.null(x$performance_table)) nrow(x$performance_table) else NA_integer_,
      has_fits = !is.null(x$fits) && length(x$fits) > 0,
      files = paths
    )
    mfile <- .write_manifest(man, dir, quiet)
    paths$manifest <- mfile
  }

  invisible(paths)
}

`%||%` <- function(a, b) if (is.null(a)) b else a