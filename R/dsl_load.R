# rmvpa DSL: reading project and analysis files, inheritance, merging with
# origins, and the resolved (fully defaulted) analysis.

`%|NA|%` <- function(x, y) if (length(x) == 0L || is.na(x[[1]])) y else x[[1]]

.dsl_require_yaml <- function() {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("The rmvpa DSL needs the 'yaml' package. Install it with install.packages(\"yaml\").",
         call. = FALSE)
  }
}

.dsl_read_yaml <- function(path) {
  .dsl_require_yaml()
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  x <- tryCatch(yaml::read_yaml(path),
                error = function(e) stop(sprintf("%s: not valid YAML (%s)", path, conditionMessage(e)),
                                         call. = FALSE))
  if (is.null(x)) list() else x
}

.dsl_abort <- function(errs, header) {
  stop(paste(c(header, paste0("  - ", errs)), collapse = "\n"), call. = FALSE)
}

#' Open an rmvpa DSL project
#'
#' Reads and validates the project file (`rmvpa.yml`) and lists the analysis
#' files under `analyses/`.
#'
#' @param path The project directory, or the path to its `rmvpa.yml`.
#' @return An object of class `rmvpa_project`.
#' @seealso [rmvpa_explain()], [rmvpa_check()], [rmvpa_run()]
#' @export
rmvpa_project <- function(path = ".") {
  file <- if (dir.exists(path)) file.path(path, "rmvpa.yml") else path
  if (!file.exists(file)) {
    stop("No rmvpa.yml found at ", normalizePath(path, mustWork = FALSE),
         ". Create one with rmvpa_new() or `rmvpa new <bids_root>`.", call. = FALSE)
  }
  file <- normalizePath(file)
  cfg <- .dsl_read_yaml(file)
  errs <- .dsl_validate(cfg, .dsl_project_node(), where = paste0(basename(file), ": "))
  if (length(errs)) .dsl_abort(errs, sprintf("%s is not a valid rmvpa project file:", file))
  root <- dirname(file)
  adir <- file.path(root, "analyses")
  afiles <- if (dir.exists(adir)) sort(list.files(adir, pattern = "\\.ya?ml$", full.names = TRUE)) else character()
  names(afiles) <- tools::file_path_sans_ext(basename(afiles))
  structure(list(root = root, file = file, config = cfg, version = cfg$rmvpa,
                 analyses = afiles),
            class = "rmvpa_project")
}

#' @export
print.rmvpa_project <- function(x, ...) {
  cat("rmvpa project (DSL ", x$version, ")\n", sep = "")
  cat("  root:     ", x$root, "\n", sep = "")
  cat("  sets:     ", paste(names(x$config$data$sets), collapse = ", "), "\n", sep = "")
  cat("  analyses: ", if (length(x$analyses)) paste(names(x$analyses), collapse = ", ") else "(none)", "\n", sep = "")
  invisible(x)
}

# ---- merging with origins --------------------------------------------------

# Merge `over` into `base`. Maps merge key by key; scalars and lists replace;
# an explicit NULL (`~` in YAML) removes the inherited key. `origins` maps
# dotted leaf paths to the layer that set them; `notes` records replacements
# and removals so `explain` can show them.
.dsl_merge <- function(base, over, layer, origins, notes, path = "") {
  for (k in names(over)) {
    p <- if (nzchar(path)) paste0(path, ".", k) else k
    v <- over[[k]]
    if (is.null(v)) {
      if (!is.null(base[[k]])) notes[[length(notes) + 1L]] <- sprintf("%s: removed by %s", p, layer)
      base[k] <- NULL
      origins <- origins[!startsWith(names(origins), p)]
      next
    }
    if (.dsl_is_map(v) && length(v) > 0L && .dsl_is_map(base[[k]])) {
      r <- .dsl_merge(base[[k]], v, layer, origins, notes, p)
      base[[k]] <- r$value; origins <- r$origins; notes <- r$notes
    } else {
      old <- base[[k]]
      if (!is.null(old) && !identical(old, v) && (length(v) > 1L || length(old) > 1L)) {
        notes[[length(notes) + 1L]] <- sprintf("%s: %s replaced the %s list %s", p, layer,
                                             origins[p] %|NA|% "inherited", .dsl_inline(old))
      }
      base[[k]] <- v
      origins <- origins[!startsWith(names(origins), paste0(p, "."))]
      origins[[p]] <- layer
      if (.dsl_is_map(v)) origins <- .dsl_tag_leaves(v, p, layer, origins)
    }
  }
  list(value = base, origins = origins, notes = notes)
}

.dsl_tag_leaves <- function(x, path, layer, origins) {
  for (k in names(x)) {
    p <- paste0(path, ".", k)
    if (.dsl_is_map(x[[k]]) && length(x[[k]])) origins <- .dsl_tag_leaves(x[[k]], p, layer, origins)
    else origins[[p]] <- layer
  }
  origins
}

.dsl_layer <- function(acc, over, layer) {
  if (is.null(over) || length(over) == 0L) return(acc)
  r <- .dsl_merge(acc$value, over, layer, acc$origins, acc$notes)
  r
}

# ---- resolving an analysis -----------------------------------------------

.dsl_analysis_file <- function(project, analysis) {
  if (file.exists(analysis) && !dir.exists(analysis)) return(normalizePath(analysis))
  if (analysis %in% names(project$analyses)) return(project$analyses[[analysis]])
  stop(sprintf("No analysis '%s' in %s/analyses%s", analysis, project$root,
               .dsl_did_you_mean(analysis, names(project$analyses))), call. = FALSE)
}

#' Resolve one analysis of a project
#'
#' Reads an analysis file, follows its `extends:` line to the project file,
#' and fills every value that was not written from the versioned defaults.
#' Each resolved value records where it came from.
#'
#' @param project An `rmvpa_project` or a project path.
#' @param analysis Analysis name (file name under `analyses/` without the
#'   extension) or a path to an analysis file.
#' @return An object of class `rmvpa_analysis_spec`.
#' @export
rmvpa_analysis <- function(project, analysis) {
  if (!inherits(project, "rmvpa_project")) project <- rmvpa_project(project)
  file <- .dsl_analysis_file(project, analysis)
  name <- tools::file_path_sans_ext(basename(file))
  where <- paste0("analyses/", basename(file), ": ")
  a <- .dsl_read_yaml(file)

  errs <- .dsl_validate(a, .dsl_analysis_node(), where = where)
  if (length(errs)) .dsl_abort(errs, sprintf("%s is not a valid analysis file:", file))

  parent <- normalizePath(file.path(dirname(file), a$extends), mustWork = FALSE)
  if (!identical(parent, project$file)) {
    .dsl_abort(sprintf("%sextends: '%s' resolves to %s, but in DSL 0.1 an analysis must extend its project file %s",
                       where, a$extends, parent, project$file),
               "Inheritance error:")
  }

  verbs <- .dsl_verbs()
  verb <- intersect(names(a), names(verbs))
  defaults <- .dsl_defaults(project$version)
  scope_kind <- names(a$scope)[1]

  acc <- list(value = list(), origins = character(), notes = list())
  dlayer <- paste0("default@", project$version)
  dflt <- list(folds = defaults$folds, run = defaults$run)
  dflt[[verb]] <- defaults$verbs[[verb]]
  if (identical(scope_kind, "searchlight")) dflt$scope <- defaults$scope
  acc <- .dsl_layer(acc, dflt, dlayer)
  acc <- .dsl_layer(acc, project$config[c("folds", "run")], "project")
  analysis_part <- a[setdiff(names(a), c("extends", "description"))]
  acc <- .dsl_layer(acc, analysis_part, "analysis")
  cfg <- acc$value
  origins <- acc$origins

  # Data: project only (analyses refer to sets by name).
  data <- utils::modifyList(defaults$data, project$config$data)
  vb <- cfg[[verb]]

  set_errs <- character()
  sets <- names(data$sets)
  set_keys <- intersect(c("set", "train_set", "test_set"), names(verbs[[verb]]$node$fields))
  if ("set" %in% set_keys && is.null(vb$set)) {
    if (length(sets) == 1L) {
      vb$set <- sets
      origins[[paste0(verb, ".set")]] <- "inferred (only set in project)"
    } else {
      set_errs <- c(set_errs, sprintf("%s%s.set is required because the project defines %d sets (%s)",
                                      where, verb, length(sets), paste(sets, collapse = ", ")))
    }
  }
  for (k in set_keys) {
    s <- vb[[k]]
    if (!is.null(s) && !s %in% sets) {
      set_errs <- c(set_errs, sprintf("%s%s.%s: no set '%s' in rmvpa.yml%s", where, verb, k, s,
                                      .dsl_did_you_mean(s, sets)))
    }
  }
  if (length(set_errs)) .dsl_abort(set_errs, "Data set error:")
  cfg[[verb]] <- vb

  if (identical(verb, "transfer") && !.dsl_transfer_method_ok(vb$method)) {
    .dsl_abort(sprintf("%stransfer.method '%s' is not naive, remap_rrr, subspace_alignment or a registered classifier%s",
                       where, vb$method,
                       .dsl_did_you_mean(vb$method, c("naive", "remap_rrr", "subspace_alignment", names(MVPAModels)))),
               "Invalid analysis:")
  }

  # Inference: what the R code supports for this verb, method and scope.
  support <- .dsl_inference_support(verb, vb, scope_kind)
  inf_user <- a$inference %||% list()
  inf <- .dsl_resolve_inference(support, inf_user, defaults$permutation, dlayer)
  cfg$inference <- inf$value
  origins <- origins[!startsWith(names(origins), "inference.")]
  origins <- c(origins, inf$origins)

  structure(list(
    name = name, file = file, project = project, verb = verb,
    family = verbs[[verb]]$family, scope_kind = scope_kind,
    config = cfg, data = data, origins = origins, notes = acc$notes,
    inference_support = support, description = a$description
  ), class = "rmvpa_analysis_spec")
}

.dsl_transfer_method_ok <- function(m) {
  m %in% c("naive", "remap_rrr", "subspace_alignment") || m %in% names(MVPAModels)
}

#' @export
print.rmvpa_analysis_spec <- function(x, ...) {
  cat(.dsl_explain_text(x), sep = "\n")
  invisible(x)
}

# ---- explain ---------------------------------------------------------------

.dsl_inline <- function(v) {
  if (is.null(v)) return("~")
  if (.dsl_is_map(v)) {
    if (!length(v)) return("{}")
    return(paste0("{ ", paste(sprintf("%s: %s", names(v), vapply(v, .dsl_inline, "")), collapse = ", "), " }"))
  }
  if (is.list(v)) v <- unlist(v)
  f <- vapply(v, .dsl_scalar, character(1))
  if (length(v) == 1L) f else paste0("[", paste(f, collapse = ", "), "]")
}

.dsl_scalar <- function(v) {
  if (is.logical(v)) return(if (isTRUE(v)) "true" else "false")
  if (is.numeric(v)) return(format(v, scientific = FALSE, trim = TRUE))
  s <- as.character(v)
  reserved <- c("true", "false", "yes", "no", "on", "off", "null", "~", "")
  if (grepl("^[A-Za-z_./][A-Za-z0-9_./-]*$", s) && !tolower(s) %in% reserved) return(s)
  paste0("\"", gsub("([\"\\\\])", "\\\\\\1", s), "\"")
}

.dsl_emit <- function(x, origins, path = "", indent = 0L) {
  out <- character()
  pad <- strrep("  ", indent)
  for (k in names(x)) {
    p <- if (nzchar(path)) paste0(path, ".", k) else k
    v <- x[[k]]
    if (.dsl_is_map(v) && length(v) > 0L && !identical(names(v), "in")) {
      out <- c(out, paste0(pad, k, ":"), .dsl_emit(v, origins, p, indent + 1L))
    } else {
      o <- origins[p] %|NA|% ""
      line <- paste0(pad, k, ": ", .dsl_inline(v))
      out <- c(out, if (nzchar(o)) sprintf("%-46s # %s", line, o) else line)
    }
  }
  out
}

.dsl_explain_text <- function(spec) {
  cfg <- spec$config
  order <- c(spec$verb, "folds", "scope", "inference", "run")
  cfg <- cfg[intersect(order, names(cfg))]
  sup <- spec$inference_support
  c(
    sprintf("# rmvpa explain: analyses/%s.yml (DSL %s)", basename(spec$file), spec$project$version),
    sprintf("# every value is annotated with its origin; this output is itself a valid analysis file"),
    "extends: ../rmvpa.yml",
    .dsl_emit(cfg, spec$origins),
    sprintf("# null: %s", sup$description),
    if (length(spec$notes)) paste0("# note: ", unlist(spec$notes))
  )
}

#' Explain a resolved analysis
#'
#' Prints the analysis with every default filled in and each value annotated
#' with its origin (`analysis`, `project`, `default@0.1` or `inferred`).
#'
#' @inheritParams rmvpa_analysis
#' @return The explanation lines, invisibly.
#' @export
rmvpa_explain <- function(project = ".", analysis) {
  spec <- rmvpa_analysis(project, analysis)
  txt <- .dsl_explain_text(spec)
  cat(txt, sep = "\n")
  invisible(txt)
}
