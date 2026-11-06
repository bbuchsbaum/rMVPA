# rmvpa DSL: check (no fitting), and the generators new, add and migrate.

# ---- check -----------------------------------------------------------------

.dsl_finding <- function(analysis, subject, level, message) {
  data.frame(analysis = analysis, subject = subject, level = level, message = message,
             stringsAsFactors = FALSE)
}

.dsl_column_keys <- function(spec) {
  vb <- spec$config[[spec$verb]]
  cols <- c(vb$target, vb$item, vb$link, vb$terms)
  if (!identical(spec$config$folds$scheme, "kfold")) cols <- c(cols, spec$config$folds$by)
  unique(cols)
}

#' Check analyses without fitting anything
#'
#' Validates every analysis file, then for each subject finds the input
#' files, reads the design tables (not the images), checks that referenced
#' columns exist, runs the fold assertions, and confirms that atlases match
#' the data grid. Also estimates the cost of the run.
#'
#' @inheritParams rmvpa_run
#' @param quiet Suppress the printed report.
#' @return An object of class `rmvpa_check`: a data frame of findings with
#'   columns analysis, subject, level (`ok`, `warning`, `error`) and message,
#'   plus an `ok` attribute.
#' @export
rmvpa_check <- function(project = ".", analyses = NULL, subjects = NULL, quiet = FALSE) {
  if (!inherits(project, "rmvpa_project")) project <- rmvpa_project(project)
  analyses <- analyses %||% names(project$analyses)
  out <- list()
  add <- function(...) out[[length(out) + 1L]] <<- .dsl_finding(...)
  passes <- 0
  for (a in analyses) {
    spec <- tryCatch(rmvpa_analysis(project, a), error = function(e) { add(a, "", "error", conditionMessage(e)); NULL })
    if (is.null(spec)) next
    subs <- tryCatch(subjects %||% .dsl_subjects(spec$data, project$root),
                     error = function(e) { add(a, "", "error", conditionMessage(e)); character() })
    inf <- spec$config$inference
    perms <- if (identical(inf$test, "permutation") &&
                 spec$inference_support$kind %in% c("searchlight_permutation", "regional_permutation")) inf$permutations else 0
    passes <- passes + length(subs) * (1 + perms)
    if (spec$inference_support$kind %in% c("model_nperm", "regional_permutation") && identical(inf$test, "permutation") &&
        inf$permutations < 1000) {
      add(a, "", "warning", sprintf("inference.permutations = %d: p-values cannot go below %.4f; use at least 1000",
                                    inf$permutations, 1 / (inf$permutations + 1)))
    }
    if (identical(spec$inference_support$kind, "searchlight_permutation") && identical(inf$test, "permutation") &&
        inf$permutations < 1000) {
      add(a, "", "warning", sprintf(paste0(
        "inference.permutations = %d: centres within one permutation are correlated, so the pooled null resolves ",
        "p only to about %.4f; use at least 1000 before reporting FDR-corrected maps"),
        inf$permutations, 1 / (inf$permutations + 1)))
    }
    if (identical(inf$test, "none")) add(a, "", "ok", paste("inference:", spec$inference_support$description))
    for (pk in .dsl_required_packages(spec$verb, spec$config[[spec$verb]])) {
      if (!requireNamespace(pk, quietly = TRUE)) {
        add(a, "", "error", sprintf("%s needs the '%s' package: install.packages(\"%s\")", spec$verb, pk, pk))
      }
    }
    for (f in .dsl_referenced_files(spec)) {
      if (!file.exists(f)) add(a, "", "error", sprintf("referenced file not found: %s", f))
    }
    sets <- .dsl_cell_sets(spec)
    for (s in subs) {
      tryCatch({
        mask_file <- .dsl_mask_file(spec$data, s, project$root)
        dd <- lapply(unique(unname(sets)), function(st) .dsl_load_design(spec$data, st, s, project$root))
        names(dd) <- unique(unname(sets))
        tr <- dd[[sets[["train"]]]]$design[dd[[sets[["train"]]]]$keep, , drop = FALSE]
        te <- if (!is.na(sets["test"])) dd[[sets[["test"]]]]$design[dd[[sets[["test"]]]]$keep, , drop = FALSE]
        for (col in .dsl_column_keys(spec)) {
          if (!col %in% names(tr)) add(a, s, "error", sprintf("train set has no column '%s'%s", col,
                                                              .dsl_did_you_mean(col, names(tr))))
          if (!is.null(te) && col %in% c(spec$config[[spec$verb]]$target, spec$config[[spec$verb]]$link) &&
              !col %in% names(te)) {
            add(a, s, "error", sprintf("test set has no column '%s'", col))
          }
        }
        if (all(.dsl_column_keys(spec) %in% names(tr))) {
          for (r in .dsl_check_assertions(spec, tr, te)) {
            add(a, s, if (isTRUE(r$pass)) "ok" else "error", sprintf("assert %s: %s", r$name, r$detail))
          }
        }
        if (identical(spec$scope_kind, "regions")) {
          atlas <- neuroim2::read_vol(.dsl_resolve_path(spec$config$scope$regions$atlas, project$root))
          mask <- neuroim2::read_vol(mask_file)
          if (!identical(dim(atlas)[1:3], dim(mask)[1:3])) {
            add(a, s, "error", sprintf("atlas grid %s does not match the mask grid %s",
                                       paste(dim(atlas)[1:3], collapse = "x"), paste(dim(mask)[1:3], collapse = "x")))
          }
        }
        add(a, s, "ok", sprintf("%d of %d train rows kept%s", nrow(tr), length(dd[[sets[["train"]]]]$keep),
                                if (!is.null(te)) sprintf("; %d test rows", nrow(te)) else ""))
      }, error = function(e) add(a, s, "error", conditionMessage(e)))
    }
  }
  res <- if (length(out)) do.call(rbind, out) else .dsl_finding(character(), character(), character(), character())
  attr(res, "ok") <- !any(res$level == "error")
  attr(res, "passes") <- passes
  class(res) <- c("rmvpa_check", "data.frame")
  if (!quiet) print(res)
  invisible(res)
}

#' @export
print.rmvpa_check <- function(x, ...) {
  df <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(df))) {
    tag <- switch(df$level[i], ok = "  ok ", warning = "WARN ", error = "ERROR")
    who <- if (nzchar(df$subject[i])) sprintf("%s sub-%s", df$analysis[i], df$subject[i]) else df$analysis[i]
    cat(sprintf("%s %s: %s\n", tag, who, df$message[i]))
  }
  n_err <- sum(df$level == "error")
  cat(sprintf("%s: %d error(s), %d warning(s). Estimated searchlight/ROI passes: %s\n",
              if (n_err) "check failed" else "check passed", n_err, sum(df$level == "warning"),
              format(attr(x, "passes"), big.mark = ",")))
  invisible(x)
}

# ---- new -------------------------------------------------------------------

.dsl_yaml_lines <- function(x, indent = 0L) {
  pad <- strrep("  ", indent)
  out <- character()
  for (k in names(x)) {
    v <- x[[k]]
    if (.dsl_is_map(v) && length(v) && !identical(names(v), "in")) {
      out <- c(out, paste0(pad, k, ":"), .dsl_yaml_lines(v, indent + 1L))
    } else {
      out <- c(out, paste0(pad, k, ": ", .dsl_inline(v)))
    }
  }
  out
}

.dsl_target_candidates <- function(d) {
  skip <- c("onset", "duration", "run", "trial", "trial_index", "volume", "sample_index")
  cols <- setdiff(names(d), skip)
  cols[vapply(cols, function(c) {
    v <- d[[c]]; k <- length(unique(v))
    (is.character(v) || is.factor(v) || is.logical(v)) && k >= 2L && k <= max(2L, nrow(d) %/% 2L)
  }, TRUE)]
}

#' Create an rmvpa DSL project from a BIDS-derivatives folder
#'
#' Scans `bids_root` for beta series (`*_betas.nii[.gz]` with design tables
#' alongside), writes `rmvpa.yml` with one set per task, and a first analysis
#' `analyses/decode_<target>.yml` (a 6 mm searchlight).
#'
#' @param bids_root Path to the BIDS-derivatives folder holding the betas.
#' @param dir Project directory to create (default: current directory).
#' @param target Design column to decode. When `NULL`, the candidates are
#'   listed; in an interactive session you are asked to pick one.
#' @param suffix BIDS suffix of the beta images.
#' @param overwrite Replace an existing `rmvpa.yml`.
#' @return The project directory, invisibly.
#' @export
rmvpa_new <- function(bids_root, dir = ".", target = NULL, suffix = "betas", overwrite = FALSE) {
  if (!dir.exists(bids_root)) stop("Not a directory: ", bids_root, call. = FALSE)
  .dsl_ensure_dir(dir)
  proj_file <- file.path(dir, "rmvpa.yml")
  if (file.exists(proj_file) && !overwrite) stop(proj_file, " already exists (use overwrite = TRUE).", call. = FALSE)
  imgs <- list.files(bids_root, pattern = paste0("_", suffix, "\\.nii(\\.gz)?$"), recursive = TRUE, full.names = TRUE)
  imgs <- imgs[grepl("(^|/)sub-[^/]+/", imgs)]
  if (!length(imgs)) stop(sprintf("No sub-*/**/*_%s.nii[.gz] files under %s", suffix, bids_root), call. = FALSE)
  ents <- lapply(imgs, function(f) .dsl_parse_bids_name(f)$entities)
  tasks <- sort(unique(vapply(ents, function(e) e$task %||% NA_character_, "")))
  tasks <- tasks[!is.na(tasks)]
  if (!length(tasks)) stop("The beta images carry no task- entity; write rmvpa.yml by hand.", call. = FALSE)
  spaces <- unique(vapply(ents, function(e) e$space %||% NA_character_, ""))
  space <- if (length(spaces) == 1L && !is.na(spaces)) spaces else NULL
  masks <- list.files(bids_root, pattern = "_desc-[A-Za-z0-9]+_mask\\.nii(\\.gz)?$", recursive = TRUE)
  mask_desc <- if (length(masks)) .dsl_parse_bids_name(masks[1])$entities$desc else "brain"

  design <- .dsl_read_table(sub("\\.nii(\\.gz)?$", ".tsv", imgs[1]))
  cands <- .dsl_target_candidates(design)
  if (is.null(target)) {
    if (!length(cands)) stop("No categorical design columns found; pass target = '<column>'.", call. = FALSE)
    if (interactive()) {
      i <- utils::menu(cands, title = "Which design column should the first analysis decode?")
      if (i == 0L) stop("No target chosen.", call. = FALSE)
      target <- cands[i]
    } else {
      stop(sprintf("Choose a target column: %s (rmvpa new %s --target <column>)",
                   paste(cands, collapse = ", "), bids_root), call. = FALSE)
    }
  }
  if (!target %in% names(design)) {
    stop(sprintf("Design tables have no column '%s'%s", target, .dsl_did_you_mean(target, names(design))), call. = FALSE)
  }

  rel <- function(p) {
    a <- normalizePath(p); b <- normalizePath(dir)
    if (startsWith(a, paste0(b, "/"))) substring(a, nchar(b) + 2L) else a
  }
  sets <- stats::setNames(lapply(tasks, function(t) list(task = t)), tasks)
  proj <- list(
    rmvpa = "0.1",
    data = c(list(bids = rel(bids_root)), if (!is.null(space)) list(space = space),
             list(mask = mask_desc, sets = sets)),
    folds = c(list(by = if ("run" %in% names(design)) "run" else "<run column>"),
              list(assert = "classes_in_every_fold")),
    run = list(workers = 1L, seed = as.integer(format(Sys.Date(), "%Y%m%d")))
  )
  writeLines(c("# rmvpa project file. Analyses live in analyses/*.yml and extend this file.",
               .dsl_yaml_lines(proj)), proj_file)
  .dsl_ensure_dir(file.path(dir, "analyses"))
  name <- paste0("decode_", gsub("[^A-Za-z0-9]+", "_", target))
  an <- list(extends = "../rmvpa.yml",
             decode = c(if (length(tasks) > 1L) list(set = tasks[1]), list(target = target)),
             scope = list(searchlight = list(radius = 6L)))
  writeLines(.dsl_yaml_lines(an), file.path(dir, "analyses", paste0(name, ".yml")))
  message(sprintf("Wrote %s and analyses/%s.yml (%d set%s: %s).", proj_file, name, length(tasks),
                  if (length(tasks) > 1L) "s" else "", paste(tasks, collapse = ", ")))
  message("next: rmvpa draft")
  invisible(dir)
}

# ---- add -------------------------------------------------------------------

.dsl_placeholders <- list(
  target = "<column to predict>", item = "<column naming the stimulus>", link = "<column pairing train and test rows>",
  train_set = "<set name>", test_set = "<set name>", r_file = "specs/build_spec.R",
  distances = "rdms/<name>.csv", contrasts = "contrasts/<name>.tsv", seed_region = "rdms/<seed>.csv",
  predictor = "rdms/<x>.csv", outcome = "rdms/<y>.csv"
)

#' Add an analysis file to a project
#'
#' Writes `analyses/<name>.yml` for `verb`, with required keys filled with
#' placeholders and the optional keys listed as comments with their defaults.
#'
#' @param project An `rmvpa_project` or a project path.
#' @param verb One of the DSL verbs, e.g. `"decode"`, `"rsa"`, `"transfer"`.
#' @param name File name for the analysis (without `.yml`).
#' @param scope `"searchlight"`, `"regions"` or `"whole_brain"`.
#' @return The path of the new file, invisibly.
#' @export
rmvpa_add <- function(project = ".", verb, name, scope = NULL) {
  if (!inherits(project, "rmvpa_project")) project <- rmvpa_project(project)
  verbs <- .dsl_verbs()
  if (!verb %in% names(verbs)) {
    stop(sprintf("Unknown verb '%s'%s. Verbs: %s", verb, .dsl_did_you_mean(verb, names(verbs)),
                 paste(names(verbs), collapse = ", ")), call. = FALSE)
  }
  if (!grepl("^[A-Za-z][A-Za-z0-9_-]*$", name)) stop("name must start with a letter and use letters, digits, _ or -.", call. = FALSE)
  f <- file.path(project$root, "analyses", paste0(name, ".yml"))
  if (file.exists(f)) stop(f, " already exists.", call. = FALSE)
  scope <- scope %||% .dsl_scope_support(verb)[1]
  node <- verbs[[verb]]$node
  defaults <- .dsl_defaults(project$version)$verbs[[verb]]
  req <- node$required
  sets <- names(project$config$data$sets)
  lines <- c("extends: ../rmvpa.yml", paste0(verb, ":"))
  if ("set" %in% names(node$fields) && length(sets) > 1L) lines <- c(lines, paste0("  set: ", sets[1]))
  for (k in req) {
    ph <- switch(node$fields[[k]]$type,
      dict = paste0("\n    <name>: ", if (k == "feature_spaces") "features/<name>.tsv" else "rdms/<name>.csv"),
      strs = " [<column>]",
      paste0(" ", .dsl_scalar(.dsl_placeholders[[k]] %||% "<value>")))
    lines <- c(lines, paste0("  ", k, ":", ph))
  }
  for (k in setdiff(names(node$fields), c(req, "set"))) {
    d <- defaults[[k]]
    lines <- c(lines, sprintf("  # %s: %s%s", k, if (!is.null(d)) .dsl_inline(d) else "",
                              if (!is.null(node$fields[[k]]$doc)) paste0("   (", node$fields[[k]]$doc, ")") else ""))
  }
  sc <- switch(scope,
    searchlight = c("scope:", "  searchlight:", "    radius: 6"),
    regions = c("scope:", "  regions:", "    atlas: <atlas_dseg.nii.gz>"),
    whole_brain = c("scope:", "  whole_brain: {}"))
  lines <- c(lines, sc)
  .dsl_ensure_dir(dirname(f))
  writeLines(lines, f)
  message(sprintf("Wrote %s. Fill in the <placeholders>, then run `rmvpa check %s`.", f, name))
  invisible(f)
}

# ---- migrate ---------------------------------------------------------------

#' Convert a config for the old CLI into a DSL project
#'
#' Reads a flat configuration written for `rmvpa-searchlight` or
#' `rmvpa-regional` and writes an equivalent project (`rmvpa.yml` with
#' explicit file paths, one subject) and one analysis file.
#'
#' @param config Path to the old YAML config.
#' @param dir Project directory to write.
#' @param mode `"searchlight"` or `"regional"`; guessed from the config when `NULL`.
#' @param name Analysis name.
#' @return The project directory, invisibly.
#' @export
rmvpa_migrate <- function(config, dir = ".", mode = NULL, name = "migrated") {
  old <- .dsl_read_yaml(config)
  need <- c("train_design", "train_data", "mask", "label_column")
  miss <- setdiff(need, names(old))
  if (length(miss)) stop("The old config lacks: ", paste(miss, collapse = ", "), call. = FALSE)
  base <- dirname(normalizePath(config))
  abs <- function(p) if (grepl("^(/|~)", p)) p else normalizePath(file.path(base, p), mustWork = FALSE)
  mode <- mode %||% (if (!is.null(old$radius) || !is.null(old$type)) "searchlight" else "regional")
  .dsl_ensure_dir(dir)
  if (file.exists(file.path(dir, "rmvpa.yml"))) stop(file.path(dir, "rmvpa.yml"), " already exists.", call. = FALSE)
  if (!is.null(old$test_data)) {
    stop("rmvpa_migrate() converts single-set configs. For a train/test config, name both sets in ",
         "rmvpa.yml (data.files with a {task} template, or BIDS discovery) and use a transfer analysis.",
         call. = FALSE)
  }
  proj <- list(rmvpa = "0.1",
               data = list(files = list(betas = abs(old$train_data), design = abs(old$train_design)),
                           subjects = "01", mask = abs(old$mask), sets = list(train = list(task = "train"))),
               folds = list(by = old$block_column %||% "run"),
               run = list(workers = as.integer(old$workers %||% 1L), seed = as.integer(old$seed %||% 1L)))
  vb <- list(target = old$label_column, method = old$model %||% "sda_notune")
  scope <- if (mode == "searchlight") {
    list(searchlight = c(list(radius = old$radius %||% 8),
                         if (!is.null(old$type)) list(sampling = old$type),
                         if (!is.null(old$niter)) list(iterations = as.integer(old$niter))))
  } else list(regions = list(atlas = abs(old$mask)))
  notes <- "# Migrated from an rmvpa-searchlight/rmvpa-regional config by rmvpa_migrate()."
  if (mode == "regional") notes <- c(notes, "# The old mask was a labelled ROI image; it serves as both the brain mask and the atlas.")
  if (!is.null(old$train_subset)) notes <- c(notes, sprintf("# TODO: old train_subset '%s' needs a data.sets.train.keep filter", old$train_subset))
  if (!is.null(old$tune_grid)) notes <- c(notes, "# NOTE: tune_grid has no DSL key in 0.1; use custom: for tuned models.")
  writeLines(c(notes, .dsl_yaml_lines(proj)), file.path(dir, "rmvpa.yml"))
  .dsl_ensure_dir(file.path(dir, "analyses"))
  an <- list(extends = "../rmvpa.yml", decode = vb, scope = scope)
  writeLines(.dsl_yaml_lines(an), file.path(dir, "analyses", paste0(name, ".yml")))
  message(sprintf("Wrote %s/rmvpa.yml and analyses/%s.yml. Review the TODO lines, then run `rmvpa check`.", dir, name))
  invisible(dir)
}
