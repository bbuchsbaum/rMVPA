# rmvpa DSL: building and running subject x analysis cells, writing outputs
# with provenance sidecars, incremental re-runs, status and the optional lock.

# ---- folds -----------------------------------------------------------------

.dsl_crossval <- function(folds, design, where) {
  if (identical(folds$scheme, "kfold")) {
    return(kfold_cross_validation(len = nrow(design), nfolds = folds$nfolds %||% 5L))
  }
  .dsl_require_columns(design, folds$by, where)
  block <- design[[folds$by]]
  if (length(unique(block)) < 2L) {
    stop(sprintf("%s folds.by '%s' has a single value; at least two are needed", where, folds$by), call. = FALSE)
  }
  switch(folds$scheme,
    leave_one_out = blocked_cross_validation(block),
    twofold = twofold_blocked_cross_validation(block, nreps = folds$repeats %||% 10L),
    bootstrap = bootstrap_blocked_cross_validation(block, nreps = folds$repeats %||% 10L)
  )
}

# Each assertion returns list(name, pass, detail). Assertions that do not apply
# to the verb are an error in check (the closed enum says where they apply).
.dsl_assertion_applies <- list(
  classes_in_every_fold = c("decode", "decode_trialwise", "transfer", "pattern"),
  no_item_across_folds = c("rsa", "rsa_vector", "rsa_contrast", "rsa_feature", "rsa_connect",
                           "rsa_mediate", "encode"),
  one_run_per_fold = NULL,
  train_test_disjoint_runs = c("transfer", "rsa_era")
)

.dsl_check_assertions <- function(spec, train_design, test_design = NULL) {
  folds <- spec$config$folds
  vb <- spec$config[[spec$verb]]
  out <- list()
  for (a in folds$assert %||% character()) {
    applies <- .dsl_assertion_applies[[a]]
    if (!is.null(applies) && !spec$verb %in% applies) {
      out[[a]] <- list(name = a, pass = FALSE,
                       detail = sprintf("does not apply to %s (applies to: %s)", spec$verb, paste(applies, collapse = ", ")))
      next
    }
    by <- train_design[[folds$by]]
    out[[a]] <- switch(a,
      classes_in_every_fold = {
        y <- train_design[[vb$target]]
        tab <- table(by, y)
        bad <- rownames(tab)[apply(tab, 1, function(r) any(r == 0))]
        list(name = a, pass = !length(bad),
             detail = if (length(bad)) sprintf("%s %s missing a class", folds$by, paste(bad, collapse = ", ")) else
               sprintf("min per class per fold: %d", min(tab)))
      },
      no_item_across_folds = {
        it <- train_design[[vb$item]]
        n <- tapply(by, it, function(b) length(unique(b)))
        bad <- names(n)[n > 1]
        list(name = a, pass = !length(bad),
             detail = if (length(bad)) sprintf("%d item(s) appear in more than one %s (e.g. %s)", length(bad),
                                                folds$by, bad[1]) else "each item in one fold")
      },
      one_run_per_fold = list(name = a, pass = TRUE,
                              detail = sprintf("%d folds, one %s each", length(unique(by)), folds$by)),
      train_test_disjoint_runs = {
        shared <- if (is.null(test_design) || !folds$by %in% names(test_design)) character() else
          intersect(unique(as.character(by)), unique(as.character(test_design[[folds$by]])))
        list(name = a, pass = !length(shared),
             detail = if (length(shared)) sprintf("train and test share %s %s", folds$by,
                                                   paste(shared, collapse = ", ")) else "disjoint")
      }
    )
  }
  unname(out)
}

# ---- cell construction -----------------------------------------------------

.dsl_cell_sets <- function(spec) {
  vb <- spec$config[[spec$verb]]
  if (!is.null(vb$train_set)) c(train = vb$train_set, test = vb$test_set) else c(train = vb$set)
}

.dsl_path_keys <- function(node, x, prefix = "") {
  out <- character()
  if (is.null(x)) return(out)
  if (node$type %in% c("path")) return(stats::setNames(x, prefix))
  if (node$type == "map") {
    for (k in intersect(names(x), names(node$fields))) {
      out <- c(out, .dsl_path_keys(node$fields[[k]], x[[k]], paste0(prefix, if (nzchar(prefix)) "." else "", k)))
    }
  } else if (node$type == "dict") {
    for (k in names(x)) out <- c(out, .dsl_path_keys(node$value, x[[k]], paste0(prefix, ".", k)))
  }
  out
}

# Files an analysis references beyond the beta series (RDMs, features, atlas).
.dsl_referenced_files <- function(spec) {
  node <- .dsl_analysis_node()
  keys <- c(.dsl_path_keys(node$fields[[spec$verb]], spec$config[[spec$verb]], spec$verb),
            .dsl_path_keys(node$fields$scope, spec$config$scope, "scope"))
  vapply(keys, .dsl_resolve_path, "", root = spec$project$root)
}

#' Build the model specification for one subject
#'
#' Loads the subject's data for an analysis and compiles it to an ordinary
#' rMVPA `model_spec`, exactly as the run would. Useful for inspecting or
#' running the compiled analysis by hand.
#'
#' @inheritParams rmvpa_analysis
#' @param subject Subject label, without the `sub-` prefix.
#' @return A list with `model_spec`, the loaded `train` (and `test`) sets,
#'   `mask`, fold `assertions` and compile `info`.
#' @export
rmvpa_compile <- function(project = ".", analysis, subject) {
  spec <- if (inherits(analysis, "rmvpa_analysis_spec")) analysis else rmvpa_analysis(project, analysis)
  .dsl_build_cell(spec, subject)
}

.dsl_build_cell <- function(spec, subject) {
  root <- spec$project$root
  data <- spec$data
  sets <- .dsl_cell_sets(spec)
  mask_file <- .dsl_mask_file(data, subject, root)
  mask <- methods::as(neuroim2::read_vol(mask_file), "LogicalNeuroVol")
  train <- .dsl_load_set(data, sets[["train"]], subject, root, mask)
  test <- if (!is.na(sets["test"])) .dsl_load_set(data, sets[["test"]], subject, root, mask) else NULL
  dataset <- mvpa_dataset(train_data = train$vec, test_data = test$vec, mask = mask)

  where <- sprintf("%s, sub-%s:", spec$name, subject)
  ctx <- list(spec = spec, vb = spec$config[[spec$verb]], cfg = spec$config, root = root,
              subject = subject, dataset = dataset, train = train, test = test,
              mask = mask, where = where)
  compiled <- .dsl_compilers()[[spec$verb]](ctx)
  list(
    model_spec = compiled$model_spec,
    info = compiled$info %||% list(),
    rsa_joint = isTRUE(compiled$rsa_joint),
    train = train, test = test, mask = mask, mask_file = mask_file,
    assertions = .dsl_check_assertions(spec, train$design, test$design),
    inputs = c(mask_file, train$files, test$files, .dsl_referenced_files(spec))
  )
}

# ---- running a scope -------------------------------------------------------

.dsl_execute <- function(spec, cell, draft = FALSE) {
  cfg <- spec$config
  ms <- cell$model_spec
  backend <- cfg$run$backend
  kind <- spec$scope_kind
  inf <- cfg$inference
  use_perm <- !draft && identical(spec$inference_support$kind, "searchlight_permutation") &&
    identical(inf$test, "permutation")

  if (kind == "searchlight") {
    sl <- cfg$scope$searchlight
    sampling <- if (draft) "randomized" else sl$sampling
    iters <- if (draft) 1L else sl$iterations
    if (use_perm) {
      if (identical(sampling, "resampled")) {
        stop("inference: the permutation engine supports searchlight sampling standard or randomized, not resampled.",
             call. = FALSE)
      }
      pc <- .dsl_permutation_control(inf, seed = cfg$run$seed, rsa_joint = cell$rsa_joint)
      return(run_permutation_searchlight(ms, radius = sl$radius, method = sampling, perm_ctrl = pc))
    }
    return(run_searchlight(ms, radius = sl$radius, method = sampling, niter = iters, backend = backend))
  }
  if (kind == "regions") {
    atlas <- neuroim2::read_vol(.dsl_resolve_path(cfg$scope$regions$atlas, spec$project$root))
    return(run_regional(ms, region_mask = atlas, backend = backend))
  }
  run_global(ms)
}

# ---- outputs ---------------------------------------------------------------

.dsl_out_dir <- function(spec, draft = FALSE) {
  file.path(spec$project$root, "derivatives", if (draft) "rmvpa-draft" else "rmvpa", spec$name)
}

.dsl_prefix <- function(spec, subject) {
  sets <- .dsl_cell_sets(spec)
  set <- spec$data$sets[[sets[["train"]]]]
  ents <- c(sub = subject, ses = set$ses, task = set$task, space = spec$data$space,
            desc = gsub("[^A-Za-z0-9]", "", spec$name))
  paste(sprintf("%s-%s", names(ents), ents), collapse = "_")
}

.dsl_stat_label <- function(x) {
  x <- gsub("[^A-Za-z0-9]", "", x)
  if (!nzchar(x)) "value" else x
}

.dsl_write_vol <- function(vol, file) {
  neuroim2::write_vol(vol, file)
  basename(file)
}

.dsl_write_table <- function(df, file) {
  utils::write.table(df, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "n/a")
  basename(file)
}

.dsl_write_outputs <- function(spec, cell, result, dir, prefix) {
  written <- character()
  maps_from <- function(res, tag = "") {
    if (is.null(res$results)) return(character())
    out <- character()
    for (m in names(res$results)) {
      v <- res$results[[m]]
      if (inherits(v, "NeuroVol")) {
        out <- c(out, .dsl_write_vol(v, file.path(dir, sprintf("%s_stat-%s%s_statmap.nii.gz", prefix,
                                                               .dsl_stat_label(m), tag))))
      }
    }
    out
  }
  if (inherits(result, "permutation_result_set") || inherits(result, "permutation_result")) {
    set <- if (inherits(result, "permutation_result")) list(result) else result
    written <- c(written, maps_from(set[[1]]$observed))
    for (r in set) {
      lab <- .dsl_stat_label(r$metric)
      if (inherits(r$p_map, "NeuroVol")) {
        written <- c(written, .dsl_write_vol(r$p_map, file.path(dir, sprintf("%s_stat-%sp_statmap.nii.gz", prefix, lab))))
      }
      if (inherits(r$p_adj_map, "NeuroVol")) {
        written <- c(written, .dsl_write_vol(r$p_adj_map, file.path(dir, sprintf("%s_stat-%spfdr_statmap.nii.gz", prefix, lab))))
      }
    }
  } else if (inherits(result, "searchlight_result")) {
    written <- c(written, maps_from(result))
  } else if (inherits(result, "regional_mvpa_result")) {
    pt <- as.data.frame(result$performance_table)
    labels <- spec$config$scope$regions$labels
    if (!is.null(labels) && "roinum" %in% names(pt)) {
      lt <- .dsl_read_table(.dsl_resolve_path(labels, spec$project$root))
      names(lt)[1] <- "roinum"
      pt <- merge(lt, pt, by = "roinum", all.y = TRUE)
    }
    written <- c(written, .dsl_write_table(pt, file.path(dir, sprintf("%s_regions.tsv", prefix))))
    if (!is.null(result$vol_results)) {
      for (m in names(result$vol_results)) {
        v <- result$vol_results[[m]]
        if (inherits(v, "NeuroVol")) {
          written <- c(written, .dsl_write_vol(v, file.path(dir, sprintf("%s_stat-%s_statmap.nii.gz", prefix,
                                                                         .dsl_stat_label(m)))))
        }
      }
    }
  } else if (inherits(result, "global_mvpa_result")) {
    perf <- result$performance_table %||% result$performance
    df <- if (is.data.frame(perf)) perf else as.data.frame(as.list(unlist(perf)))
    written <- c(written, .dsl_write_table(df, file.path(dir, sprintf("%s_global.tsv", prefix))))
  } else {
    saveRDS(result, file.path(dir, sprintf("%s_result.rds", prefix)))
    written <- c(written, sprintf("%s_result.rds", prefix))
  }
  written
}

.dsl_chance <- function(spec, cell) {
  vb <- spec$config[[spec$verb]]
  if (!spec$verb %in% c("decode", "transfer", "decode_trialwise", "pattern")) return(NULL)
  y <- cell$train$design[[vb$target]]
  if (is.numeric(y) && !identical(vb$target_type, "classification")) return(NULL)
  k <- length(unique(y))
  list(accuracy = 1 / k, auc = 0.5, classes = k)
}

.dsl_require_json <- function() {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Writing rmvpa sidecars needs the 'jsonlite' package.", call. = FALSE)
  }
}

.dsl_hash <- function(x) digest::digest(x, algo = "sha256")

.dsl_spec_hash <- function(spec) {
  sets <- .dsl_cell_sets(spec)
  .dsl_hash(list(version = spec$project$version, verb = spec$verb, config = spec$config,
                 data = spec$data[setdiff(names(spec$data), c("sets", "subjects"))],
                 sets = spec$data$sets[unname(sets)]))
}

.dsl_cell_hash <- function(spec, inputs, draft) {
  md5 <- unname(tools::md5sum(unique(inputs)))
  .dsl_hash(list(spec = .dsl_spec_hash(spec), inputs = md5, draft = draft,
                 rmvpa = as.character(utils::packageVersion("rMVPA"))))
}

.dsl_cell_record_file <- function(spec, subject, draft) {
  file.path(.dsl_out_dir(spec, draft), paste0("sub-", subject), ".rmvpa-cell.json")
}

# ---- lock ------------------------------------------------------------------

.dsl_lock_file <- function(project) file.path(project$root, "rmvpa.lock")

.dsl_read_lock <- function(project) {
  f <- .dsl_lock_file(project)
  if (!file.exists(f)) return(list())
  .dsl_require_json()
  lines <- readLines(f, warn = FALSE)
  lapply(lines[nzchar(lines)], jsonlite::fromJSON)
}

.dsl_lock_status <- function(spec) {
  entries <- Filter(function(e) identical(e$analysis, spec$name), .dsl_read_lock(spec$project))
  if (!length(entries)) return("none")
  if (identical(entries[[length(entries)]]$spec_hash, .dsl_spec_hash(spec))) "match" else "mismatch"
}

#' Record a preregistration lock
#'
#' Appends the hash of each resolved analysis to `rmvpa.lock` in the project
#' root. The file is append-only. Later runs stamp each output with whether
#' its analysis still matches the most recent lock entry. Locking is optional
#' and does not affect re-runs or caching.
#'
#' @param project An `rmvpa_project` or a project path.
#' @param analyses Analysis names; all analyses when `NULL`.
#' @param note Optional free-text note stored with each entry.
#' @return The appended entries, invisibly.
#' @export
rmvpa_lock <- function(project = ".", analyses = NULL, note = NULL) {
  if (!inherits(project, "rmvpa_project")) project <- rmvpa_project(project)
  .dsl_require_json()
  analyses <- analyses %||% names(project$analyses)
  entries <- lapply(analyses, function(a) {
    spec <- rmvpa_analysis(project, a)
    list(analysis = spec$name, spec_hash = .dsl_spec_hash(spec), dsl = project$version,
         locked_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), note = note %||% "")
  })
  con <- file(.dsl_lock_file(project), open = "a")
  on.exit(close(con))
  for (e in entries) writeLines(as.character(jsonlite::toJSON(e, auto_unbox = TRUE)), con)
  message(sprintf("Locked %d analysis spec(s) in %s", length(entries), .dsl_lock_file(project)))
  invisible(entries)
}

# ---- run -------------------------------------------------------------------

.dsl_with_workers <- function(workers, expr) {
  if (workers > 1L) {
    old <- future::plan(future::multisession, workers = workers)
    on.exit(future::plan(old), add = TRUE)
  }
  force(expr)
}

.dsl_run_cell <- function(spec, subject, draft = FALSE, force = FALSE) {
  rec_file <- .dsl_cell_record_file(spec, subject, draft)
  dir <- dirname(rec_file)
  cell <- .dsl_build_cell(spec, subject)
  hash <- .dsl_cell_hash(spec, cell$inputs, draft)
  if (!force && file.exists(rec_file)) {
    rec <- jsonlite::fromJSON(rec_file)
    if (identical(rec$cell_hash, hash) && identical(rec$status, "done")) {
      return(list(subject = subject, status = "up to date", dir = dir))
    }
  }
  failed <- Filter(function(a) !isTRUE(a$pass), cell$assertions)
  if (length(failed)) {
    stop(sprintf("%s, sub-%s: fold assertion(s) failed: %s", spec$name, subject,
                 paste(vapply(failed, function(a) sprintf("%s (%s)", a$name, a$detail), ""), collapse = "; ")),
         call. = FALSE)
  }
  .dsl_ensure_dir(dir)
  old <- list.files(dir, full.names = TRUE)
  unlink(old[!grepl("\\.rmvpa-cell\\.json$", old)])

  set.seed(spec$config$run$seed)
  t0 <- Sys.time()
  result <- .dsl_with_workers(spec$config$run$workers, .dsl_execute(spec, cell, draft = draft))
  prefix <- .dsl_prefix(spec, subject)
  written <- .dsl_write_outputs(spec, cell, result, dir, prefix)

  inf <- spec$config$inference
  sidecar <- list(
    rmvpa_dsl = spec$project$version,
    rmvpa_version = as.character(utils::packageVersion("rMVPA")),
    analysis = spec$name, subject = subject, verb = spec$verb,
    draft = draft,
    model_class = class(cell$model_spec)[1],
    scope = spec$config$scope,
    resolved = spec$config[c(spec$verb, "folds")],
    inference = c(if (draft) list(test = "none (draft)") else inf,
                  list(null = spec$inference_support$description,
                       p_floor = if (!draft && identical(inf$test, "permutation")) 1 / (inf$permutations + 1) else NULL)),
    chance = .dsl_chance(spec, cell),
    folds = list(by = spec$config$folds$by, scheme = spec$config$folds$scheme,
                 assertions = cell$assertions),
    rows = list(train = list(before_keep = cell$train$n_rows, kept = cell$train$n_kept),
                test = if (!is.null(cell$test)) list(before_keep = cell$test$n_rows, kept = cell$test$n_kept)),
    compile = cell$info,
    lock = .dsl_lock_status(spec),
    spec_hash = .dsl_spec_hash(spec),
    inputs = stats::setNames(as.list(unname(tools::md5sum(unique(cell$inputs)))),
                             vapply(unique(cell$inputs), function(f) {
                               rel <- sub(paste0("^", normalizePath(spec$project$root), "/?"), "",
                                          normalizePath(f, mustWork = FALSE))
                               rel
                             }, "")),
    outputs = written,
    started = format(t0, "%Y-%m-%dT%H:%M:%S%z"),
    seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  )
  if (!is.null(spec$inference_support$missing) && !draft) {
    sidecar$inference$descriptive_only <- TRUE
  }
  jsonlite::write_json(sidecar, file.path(dir, sprintf("%s_rmvpa.json", prefix)),
                       auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
  jsonlite::write_json(list(cell_hash = hash, status = "done", spec_hash = sidecar$spec_hash,
                            finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
                       rec_file, auto_unbox = TRUE)
  list(subject = subject, status = "ran", dir = dir, outputs = written)
}

.dsl_ensure_dir <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

#' Run analyses of an rmvpa DSL project
#'
#' Checks each analysis, then runs every subject x analysis cell whose
#' resolved spec or input files changed since its last successful run.
#' Outputs go to `derivatives/rmvpa/<analysis>/sub-<label>/` with BIDS
#' entity names and a JSON sidecar per subject.
#'
#' @param project An `rmvpa_project` or a project path.
#' @param analyses Analysis names; all analyses when `NULL`.
#' @param subjects Subject labels; all subjects when `NULL`.
#' @param draft If `TRUE`, run the first subject only, with a randomized
#'   searchlight (one iteration) and no permutations, into
#'   `derivatives/rmvpa-draft/`. Outputs are stamped as drafts.
#' @param force Re-run cells even when they are up to date.
#' @return A data frame with one row per cell (analysis, subject, status, dir),
#'   invisibly.
#' @export
rmvpa_run <- function(project = ".", analyses = NULL, subjects = NULL, draft = FALSE, force = FALSE) {
  if (!inherits(project, "rmvpa_project")) project <- rmvpa_project(project)
  .dsl_require_json()
  analyses <- analyses %||% names(project$analyses)
  if (!length(analyses)) stop("The project has no analyses. Add one with `rmvpa add <verb> <name>`.", call. = FALSE)
  rows <- list()
  for (a in analyses) {
    spec <- rmvpa_analysis(project, a)
    subs <- subjects %||% .dsl_subjects(spec$data, project$root)
    if (draft) subs <- subs[1]
    for (s in subs) {
      message(sprintf("[%s] sub-%s ...", spec$name, s))
      r <- .dsl_run_cell(spec, s, draft = draft, force = force)
      message(sprintf("[%s] sub-%s %s", spec$name, s, r$status))
      rows[[length(rows) + 1L]] <- data.frame(analysis = spec$name, subject = s, status = r$status,
                                              dir = r$dir, stringsAsFactors = FALSE)
    }
  }
  invisible(do.call(rbind, rows))
}

#' Status of every subject x analysis cell
#'
#' Reports whether each cell is `done`, `stale` (its resolved spec or inputs
#' changed since the last run), or `missing`, and whether its analysis matches
#' the most recent lock entry. Reads no neuroimaging data beyond hashing files.
#'
#' @inheritParams rmvpa_run
#' @return A data frame with columns analysis, subject, status, lock.
#' @export
rmvpa_status <- function(project = ".", analyses = NULL) {
  if (!inherits(project, "rmvpa_project")) project <- rmvpa_project(project)
  .dsl_require_json()
  analyses <- analyses %||% names(project$analyses)
  rows <- list()
  for (a in analyses) {
    spec <- rmvpa_analysis(project, a)
    lock <- .dsl_lock_status(spec)
    for (s in .dsl_subjects(spec$data, project$root)) {
      rec_file <- .dsl_cell_record_file(spec, s, FALSE)
      status <- if (!file.exists(rec_file)) "missing" else {
        rec <- jsonlite::fromJSON(rec_file)
        inputs <- tryCatch({
          sets <- .dsl_cell_sets(spec)
          f <- unlist(lapply(unique(unname(sets)), function(st) {
            x <- .dsl_set_files(spec$data, st, s, project$root); c(x$betas, x$designs)
          }))
          c(.dsl_mask_file(spec$data, s, project$root), f, .dsl_referenced_files(spec))
        }, error = function(e) NULL)
        if (is.null(inputs)) "stale" else if (identical(rec$cell_hash, .dsl_cell_hash(spec, inputs, FALSE))) "done" else "stale"
      }
      rows[[length(rows) + 1L]] <- data.frame(analysis = spec$name, subject = s, status = status, lock = lock,
                                              stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}
