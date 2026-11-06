# rmvpa DSL: finding and loading beta series.
#
# BIDS convention (data.bids): for subject XX and a set with entities
# {task, ses, ...}, the beta images are files under <bids>/sub-XX/ whose names
# carry those entities, data.space as `space-`, and end in
# `_<suffix>.nii[.gz]` (suffix defaults to `betas`). Each image has a design
# table with the same name ending in `_<suffix>.tsv`, one row per volume.
# Several matching files (one per run) are concatenated in run order, and a
# `run` column is added from the run- entity when the table lacks one.
#
# Explicit templates (data.files) replace discovery: {subject} and the set's
# entities are substituted, and the result may contain * wildcards.

.dsl_parse_bids_name <- function(fname) {
  base <- basename(fname)
  ext <- regmatches(base, regexpr("\\.(nii\\.gz|nii|tsv)$", base))
  stem <- sub("\\.(nii\\.gz|nii|tsv)$", "", base)
  parts <- strsplit(stem, "_", fixed = TRUE)[[1]]
  kv <- parts[grepl("-", parts, fixed = TRUE)]
  ents <- stats::setNames(sub("^[^-]*-", "", kv), sub("-.*$", "", kv))
  suffix <- if (length(parts) && !grepl("-", parts[length(parts)], fixed = TRUE)) parts[length(parts)] else ""
  list(entities = as.list(ents), suffix = suffix, ext = ext)
}

.dsl_resolve_path <- function(p, root) {
  if (grepl("^(/|~|[A-Za-z]:)", p)) path.expand(p) else file.path(root, p)
}

.dsl_subjects <- function(data, root) {
  if (!is.null(data$subjects)) return(as.character(data$subjects))
  if (!is.null(data$bids)) {
    d <- .dsl_resolve_path(data$bids, root)
    if (!dir.exists(d)) stop("data.bids: directory not found: ", d, call. = FALSE)
    dirs <- basename(list.dirs(d, recursive = FALSE))
    subs <- sub("^sub-", "", dirs[grepl("^sub-", dirs)])
    if (!length(subs)) stop("data.bids: no sub-* folders in ", d, call. = FALSE)
    return(sort(subs))
  }
  stop("data.subjects is required when data.files templates are used.", call. = FALSE)
}

.dsl_set_entities_of <- function(set) set[intersect(names(set), .dsl_set_entities)]

.dsl_fill_template <- function(tpl, subject, set, root) {
  vals <- c(list(subject = subject), .dsl_set_entities_of(set))
  out <- tpl
  for (k in names(vals)) out <- gsub(paste0("{", k, "}"), vals[[k]], out, fixed = TRUE)
  if (grepl("\\{[A-Za-z_]+\\}", out)) {
    stop(sprintf("data.files: template '%s' has a placeholder the set does not define", tpl), call. = FALSE)
  }
  .dsl_resolve_path(out, root)
}

.dsl_set_files <- function(data, set_name, subject, root) {
  set <- data$sets[[set_name]]
  if (!is.null(data$files)) {
    betas <- sort(Sys.glob(.dsl_fill_template(data$files$betas, subject, set, root)))
    designs <- sort(Sys.glob(.dsl_fill_template(data$files$design, subject, set, root)))
    if (!length(betas)) stop(sprintf("set '%s', sub-%s: no file matches data.files.betas", set_name, subject), call. = FALSE)
    if (length(designs) != length(betas)) {
      stop(sprintf("set '%s', sub-%s: %d beta files but %d design tables", set_name, subject,
                   length(betas), length(designs)), call. = FALSE)
    }
    runs <- vapply(betas, function(f) .dsl_parse_bids_name(f)$entities$run %||% NA_character_, "")
    return(list(betas = betas, designs = designs, runs = unname(runs)))
  }

  sdir <- file.path(.dsl_resolve_path(data$bids, root), paste0("sub-", subject))
  cand <- list.files(sdir, pattern = "\\.nii(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  want <- .dsl_set_entities_of(set)
  if (!is.null(data$space)) want$space <- data$space
  hit <- vapply(cand, function(f) {
    p <- .dsl_parse_bids_name(f)
    identical(p$suffix, data$suffix) && identical(p$entities$sub, subject) &&
      all(vapply(names(want), function(k) identical(p$entities[[k]], want[[k]]), TRUE))
  }, TRUE)
  betas <- cand[hit]
  if (!length(betas)) {
    stop(sprintf("set '%s', sub-%s: no *_%s.nii[.gz] under %s with %s", set_name, subject, data$suffix, sdir,
                 paste(sprintf("%s-%s", names(want), unlist(want)), collapse = ", ")), call. = FALSE)
  }
  runs <- vapply(betas, function(f) .dsl_parse_bids_name(f)$entities$run %||% NA_character_, "")
  ord <- order(suppressWarnings(as.numeric(runs)), runs, betas)
  betas <- betas[ord]; runs <- unname(runs[ord])
  designs <- sub("\\.nii(\\.gz)?$", ".tsv", betas)
  missing <- designs[!file.exists(designs)]
  if (length(missing)) stop("Design table not found next to its beta image: ", missing[1], call. = FALSE)
  list(betas = unname(betas), designs = unname(designs), runs = runs)
}

.dsl_mask_file <- function(data, subject, root) {
  m <- data$mask
  if (grepl("[/{]|\\.nii", m)) {
    f <- Sys.glob(.dsl_fill_template(m, subject, list(), root))
    if (length(f) != 1L) stop(sprintf("data.mask: '%s' matches %d files for sub-%s", m, length(f), subject), call. = FALSE)
    return(f)
  }
  if (is.null(data$bids)) stop("data.mask must be a path or template when data.files is used.", call. = FALSE)
  broot <- .dsl_resolve_path(data$bids, root)
  pick <- function(dir, recursive) {
    cand <- list.files(dir, pattern = "_mask\\.nii(\\.gz)?$", recursive = recursive, full.names = TRUE)
    cand[vapply(cand, function(f) {
      p <- .dsl_parse_bids_name(f)
      identical(p$entities$desc, m) &&
        (is.null(data$space) || identical(p$entities$space, data$space)) &&
        (is.null(p$entities$sub) || identical(p$entities$sub, subject))
    }, TRUE)]
  }
  f <- pick(file.path(broot, paste0("sub-", subject)), TRUE)
  if (!length(f)) f <- pick(broot, FALSE)
  if (length(f) != 1L) {
    stop(sprintf("data.mask: found %d files *_desc-%s_mask.nii[.gz]%s for sub-%s under %s",
                 length(f), m, if (!is.null(data$space)) paste0(" (space-", data$space, ")") else "",
                 subject, broot), call. = FALSE)
  }
  f
}

.dsl_read_table <- function(f) {
  sep <- if (grepl("\\.csv$", f, ignore.case = TRUE)) "," else "\t"
  utils::read.table(f, header = TRUE, sep = sep, quote = "\"", comment.char = "",
                    stringsAsFactors = FALSE, check.names = FALSE)
}

.dsl_apply_keep <- function(design, keep, where) {
  if (is.null(keep)) return(rep(TRUE, nrow(design)))
  ok <- rep(TRUE, nrow(design))
  for (col in names(keep)) {
    if (!col %in% names(design)) {
      stop(sprintf("%s keep: design has no column '%s'%s", where, col, .dsl_did_you_mean(col, names(design))),
           call. = FALSE)
    }
    rule <- keep[[col]]
    vals <- if (.dsl_is_map(rule)) rule[["in"]] else rule
    ok <- ok & as.character(design[[col]]) %in% as.character(vals)
  }
  ok
}

# Design table of one set for one subject, concatenated across runs, with the
# `keep` mask. Reads no image data.
.dsl_load_design <- function(data, set_name, subject, root, files = NULL) {
  files <- files %||% .dsl_set_files(data, set_name, subject, root)
  where <- sprintf("set '%s', sub-%s:", set_name, subject)
  designs <- lapply(seq_along(files$designs), function(i) {
    d <- .dsl_read_table(files$designs[i])
    if (!"run" %in% names(d) && !is.na(files$runs[i])) d$run <- files$runs[i]
    d
  })
  design <- do.call(rbind, designs)
  keep <- .dsl_apply_keep(design, data$sets[[set_name]]$keep, where)
  if (!any(keep)) stop(sprintf("%s keep filters remove every row", where), call. = FALSE)
  list(design = design, keep = keep, per_file = vapply(designs, nrow, 1L), files = files, where = where)
}

# Load one set for one subject: concatenated beta series (masked) and its
# design table after `keep` filtering.
.dsl_load_set <- function(data, set_name, subject, root, mask) {
  dd <- .dsl_load_design(data, set_name, subject, root)
  files <- dd$files
  vecs <- lapply(seq_along(files$betas), function(i) {
    v <- neuroim2::read_vec(files$betas[i], mask = mask)
    nv <- dim(v)[4]
    if (nv != dd$per_file[i]) {
      stop(sprintf("%s %s has %d volumes but its design table has %d rows", dd$where,
                   basename(files$betas[i]), nv, dd$per_file[i]), call. = FALSE)
    }
    v
  })
  vec <- if (length(vecs) == 1L) vecs[[1]] else do.call(neuroim2::concat, vecs)
  design <- dd$design
  keep <- dd$keep
  if (!all(keep)) {
    vec <- neuroim2::sub_vector(vec, which(keep))
    design <- design[keep, , drop = FALSE]
  }
  rownames(design) <- NULL
  list(vec = vec, design = design, files = c(files$betas, files$designs),
       n_rows = length(keep), n_kept = sum(keep))
}

.dsl_require_columns <- function(design, cols, where) {
  for (col in cols) {
    if (!col %in% names(design)) {
      stop(sprintf("%s design has no column '%s'%s", where, col, .dsl_did_you_mean(col, names(design))),
           call. = FALSE)
    }
  }
}

# Square item x item table: first column holds item names unless the table
# already has row names matching its column names.
.dsl_read_rdm <- function(path, root, items, what) {
  f <- .dsl_resolve_path(path, root)
  if (!file.exists(f)) stop(sprintf("%s: file not found: %s", what, f), call. = FALSE)
  t <- .dsl_read_table(f)
  rn <- as.character(t[[1]])
  m <- as.matrix(t[, -1, drop = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- rn
  if (!identical(sort(rownames(m)), sort(colnames(m)))) {
    stop(sprintf("%s: %s must be square with matching row and column item names", what, f), call. = FALSE)
  }
  miss <- setdiff(items, rownames(m))
  if (length(miss)) {
    stop(sprintf("%s: %s has no row for item(s) %s", what, f, paste(utils::head(miss, 5), collapse = ", ")),
         call. = FALSE)
  }
  m[items, items, drop = FALSE]
}

# Item x feature table: first column holds item names.
.dsl_read_features <- function(path, root, items, what) {
  f <- .dsl_resolve_path(path, root)
  if (!file.exists(f)) stop(sprintf("%s: file not found: %s", what, f), call. = FALSE)
  t <- .dsl_read_table(f)
  rn <- as.character(t[[1]])
  m <- as.matrix(t[, -1, drop = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- rn
  miss <- setdiff(items, rn)
  if (length(miss)) {
    stop(sprintf("%s: %s has no row for item(s) %s", what, f, paste(utils::head(miss, 5), collapse = ", ")),
         call. = FALSE)
  }
  m
}
