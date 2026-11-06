# Images the studio draws under the results, found in the project or fetched
# on request, so a user never has to name a file.
#
# Backgrounds (the grey anatomy under every slice):
#   mask                 the analysis mask (always there)
#   file:<path>          an anatomical T1w in the data space, found next to the data
#   templateflow:<tpl>   the TemplateFlow template for the data space, fetched once
# Cortical surfaces for regional results:
#   schaefer-<n>-<k>     CBIG's own fsaverage5 annotation of a Schaefer-2018 atlas (exact)
#   project-<md5>        the volumetric atlas sampled at fsaverage5 mid-thickness
#                        (approximate; only for atlases in an MNI space)
# Downloads land in tools::R_user_dir("rmvpadsl", "cache") (or $RMVPADSL_CACHE)
# and happen only when the user asks: the studio's Fetch buttons, or
# rmvpa_export_studio(fetch = TRUE).

.studio_cache_dir <- function(...) {
  d <- file.path(Sys.getenv("RMVPADSL_CACHE", tools::R_user_dir("rmvpadsl", "cache")), ...)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

.studio_download <- function(url, dest) {
  if (file.exists(dest) && file.size(dest) > 0) return(dest)
  tmp <- paste0(dest, ".part")
  ok <- tryCatch(utils::download.file(url, tmp, mode = "wb", quiet = TRUE) == 0L, error = function(e) FALSE,
                 warning = function(w) FALSE)
  if (!ok || !file.exists(tmp) || file.size(tmp) == 0) {
    unlink(tmp)
    stop(sprintf("could not download %s (offline, or a proxy refused it)", url), call. = FALSE)
  }
  file.rename(tmp, dest)
  dest
}

# ---- backgrounds ----------------------------------------------------------------

# TemplateFlow images at 2 mm, brain-extracted where TemplateFlow has one.
.studio_tf_templates <- list(
  MNI152NLin2009cAsym = list(file = "tpl-MNI152NLin2009cAsym_res-02_desc-brain_T1w.nii.gz", size = "1.8 MB"),
  MNI152NLin6Asym = list(file = "tpl-MNI152NLin6Asym_res-02_desc-brain_T1w.nii.gz", size = "1.4 MB"),
  MNI152Lin = list(file = "tpl-MNI152Lin_res-02_T1w.nii.gz", size = "1.4 MB")
)

.studio_tf_file <- function(space) file.path(.studio_cache_dir("templateflow"), .studio_tf_templates[[space]]$file)

# T1w images in the data space for the first subject, in the data's BIDS
# folder or any other derivatives pipeline (fMRIPrep's anat/ outputs).
.studio_find_anat <- function(project, data, subject) {
  roots <- unique(normalizePath(c(if (!is.null(data$bids)) .dsl_resolve_path(data$bids, project$root),
                                  list.dirs(file.path(project$root, "derivatives"), recursive = FALSE)),
                                mustWork = FALSE))
  roots <- roots[dir.exists(roots) & basename(roots) != "rmvpa" & basename(roots) != "rmvpa-draft"]
  sdirs <- file.path(roots, paste0("sub-", subject))
  anat <- unlist(lapply(sdirs, function(s) c(file.path(s, "anat"), Sys.glob(file.path(s, "ses-*", "anat")))))
  f <- unlist(lapply(anat[dir.exists(anat)], list.files, pattern = "_T1w\\.nii(\\.gz)?$", full.names = TRUE))
  if (!length(f)) return(character())
  ents <- lapply(f, function(x) .dsl_parse_bids_name(x)$entities)
  same_space <- vapply(ents, function(e) {
    if (is.null(data$space)) is.null(e$space) else identical(e$space, data$space)
  }, TRUE)
  f <- f[same_space]; ents <- ents[same_space]
  if (!length(f)) return(character())
  rank <- vapply(ents, function(e) if (identical(e$desc, "preproc")) 0L else if (is.null(e$desc)) 1L else 2L, 1L)
  f[order(rank, f)]
}

.studio_backgrounds <- function(project) {
  data <- utils::modifyList(.dsl_defaults(project$version)$data, project$config$data)
  s1 <- .dsl_subjects(data, project$root)[1]
  out <- list(list(id = "mask", label = "Analysis mask", source = "the brain mask the analysis used", status = "ready"))
  for (f in .studio_find_anat(project, data, s1)) {
    rel <- sub(paste0("^", normalizePath(project$root), "/?"), "", normalizePath(f))
    out[[length(out) + 1L]] <- list(id = paste0("file:", rel), label = sprintf("sub-%s T1w", s1),
                                    source = rel, status = "ready", path = f)
  }
  sp <- data$space
  if (!is.null(sp) && sp %in% names(.studio_tf_templates)) {
    f <- .studio_tf_file(sp)
    out[[length(out) + 1L]] <- list(id = paste0("templateflow:", sp), label = sprintf("%s template", sp),
                                    source = "TemplateFlow", status = if (file.exists(f)) "ready" else "fetch",
                                    size = .studio_tf_templates[[sp]]$size, path = f)
  }
  # Default: a subject's own anatomy, else a cached template, else the mask.
  ready <- Filter(function(b) identical(b$status, "ready"), out)
  pick <- Find(function(b) startsWith(b$id, "file:"), ready) %||% Find(function(b) startsWith(b$id, "templateflow:"), ready) %||% out[[1]]
  for (i in seq_along(out)) out[[i]]$default <- identical(out[[i]]$id, pick$id)
  out
}

.studio_fetch_background <- function(project, id) {
  b <- Find(function(x) identical(x$id, id), .studio_backgrounds(project))
  if (is.null(b)) stop("no background '", id, "'", call. = FALSE)
  if (identical(b$status, "fetch")) {
    sp <- sub("^templateflow:", "", id)
    tf <- .studio_tf_templates[[sp]]
    .studio_download(sprintf("https://templateflow.s3.amazonaws.com/tpl-%s/%s", sp, tf$file), b$path)
  }
  invisible(id)
}

# ---- surfaces -------------------------------------------------------------------

.studio_fs5_files <- c(pial = "pial.surf.gii", white = "white.surf.gii", sulc = "sulc.shape.gii")

.studio_fs5_paths <- function() {
  d <- .studio_cache_dir("fsaverage5")
  unlist(lapply(c("L", "R"), function(H) stats::setNames(
    file.path(d, sprintf("tpl-fsaverage_hemi-%s_den-10k_%s", H, .studio_fs5_files)),
    paste0(H, "_", names(.studio_fs5_files)))))
}

.studio_schaefer_of <- function(atlas_file, lt) {
  m <- regmatches(basename(atlas_file), regexec("Schaefer2018_([0-9]+)Parcels_(7|17)Networks", basename(atlas_file)))[[1]]
  if (length(m) == 3L) return(list(n = as.integer(m[2]), k = as.integer(m[3])))
  nmcol <- if (!is.null(lt)) c(intersect(c("name", "label"), names(lt)), names(lt)[2])[1] else NA
  nm <- if (!is.na(nmcol)) as.character(lt[[nmcol]]) else character()
  hit <- regmatches(nm, regexec("^(7|17)Networks_[LR]H_", nm))
  k <- unique(vapply(hit, function(h) if (length(h) == 2L) h[2] else NA_character_, ""))
  if (length(nm) %in% c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000) && length(k) == 1L && !is.na(k)) {
    return(list(n = length(nm), k = as.integer(k)))
  }
  NULL
}

.studio_annot_paths <- function(sc) {
  d <- .studio_cache_dir("schaefer2018")
  stats::setNames(file.path(d, sprintf("%s.Schaefer2018_%dParcels_%dNetworks_order.annot", c("lh", "rh"), sc$n, sc$k)),
                  c("lh", "rh"))
}

# What the cortex view can offer for one regional analysis.
.studio_surface_option <- function(project, spec, provided = list()) {
  rg <- spec$config$scope$regions
  afile <- .dsl_resolve_path(rg$atlas, project$root)
  md5 <- unname(tools::md5sum(afile))
  given <- Find(function(s) identical(s$atlas_md5, md5), provided)
  if (!is.null(given)) {
    return(list(id = paste0("given-", md5), template = given$template, exact = TRUE, status = "ready",
                label = sprintf("%s (from the surface file given to the studio)", given$template)))
  }
  lt <- if (!is.null(rg$labels)) .dsl_read_table(.dsl_resolve_path(rg$labels, project$root)) else NULL
  have_fs5 <- all(file.exists(.studio_fs5_paths()))
  sc <- .studio_schaefer_of(afile, lt)
  if (!is.null(sc)) {
    ready <- have_fs5 && all(file.exists(.studio_annot_paths(sc)))
    return(list(id = sprintf("schaefer-%d-%d", sc$n, sc$k), template = "fsaverage5", exact = TRUE,
                status = if (ready) "ready" else "fetch", size = "2.9 MB",
                label = sprintf("Schaefer-2018 %d parcels, %d networks on fsaverage5", sc$n, sc$k),
                source = "fsaverage5 pial and sulc from TemplateFlow; parcel labels from CBIG's own fsaverage5 annotation",
                atlas_md5 = md5))
  }
  space <- spec$data$space %||% ""
  if (grepl("^MNI", space)) {
    ready <- have_fs5 && file.exists(file.path(.studio_cache_dir("projected"), paste0(md5, ".rds")))
    return(list(id = paste0("project-", md5), template = "fsaverage5", exact = FALSE,
                status = if (ready) "ready" else "fetch", size = if (have_fs5) "computed here" else "2.6 MB",
                label = "Approximate projection onto fsaverage5",
                source = paste("fsaverage5 meshes from TemplateFlow; each vertex takes the atlas label at its",
                               "mid-thickness position (nearest labelled voxel within 2 voxels). fsaverage",
                               "coordinates are close to, not the same as, MNI152, so borders can shift by a few mm."),
                atlas_md5 = md5))
  }
  list(id = NULL, status = "unavailable", exact = FALSE,
       reason = sprintf("The atlas is in space '%s'. The studio knows how to put only Schaefer-2018 or MNI-space atlases on a cortical surface.",
                        if (nzchar(space)) space else "unknown"))
}

.studio_read_fs5 <- function() {
  p <- .studio_fs5_paths()
  hemi <- function(H) {
    pial <- gifti::read_gifti(p[[paste0(H, "_pial")]])$data
    white <- gifti::read_gifti(p[[paste0(H, "_white")]])$data
    list(pos = pial$pointset, faces = pial$triangle, white = white$pointset,
         sulc = as.numeric(gifti::read_gifti(p[[paste0(H, "_sulc")]])$data[[1]]))
  }
  list(left = hemi("L"), right = hemi("R"))
}

.studio_fetch_fs5 <- function() {
  p <- .studio_fs5_paths()
  for (f in p) .studio_download(paste0("https://templateflow.s3.amazonaws.com/tpl-fsaverage/", basename(f)), f)
  invisible(p)
}

# Vertex labels by sampling a volumetric atlas at mid-thickness: the voxel
# under the vertex, else the most common label among voxels within 2. A
# hemisphere only takes parcels whose centre lies on its side (or within 4 mm
# of the midline), so medial vertices do not borrow the other side's parcels.
.studio_project_atlas <- function(atlas, mid, side = c("left", "right")) {
  side <- match.arg(side)
  sp <- neuroim2::space(atlas)
  a <- as.array(atlas)
  if (length(dim(a)) > 3L) a <- a[, , , 1]
  d <- dim(a)[1:3]
  g <- round(neuroim2::coord_to_grid(sp, mid))
  ids <- sort(unique(as.integer(round(a[a > 0]))))
  vox <- which(round(a) > 0, arr.ind = TRUE)
  cx <- tapply(neuroim2::grid_to_coord(sp, vox)[, 1], as.integer(round(a[vox])), mean)
  allowed <- integer(max(ids, 0L) + 1L)
  allowed[as.integer(names(cx)) + 1L] <- as.integer(if (side == "left") cx <= 4 else cx >= -4)
  look <- function(gg) {
    ok <- rowSums(gg >= 1) == 3 & gg[, 1] <= d[1] & gg[, 2] <= d[2] & gg[, 3] <= d[3]
    v <- integer(nrow(gg)); v[ok] <- as.integer(round(a[gg[ok, , drop = FALSE]]))
    v[v > 0][allowed[v[v > 0] + 1L] == 0L] <- 0L
    v
  }
  lab <- look(g)
  offs <- as.matrix(expand.grid(-2:2, -2:2, -2:2)); offs <- offs[order(rowSums(offs^2)), , drop = FALSE][-1, , drop = FALSE]
  todo <- which(lab == 0L)
  if (length(todo)) {
    votes <- vapply(seq_len(nrow(offs)), function(o) look(sweep(g[todo, , drop = FALSE], 2, offs[o, ], "+")), integer(length(todo)))
    if (is.null(dim(votes))) votes <- matrix(votes, nrow = length(todo))
    lab[todo] <- apply(votes, 1, function(v) { v <- v[v > 0]; if (length(v)) as.integer(names(which.max(table(v)))) else 0L })
  }
  lab
}

.studio_surface_payload <- function(opt, fs5, labels) {
  hemi <- function(h, l) list(n = nrow(fs5[[h]]$pos), positions = .studio_f32(t(fs5[[h]]$pos)),
                              faces = .studio_i32(t(fs5[[h]]$faces)), sulc = .studio_f32(fs5[[h]]$sulc),
                              labels = .studio_i32(l))
  list(id = opt$id, template = "fsaverage5", mesh = "pial", exact = opt$exact, source = opt$source,
       atlas_md5 = opt$atlas_md5, left = hemi("left", labels$left), right = hemi("right", labels$right))
}

# Builds (fetching what is missing when `fetch`) the surface for an option.
.studio_surface_build <- function(project, spec, opt, fetch = FALSE) {
  if (identical(opt$status, "unavailable")) stop(opt$reason, call. = FALSE)
  if (fetch) .studio_fetch_fs5() else if (!all(file.exists(.studio_fs5_paths()))) stop("fsaverage5 is not fetched yet", call. = FALSE)
  fs5 <- .studio_read_fs5()
  rg <- spec$config$scope$regions
  afile <- .dsl_resolve_path(rg$atlas, project$root)
  if (startsWith(opt$id, "schaefer-")) {
    lt <- if (!is.null(rg$labels)) .dsl_read_table(.dsl_resolve_path(rg$labels, project$root)) else NULL
    sc <- .studio_schaefer_of(afile, lt)
    ap <- .studio_annot_paths(sc)
    if (fetch) {
      base <- "https://raw.githubusercontent.com/ThomasYeoLab/CBIG/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage5/label/"
      for (f in ap) .studio_download(paste0(base, basename(f)), f)
    }
    an <- lapply(ap, freesurferformats::read.fs.annot)
    drop <- c("Background+FreeSurfer_Defined_Medial_Wall", "Unknown", "???", "unknown")
    full <- c(setdiff(an$lh$colortable$struct_names, drop), setdiff(an$rh$colortable$struct_names, drop))
    # Parcel ids: the label table's own ids when it names the parcels, else
    # Schaefer's order (LH then RH, as in the volumetric atlas).
    ids_of <- function(nm) {
      nmcol <- if (!is.null(lt)) intersect(c("name", "label"), names(lt))[1] else NA
      if (!is.na(nmcol) && all(full %in% lt[[nmcol]])) lt[[1]][match(nm, lt[[nmcol]])] else match(nm, full)
    }
    labels <- lapply(an, function(a) { v <- ids_of(a$label_names); v[is.na(v)] <- 0L; as.integer(v) })
    names(labels) <- c("left", "right")
  } else {
    cached <- file.path(.studio_cache_dir("projected"), paste0(opt$atlas_md5, ".rds"))
    labels <- if (file.exists(cached)) readRDS(cached) else {
      atlas <- neuroim2::read_vol(afile)
      l <- lapply(c(left = "left", right = "right"), function(h) .studio_project_atlas(atlas, (fs5[[h]]$pos + fs5[[h]]$white) / 2, h))
      saveRDS(l, cached)
      l
    }
  }
  .studio_surface_payload(opt, fs5, labels)
}
