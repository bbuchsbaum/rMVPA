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
#   native-<sub>-<md5>   a native-space atlas sampled on that subject's own
#                        fMRIPrep surfaces (same space, no template warp)
# Anatomy and native surfaces are located with bidser when it is installed
# (bids_project(fmriprep = TRUE): search_files(), surface_files()), else by a
# direct scan of the derivatives folders.
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

# Native spaces: no space entity, or T1w/fsnative. fMRIPrep writes native
# outputs without a space entity (or space-T1w) and native surfaces in T1w
# coordinates.
.studio_native_space <- function(space) is.null(space) || !nzchar(space) || space %in% c("T1w", "fsnative", "native")

.studio_space_matches <- function(entity_space, data_space) {
  if (.studio_native_space(data_space)) is.null(entity_space) || entity_space %in% c("T1w", "fsnative")
  else identical(entity_space, data_space)
}

# BIDS roots that may hold fMRIPrep derivatives: the project itself, and the
# dataset above the data's derivatives folder.
.studio_bids_roots <- function(project, data) {
  r <- project$root
  if (!is.null(data$bids)) {
    b <- normalizePath(.dsl_resolve_path(data$bids, project$root), mustWork = FALSE)
    parts <- strsplit(b, .Platform$file.sep, fixed = TRUE)[[1]]
    at <- which(parts == "derivatives")
    if (length(at)) r <- c(r, paste(parts[seq_len(at[length(at)] - 1L)], collapse = .Platform$file.sep))
  }
  unique(normalizePath(r[dir.exists(r)], mustWork = FALSE))
}

.studio_bidser_project <- function(root) {
  if (!requireNamespace("bidser", quietly = TRUE)) return(NULL)
  tryCatch(suppressMessages(suppressWarnings(
    bidser::bids_project(root, fmriprep = TRUE, strict_participants = FALSE))), error = function(e) NULL)
}

# Anatomical-side files of one subject (T1w images, surfaces, shape maps):
# through bidser when it can read the dataset, else a direct scan of
# derivatives/*/sub-<s>/[ses-*/]anat.
.studio_subject_files <- function(project, data, subject, regex) {
  roots <- .studio_bids_roots(project, data)
  for (root in roots) {
    bp <- .studio_bidser_project(root)
    if (is.null(bp)) next
    f <- tryCatch(bidser::search_files(bp, regex = regex, subid = subject, full_path = TRUE, strict = FALSE),
                  error = function(e) NULL)
    if (length(f)) return(unique(normalizePath(f)))
  }
  dirs <- unique(normalizePath(c(if (!is.null(data$bids)) .dsl_resolve_path(data$bids, project$root),
                                 unlist(lapply(roots, function(r) list.dirs(file.path(r, "derivatives"), recursive = FALSE)))),
                               mustWork = FALSE))
  dirs <- dirs[dir.exists(dirs) & !basename(dirs) %in% c("rmvpa", "rmvpa-draft")]
  sdirs <- file.path(dirs, paste0("sub-", subject))
  anat <- unlist(lapply(sdirs, function(s) c(file.path(s, "anat"), Sys.glob(file.path(s, "ses-*", "anat")))))
  unique(unlist(lapply(anat[dir.exists(anat)], list.files, pattern = regex, full.names = TRUE)))
}

# T1w images in the data space for the first subject.
.studio_find_anat <- function(project, data, subject) {
  f <- .studio_subject_files(project, data, subject, "_T1w\\.nii(\\.gz)?$")
  if (!length(f)) return(character())
  ents <- lapply(f, function(x) .dsl_parse_bids_name(x)$entities)
  same_space <- vapply(ents, function(e) .studio_space_matches(e$space, data$space), TRUE)
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
  if (.studio_native_space(space)) {
    s1 <- .dsl_subjects(spec$data, project$root)[1]
    ns <- .studio_native_surfaces(project, spec$data, s1)
    if (!is.null(ns)) {
      cached <- file.path(.studio_cache_dir("projected"), sprintf("native-%s-%s.rds", s1, md5))
      return(list(id = sprintf("native-%s-%s", s1, md5), template = sprintf("sub-%s", s1), exact = FALSE, tag = "sampled",
                  status = if (file.exists(cached)) "ready" else "fetch", size = "computed here",
                  label = sprintf("sub-%s's own cortical surface", s1),
                  source = paste0("sub-", s1, "'s pial and white surfaces from fMRIPrep (", ns$found_by, "). Each vertex takes ",
                                  "the atlas label at its mid-thickness position; atlas and surface share the subject's ",
                                  "T1w space, so no template warp is involved, but labels near parcel borders follow ",
                                  "the volume's voxels. Values are the group means drawn on this one subject's anatomy."),
                  atlas_md5 = md5, subject = s1, files = ns))
    }
    return(list(id = NULL, status = "unavailable", exact = FALSE,
                reason = sprintf(paste("The atlas is in the native (T1w) space and no fMRIPrep surfaces were found for sub-%s",
                                       "(sub-%s_hemi-L_pial.surf.gii and _white.surf.gii in an anat/ folder)."), s1, s1)))
  }
  if (grepl("^MNI", space)) {
    ready <- have_fs5 && file.exists(file.path(.studio_cache_dir("projected"), paste0(md5, ".rds")))
    return(list(id = paste0("project-", md5), template = "fsaverage5", exact = FALSE, tag = "approximate projection",
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

# fMRIPrep surfaces of one subject: pial and white per hemisphere (required),
# sulc (optional; shading is flat without it).
.studio_native_surfaces <- function(project, data, subject) {
  f <- .studio_subject_files(project, data, subject, "_hemi-[LR]_.*(pial|white|smoothwm|midthickness)\\.surf\\.gii$|_hemi-[LR]_.*sulc\\.shape\\.gii$")
  if (!length(f)) return(NULL)
  pick <- function(H, kind) {
    hit <- f[grepl(sprintf("_hemi-%s_", H), basename(f)) & grepl(sprintf("_%s\\.(surf|shape)\\.gii$", kind), basename(f))]
    if (length(hit)) hit[1] else NULL
  }
  out <- list()
  for (H in c("L", "R")) {
    out[[paste0(H, "_pial")]] <- pick(H, "pial")
    out[[paste0(H, "_white")]] <- pick(H, "white") %||% pick(H, "smoothwm")
    out[[paste0(H, "_sulc")]] <- pick(H, "sulc")
  }
  if (any(vapply(c("L_pial", "L_white", "R_pial", "R_white"), function(k) is.null(out[[k]]), TRUE))) return(NULL)
  out$found_by <- if (requireNamespace("bidser", quietly = TRUE)) "located with bidser" else "found in the derivatives folders"
  out
}

.studio_read_native <- function(files) {
  hemi <- function(H) {
    pial <- gifti::read_gifti(files[[paste0(H, "_pial")]])$data
    white <- gifti::read_gifti(files[[paste0(H, "_white")]])$data
    sulc <- if (!is.null(files[[paste0(H, "_sulc")]])) as.numeric(gifti::read_gifti(files[[paste0(H, "_sulc")]])$data[[1]]) else rep(0, nrow(pial$pointset))
    list(pos = pial$pointset, faces = pial$triangle, white = white$pointset, sulc = sulc)
  }
  list(left = hemi("L"), right = hemi("R"))
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
.studio_project_atlas <- function(atlas, mid, side = c("left", "right"), midline = 0) {
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
  allowed[as.integer(names(cx)) + 1L] <- as.integer(if (side == "left") cx <= midline + 4 else cx >= midline - 4)
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
  list(id = opt$id, template = "fsaverage5", mesh = "pial", exact = opt$exact, tag = opt$tag, source = opt$source,
       atlas_md5 = opt$atlas_md5, left = hemi("left", labels$left), right = hemi("right", labels$right))
}

# Builds (fetching what is missing when `fetch`) the surface for an option.
.studio_surface_build <- function(project, spec, opt, fetch = FALSE) {
  if (identical(opt$status, "unavailable")) stop(opt$reason, call. = FALSE)
  rg <- spec$config$scope$regions
  afile <- .dsl_resolve_path(rg$atlas, project$root)
  if (startsWith(opt$id, "native-")) {
    mesh <- .studio_read_native(opt$files)
    cached <- file.path(.studio_cache_dir("projected"), sprintf("native-%s-%s.rds", opt$subject, opt$atlas_md5))
    labels <- if (file.exists(cached)) readRDS(cached) else {
      atlas <- neuroim2::read_vol(afile)
      # Midline from the subject's own surfaces (scanner x need not be 0 there).
      midline <- (max(mesh$left$pos[, 1]) + min(mesh$right$pos[, 1])) / 2
      l <- lapply(c(left = "left", right = "right"), function(h)
        .studio_project_atlas(atlas, (mesh[[h]]$pos + mesh[[h]]$white) / 2, h, midline = midline))
      saveRDS(l, cached)
      l
    }
    out <- .studio_surface_payload(opt, mesh, labels)
    out$template <- opt$template
    return(out)
  }
  if (fetch) .studio_fetch_fs5() else if (!all(file.exists(.studio_fs5_paths()))) stop("fsaverage5 is not fetched yet", call. = FALSE)
  fs5 <- .studio_read_fs5()
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
