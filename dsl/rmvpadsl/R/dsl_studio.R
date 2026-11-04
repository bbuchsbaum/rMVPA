# The studio: the React app in inst/app, served on the local machine by
# rmvpa_studio() or written with its data to a folder by rmvpa_export_studio().
#
# Both read the same payloads, built here from a project and its outputs on
# disk (the app's src/api.js documents how it finds them):
#   project   analyses, explain rows, check findings, cell status
#   template  the display grid: dims, affine and a grey background
#   results/<analysis>   per-subject maps or region tables, summarised
#   surface/<template>   optional cortical meshes with per-vertex parcel ids
#
# Nothing here fits a model. Results come from rmvpa_run()'s outputs only.

# ---- small helpers ------------------------------------------------------------

.studio_b64 <- function(x, size = 4L) jsonlite::base64_enc(writeBin(x, raw(), size = size))
.studio_f32 <- function(x) { x <- as.numeric(x); x[!is.finite(x)] <- NaN; .studio_b64(x, 4L) }
.studio_i32 <- function(x) .studio_b64(as.integer(x), 4L)
.studio_json <- function(x) jsonlite::toJSON(x, auto_unbox = TRUE, digits = NA, null = "null", na = "null")
.studio_round <- function(x, d = 5) if (is.null(x)) NULL else signif(x, d)

.studio_app_dir <- function() {
  d <- system.file("app", package = "rmvpadsl")
  if (!nzchar(d) || !file.exists(file.path(d, "index.html"))) {
    stop("The studio app is not built into this installation of rmvpadsl. From the repository, run ",
         "`npm ci && npm run build` in dsl/app, then reinstall dsl/rmvpadsl.", call. = FALSE)
  }
  d
}

.studio_tilde <- function(path) {
  home <- normalizePath("~", mustWork = FALSE)
  if (startsWith(path, home)) paste0("~", substring(path, nchar(home) + 1L)) else path
}

# Nearest-neighbour lookup of `vol` at the voxel centres of `space`, through
# world coordinates, so atlases and backgrounds on other grids line up.
.studio_on_grid <- function(vol, space) {
  if (identical(dim(vol)[1:3], dim(space)[1:3]) &&
      isTRUE(all.equal(neuroim2::trans(neuroim2::space(vol)), neuroim2::trans(space), tolerance = 1e-4))) {
    a <- as.array(vol)
    if (length(dim(a)) > 3L) a <- a[, , , 1]
    dim(a) <- dim(space)[1:3]
    return(a)
  }
  d <- dim(space)[1:3]
  g <- as.matrix(expand.grid(seq_len(d[1]), seq_len(d[2]), seq_len(d[3])))
  world <- neuroim2::grid_to_coord(space, g)
  sg <- round(neuroim2::coord_to_grid(neuroim2::space(vol), world))
  sd <- dim(vol)[1:3]
  ok <- rowSums(sg >= 1) == 3 & sg[, 1] <= sd[1] & sg[, 2] <= sd[2] & sg[, 3] <= sd[3]
  out <- array(NA_real_, d)
  out[ok] <- as.array(vol)[sg[ok, , drop = FALSE]]
  out
}

# ---- project ------------------------------------------------------------------

.studio_explain_rows <- function(spec) {
  rows <- list()
  walk <- function(x, path, depth) {
    for (k in names(x)) {
      p <- if (nzchar(path)) paste0(path, ".", k) else k
      v <- x[[k]]
      if (.dsl_is_map(v) && length(v)) {
        rows[[length(rows) + 1L]] <<- list(path = p, key = k, depth = depth, value = NULL, origin = NULL)
        walk(v, p, depth + 1L)
      } else {
        o <- unname(spec$origins[p])
        rows[[length(rows) + 1L]] <<- list(path = p, key = k, depth = depth, value = .dsl_inline(v),
                                           origin = if (length(o) && !is.na(o)) o else NULL)
      }
    }
  }
  walk(spec$config[intersect(c(spec$verb, "folds", "scope", "inference", "run"), names(spec$config))], "", 0L)
  rows
}

.studio_project <- function(project) {
  root <- project$root
  data <- utils::modifyList(.dsl_defaults(project$version)$data, project$config$data)
  subjects <- tryCatch(.dsl_subjects(data, root), error = function(e) character())
  chk <- as.data.frame(unclass(rmvpa_check(project, quiet = TRUE)), stringsAsFactors = FALSE)
  analyses <- lapply(names(project$analyses), function(a) {
    file <- readLines(project$analyses[[a]], warn = FALSE)
    spec <- tryCatch(rmvpa_analysis(project, a), error = function(e) e)
    if (inherits(spec, "error")) {
      return(list(name = a, verb = "invalid", scope = "", file = paste(file, collapse = "\n"), explain = list(),
                  notes = conditionMessage(spec), inference = list(kind = "none", description = conditionMessage(spec)),
                  lock = "none", cells = lapply(subjects, function(s) list(subject = s, status = "blocked"))))
    }
    st <- tryCatch(rmvpa_status(project, a), error = function(e) NULL)
    cells <- if (is.null(st)) lapply(subjects, function(s) list(subject = s, status = "blocked"))
             else lapply(seq_len(nrow(st)), function(i) list(subject = st$subject[i], status = st$status[i]))
    inf <- spec$config$inference
    list(name = a, verb = spec$verb, family = spec$family, scope = spec$scope_kind,
         scope_detail = spec$config$scope[[1]], file = paste(file, collapse = "\n"),
         explain = .studio_explain_rows(spec), notes = unlist(spec$notes),
         inference = list(kind = spec$inference_support$kind, description = spec$inference_support$description,
                          test = inf$test, permutations = inf$permutations),
         lock = .dsl_lock_status(spec), cells = cells)
  })
  sets <- lapply(names(data$sets), function(s) {
    st <- data$sets[[s]]
    first <- tryCatch(.dsl_set_files(data, s, subjects[1], root), error = function(e) NULL)
    rows <- if (is.null(first)) NULL else sum(vapply(first$designs, function(f) nrow(.dsl_read_table(f)), 1L))
    list(name = s, entities = st[setdiff(names(st), "keep")], keep = st$keep, subjects = length(subjects),
         runs = if (is.null(first)) NULL else length(first$betas), rows_per_subject = rows)
  })
  list(name = basename(root), root = .studio_tilde(root), version = project$version,
       data = project$config$data[intersect(c("bids", "space", "mask"), names(project$config$data))],
       sets = sets, folds = project$config$folds, run = project$config$run, subjects = subjects,
       project_file = readLines(project$file, warn = FALSE), analyses = analyses, check = chk)
}

# ---- template -----------------------------------------------------------------

# The display grid is the first subject's mask. The grey background is
# `background` (any image; resampled to that grid) or else the mask itself.
.studio_grid <- function(project, background = NULL) {
  data <- utils::modifyList(.dsl_defaults(project$version)$data, project$config$data)
  s1 <- .dsl_subjects(data, project$root)[1]
  mask <- neuroim2::read_vol(.dsl_mask_file(data, s1, project$root))
  sp <- neuroim2::space(mask)
  m <- as.array(mask) > 0
  bg <- if (is.null(background)) {
    ifelse(m, 150, 0)
  } else {
    b <- .studio_on_grid(neuroim2::read_vol(.dsl_resolve_path(background, project$root)), sp)
    b[!is.finite(b)] <- 0
    hi <- stats::quantile(b[b > 0], 0.995, names = FALSE)
    if (!is.finite(hi) || hi <= 0) hi <- max(b, 1)
    pmin(255, round(255 * pmax(b, 0) / hi))
  }
  list(space = sp, mask = m, bg = bg)
}

.studio_template <- function(project, grid) {
  A <- neuroim2::trans(grid$space)
  list(dims = dim(grid$space)[1:3], spacing = neuroim2::spacing(grid$space)[1:3],
       affine = lapply(1:4, function(i) unname(A[i, ])), space = project$config$data$space,
       template = jsonlite::base64_enc(as.raw(as.integer(grid$bg))),
       mask = jsonlite::base64_enc(as.raw(as.integer(grid$mask))))
}

# ---- results ------------------------------------------------------------------

.studio_cells <- function(project, spec) {
  st <- rmvpa_status(project, spec$name)
  out <- .dsl_out_dir(spec)
  lapply(seq_len(nrow(st)), function(i) {
    dir <- file.path(out, paste0("sub-", st$subject[i]))
    side <- list.files(dir, pattern = "_rmvpa\\.json$", full.names = TRUE)
    list(id = st$subject[i], status = st$status[i], dir = dir,
         sidecar = if (length(side)) jsonlite::fromJSON(side[1], simplifyVector = FALSE) else NULL)
  })
}

# Observed-metric labels among a cell's statmaps: every stat-<x> label that is
# not another label's p or pfdr map.
.studio_metric_labels <- function(files) {
  labs <- unique(sub("^.*_stat-([A-Za-z0-9]+)_statmap\\.nii(\\.gz)?$", "\\1", grep("_stat-[A-Za-z0-9]+_statmap\\.nii", files, value = TRUE)))
  labs[!vapply(labs, function(l) any(l %in% c(paste0(labs, "p"), paste0(labs, "pfdr"))), TRUE)]
}

.studio_stat_file <- function(dir, label, suffix = "_statmap.nii") {
  f <- list.files(dir, pattern = paste0("_stat-", label, suffix), full.names = TRUE)
  if (length(f)) f[1] else NULL
}

.studio_null_hist <- function(file, edges) {
  if (is.null(file)) return(NULL)
  q <- .dsl_read_table(file)
  v <- q$value; p <- q$quantile
  keep <- !duplicated(v)
  F <- if (sum(keep) >= 2) stats::approx(v[keep], p[keep], xout = edges, yleft = 0, yright = 1, ties = "ordered")$y
       else as.numeric(edges >= v[1])
  round(diff(F) * q$n[1])
}

.studio_hist <- function(v, edges) {
  v <- v[is.finite(v)]
  as.integer(table(cut(pmin(pmax(v, edges[1]), edges[length(edges)]), edges, include.lowest = TRUE)))
}

.studio_searchlight <- function(spec, cells, grid, k) {
  have <- Filter(function(c) c$status %in% c("done", "stale") && !is.null(c$sidecar), cells)
  if (!length(have)) return(NULL)
  files <- list.files(have[[1]]$dir)
  labels <- .studio_metric_labels(files)
  if (!length(labels)) return(NULL)
  chance_tab <- have[[1]]$sidecar$chance
  metric <- if ("Accuracy" %in% labels) "Accuracy" else labels[1]
  chance <- chance_tab[[metric]] %||% 0
  m <- grid$mask
  read <- function(f) if (is.null(f)) NULL else .studio_on_grid(neuroim2::read_vol(f), grid$space)
  per <- lapply(have, function(c) {
    obs <- read(.studio_stat_file(c$dir, metric)); obs[!m] <- NA
    pf <- read(.studio_stat_file(c$dir, paste0(metric, "pfdr")))
    list(cell = c, obs = obs, sig = if (is.null(pf)) NULL else (!is.na(pf) & pf < 0.05 & m),
         null = .studio_stat_file(c$dir, metric, "_null\\.tsv$"))
  })
  done <- vapply(per, function(p) identical(p$cell$status, "done"), TRUE)
  use <- if (any(done)) per[done] else per
  group <- array(rowMeans(do.call(cbind, lapply(use, function(p) as.vector(p$obs))), na.rm = TRUE), dim(m))
  group[!m | !is.finite(group)] <- NA
  tested <- all(vapply(per, function(p) !is.null(p$sig), TRUE))
  count <- function(ps) { x <- Reduce(`+`, lapply(ps, function(p) p$sig + 0L), array(0L, dim(m))); x[!m] <- NA; x }
  peak <- which.max(group)
  pk <- arrayInd(peak, dim(group))
  vals <- unlist(lapply(per, function(p) p$obs[m]))
  vals <- vals[is.finite(vals)]
  lo <- min(chance, stats::quantile(vals, 0.001, names = FALSE)); hi <- max(stats::quantile(vals, 0.999, names = FALSE), chance)
  pad <- 0.05 * (hi - lo + 1e-9); lo <- lo - pad; hi <- hi + pad
  if (metric == "Accuracy") { lo <- max(lo, 0); hi <- min(hi, 1) }
  edges <- seq(lo, hi, length.out = 101L)
  sig_done <- if (tested) lapply(per[done], `[[`, "sig") else list()
  subjects <- lapply(seq_along(per), function(i) {
    p <- per[[i]]; c <- p$cell
    sl <- p$obs[, , pk[3]]
    loso <- NULL
    if (tested) {
      others <- per[done & seq_along(per) != i]
      if (length(others) >= k) {
        cnt <- Reduce(`+`, lapply(others, function(o) o$sig + 0L))
        roi <- cnt >= k & m
        if (any(roi)) loso <- round(mean(p$obs[roi], na.rm = TRUE), 4)
      }
    }
    surv <- if (tested) p$obs[p$sig] else numeric()
    list(id = c$id, status = c$status, loso_roi = loso, mean = round(mean(p$obs[m], na.rm = TRUE), 4),
         seconds = c$sidecar$seconds, kept = c$sidecar$rows$train$kept,
         slice = .studio_f32(sl),
         fdr = if (tested) list(threshold = if (length(surv)) round(min(surv, na.rm = TRUE), 4) else NULL,
                                surviving = sum(p$sig)) else NULL,
         centres_hist = .studio_hist(p$obs[m], edges),
         null_hist = .studio_null_hist(p$null, edges))
  })
  inf <- have[[1]]$sidecar$inference
  maps <- list(group_mean = .studio_f32(group))
  if (tested) {
    maps$prevalence <- .studio_f32(count(per[done]))
    maps$prevalence_all <- .studio_f32(count(per))
  }
  perms <- if (identical(inf$test, "permutation")) inf$permutations else NULL
  list(analysis = spec$name, metric = metric, metrics = labels, chance = chance, classes = chance_tab$classes,
       null = list(kind = if (tested) spec$inference_support$kind else "none", permutations = perms,
                   p_resolution = if (!is.null(perms)) round(1 / (perms + 1), 4) else NULL,
                   adjust = if (identical(inf$null_method, "adjusted")) "binned by searchlight size (number of voxels), so edge searchlights are compared with equally small ones" else "one pooled null for all centres",
                   description = inf$null %||% spec$inference_support$description),
       domain = c(lo, hi), top = signif(stats::quantile(group[m], 0.995, na.rm = TRUE, names = FALSE), 3),
       peak_z = pk[3] - 1L, peak = list(ijk = as.integer(pk) - 1L, value = round(group[peak], 4)),
       group = list(roi_rule = "centres where at least k of the other subjects survive their own FDR (q < .05)", roi_k = k,
                    note = "no group-level test is computed; these are per-subject values"),
       subjects = subjects, maps = maps, hist_edges = signif(edges, 6), sidecar = have[[1]]$sidecar)
}

.studio_regional <- function(spec, cells, grid, surfaces) {
  have <- Filter(function(c) identical(c$status, "done") && !is.null(c$sidecar), cells)
  if (!length(have)) return(NULL)
  tabs <- lapply(have, function(c) {
    f <- list.files(c$dir, pattern = "_regions\\.tsv$", full.names = TRUE)
    if (length(f)) .dsl_read_table(f[1]) else NULL
  })
  ok <- !vapply(tabs, is.null, TRUE)
  have <- have[ok]; tabs <- tabs[ok]
  if (!length(have)) return(NULL)
  t1 <- tabs[[1]]
  chance_tab <- have[[1]]$sidecar$chance
  num <- names(t1)[vapply(t1, is.numeric, TRUE) & names(t1) != "roinum"]
  metric <- if ("Accuracy" %in% num) "Accuracy" else num[1]
  if (is.na(metric)) return(NULL)
  chance <- chance_tab[[metric]] %||% 0
  ids <- sort(unique(unlist(lapply(tabs, `[[`, "roinum"))))
  rg <- spec$config$scope$regions
  afile <- .dsl_resolve_path(rg$atlas, spec$project$root)
  lt <- if (!is.null(rg$labels)) .dsl_read_table(.dsl_resolve_path(rg$labels, spec$project$root)) else NULL
  col <- function(nm) if (!is.null(lt) && nm %in% names(lt)) lt[[nm]] else NULL
  lname <- col("name") %||% col("label")
  hemi_guess <- function(n) {
    if (grepl("(^|[_ ])(LH|lh|L)([_ ]|$)|^left|_left$|^Left", n)) "left"
    else if (grepl("(^|[_ ])(RH|rh|R)([_ ]|$)|^right|_right$|^Right", n)) "right" else NULL
  }
  regions <- lapply(ids, function(r) {
    li <- if (!is.null(lt)) match(r, lt[[1]]) else NA
    full <- if (!is.na(li) && !is.null(lname)) as.character(lname[li]) else paste("region", r)
    hemi <- if (!is.na(li) && !is.null(col("hemisphere"))) as.character(col("hemisphere")[li]) else hemi_guess(full)
    short <- if (!is.na(li) && !is.null(col("label"))) as.character(col("label")[li])
             else trimws(sub("^(LH|RH|Left|Right|left|right|L|R)[_ ]+", "", sub("^7Networks_[LR]H_", "", full)))
    net <- if (!is.na(li) && !is.null(col("network"))) as.character(col("network")[li])
           else if (grepl("^[0-9]+Networks_[LR]H_", full)) sub("^[0-9]+Networks_[LR]H_([A-Za-z]+)_.*$", "\\1", full) else NULL
    list(index = r, name = short, full = full, hemisphere = hemi, network = net,
         values = vapply(tabs, function(t) { v <- t[[metric]][match(r, t$roinum)]; if (length(v) && is.finite(v)) round(v, 4) else NA_real_ }, 1))
  })
  lab <- .studio_on_grid(neuroim2::read_vol(afile), grid$space)
  lab[!is.finite(lab)] <- 0
  md5 <- unname(tools::md5sum(afile))
  surf <- Filter(function(s) identical(s$atlas_md5, md5), surfaces)
  inf <- have[[1]]$sidecar$inference
  list(analysis = spec$name, metric = metric, chance = chance, classes = chance_tab$classes,
       chance_source = if (!is.null(chance_tab$classes) && identical(metric, "Accuracy")) sprintf("1 / %d classes", chance_tab$classes) else NULL,
       subjects = vapply(have, `[[`, "", "id"), regions = regions,
       atlas = list(name = basename(afile), short = sub("\\.nii(\\.gz)?$", "", basename(afile)), file = rg$atlas, md5 = md5,
                    surface = if (length(surf)) surf[[1]]$template else NULL),
       labels = .studio_i32(round(lab)),
       null = list(kind = if (identical(inf$test, "permutation")) spec$inference_support$kind else "none",
                   description = inf$null %||% spec$inference_support$description),
       sidecar = have[[1]]$sidecar)
}

.studio_results <- function(project, name, grid, k = 3L, surfaces = list()) {
  if (!name %in% names(project$analyses)) stop("no analysis named '", name, "'", call. = FALSE)
  spec <- rmvpa_analysis(project, name)
  cells <- .studio_cells(project, spec)
  out <- switch(spec$scope_kind,
                searchlight = .studio_searchlight(spec, cells, grid, k),
                regions = .studio_regional(spec, cells, grid, surfaces),
                list(analysis = name, note = sprintf("The studio has no view for %s results yet; see %s.",
                                                     spec$scope_kind, .studio_tilde(.dsl_out_dir(spec)))))
  out %||% list(missing = TRUE)
}

# A surface file is JSON as written by dsl/app/fixture/make-regions.R:
# template, left/right {positions, faces, sulc, labels}, and atlas_md5 naming
# the volumetric atlas whose parcel ids the vertex labels use.
.studio_read_surfaces <- function(paths) {
  out <- lapply(paths, function(p) {
    s <- jsonlite::fromJSON(p, simplifyVector = FALSE)
    if (is.null(s$template) || is.null(s$left) || is.null(s$right)) stop("not a studio surface file: ", p, call. = FALSE)
    s
  })
  stats::setNames(out, vapply(out, `[[`, "", "template"))
}

# Everything the app can ask for, keyed as in src/api.js.
.studio_source <- function(project, background = NULL, surfaces = NULL, k = 3L) {
  if (!inherits(project, "rmvpa_project")) project <- rmvpa_project(project)
  .dsl_require_json()
  surf <- .studio_read_surfaces(surfaces %||% character())
  grid <- NULL
  get_grid <- function() { if (is.null(grid)) grid <<- .studio_grid(project, background); grid }
  cache <- new.env(parent = emptyenv())
  stamp <- function(name) {
    spec <- rmvpa_analysis(project, name)
    f <- list.files(.dsl_out_dir(spec), recursive = TRUE, all.files = TRUE, full.names = TRUE)
    paste(c(f, format(file.mtime(f), "%s")), collapse = "|")
  }
  list(
    project = function() .studio_project(rmvpa_project(project$file)),
    template = function() .studio_template(project, get_grid()),
    results = function(name) {
      key <- digest::digest(c(name, stamp(name)))
      if (is.null(cache[[key]])) cache[[key]] <- .studio_results(project, name, get_grid(), k, surf)
      cache[[key]]
    },
    surface = function(template) {
      if (is.null(surf[[template]])) stop("no surface '", template, "' was given to the studio", call. = FALSE)
      surf[[template]]
    },
    analyses = function() names(rmvpa_project(project$file)$analyses)
  )
}

# ---- server -------------------------------------------------------------------

.studio_token <- function() {
  raw <- tryCatch({ con <- file("/dev/urandom", "rb", raw = TRUE); on.exit(close(con)); readBin(con, "raw", 24L) },
                  error = function(e) NULL)
  if (length(raw) != 24L) raw <- as.raw(sample.int(256L, 24L, replace = TRUE) - 1L)
  paste(sprintf("%02x", as.integer(raw)), collapse = "")
}

.studio_mime <- function(path) {
  switch(tolower(tools::file_ext(path)),
         html = "text/html; charset=utf-8", js = "text/javascript; charset=utf-8", css = "text/css; charset=utf-8",
         json = "application/json", svg = "image/svg+xml", png = "image/png", ico = "image/x-icon",
         woff2 = "font/woff2", "application/octet-stream")
}

.studio_headers <- function(type, extra = list()) {
  c(list(`Content-Type` = type, `Cache-Control` = "no-store", `X-Content-Type-Options` = "nosniff",
         `Referrer-Policy` = "no-referrer", `Content-Security-Policy` = "frame-ancestors 'none'",
         `X-Frame-Options` = "DENY"), extra)
}

.studio_reply <- function(status, body, type = "text/plain; charset=utf-8") {
  list(status = status, headers = .studio_headers(type), body = body)
}

.studio_query <- function(qs) {
  qs <- sub("^\\?", "", qs %||% "")
  if (!nzchar(qs)) return(list())
  kv <- strsplit(strsplit(qs, "&", fixed = TRUE)[[1]], "=", fixed = TRUE)
  stats::setNames(lapply(kv, function(x) utils::URLdecode(if (length(x) > 1) x[2] else "")),
                  vapply(kv, function(x) utils::URLdecode(x[1]), ""))
}

.studio_const_eq <- function(a, b) {
  is.character(a) && length(a) == 1L && !is.na(a) && nchar(a) == nchar(b) &&
    sum(utf8ToInt(a) != utf8ToInt(b)) == 0L
}

# The request handler, separate from httpuv so it can be tested directly.
# Every request must name this server in Host (DNS-rebinding guard) and, when
# it carries an Origin, come from this server's own origin. The page needs
# ?token= and the API needs the X-Rmvpa-Token header: another site in the
# browser can neither read the token nor set that header without a CORS
# preflight, which this server never grants.
.studio_handler <- function(source, app_dir, token, port) {
  hosts <- paste0(c("127.0.0.1", "localhost"), ":", port)
  origins <- paste0("http://", hosts)
  index <- paste(readLines(file.path(app_dir, "index.html"), warn = FALSE), collapse = "\n")
  index <- sub("<head>", sprintf('<head>\n    <meta name="rmvpa-token" content="%s" />', token), index, fixed = TRUE)
  app_root <- normalizePath(app_dir)
  function(req) {
    path <- req$PATH_INFO %||% "/"
    if (!(req$HTTP_HOST %||% "") %in% hosts) return(.studio_reply(403L, "Forbidden: unexpected Host header."))
    origin <- req$HTTP_ORIGIN
    if (!is.null(origin) && !origin %in% origins) return(.studio_reply(403L, "Forbidden: cross-origin request."))
    if (!req$REQUEST_METHOD %in% c("GET", "HEAD")) return(.studio_reply(405L, "Method not allowed."))
    q <- .studio_query(req$QUERY_STRING)
    if (path %in% c("/", "/index.html")) {
      if (!.studio_const_eq(q$token, token)) {
        return(.studio_reply(403L, "Open the studio with the address rmvpa_studio() printed (it carries a session token)."))
      }
      return(.studio_reply(200L, index, "text/html; charset=utf-8"))
    }
    if (startsWith(path, "/api/")) {
      if (!.studio_const_eq(req$HTTP_X_RMVPA_TOKEN, token)) return(.studio_reply(403L, "Forbidden: missing session token."))
      res <- tryCatch(switch(sub("^/api/", "", path),
                             project = source$project(),
                             template = source$template(),
                             results = source$results(q$analysis %||% ""),
                             surface = source$surface(q$template %||% ""),
                             structure(list(), class = "studio_404")),
                      error = function(e) structure(list(message = conditionMessage(e)), class = "studio_error"))
      if (inherits(res, "studio_404")) return(.studio_reply(404L, "Not found."))
      if (inherits(res, "studio_error")) {
        return(.studio_reply(500L, as.character(.studio_json(list(error = res$message))), "application/json"))
      }
      return(.studio_reply(200L, as.character(.studio_json(res)), "application/json"))
    }
    f <- normalizePath(file.path(app_root, sub("^/", "", utils::URLdecode(path))), mustWork = FALSE)
    if (!startsWith(f, paste0(app_root, .Platform$file.sep)) || !file.exists(f) || dir.exists(f)) {
      return(.studio_reply(404L, "Not found."))
    }
    list(status = 200L, headers = .studio_headers(.studio_mime(f)),
         body = readBin(f, "raw", file.info(f)$size))
  }
}

#' Open the studio for a project
#'
#' Serves the studio (a browser app bundled with rmvpadsl) on this machine and
#' opens it. The studio shows the project's analyses, what `rmvpa explain` and
#' `rmvpa check` report, the status of every subject x analysis cell, and the
#' results that [rmvpa_run()] wrote. It reads outputs; it never fits a model.
#'
#' The server listens on 127.0.0.1 only. The address it prints carries a
#' session token; the page sends it back with every data request, and requests
#' with another Host or Origin are refused, so other web pages open in the same
#' browser cannot read the project.
#'
#' @param project An `rmvpa_project` or a project path.
#' @param port Port to listen on; a free one by default.
#' @param launch Open the studio in the browser.
#' @param background Optional image used as the grey background of every
#'   slice (an anatomical in the data space, say); resampled to the grid of the
#'   first subject's mask. Without it the mask itself is drawn.
#' @param surfaces Optional studio surface files (JSON: fsaverage meshes with
#'   per-vertex parcel ids, as written by `dsl/app/fixture/make-regions.R`).
#'   Regional results on an atlas whose md5 matches a file's `atlas_md5` are
#'   also drawn on that cortical surface.
#' @param k Subjects (other than the one shown) that must survive their own FDR
#'   for a centre to enter a subject's leave-one-subject-out region.
#' @param block Keep serving until interrupted. The default does so outside
#'   interactive sessions (the `rmvpa studio` command).
#' @return An `rmvpa_studio` object (invisibly) with the `url`; stop it with
#'   [rmvpa_studio_stop()].
#' @export
rmvpa_studio <- function(project = ".", port = NULL, launch = interactive(), background = NULL,
                         surfaces = NULL, k = 3L, block = !interactive()) {
  if (!requireNamespace("httpuv", quietly = TRUE)) {
    stop("The studio server needs the 'httpuv' package: install.packages(\"httpuv\")", call. = FALSE)
  }
  if (!inherits(project, "rmvpa_project")) project <- rmvpa_project(project)
  app <- .studio_app_dir()
  port <- as.integer(port %||% httpuv::randomPort(min = 7700L, max = 7999L))
  token <- .studio_token()
  handler <- .studio_handler(.studio_source(project, background, surfaces, k), app, token, port)
  server <- httpuv::startServer("127.0.0.1", port, list(call = handler))
  url <- sprintf("http://127.0.0.1:%d/?token=%s", port, token)
  message("rmvpa studio for ", .studio_tilde(project$root), " at\n  ", url)
  out <- structure(list(url = url, port = port, server = server, project = project$root), class = "rmvpa_studio")
  if (launch) utils::browseURL(url)
  if (block) {
    message("Press Ctrl-C to stop.")
    on.exit(httpuv::stopServer(server))
    repeat httpuv::service(250)
  }
  invisible(out)
}

#' @rdname rmvpa_studio
#' @param x An `rmvpa_studio` object.
#' @export
rmvpa_studio_stop <- function(x) {
  if (inherits(x, "rmvpa_studio") && !is.null(x$server)) httpuv::stopServer(x$server)
  invisible(NULL)
}

#' @export
print.rmvpa_studio <- function(x, ...) {
  cat("rmvpa studio for", .studio_tilde(x$project), "\n ", x$url, "\n")
  invisible(x)
}

# ---- static export ------------------------------------------------------------

#' Write the studio and a project's results to a folder
#'
#' Copies the studio app and writes every payload it reads (project, display
#' template, results of each analysis that has any, surfaces) next to it as
#' small script files. The folder needs no R and no server: open its
#' `index.html` from disk, or put it on any static web host. It is a snapshot,
#' stamped with the time it was written; re-export after new runs.
#'
#' The export contains result summaries and maps, the analysis files and the
#' project's paths, but not the beta images.
#'
#' @inheritParams rmvpa_studio
#' @param dir Folder to write; created if needed.
#' @param analyses Analyses to include; all by default.
#' @param overwrite Replace `dir` if it already holds files.
#' @return The path of the exported `index.html`, invisibly.
#' @export
rmvpa_export_studio <- function(project = ".", dir, analyses = NULL, background = NULL, surfaces = NULL,
                                k = 3L, overwrite = FALSE) {
  if (!inherits(project, "rmvpa_project")) project <- rmvpa_project(project)
  app <- .studio_app_dir()
  if (dir.exists(dir) && length(list.files(dir, all.files = TRUE, no.. = TRUE))) {
    if (!overwrite) stop(dir, " is not empty; pass overwrite = TRUE to replace it.", call. = FALSE)
    unlink(dir, recursive = TRUE)
  }
  dir.create(file.path(dir, "data"), recursive = TRUE, showWarnings = FALSE)
  file.copy(list.files(app, full.names = TRUE), dir, recursive = TRUE)

  src <- .studio_source(project, background, surfaces, k)
  payload <- list()
  proj <- src$project()
  proj$exported <- format(Sys.time(), "%Y-%m-%d %H:%M %Z")
  payload$project <- proj
  payload$template <- src$template()
  analyses <- analyses %||% src$analyses()
  for (a in analyses) {
    r <- tryCatch(src$results(a), error = function(e) { warning(a, ": ", conditionMessage(e), call. = FALSE); NULL })
    if (!is.null(r) && !isTRUE(r$missing)) payload[[paste0("results/", a)]] <- r
    tpl <- r$atlas$surface
    if (!is.null(tpl) && is.null(payload[[paste0("surface/", tpl)]])) payload[[paste0("surface/", tpl)]] <- src$surface(tpl)
  }
  tags <- character()
  for (key in names(payload)) {
    f <- sprintf("data/%s.js", gsub("[^A-Za-z0-9_-]", "_", key))
    writeLines(sprintf("(window.__RMVPA_DATA__ = window.__RMVPA_DATA__ || {})[%s] = %s;",
                       .studio_json(key), .studio_json(payload[[key]])), file.path(dir, f), useBytes = TRUE)
    tags <- c(tags, sprintf('    <script src="%s"></script>', f))
  }
  index <- file.path(dir, "index.html")
  html <- readLines(index, warn = FALSE)
  at <- grep("<script defer", html, fixed = TRUE)[1]
  if (is.na(at)) stop("unexpected studio index.html: no deferred app script", call. = FALSE)
  html <- append(html, tags, after = at - 1L)
  writeLines(html, index)
  message(sprintf("Wrote the studio for %s (%d analyses with results) to %s", basename(project$root),
                  sum(startsWith(names(payload), "results/")), normalizePath(dir)))
  invisible(normalizePath(index))
}
