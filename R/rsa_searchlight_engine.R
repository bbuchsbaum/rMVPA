# Exact searchlight engine for rsa_model.
#
# An RSA searchlight's per-sphere work is the RDM (n^2 p) and a small
# regression against the precomputed model design. Summing over spheres does
# not reduce that work, and summed RDMs would differ from the per-sphere ones
# in the last bit (enough to reorder near-tied ranks under Spearman). The
# general path's cost is the per-sphere iteration around it: ROI objects,
# filtering, result tibbles and merging. This engine extracts the data once
# and, for each sphere, calls train_model.rsa_model() on the same columns the
# general path uses (the filter_roi() rule: no missing values and non-zero
# range, plus the centre voxel; in get_searchlight() order). Results are
# therefore bit-identical, and every RSA variant train_model.rsa_model()
# supports (regression types, include masks, between-set pairs, item
# permutations) is covered without reimplementation.

#' @keywords internal
#' @noRd
.is_rsa_fast_path <- function(model_spec, method) {
  inherits(model_spec, "rsa_model") &&
    identical(method, "standard") &&
    inherits(model_spec$dataset, "mvpa_image_dataset") &&
    !inherits(model_spec$dataset, "mvpa_multibasis_image_dataset") &&
    !isTRUE(model_spec$return_fingerprint)
}

#' Label-independent setup of the RSA engine: data, the filter_roi() voxel
#' rule, and each sphere's columns (optionally only for `centers`).
#' @keywords internal
#' @noRd
.rsa_engine_prepare <- function(model_spec, radius, centers = NULL) {
  ds <- model_spec$dataset
  mask_indices <- ds$mask_indices %||% compute_mask_indices(ds$mask)
  sl <- get_searchlight(ds, "standard", radius)
  sp <- neuroim2::space(ds$mask)
  x_all <- as.matrix(neuroim2::series(ds$train_data, mask_indices))
  if (!all(is.finite(x_all))) {
    stop(.aggregate_ineligible("data contain missing or non-finite values"))
  }
  # filter_roi(): keep voxels whose range across observations is non-zero.
  span <- matrixStats::colMaxs(x_all) - matrixStats::colMins(x_all)
  keep_voxel <- is.finite(span) & span > 0
  # The searchlight list builds each sphere on access; build them once.
  rois <- as.list(sl)
  ids <- vapply(rois, function(w) as.integer(w@parent_index), 1L)
  sel <- if (is.null(centers)) seq_along(rois) else which(ids %in% centers)
  # Grid index -> mask column, built once (a per-sphere match() against the
  # mask indices rehashes them for every centre). Out-of-mask voxels map to 0.
  dims <- dim(sp)[1:3]
  col_of <- integer(prod(dims))
  col_of[mask_indices] <- seq_along(mask_indices)
  cols <- lapply(sel, function(i) {
    co <- rois[[i]]@coords
    cc <- col_of[co[, 1] + (co[, 2] - 1L) * dims[1] + (co[, 3] - 1L) * (dims[1] * dims[2])]
    cc <- cc[cc > 0L]
    cc[keep_voxel[cc] | mask_indices[cc] == ids[i]]
  })
  list(ds = ds, x_all = x_all, mask_indices = mask_indices, centers = ids[sel], cols = cols,
       cache = new.env(parent = emptyenv()))
}

#' Direct RDM plan for the common rsa_model variants
#'
#' train_model.rsa_model() spends most of a small sphere's time on argument
#' handling, copies and sweep(). For within-set designs without pattern
#' centring, a correlation distance and the cached correlation or plain-lm
#' kernel, the same arithmetic can be done directly: the row order
#' (row_idx_a, then item_perm), the lower-triangle positions (after
#' `include`) and the regression function are resolved once per call.
#' Returns NULL when the spec needs the general path, or when train_model
#' would stop (it then still does, per sphere).
#' @keywords internal
#' @noRd
.rsa_lean_plan <- function(spec, n_rows) {
  design <- spec$design
  if (!identical(design$pair_kind %||% "within", "within")) return(NULL)
  if (!identical(spec$pattern_center %||% "none", "none")) return(NULL)
  if (!(spec$distmethod %in% c("pearson", "spearman"))) return(NULL)
  kernel <- spec$.fast_kernel
  regfun <- switch(spec$regtype,
    pearson = ,
    spearman = if (!is.null(kernel$cor)) run_cor_fast,
    lm = if (!is.null(kernel$lm) && length(spec$nneg) == 0L && !isTRUE(spec$semipartial)) run_lm_fast,
    NULL
  )
  if (is.null(regfun)) return(NULL)
  base_rows <- design$row_idx_a %||% seq_len(n_rows)
  if (length(base_rows) < 2L || min(base_rows) < 1L || max(base_rows) > n_rows) return(NULL)
  perm <- design$item_perm
  if (!is.null(perm)) {
    perm <- tryCatch(.rsa_apply_item_perm(seq_along(base_rows), perm), error = function(e) NULL)
    if (is.null(perm)) return(NULL)
  }
  K <- length(base_rows)
  pairs <- .rdm_pair_indices(K)
  list(
    base_rows = base_rows,
    perm = perm,
    K = K,
    pairs = pairs,
    lin = (pairs$j - 1L) * K + pairs$i,
    include = design$include,
    spearman = identical(spec$distmethod, "spearman"),
    regfun = regfun
  )
}

#' One sphere's correlation-distance RDM (lower triangle, all pairs)
#'
#' Same operations as .rdm_vector_correlation(), so identical values. NULL
#' for spheres the general path must handle: a single voxel (the Spearman
#' path fails there) or a constant pattern (it warns there).
#' @keywords internal
#' @noRd
.rsa_lean_rdm <- function(X, lin, spearman) {
  if (ncol(X) <= 1L) return(NULL)
  if (spearman) X <- matrixStats::rowRanks(X, ties.method = "average")
  X <- X - rowMeans(X)
  norms <- sqrt(rowSums(X^2))
  if (any(norms == 0)) return(NULL)
  X <- X / norms
  1 - tcrossprod(X)[lin]
}

#' Each permuted pair's position in the unpermuted lower triangle
#'
#' Permuted row a is unpermuted row perm[a], so pair (i, j) reads base pair
#' (perm[i], perm[j]); the correlation matrix is exactly symmetric.
#' @keywords internal
#' @noRd
.rsa_perm_pair_pos <- function(plan) {
  a <- plan$perm[plan$pairs$i]
  b <- plan$perm[plan$pairs$j]
  hi <- pmax(a, b)
  lo <- pmin(a, b)
  ((lo - 1L) * (2L * plan$K - lo)) %/% 2L + (hi - lo)
}

#' Unpermuted RDMs for every sphere, cached on the prep for permutations
#'
#' A column of NA marks a sphere the general path handles. Built only when it
#' fits the `rMVPA.rsa_perm_cache_bytes` budget (default 512 MiB); NULL
#' otherwise.
#' @keywords internal
#' @noRd
.rsa_perm_rdm_cache <- function(prep, plan, distmethod) {
  key <- list(rows = plan$base_rows, distmethod = distmethod)
  if (identical(prep$cache$rdm_key, key)) return(prep$cache$rdm)
  if (!.rsa_perm_cache_fits(length(plan$lin), length(prep$cols))) return(NULL)
  n_pairs <- length(plan$lin)
  x <- prep$x_all[plan$base_rows, , drop = FALSE]
  rdm <- matrix(NA_real_, n_pairs, length(prep$cols))
  for (s in seq_along(prep$cols)) {
    cols <- prep$cols[[s]]
    if (length(cols) < 1L) next
    d <- .rsa_lean_rdm(x[, cols, drop = FALSE], plan$lin, plan$spearman)
    if (!is.null(d)) rdm[, s] <- d
  }
  prep$cache$rdm_key <- key
  prep$cache$rdm <- rdm
  rdm
}

#' @keywords internal
#' @noRd
.rsa_perm_cache_fits <- function(n_pairs, n_spheres) {
  budget <- getOption("rMVPA.rsa_perm_cache_bytes", 512 * 1024^2)
  is.numeric(budget) && length(budget) == 1L && !is.na(budget) &&
    as.double(n_pairs) * n_spheres * 8 <= budget
}

#' Score every prepared sphere; a list of named numeric vectors (or NULL)
#' @keywords internal
#' @noRd
.rsa_engine_outputs <- function(prep, spec, use_cache = TRUE) {
  general <- function(cols) {
    tryCatch(
      train_model(spec, prep$x_all[, cols, drop = FALSE], y = NULL,
                  indices = prep$mask_indices[cols]),
      error = function(e) NULL
    )
  }
  plan <- if (.rsa_fast_kernel_enabled()) .rsa_lean_plan(spec, nrow(prep$x_all))
  rdm <- NULL
  if (!is.null(plan) && !is.null(plan$perm) && isTRUE(use_cache)) {
    rdm <- .rsa_perm_rdm_cache(prep, plan, spec$distmethod)
  }
  if (!is.null(rdm)) {
    pos <- .rsa_perm_pair_pos(plan)
    if (!is.null(plan$include)) pos <- pos[plan$include]
  } else if (!is.null(plan)) {
    rows <- plan$base_rows
    if (!is.null(plan$perm)) rows <- rows[plan$perm]
    x <- if (identical(rows, seq_len(nrow(prep$x_all)))) prep$x_all else prep$x_all[rows, , drop = FALSE]
    lin <- plan$lin
  }
  lapply(seq_along(prep$cols), function(i) {
    cols <- prep$cols[[i]]
    if (length(cols) < 1L) return(NULL)
    res <- NULL
    if (!is.null(rdm)) {
      if (!is.na(rdm[1L, i])) res <- tryCatch(plan$regfun(rdm[pos, i], spec), error = function(e) NULL)
      else res <- general(cols)
    } else if (!is.null(plan)) {
      d <- .rsa_lean_rdm(x[, cols, drop = FALSE], lin, plan$spearman)
      if (is.null(d)) {
        res <- general(cols)
      } else {
        if (!is.null(plan$include)) d <- d[plan$include]
        res <- tryCatch(plan$regfun(d, spec), error = function(e) NULL)
      }
    } else {
      res <- general(cols)
    }
    if (!is.numeric(res) || length(res) == 0L) return(NULL)
    attributes(res)[["fingerprint"]] <- NULL
    res
  })
}

#' Split a prep into spatially contiguous chunks carrying only their columns
#' @keywords internal
#' @noRd
.rsa_engine_chunks <- function(prep, chunk_size) {
  idx <- split(seq_along(prep$cols), ceiling(seq_along(prep$cols) / chunk_size))
  lapply(idx, function(g) {
    used <- sort(unique(unlist(prep$cols[g], use.names = FALSE)))
    local_col <- integer(ncol(prep$x_all))
    local_col[used] <- seq_along(used)
    list(x_all = prep$x_all[, used, drop = FALSE],
         mask_indices = prep$mask_indices[used],
         cols = lapply(prep$cols[g], function(cc) local_col[cc]))
  })
}

#' Run the RSA kernel on every prepared sphere
#'
#' `spec` carries the design, so permuted designs (item_perm) are honoured.
#' Returns a centre x output matrix (row names are centre ids; NA rows for
#' spheres that fail), as the per-ROI path would produce; attribute `scored`
#' indexes the rows that returned a result. With several
#' future workers (and no permutation RDM cache to reuse), contiguous chunks
#' of spheres run in parallel, each shipping only its own data columns and a
#' dataset-free spec; set `options(rMVPA.rsa_fast_parallel = FALSE)` to
#' disable.
#' @keywords internal
#' @noRd
.rsa_engine_score <- function(prep, spec) {
  n <- length(prep$cols)
  nworkers <- future::nbrOfWorkers()
  # Permutations reuse the cached RDMs serially when they fit; reindexing is
  # cheaper than shipping them to workers.
  if (!is.null(spec$design$item_perm)) {
    K <- length(spec$design$row_idx_a %||% seq_len(nrow(prep$x_all)))
    cached <- .rsa_perm_cache_fits(K * (K - 1) / 2, n)
  } else {
    cached <- FALSE
  }
  parallel <- isTRUE(getOption("rMVPA.rsa_fast_parallel", TRUE)) &&
    nworkers > 1L && n >= 256L && !cached
  if (parallel) {
    chunks <- .rsa_engine_chunks(prep, .searchlight_chunk_size(n, nworkers, min_chunk = 64L))
    outputs <- unlist(
      future.apply::future_lapply(chunks, .rsa_engine_outputs, spec = as_worker_spec(spec),
                                  use_cache = FALSE, future.seed = FALSE),
      recursive = FALSE, use.names = FALSE
    )
  } else {
    outputs <- .rsa_engine_outputs(prep, spec)
  }
  good <- which(!vapply(outputs, is.null, logical(1)))
  if (length(good) == 0L) return(NULL)
  metric_names <- names(outputs[[good[1]]])
  perf <- matrix(NA_real_, length(outputs), length(metric_names),
                 dimnames = list(as.character(prep$centers), metric_names))
  for (g in good) perf[g, ] <- outputs[[g]][metric_names]
  # Spheres that returned a result, even an all-NA one (e.g. a constant
  # pattern row); the general path writes those centres as NA, not 0.
  attr(perf, "scored") <- good
  perf
}

#' @keywords internal
#' @noRd
run_searchlight_rsa_fast <- function(model_spec, radius, verbose = FALSE, ...) {
  ds <- model_spec$dataset
  prep <- .rsa_engine_prepare(model_spec, radius)
  if (length(prep$centers) == 0L) return(empty_searchlight_result(ds))
  perf <- .rsa_engine_score(prep, model_spec)
  if (is.null(perf)) return(empty_searchlight_result(ds))
  good <- attr(perf, "scored")
  out <- wrap_out(perf[good, , drop = FALSE], ds, ids = prep$centers[good])
  attr(out, "bad_results") <- tibble::tibble()
  out
}

#' @keywords internal
#' @noRd
.run_searchlight_engine.rsa_model <- function(model_spec, radius, method,
                                              engine = "auto",
                                              niter = 4L,
                                              combiner = "average",
                                              drop_probs = FALSE,
                                              fail_fast = FALSE,
                                              backend = c("default", "shard", "auto"),
                                              incremental = TRUE,
                                              gamma = NULL,
                                              verbose = FALSE,
                                              ...) {
  requested <- .match_searchlight_engine(engine)
  backend <- match.arg(backend)
  if (identical(requested, "legacy") || !(requested %in% c("auto", "rsa_fast"))) {
    return(list(handled = FALSE, result = NULL, engine = "legacy"))
  }
  eligible <- .is_rsa_fast_path(model_spec, method) &&
    backend %in% c("default", "auto") && !isTRUE(fail_fast)
  if (!isTRUE(eligible)) {
    if (identical(requested, "rsa_fast")) {
      stop("Requested searchlight engine 'rsa_fast' is not eligible for this analysis.",
           call. = FALSE)
    }
    return(list(handled = FALSE, result = NULL, engine = "legacy"))
  }
  res <- tryCatch(
    run_searchlight_rsa_fast(model_spec, radius = radius, verbose = verbose),
    rmvpa_engine_ineligible = function(e) {
      if (identical(requested, "rsa_fast")) stop(conditionMessage(e), call. = FALSE)
      futile.logger::flog.info("searchlight engine 'rsa_fast' not applicable (%s)",
                               conditionMessage(e))
      NULL
    }
  )
  if (is.null(res)) return(list(handled = FALSE, result = NULL, engine = "legacy"))
  list(handled = TRUE, result = res, engine = "rsa_fast")
}
