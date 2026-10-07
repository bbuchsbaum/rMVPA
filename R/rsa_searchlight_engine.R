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
  ids <- vapply(sl, function(w) as.integer(w@parent_index), 1L)
  sel <- if (is.null(centers)) seq_along(sl) else which(ids %in% centers)
  cols <- lapply(sel, function(i) {
    cc <- match(as.integer(neuroim2::grid_to_index(sp, sl[[i]]@coords)), mask_indices)
    cc <- cc[!is.na(cc)]
    cc[keep_voxel[cc] | mask_indices[cc] == ids[i]]
  })
  list(ds = ds, x_all = x_all, mask_indices = mask_indices, centers = ids[sel], cols = cols)
}

#' Run train_model.rsa_model() on every prepared sphere
#'
#' `spec` carries the design, so permuted designs (item_perm) are honoured.
#' Returns a centre x output matrix (row names are centre ids; NA rows for
#' spheres that fail), as the per-ROI path would produce.
#' @keywords internal
#' @noRd
.rsa_engine_score <- function(prep, spec) {
  outputs <- lapply(seq_along(prep$cols), function(i) {
    cols <- prep$cols[[i]]
    if (length(cols) < 1L) return(NULL)
    res <- tryCatch(
      train_model(spec, prep$x_all[, cols, drop = FALSE], y = NULL,
                  indices = prep$mask_indices[cols]),
      error = function(e) NULL
    )
    if (!is.numeric(res) || length(res) == 0L) return(NULL)
    attributes(res)[["fingerprint"]] <- NULL
    res
  })
  good <- which(!vapply(outputs, is.null, logical(1)))
  if (length(good) == 0L) return(NULL)
  metric_names <- names(outputs[[good[1]]])
  perf <- matrix(NA_real_, length(outputs), length(metric_names),
                 dimnames = list(as.character(prep$centers), metric_names))
  for (g in good) perf[g, ] <- outputs[[g]][metric_names]
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
  good <- rowSums(!is.na(perf)) > 0L
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
