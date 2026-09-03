# Predictive leave-one-band-out attribution -------------------------------

.bra_validate_delta_sets <- function(delta_sets, group_names) {
  if (is.null(delta_sets)) return(character())
  if (!is.character(delta_sets) || !length(delta_sets) || anyNA(delta_sets) ||
      any(!nzchar(delta_sets)) || anyDuplicated(delta_sets)) {
    stop("delta_sets must be NULL, 'all', or unique feature-band names.",
         call. = FALSE)
  }
  if ("all" %in% delta_sets) {
    if (length(delta_sets) != 1L) {
      stop("delta_sets = 'all' cannot be combined with named bands.",
           call. = FALSE)
    }
    delta_sets <- group_names
  }
  unknown <- setdiff(delta_sets, group_names)
  if (length(unknown)) {
    stop("delta_sets contains unknown feature bands: ",
         paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  if (length(group_names) < 2L && length(delta_sets)) {
    stop("Predictive leave-one-band-out delta R2 cannot drop the only feature band.",
         call. = FALSE)
  }
  delta_sets
}

.bra_project_theta <- function(theta, keep_groups, context) {
  theta <- theta[, keep_groups, drop = FALSE]
  mass <- rowSums(theta)
  valid <- is.finite(mass) & mass > sqrt(.Machine$double.eps)
  if (!any(valid)) {
    stop(context, " leaves no positive weight on any retained feature band.",
         call. = FALSE)
  }
  theta <- theta[valid, , drop = FALSE]
  theta <- theta / rowSums(theta)
  list(theta = theta, valid = valid)
}

.bra_reduced_candidates <- function(candidates, excluded_band) {
  group_names <- candidates$group_names
  if (!excluded_band %in% group_names) {
    stop("excluded_band must name one candidate feature band.", call. = FALSE)
  }
  keep_groups <- setdiff(group_names, excluded_band)
  if (!length(keep_groups)) {
    stop("Predictive leave-one-band-out delta R2 cannot drop every feature band.",
         call. = FALSE)
  }
  projected <- .bra_project_theta(
    candidates$theta, keep_groups,
    paste0("Dropping band '", excluded_band, "'")
  )
  alpha <- candidates$alpha[projected$valid]
  ids <- .brt_candidate_ids(alpha, projected$theta, keep_groups)
  unique_rows <- !duplicated(ids)
  out <- .brt_canonical_candidates(
    alpha[unique_rows], projected$theta[unique_rows, , drop = FALSE],
    keep_groups
  )
  attr(out, "generation") <- list(
    method = "project-and-retune",
    excluded_band = excluded_band,
    source_candidate_count = length(candidates$alpha),
    invalid_zero_mass_count = sum(!projected$valid),
    duplicate_projection_count = sum(!unique_rows)
  )
  out
}

.bra_reduced_fixed_theta <- function(fixed_theta, group_names, excluded_band) {
  if (is.null(fixed_theta)) return(NULL)
  theta <- .brt_theta_matrix(fixed_theta, group_names, "fixed_theta")
  projected <- .bra_project_theta(
    theta, setdiff(group_names, excluded_band),
    paste0("fixed_theta after dropping band '", excluded_band, "'")
  )
  drop(projected$theta[1L, ])
}

.bra_cache_counters <- function(cache) {
  if (is.null(cache)) {
    return(c(decompositions = 0L, hits = 0L, kernel_builds = 0L,
             kernel_hits = 0L, evictions = 0L, oversize = 0L))
  }
  c(decompositions = cache$decomposition_count, hits = cache$hit_count,
    kernel_builds = cache$kernel_build_count,
    kernel_hits = cache$kernel_hit_count,
    evictions = cache$eviction_count,
    oversize = cache$oversize_count)
}

.bra_work_row <- function(model_id,
                          excluded_band,
                          chunk,
                          n_features,
                          n_responses,
                          n_candidates,
                          nested,
                          before,
                          after) {
  inner_candidate_fit_calls <- sum(vapply(nested$outer_results, function(x) {
    dim(x$inner_fold_scores)[[1L]] * dim(x$inner_fold_scores)[[2L]]
  }, integer(1)))
  outer_refit_calls <- sum(vapply(
    nested$outer_results, function(x) length(x$outer_fits), integer(1)
  ))
  data.frame(
    model_id = model_id,
    excluded_band = if (is.null(excluded_band)) NA_character_ else excluded_band,
    chunk = as.integer(chunk),
    n_features = as.integer(n_features),
    n_responses = as.integer(n_responses),
    n_candidates = as.integer(n_candidates),
    inner_candidate_fit_calls = as.integer(inner_candidate_fit_calls),
    outer_refit_calls = as.integer(outer_refit_calls),
    total_fit_calls = as.integer(inner_candidate_fit_calls + outer_refit_calls),
    decompositions_added = as.integer(after[["decompositions"]] -
                                        before[["decompositions"]]),
    cache_hits_added = as.integer(after[["hits"]] - before[["hits"]]),
    band_kernel_builds_added = as.integer(after[["kernel_builds"]] -
                                            before[["kernel_builds"]]),
    band_kernel_hits_added = as.integer(after[["kernel_hits"]] -
                                          before[["kernel_hits"]]),
    cache_evictions_added = as.integer(after[["evictions"]] -
                                         before[["evictions"]]),
    cache_oversize_added = as.integer(after[["oversize"]] -
                                        before[["oversize"]]),
    stringsAsFactors = FALSE
  )
}
