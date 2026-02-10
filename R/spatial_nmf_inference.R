#' Component-level Inference for Spatial NMF
#'
#' Tests whether subjects express NMF components differently between groups,
#' using the component loadings \eqn{W} from a spatial NMF fit. This is a
#' component-wise analysis (not voxelwise): each component gets a statistic and
#' a permutation p-value.
#'
#' For two-group inference, the function fits a linear model for each component
#' (group plus optional covariates) and uses label permutations to build a null
#' distribution of the component-level t-statistics. For one-group inference,
#' the user supplies a null distribution of \eqn{W} (e.g., from permuted fits)
#' and the observed component means are compared to that null.
#'
#' Use this when you want to identify which components show reliable group
#' differences in expression, rather than testing individual voxels.
#'
#' @param fit A spatial_nmf_fit object containing a W matrix.
#' @param W Optional n x k matrix of loadings (overrides fit$W).
#' @param groups Factor or vector of group labels (length n); required for two-group tests.
#' @param covariates Optional data frame of covariates (n rows).
#' @param test One of "two_group" or "one_group".
#' @param nperm Number of permutations (ignored for one_group if null_W provided).
#' @param correction One of "maxT" or "none".
#' @param alternative Alternative hypothesis: "two.sided" or "greater". For two_group
#'   tests, "greater" evaluates positive differences for the second factor level
#'   relative to the first.
#' @param null_W Null distribution of W for one-group inference (list or 3D array).
#' @param alpha Significance threshold for counting significant components.
#' @param seed Optional RNG seed for permutations.
#' @param return_perm Logical; return permutation statistics.
#' @param parallel Logical; use future_lapply for permutations (requires future.apply).
#' @param future_seed Optional seed control for future.apply (passed to future_lapply).
#' @param progress Logical; report progress via progressr (works with parallel futures).
#'
#' @return A list with a component-level results table and summary stats.
#'
#' @details
#' \itemize{
#'   \item \code{stat} is the component-level test statistic (t-statistic for
#'   two-group tests; mean loading for one-group tests).
#'   \item \code{p_unc} is the uncorrected permutation p-value.
#'   \item \code{p_fwer} applies maxT family-wise error correction across components.
#'   \item \code{p_global} tests whether any component is significant.
#' }
#' @examples
#' \dontrun{
#'   # Requires completed spatial NMF fit
#'   W <- matrix(rnorm(20*5), 20, 5)
#'   groups <- factor(rep(c("A","B"), each=10))
#'   result <- spatial_nmf_component_test(W=W, groups=groups, nperm=100)
#' }
#' @export
spatial_nmf_component_test <- function(fit = NULL,
                                       W = NULL,
                                       groups = NULL,
                                       covariates = NULL,
                                       test = c("two_group", "one_group"),
                                       nperm = 1000,
                                       correction = c("maxT", "none"),
                                       alternative = c("two.sided", "greater"),
                                       null_W = NULL,
                                       alpha = 0.05,
                                       seed = NULL,
                                       return_perm = FALSE,
                                       parallel = FALSE,
                                       future_seed = TRUE,
                                       progress = FALSE) {
  test <- match.arg(test)
  correction <- match.arg(correction)
  alternative <- match.arg(alternative)

  if (!is.null(fit)) {
    W <- fit$W
  }
  if (is.null(W)) {
    stop("Provide either fit or W.")
  }
  W <- as.matrix(W)
  n <- nrow(W)
  k <- ncol(W)

  if (!is.null(seed)) set.seed(seed)

  if (test == "two_group") {
    if (is.null(groups)) stop("groups is required for two_group inference.")
    groups <- factor(groups)
    if (length(groups) != n) stop("groups length must match nrow(W).")
    if (nlevels(groups) != 2) stop("two_group inference requires exactly 2 groups.")
    if (!is.numeric(nperm) || length(nperm) != 1L || nperm < 1) {
      stop("nperm must be a positive integer.")
    }
    nperm <- as.integer(nperm)

    futile.logger::flog.info(
      "Spatial NMF component test: %d permutations (%s)",
      nperm,
      if (isTRUE(parallel)) "parallel" else "serial"
    )

    design <- .component_design(groups, covariates)
    stat_obs <- .component_t_stats(W, design, coef_index = 2L)
    perm_stats <- .permute_component_stats(
      nperm = nperm,
      groups = groups,
      covariates = covariates,
      W = W,
      coef_index = 2L,
      parallel = parallel,
      future_seed = future_seed,
      progress = progress
    )

    futile.logger::flog.info("Spatial NMF component test: permutations complete.")

    two_sided <- alternative == "two.sided"
    p_unc <- .perm_pvals(stat_obs, perm_stats, two_sided = two_sided)

    if (correction == "maxT") {
      if (two_sided) {
        max_stat <- apply(abs(perm_stats), 1, max)
        p_fwer <- (colSums(outer(max_stat, abs(stat_obs), ">=")) + 1) / (nperm + 1)
      } else {
        max_stat <- apply(perm_stats, 1, max)
        p_fwer <- (colSums(outer(max_stat, stat_obs, ">=")) + 1) / (nperm + 1)
      }
    } else {
      p_fwer <- p_unc
    }

    if (two_sided) {
      p_global <- (sum(apply(abs(perm_stats), 1, max) >= max(abs(stat_obs))) + 1) / (nperm + 1)
    } else {
      p_global <- (sum(apply(perm_stats, 1, max) >= max(stat_obs)) + 1) / (nperm + 1)
    }

    group_means <- vapply(levels(groups), function(g) {
      colMeans(W[groups == g, , drop = FALSE])
    }, numeric(k))
    group_means <- matrix(group_means, nrow = k)

    res_table <- data.frame(
      component = seq_len(k),
      stat = stat_obs,
      mean_group1 = group_means[, 1],
      mean_group2 = group_means[, 2],
      p_unc = p_unc,
      p_fwer = p_fwer,
      stringsAsFactors = FALSE
    )
  } else {
    if (is.null(null_W)) {
      stop("null_W is required for one_group inference.")
    }
    null_list <- .as_null_W_list(null_W, n, k)
    nperm <- length(null_list)

    stat_obs <- colMeans(W)
    perm_stats <- vapply(null_list, function(Wp) colMeans(Wp), numeric(k))
    perm_stats <- t(perm_stats)

    if (alternative == "two.sided") {
      p_unc <- .perm_pvals(stat_obs, perm_stats, two_sided = TRUE)
      obs_max <- max(abs(stat_obs))
      max_stat <- apply(abs(perm_stats), 1, max)
    } else {
      p_unc <- (colSums(perm_stats >= matrix(stat_obs, nrow = nperm, ncol = k, byrow = TRUE)) + 1) / (nperm + 1)
      obs_max <- max(stat_obs)
      max_stat <- apply(perm_stats, 1, max)
    }

    if (correction == "maxT") {
      p_fwer <- (colSums(outer(max_stat, if (alternative == "two.sided") abs(stat_obs) else stat_obs, ">=")) + 1) /
        (nperm + 1)
    } else {
      p_fwer <- p_unc
    }

    p_global <- (sum(max_stat >= obs_max) + 1) / (nperm + 1)

    res_table <- data.frame(
      component = seq_len(k),
      stat = stat_obs,
      p_unc = p_unc,
      p_fwer = p_fwer,
      stringsAsFactors = FALSE
    )
  }

  list(
    table = res_table,
    alpha = alpha,
    n_sig = sum(res_table$p_fwer < alpha, na.rm = TRUE),
    p_global = p_global,
    test = test,
    correction = correction,
    alternative = alternative,
    nperm = nperm,
    perm_stats = if (isTRUE(return_perm)) perm_stats else NULL
  )
}

#' Global Cross-validated Group Test for Spatial NMF
#'
#' Tests whether groups are distinguishable in component space, without
#' attributing the effect to any single component. This is a multivariate,
#' cross-validated test: for each fold, NMF is fit on the training data and
#' component loadings for held-out subjects are used to train/test a classifier.
#' Permuting group labels yields a null distribution for the chosen metric
#' (AUC or accuracy).
#'
#' Use this when you want a single omnibus test of group separability in the
#' learned component space.
#'
#' @param x Optional spatial_nmf_maps_result with return_data=TRUE.
#' @param X Optional subject-by-voxel matrix (overrides x$data).
#' @param groups Factor or vector of group labels (length n).
#' @param k Number of components (defaults to x$fit$k if available).
#' @param lambda Spatial regularization strength (defaults to x$fit$lambda if available).
#' @param graph Optional graph Laplacian list (required if lambda > 0). If `graph$A`
#'   is weighted, set `graph$weighted=TRUE` to preserve weights; otherwise edges are binarized.
#' @param neighbors Neighborhood size for volumetric adjacency (6/18/26).
#' @param nfolds Number of cross-validation folds.
#' @param folds Optional fold specification (vector of fold IDs or list of test indices).
#'   If NULL, folds are stratified by `groups`.
#' @param nperm Number of label permutations.
#' @param permute Permutation strategy: "labels" permutes labels with fixed folds
#'   and fixed NMF projections; "full" re-derives folds based on permuted labels
#'   and refits NMF within each fold. If `folds` is supplied, "full" falls back to
#'   "labels" (fixed folds).
#' @param metric Performance metric: "auc" or "accuracy".
#' @param classifier Classifier: "glm", "lda", or "centroid".
#' @param scale Logical; z-score W within each fold.
#' @param positive Optional positive class label (defaults to second factor level).
#' @param seed Optional RNG seed.
#' @param project_args List of arguments passed to spatial_nmf_project.
#' @param return_perm Logical; return permutation statistics.
#' @param return_cv Logical; return cross-validated predictions and fold IDs.
#' @param parallel Logical; use future_lapply for permutations (requires future.apply).
#' @param future_seed Optional seed control for future.apply (passed to future_lapply).
#' @param progress Logical; report progress via progressr (works with parallel futures).
#' @param ... Additional arguments passed to spatial_nmf_fit.
#'
#' @return A list with the observed statistic, permutation p-value, and metadata.
#'
#' @details
#' \itemize{
#'   \item \code{stat} is the cross-validated performance (AUC or accuracy).
#'   \item \code{p_value} is the permutation p-value under the chosen strategy.
#'   \item \code{permute = "labels"} keeps folds fixed; \code{"full"} re-derives
#'   folds and refits NMF for each permutation (more expensive).
#' }
#' @examples
#' \dontrun{
#'   result <- spatial_nmf_global_test(
#'     matrix(rnorm(100*10), 100, 10),
#'     k = 3, nperm = 99
#'   )
#' }
#' @export
spatial_nmf_global_test <- function(x = NULL,
                                    X = NULL,
                                    groups = NULL,
                                    k = NULL,
                                    lambda = NULL,
                                    graph = NULL,
                                    neighbors = 6,
                                    nfolds = 5,
                                    folds = NULL,
                                    nperm = 1000,
                                    permute = c("labels", "full"),
                                    metric = c("auc", "accuracy"),
                                    classifier = c("glm", "lda", "centroid"),
                                    scale = TRUE,
                                    positive = NULL,
                                    seed = NULL,
                                    project_args = list(max_iter = 100, check_every = 5),
                                    return_perm = FALSE,
                                    return_cv = FALSE,
                                    parallel = FALSE,
                                    future_seed = TRUE,
                                    progress = FALSE,
                                    ...) {
  metric <- match.arg(metric)
  classifier <- match.arg(classifier)
  permute <- match.arg(permute)

  if (!is.null(seed)) set.seed(seed)

  if (!is.null(x)) {
    if (!inherits(x, "spatial_nmf_maps_result")) {
      stop("x must be a spatial_nmf_maps_result.")
    }
    if (is.null(X) && !is.null(x$data)) X <- x$data
    if (is.null(groups) && !is.null(x$groups)) groups <- x$groups
    if (is.null(k) && !is.null(x$fit)) k <- x$fit$k
    if (missing(lambda) && !is.null(x$fit)) lambda <- x$fit$lambda
  }

  if (is.null(X)) {
    stop("X is required (pass spatial_nmf_maps with return_data=TRUE or supply X).")
  }
  X <- as.matrix(X)
  if (any(!is.finite(X))) stop("X contains non-finite values.")
  if (any(X < 0)) stop("X must be non-negative.")

  if (is.null(groups)) stop("groups is required.")
  label_info <- .prepare_binary_groups(groups, positive)
  groups <- label_info$groups
  positive <- label_info$positive
  n <- nrow(X)
  if (length(groups) != n) stop("groups length must match nrow(X).")

  if (is.null(k)) stop("k is required.")
  if (!is.numeric(k) || length(k) != 1L || k < 1) stop("k must be a positive integer.")
  k <- as.integer(k)
  if (is.null(lambda)) lambda <- 0
  if (!is.numeric(nperm) || length(nperm) != 1L || nperm < 1) {
    stop("nperm must be a positive integer.")
  }
  nperm <- as.integer(nperm)

  if (!is.list(project_args)) stop("project_args must be a list.")

  if (lambda > 0 && is.null(graph)) {
    if (!is.null(x) && inherits(x, "spatial_nmf_maps_result")) {
      if (identical(x$map_type, "volume") && !is.null(x$mask) && !is.null(x$dims)) {
        if (length(x$dims) == 2L && neighbors == 6) {
          neighbors <- 4
        }
        A <- build_voxel_adjacency(x$mask, dims = x$dims, neighbors = neighbors)
        graph <- build_graph_laplacian(A)
      } else {
        stop("graph is required for lambda > 0 when using surface inputs.")
      }
    } else {
      stop("graph is required for lambda > 0.")
    }
  }

  fit_args_base <- list(...)
  fold_info <- .build_fold_data(
    X = X,
    groups = groups,
    k = k,
    lambda = lambda,
    graph = graph,
    folds = folds,
    nfolds = nfolds,
    project_args = project_args,
    fit_args = fit_args_base,
    scale = scale
  )
  fold_data <- fold_info$fold_data
  fold_id <- fold_info$fold_id
  fold_list <- fold_info$fold_list

  cv_pred <- .cv_predict(fold_data, groups, classifier, positive)
  stat_obs <- .metric_score(groups, cv_pred, metric, positive)

  futile.logger::flog.info(
    "Spatial NMF global test: %d permutations (%s, %s)",
    nperm,
    permute,
    if (isTRUE(parallel)) "parallel" else "serial"
  )

  if (permute == "full" && !is.null(folds)) {
    warning("permute='full' requested but folds supplied; using fixed folds (equivalent to permute='labels').")
    permute <- "labels"
  }
  perm_stats <- .global_perm_stats(
    nperm = nperm,
    permute = permute,
    X = X,
    groups = groups,
    fold_data = fold_data,
    folds = fold_list,
    nfolds = nfolds,
    k = k,
    lambda = lambda,
    graph = graph,
    neighbors = neighbors,
    metric = metric,
    classifier = classifier,
    scale = scale,
    positive = positive,
    project_args = project_args,
    fit_args = fit_args_base,
    parallel = parallel,
    future_seed = future_seed,
    progress = progress
  )

  futile.logger::flog.info("Spatial NMF global test: permutations complete.")

  p_value <- (sum(perm_stats >= stat_obs, na.rm = TRUE) + 1) / (nperm + 1)

  list(
    stat = stat_obs,
    p_value = p_value,
    metric = metric,
    classifier = classifier,
    nperm = nperm,
    nfolds = length(fold_list),
    positive = positive,
    perm_stats = if (isTRUE(return_perm)) perm_stats else NULL,
    cv = if (isTRUE(return_cv)) list(pred = cv_pred, fold_id = fold_id) else NULL
  )
}

#' Bootstrap Stability for Spatial NMF Components
#'
#' Quantifies how stable the learned component maps are to resampling subjects.
#' For each bootstrap (or subsample), the NMF is re-fit, components are matched
#' to the reference solution, and summary maps are accumulated.
#'
#' Use this to assess whether components are reproducible and which voxels are
#' consistently among the strongest loadings.
#'
#' @param x A spatial_nmf_maps_result or spatial_nmf_fit object.
#' @param X Optional data matrix (n x p) if not included in x.
#' @param fit Optional spatial_nmf_fit object (if x is not provided).
#' @param graph Optional graph Laplacian list (required if lambda > 0).
#' @param lambda Spatial regularization strength (defaults to fit$lambda).
#' @param n_boot Number of bootstrap samples.
#' @param sample One of "bootstrap" or "subsample".
#' @param sample_frac Fraction of subjects to sample.
#' @param init NMF initialization for bootstrap fits.
#' @param fast Logical; use faster defaults for each bootstrap fit (e.g., random init,
#'   fewer iterations). You can still override specific optimization settings via `...`.
#' @param normalize Component normalization ("H" rescales rows to sum 1).
#' @param similarity Similarity measure for component matching ("cosine" or "cor").
#' @param match Matching strategy (currently "greedy").
#' @param top_frac Fraction of top voxels used to compute selection frequency.
#' @param seed Optional RNG seed.
#' @param return_maps Logical; return stability maps as NeuroVol/NeuroSurface.
#' @param parallel Logical; use future_lapply for bootstrap resamples (requires future.apply).
#' @param future_seed Optional seed control for future.apply (passed to future_lapply).
#' @param progress Logical; report progress via progressr (works with parallel futures).
#' @param ... Additional arguments passed to spatial_nmf_fit.
#'
#' @return A list with stability summaries and optional maps.
#'
#' @details
#' \itemize{
#'   \item \code{mean}, \code{sd}, \code{cv}: bootstrap mean/SD/CV of component maps.
#'   \item \code{selection}: frequency with which a voxel appears in the top fraction.
#'   \item \code{component_similarity}: average similarity to the reference components.
#' }
#' @examples
#' \dontrun{
#'   stab <- spatial_nmf_stability(
#'     matrix(rnorm(100*10), 100, 10),
#'     k = 3, nruns = 5
#'   )
#' }
#' @export
spatial_nmf_stability <- function(x = NULL,
                                  X = NULL,
                                  fit = NULL,
                                  graph = NULL,
                                  lambda = NULL,
                                  n_boot = 200,
                                  sample = c("bootstrap", "subsample"),
                                  sample_frac = 1,
                                  init = c("nndsvd", "random"),
                                  fast = FALSE,
                                  normalize = c("H", "none"),
                                  similarity = c("cosine", "cor"),
                                  match = c("greedy"),
                                  top_frac = 0.1,
                                  seed = NULL,
                                  return_maps = FALSE,
                                  parallel = FALSE,
                                  future_seed = TRUE,
                                  progress = FALSE,
                                  ...) {
  sample <- match.arg(sample)
  init <- match.arg(init)
  normalize <- match.arg(normalize)
  similarity <- match.arg(similarity)
  match <- match.arg(match)

  mask_idx <- NULL
  map_type <- NULL
  mask <- NULL
  dims <- NULL
  ref_map <- NULL
  full_length <- NULL

  if (!is.null(x) && inherits(x, "spatial_nmf_maps_result")) {
    fit <- x$fit
    if (is.null(X)) X <- x$data
    mask_idx <- x$mask_indices
    map_type <- x$map_type
    mask <- x$mask
    dims <- x$dims
    ref_map <- x$ref_map
    full_length <- x$full_length
  } else if (!is.null(x) && inherits(x, "spatial_nmf_fit")) {
    fit <- x
  }

  if (is.null(fit)) stop("fit is required.")
  if (is.null(X)) stop("X is required (pass return_data=TRUE in spatial_nmf_maps).")

  X <- as.matrix(X)
  n <- nrow(X)
  k <- nrow(fit$H)
  p <- ncol(fit$H)

  if (is.null(lambda)) lambda <- fit$lambda
  if (lambda > 0 && is.null(graph)) {
    stop("graph is required for stability when lambda > 0.")
  }
  if (!is.numeric(sample_frac) || sample_frac <= 0 || sample_frac > 1) {
    stop("sample_frac must be in (0, 1].")
  }
  if (!is.numeric(top_frac) || top_frac <= 0 || top_frac > 1) {
    stop("top_frac must be in (0, 1].")
  }

  fit_args <- list(...)
  if (isTRUE(fast)) {
    if (missing(init)) init <- "random"
    if (is.null(fit_args$max_iter)) fit_args$max_iter <- 100
    if (is.null(fit_args$min_iter)) fit_args$min_iter <- 3
    if (is.null(fit_args$tol)) fit_args$tol <- 1e-3
    if (is.null(fit_args$check_every)) fit_args$check_every <- 5
  }

  if (!is.null(seed)) set.seed(seed)

  H_ref <- fit$H
  if (normalize == "H") {
    H_ref <- .normalize_components(H_ref)
  }

  futile.logger::flog.info(
    "Spatial NMF stability: %d bootstrap samples (%s)",
    n_boot,
    if (isTRUE(parallel)) "parallel" else "serial"
  )

  boot_res <- .stability_bootstrap(
    n_boot = n_boot,
    n = n,
    k = k,
    p = p,
    X = X,
    graph = graph,
    lambda = lambda,
    init = init,
    sample = sample,
    sample_frac = sample_frac,
    normalize = normalize,
    similarity = similarity,
    match = match,
    H_ref = H_ref,
    top_frac = top_frac,
    parallel = parallel,
    future_seed = future_seed,
    progress = progress,
    fit_args = fit_args
  )

  futile.logger::flog.info("Spatial NMF stability: bootstrap complete.")

  mean_acc <- boot_res$mean
  m2_acc <- boot_res$m2
  sel_counts <- boot_res$sel_counts
  sim_sum <- boot_res$sim_sum

  sd_acc <- sqrt(m2_acc / pmax(n_boot - 1, 1))
  cv_acc <- sd_acc / pmax(mean_acc, .Machine$double.eps)
  selection_freq <- sel_counts / n_boot
  comp_similarity <- sim_sum / n_boot

  out <- list(
    mean = mean_acc,
    sd = sd_acc,
    cv = cv_acc,
    selection = selection_freq,
    component_similarity = comp_similarity,
    n_boot = n_boot,
    k = k,
    p = p
  )

  if (isTRUE(return_maps)) {
    if (is.null(mask_idx) || is.null(map_type)) {
      stop("return_maps requires spatial metadata (use spatial_nmf_maps result).")
    }
    out$maps <- list(
      mean = .components_to_maps(mean_acc, mask_idx, map_type, mask, dims, full_length, ref_map),
      sd = .components_to_maps(sd_acc, mask_idx, map_type, mask, dims, full_length, ref_map),
      cv = .components_to_maps(cv_acc, mask_idx, map_type, mask, dims, full_length, ref_map),
      selection = .components_to_maps(selection_freq, mask_idx, map_type, mask, dims, full_length, ref_map)
    )
  }

  out
}

#' Voxelwise Statistics from Spatial NMF Stability
#'
#' Convenience function that converts stability summaries into voxelwise z- and
#' p-value maps (based on bootstrap mean/SD).
#'
#' @param x Optional spatial_nmf_maps_result containing stability results.
#' @param stability Optional spatial_nmf_stability result (overrides x$stability).
#' @param map_type Map type ("volume" or "surface") when providing stability without maps.
#' @param mask Mask object for volumetric maps or surface mask (optional for surface).
#' @param dims Optional spatial dimensions for volumetric masks.
#' @param mask_indices Optional mask indices used to vectorize maps.
#' @param full_length Optional full map length for surface/vectorized outputs.
#' @param ref_map Optional reference map to copy metadata from.
#'
#' @return A list with `z` and `p` component maps.
#' @examples
#' \dontrun{
#'   # stats <- spatial_nmf_voxelwise_stats(nmf_result, design_matrix)
#' }
#' @importFrom neuroim2 values
#' @export
spatial_nmf_voxelwise_stats <- function(x = NULL,
                                        stability = NULL,
                                        map_type = NULL,
                                        mask = NULL,
                                        dims = NULL,
                                        mask_indices = NULL,
                                        full_length = NULL,
                                        ref_map = NULL) {
  if (is.null(stability) && !is.null(x)) {
    if (!inherits(x, "spatial_nmf_maps_result")) {
      stop("x must be a spatial_nmf_maps_result when provided.")
    }
    stability <- x$stability
    if (is.null(map_type)) map_type <- x$map_type
    if (is.null(mask)) mask <- x$mask
    if (is.null(dims)) dims <- x$dims
    if (is.null(mask_indices)) mask_indices <- x$mask_indices
    if (is.null(full_length)) full_length <- x$full_length
    if (is.null(ref_map)) ref_map <- x$ref_map
  }
  if (is.null(stability)) stop("stability is required.")

  mean_maps <- NULL
  sd_maps <- NULL
  if (!is.null(stability$maps)) {
    mean_maps <- stability$maps$mean
    sd_maps <- stability$maps$sd
  }

  if (is.null(mean_maps) || is.null(sd_maps)) {
    if (is.null(map_type) || is.null(mask_indices) || is.null(ref_map)) {
      stop("stability maps are missing; supply spatial metadata or use spatial_nmf_maps with stability$return_maps=TRUE.")
    }
    mean_maps <- .components_to_maps(
      stability$mean,
      mask_idx = mask_indices,
      map_type = map_type,
      mask = mask,
      dims = dims,
      full_length = full_length,
      ref_map = ref_map
    )
    sd_maps <- .components_to_maps(
      stability$sd,
      mask_idx = mask_indices,
      map_type = map_type,
      mask = mask,
      dims = dims,
      full_length = full_length,
      ref_map = ref_map
    )
  }

  k <- length(mean_maps)
  z_maps <- vector("list", k)
  p_maps <- vector("list", k)
  for (i in seq_len(k)) {
    z_vals <- neuroim2::values(mean_maps[[i]]) /
      pmax(neuroim2::values(sd_maps[[i]]), .Machine$double.eps)
    p_vals <- 2 * stats::pnorm(-abs(z_vals))
    ref <- mean_maps[[i]]
    if (inherits(ref, "SparseNeuroVol")) {
      z_maps[[i]] <- neuroim2::SparseNeuroVol(
        data = as.numeric(z_vals), space = neuroim2::space(ref),
        indices = as.integer(neuroim2::indices(ref)))
      p_maps[[i]] <- neuroim2::SparseNeuroVol(
        data = as.numeric(p_vals), space = neuroim2::space(ref),
        indices = as.integer(neuroim2::indices(ref)))
    } else if (inherits(ref, "NeuroSurface")) {
      z_maps[[i]] <- neurosurf::NeuroSurface(
        geometry = neurosurf::geometry(ref),
        indices = ref@indices, data = as.numeric(z_vals))
      p_maps[[i]] <- neurosurf::NeuroSurface(
        geometry = neurosurf::geometry(ref),
        indices = ref@indices, data = as.numeric(p_vals))
    } else {
      z_maps[[i]] <- neuroim2::NeuroVol(as.numeric(z_vals), neuroim2::space(ref))
      p_maps[[i]] <- neuroim2::NeuroVol(as.numeric(p_vals), neuroim2::space(ref))
    }
  }

  list(z = z_maps, p = p_maps)
}

.permute_component_stats <- function(nperm,
                                     groups,
                                     covariates,
                                     W,
                                     coef_index,
                                     parallel = FALSE,
                                     future_seed = TRUE,
                                     progress = FALSE) {
  iter_fun <- function(i) {
    perm_groups <- sample(groups)
    perm_design <- .component_design(perm_groups, covariates)
    .component_t_stats(W, perm_design, coef_index = coef_index)
  }

  if (isTRUE(parallel)) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("parallel=TRUE requires the future.apply package.")
    }
    if (isTRUE(progress)) {
      if (!requireNamespace("progressr", quietly = TRUE)) {
        stop("progress=TRUE requires the progressr package.")
      }
      return(progressr::with_progress({
        p <- progressr::progressor(steps = nperm)
        perm_list <- future.apply::future_lapply(
          seq_len(nperm),
          function(i) {
            res <- iter_fun(i)
            p()
            res
          },
          future.seed = future_seed
        )
        do.call(rbind, perm_list)
      }))
    }
    perm_list <- future.apply::future_lapply(seq_len(nperm), iter_fun, future.seed = future_seed)
    return(do.call(rbind, perm_list))
  }

  perm_stats <- matrix(NA_real_, nrow = nperm, ncol = ncol(W))
  if (isTRUE(progress)) {
    if (!requireNamespace("progressr", quietly = TRUE)) {
      stop("progress=TRUE requires the progressr package.")
    }
    return(progressr::with_progress({
      p <- progressr::progressor(steps = nperm)
      for (i in seq_len(nperm)) {
        perm_stats[i, ] <- iter_fun(i)
        p()
      }
      perm_stats
    }))
  }
  step <- .progress_step(nperm)
  for (i in seq_len(nperm)) {
    perm_stats[i, ] <- iter_fun(i)
    .maybe_log_progress(i, nperm, step, "Component permutations")
  }
  perm_stats
}

.global_perm_stats <- function(nperm,
                               permute,
                               X,
                               groups,
                               fold_data,
                               folds,
                               nfolds,
                               k,
                               lambda,
                               graph,
                               neighbors,
                               metric,
                               classifier,
                               scale,
                               positive,
                               project_args,
                               fit_args,
                               parallel = FALSE,
                               future_seed = TRUE,
                               progress = FALSE) {
  iter_fun <- function(i) {
    perm_groups <- sample(groups)
    if (permute == "labels") {
      perm_pred <- .cv_predict(fold_data, perm_groups, classifier, positive)
      .metric_score(perm_groups, perm_pred, metric, positive)
    } else {
      perm_fold <- suppressWarnings(.build_fold_data(
        X = X,
        groups = perm_groups,
        k = k,
        lambda = lambda,
        graph = graph,
        folds = NULL,
        nfolds = nfolds,
        project_args = project_args,
        fit_args = fit_args,
        scale = scale
      ))
      perm_pred <- .cv_predict(perm_fold$fold_data, perm_groups, classifier, positive)
      .metric_score(perm_groups, perm_pred, metric, positive)
    }
  }

  if (isTRUE(parallel)) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("parallel=TRUE requires the future.apply package.")
    }
    if (isTRUE(progress)) {
      if (!requireNamespace("progressr", quietly = TRUE)) {
        stop("progress=TRUE requires the progressr package.")
      }
      return(progressr::with_progress({
        p <- progressr::progressor(steps = nperm)
        unlist(future.apply::future_lapply(
          seq_len(nperm),
          function(i) {
            res <- iter_fun(i)
            p()
            res
          },
          future.seed = future_seed
        ))
      }))
    }
    return(unlist(future.apply::future_lapply(seq_len(nperm), iter_fun, future.seed = future_seed)))
  }

  perm_stats <- numeric(nperm)
  if (isTRUE(progress)) {
    if (!requireNamespace("progressr", quietly = TRUE)) {
      stop("progress=TRUE requires the progressr package.")
    }
    return(progressr::with_progress({
      p <- progressr::progressor(steps = nperm)
      for (i in seq_len(nperm)) {
        perm_stats[i] <- iter_fun(i)
        p()
      }
      perm_stats
    }))
  }
  step <- .progress_step(nperm)
  for (i in seq_len(nperm)) {
    perm_stats[i] <- iter_fun(i)
    .maybe_log_progress(i, nperm, step, "Global permutations")
  }
  perm_stats
}

.stability_bootstrap <- function(n_boot,
                                 n,
                                 k,
                                 p,
                                 X,
                                 graph,
                                 lambda,
                                 init,
                                 sample,
                                 sample_frac,
                                 normalize,
                                 similarity,
                                 match,
                                 H_ref,
                                 top_frac,
                                 parallel = FALSE,
                                 future_seed = TRUE,
                                 progress = FALSE,
                                 fit_args = list()) {
  iter_fun <- function(i) {
    if (sample == "bootstrap") {
      idx <- sample.int(n, size = ceiling(n * sample_frac), replace = TRUE)
    } else {
      idx <- sample.int(n, size = ceiling(n * sample_frac), replace = FALSE)
    }
    fit_call <- c(
      list(
        X = X[idx, , drop = FALSE],
        k = k,
        graph = graph,
        lambda = lambda,
        init = init
      ),
      fit_args
    )
    fit_b <- do.call(spatial_nmf_fit, fit_call)
    H_b <- fit_b$H
    sim <- .component_similarity(H_ref, H_b, similarity)
    order <- .match_components(sim, match)
    H_b <- H_b[order, , drop = FALSE]
    if (normalize == "H") {
      H_b <- .normalize_components(H_b)
    }
    top_idx <- NULL
    if (!is.null(top_frac) && top_frac > 0) {
      top_n <- max(1L, ceiling(ncol(H_b) * top_frac))
      top_idx <- lapply(seq_len(nrow(H_b)), function(j) {
        order(H_b[j, ], decreasing = TRUE)[seq_len(top_n)]
      })
    }
    list(H = H_b, sim = sim[cbind(seq_len(k), order)], top_idx = top_idx)
  }

  if (isTRUE(parallel)) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("parallel=TRUE requires the future.apply package.")
    }
    if (isTRUE(progress)) {
      if (!requireNamespace("progressr", quietly = TRUE)) {
        stop("progress=TRUE requires the progressr package.")
      }
      res_list <- progressr::with_progress({
        p <- progressr::progressor(steps = n_boot)
        future.apply::future_lapply(
          seq_len(n_boot),
          function(i) {
            res <- iter_fun(i)
            p()
            res
          },
          future.seed = future_seed
        )
      })
    } else {
      res_list <- future.apply::future_lapply(seq_len(n_boot), iter_fun, future.seed = future_seed)
    }
  } else {
    res_list <- vector("list", n_boot)
    if (isTRUE(progress)) {
      if (!requireNamespace("progressr", quietly = TRUE)) {
        stop("progress=TRUE requires the progressr package.")
      }
      progressr::with_progress({
        p <- progressr::progressor(steps = n_boot)
        for (b in seq_len(n_boot)) {
          res_list[[b]] <- iter_fun(b)
          p()
        }
      })
    } else {
      step <- .progress_step(n_boot)
      for (b in seq_len(n_boot)) {
        res_list[[b]] <- iter_fun(b)
        .maybe_log_progress(b, n_boot, step, "Stability bootstrap")
      }
    }
  }

  mean_acc <- matrix(0, nrow = k, ncol = p)
  m2_acc <- matrix(0, nrow = k, ncol = p)
  sel_counts <- matrix(0, nrow = k, ncol = p)
  sim_sum <- numeric(k)

  for (b in seq_len(n_boot)) {
    H_b <- res_list[[b]]$H
    sim_sum <- sim_sum + res_list[[b]]$sim
    delta <- H_b - mean_acc
    mean_acc <- mean_acc + delta / b
    m2_acc <- m2_acc + delta * (H_b - mean_acc)
    if (!is.null(res_list[[b]]$top_idx)) {
      for (j in seq_len(k)) {
        sel_counts[j, res_list[[b]]$top_idx[[j]]] <- sel_counts[j, res_list[[b]]$top_idx[[j]]] + 1
      }
    }
  }

  list(mean = mean_acc, m2 = m2_acc, sel_counts = sel_counts, sim_sum = sim_sum)
}

.progress_step <- function(n, n_updates = 10L) {
  if (!is.numeric(n) || length(n) != 1L || n < 1) return(1L)
  step <- floor(n / n_updates)
  if (step < 1L) 1L else as.integer(step)
}

.maybe_log_progress <- function(i, n, step, label) {
  if (i %% step == 0 || i == n) {
    futile.logger::flog.debug("%s %d/%d", label, i, n)
  }
  invisible(NULL)
}

.component_design <- function(groups, covariates) {
  data <- data.frame(group = groups)
  if (!is.null(covariates)) {
    if (!is.data.frame(covariates)) covariates <- as.data.frame(covariates)
    if (nrow(covariates) != nrow(data)) {
      stop("covariates must have the same number of rows as W.")
    }
    data <- cbind(data, covariates)
  }
  stats::model.matrix(~ ., data = data)
}

.component_t_stats <- function(W, design, coef_index) {
  if (!is.numeric(coef_index) || length(coef_index) != 1L) {
    stop("coef_index must be a single integer.")
  }
  coef_index <- as.integer(coef_index)
  if (coef_index < 1L || coef_index > ncol(design)) {
    stop("coef_index is out of bounds for design matrix.")
  }
  qrX <- qr(design)
  coef <- qr.coef(qrX, W)
  fitted <- design %*% coef
  resid <- W - fitted
  df <- nrow(W) - qrX$rank
  if (df <= 0) return(rep(NA_real_, ncol(W)))
  sigma2 <- colSums(resid^2) / df
  xtx_inv <- tryCatch(
    {
      R <- qr.R(qrX)
      chol2inv(R)
    },
    error = function(e) MASS::ginv(crossprod(design))
  )
  if (coef_index > nrow(xtx_inv) || coef_index > ncol(xtx_inv)) {
    return(rep(NA_real_, ncol(W)))
  }
  se <- sqrt(xtx_inv[coef_index, coef_index] * sigma2)
  out <- coef[coef_index, ] / se
  out[!is.finite(out)] <- NA_real_
  out
}

.perm_pvals <- function(stat_obs, perm_stats, two_sided = TRUE) {
  nperm <- nrow(perm_stats)
  k <- length(stat_obs)
  if (two_sided) {
    cmp <- abs(perm_stats) >= matrix(abs(stat_obs), nrow = nperm, ncol = k, byrow = TRUE)
  } else {
    cmp <- perm_stats >= matrix(stat_obs, nrow = nperm, ncol = k, byrow = TRUE)
  }
  (colSums(cmp) + 1) / (nperm + 1)
}

.as_null_W_list <- function(null_W, n, k) {
  if (is.list(null_W)) {
    if (!all(vapply(null_W, function(mat) is.matrix(mat) && all(dim(mat) == c(n, k)), logical(1)))) {
      stop("Each null_W entry must be an n x k matrix.")
    }
    return(null_W)
  }

  if (is.array(null_W) && length(dim(null_W)) == 3L) {
    dims <- dim(null_W)
    if (dims[2] == n && dims[3] == k) {
      return(lapply(seq_len(dims[1]), function(i) null_W[i, , ]))
    }
    if (dims[1] == n && dims[2] == k) {
      return(lapply(seq_len(dims[3]), function(i) null_W[, , i]))
    }
  }

  stop("null_W must be a list of matrices or a 3D array matching n and k.")
}

.normalize_components <- function(H) {
  scale <- pmax(rowSums(H), .Machine$double.eps)
  H / scale
}

.component_similarity <- function(H_ref, H_boot, method) {
  if (method == "cosine") {
    ref_norm <- H_ref / sqrt(pmax(rowSums(H_ref^2), .Machine$double.eps))
    boot_norm <- H_boot / sqrt(pmax(rowSums(H_boot^2), .Machine$double.eps))
    ref_norm %*% t(boot_norm)
  } else {
    sim <- stats::cor(t(H_ref), t(H_boot))
    sim[is.na(sim)] <- 0
    sim
  }
}

.match_components <- function(sim_mat, method) {
  if (method != "greedy") stop("Only greedy matching is implemented.")
  k <- nrow(sim_mat)
  pairs <- expand.grid(ref = seq_len(k), boot = seq_len(k))
  pairs$sim <- as.vector(sim_mat)
  pairs <- pairs[order(-pairs$sim), , drop = FALSE]

  ref_used <- rep(FALSE, k)
  boot_used <- rep(FALSE, k)
  order <- integer(k)

  for (i in seq_len(nrow(pairs))) {
    r <- pairs$ref[i]
    b <- pairs$boot[i]
    if (!ref_used[r] && !boot_used[b]) {
      order[r] <- b
      ref_used[r] <- TRUE
      boot_used[b] <- TRUE
    }
    if (all(ref_used)) break
  }

  if (any(order == 0)) {
    remaining_ref <- which(order == 0)
    remaining_boot <- setdiff(seq_len(k), order)
    order[remaining_ref] <- remaining_boot
  }

  order
}

.prepare_binary_groups <- function(groups, positive) {
  groups <- factor(groups)
  if (nlevels(groups) != 2L) {
    stop("groups must have exactly 2 levels.")
  }
  if (is.null(positive)) {
    positive <- levels(groups)[2]
  }
  if (!positive %in% levels(groups)) {
    stop("positive must match a level in groups.")
  }
  list(groups = groups, positive = positive)
}

.coerce_folds <- function(folds, groups, nfolds) {
  n <- length(groups)
  if (is.null(folds)) {
    fold_id <- .stratified_folds(groups, nfolds)
  } else if (is.list(folds)) {
    fold_list <- lapply(folds, function(idx) as.integer(idx))
    all_idx <- sort(unlist(fold_list, use.names = FALSE))
    if (!identical(all_idx, seq_len(n))) {
      stop("folds must partition all indices 1..n without overlap.")
    }
    fold_id <- integer(n)
    for (i in seq_along(fold_list)) {
      fold_id[fold_list[[i]]] <- i
    }
  } else if (is.numeric(folds) && length(folds) == n) {
    fold_id <- as.integer(folds)
    if (any(is.na(fold_id))) stop("folds contains NA values.")
  } else {
    stop("folds must be NULL, a fold ID vector, or a list of test indices.")
  }

  fold_list <- split(seq_len(n), fold_id)
  fold_list <- fold_list[order(as.integer(names(fold_list)))]
  list(folds = fold_list, fold_id = fold_id)
}

.stratified_folds <- function(groups, nfolds) {
  groups <- factor(groups)
  if (!is.numeric(nfolds) || length(nfolds) != 1L) {
    stop("nfolds must be a single integer.")
  }
  nfolds <- as.integer(nfolds)
  counts <- table(groups)
  min_count <- min(counts)
  if (nfolds > min_count) {
    warning("nfolds reduced to ", min_count, " to preserve stratification.")
    nfolds <- min_count
  }
  if (nfolds < 2L) stop("nfolds must be at least 2.")

  fold_id <- integer(length(groups))
  for (lvl in levels(groups)) {
    idx <- which(groups == lvl)
    idx <- sample(idx)
    fold_id[idx] <- rep(seq_len(nfolds), length.out = length(idx))
  }
  fold_id
}

.scale_train_test <- function(W_train, W_test) {
  mu <- colMeans(W_train)
  sigma <- apply(W_train, 2, stats::sd)
  sigma[!is.finite(sigma) | sigma == 0] <- 1
  train <- sweep(W_train, 2, mu, "-")
  train <- sweep(train, 2, sigma, "/")
  test <- sweep(W_test, 2, mu, "-")
  test <- sweep(test, 2, sigma, "/")
  list(train = train, test = test)
}

.build_fold_data <- function(X, groups, k, lambda, graph, folds, nfolds, project_args, fit_args, scale) {
  n <- nrow(X)
  fold_info <- .coerce_folds(folds, groups, nfolds)
  fold_list <- fold_info$folds
  fold_id <- fold_info$fold_id

  min_train <- min(vapply(fold_list, function(idx) n - length(idx), integer(1)))
  if (k > min_train) stop("k must be <= min training size across folds.")
  if (k > ncol(X)) stop("k must be <= ncol(X).")
  if (any(table(groups) < 2)) stop("Each group must have at least 2 samples.")

  fold_data <- vector("list", length(fold_list))
  for (i in seq_along(fold_list)) {
    test_idx <- fold_list[[i]]
    train_idx <- setdiff(seq_len(n), test_idx)
    if (length(unique(groups[train_idx])) < 2L) {
      stop("Each fold must have both groups represented in training.")
    }

    fit_args$X <- X[train_idx, , drop = FALSE]
    fit_args$k <- k
    fit_args$graph <- graph
    fit_args$lambda <- lambda

    fit <- do.call(spatial_nmf_fit, fit_args)

    proj_args <- project_args
    proj_args$X <- X[test_idx, , drop = FALSE]
    proj_args$H <- fit$H
    proj <- do.call(spatial_nmf_project, proj_args)

    W_train <- fit$W
    W_test <- proj$W
    if (isTRUE(scale)) {
      scaled <- .scale_train_test(W_train, W_test)
      W_train <- scaled$train
      W_test <- scaled$test
    }

    fold_data[[i]] <- list(
      train_idx = train_idx,
      test_idx = test_idx,
      W_train = W_train,
      W_test = W_test
    )
  }

  list(fold_data = fold_data, fold_id = fold_id, fold_list = fold_list)
}

.cv_predict <- function(fold_data, groups, classifier, positive) {
  n <- length(groups)
  preds <- rep(NA_real_, n)
  for (i in seq_along(fold_data)) {
    info <- fold_data[[i]]
    model <- .fit_binary_classifier(info$W_train, groups[info$train_idx], classifier, positive)
    preds[info$test_idx] <- .predict_binary_classifier(model, info$W_test, positive)
  }
  if (any(!is.finite(preds))) stop("Non-finite predictions in CV.")
  preds
}

.fit_binary_classifier <- function(W_train, groups, classifier, positive) {
  df <- as.data.frame(W_train)
  if (is.null(colnames(df))) colnames(df) <- paste0("V", seq_len(ncol(df)))

  if (classifier == "centroid") {
    mu_pos <- colMeans(W_train[groups == positive, , drop = FALSE])
    mu_neg <- colMeans(W_train[groups != positive, , drop = FALSE])
    return(list(type = "centroid", mu_pos = mu_pos, mu_neg = mu_neg))
  }

  if (classifier == "lda") {
    fit <- tryCatch(
      MASS::lda(x = df, grouping = factor(groups)),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(.fit_binary_classifier(W_train, groups, "centroid", positive))
    }
    return(list(type = "lda", fit = fit))
  }

  y <- as.integer(groups == positive)
  df$y <- y
  fit <- suppressWarnings(
    tryCatch(
      stats::glm(y ~ ., data = df, family = stats::binomial()),
      error = function(e) NULL
    )
  )
  if (is.null(fit) || !isTRUE(fit$converged)) {
    return(.fit_binary_classifier(W_train, groups, "lda", positive))
  }
  list(type = "glm", fit = fit)
}

.predict_binary_classifier <- function(model, W_test, positive) {
  if (model$type == "centroid") {
    score <- as.vector(W_test %*% (model$mu_pos - model$mu_neg))
    return(stats::plogis(score))
  }

  df <- as.data.frame(W_test)
  if (is.null(colnames(df))) colnames(df) <- paste0("V", seq_len(ncol(df)))

  if (model$type == "lda") {
    pred <- predict(model$fit, df)
    probs <- pred$posterior[, positive]
    return(as.numeric(probs))
  }

  as.numeric(stats::predict(model$fit, newdata = df, type = "response"))
}

.metric_score <- function(groups, preds, metric, positive) {
  y <- as.integer(groups == positive)
  if (metric == "accuracy") {
    mean((preds >= 0.5) == y)
  } else {
    .auc_score(y, preds)
  }
}

.auc_score <- function(y, score) {
  y <- as.integer(y)
  pos <- y == 1L
  n_pos <- sum(pos)
  n_neg <- length(y) - n_pos
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[pos]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}
