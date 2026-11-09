# Decoding a structured feature space from local activity patterns -----------

#' Feature decoding model: predict a feature space from activity patterns
#'
#' The pattern-level complement of [banded_ridge_model()]. Within each region
#' or searchlight, the feature columns of a [feature_sets_design()] are
#' predicted from the local activity pattern by cross-validated ridge
#' regression, and the held-out predictions are summarised over the declared
#' structure of the feature space: every set, every family of sets, and every
#' contrast (see [feature_sets()]). Family and contrast summaries are computed
#' from the same held-out predictions, so they need no additional fits.
#'
#' The model shares its cross-validation, temporal purge, nuisance handling and
#' ridge solvers with [banded_ridge_model()], with the roles of the two
#' matrices exchanged. It implements [fit_roi()], so [run_searchlight()] and
#' [run_regional()] work unchanged.
#'
#' @section Metrics:
#' One value per region or searchlight centre:
#' \describe{
#'   \item{`all_r2`, `all_cor`}{Mean held-out R2 and correlation over all
#'     feature columns.}
#'   \item{`<set>_r2`, `<set>_cor`}{The same, over the columns of one set or
#'     family. Every feature column counts equally.}
#'   \item{`<contrast>_r2`, `<contrast>_cor`}{Accuracy for the contrast series,
#'     i.e. observed and predicted features projected onto the contrast.}
#'   \item{`<set>_partial_r2`, `<set>_partial_cor`}{With `partial = TRUE`:
#'     accuracy for a set or family after every other feature column has been
#'     regressed out of both the pattern and the targets, with coefficients
#'     fitted on each training split. This is the pattern-level analogue of
#'     unique variance.}
#'   \item{`<metric>_p`}{With `n_perm > 0`: one-sided p-value against a null in
#'     which the predictions are shifted in time against the observations.}
#' }
#'
#' @param dataset An `mvpa_dataset`; rows of its training data are the
#'   observations (for example volumes of one run).
#' @param design A [feature_sets_design()] whose training features are the
#'   decoding targets, with the same number of rows as the dataset.
#' @param outer_crossval,tune_crossval,purge Cross-validation, as in
#'   [banded_ridge_model()]: an integer fold count, a `cross_validation`
#'   object, or explicit folds; `tune_crossval` sets the inner folds used to
#'   choose the ridge penalty; `purge` drops training rows within that many
#'   rows of a held-out row.
#' @param alphas Candidate ridge penalties. Each feature column selects its
#'   own penalty inside each outer-training set.
#' @param nuisance,nuisance_fit Optional nuisance matrix (rows = observations)
#'   regressed out of patterns and targets inside every split, and how held-out
#'   rows are residualized; see [banded_ridge_model()].
#' @param partial If `TRUE`, also report partial accuracies (see Metrics).
#' @param n_perm Number of time-shift permutations; `0` (default) disables the
#'   null.
#' @param perm_strategy,perm_block,min_shift Permutation scheme:
#'   `"circular_shift"` (default) or `"block_shuffle"`, the block length for the
#'   latter, and the smallest circular shift in rows.
#' @param min_voxels Regions with fewer voxels are returned as errors.
#' @param seed Integer seed for fold construction and permutations.
#' @param ... Passed to the model specification.
#'
#' @return A model specification of class `feature_decode_model`.
#' @seealso [banded_ridge_model()], [feature_sets()], [run_searchlight()]
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 80
#' ds <- gen_sample_dataset(c(5, 5, 5), nobs = n, blocks = 4)
#' X <- matrix(rnorm(n * 3), n, 3, dimnames = list(NULL, c("m1", "p1", "ac")))
#' fs <- feature_sets(X, by_set(c("manner", "place", "acoustic")),
#'                    families = list(phonology = c("manner", "place")))
#' des <- feature_sets_design(fs, block_var_train = rep(1:4, each = 20),
#'                            time_series = TRUE)
#' model <- feature_decode_model(ds$dataset, des, outer_crossval = 4)
#' res <- run_searchlight(model, radius = 4, method = "standard")
#' }
#' @export
feature_decode_model <- function(dataset,
                                 design,
                                 outer_crossval = NULL,
                                 tune_crossval = NULL,
                                 alphas = 10^seq(-1, 5),
                                 purge = 0L,
                                 nuisance = NULL,
                                 nuisance_fit = c("within", "train"),
                                 partial = FALSE,
                                 n_perm = 0L,
                                 perm_strategy = c("circular_shift", "block_shuffle"),
                                 perm_block = NULL,
                                 min_shift = 1L,
                                 min_voxels = 2L,
                                 seed = 1L,
                                 ...) {
  perm_strategy <- match.arg(perm_strategy)
  nuisance_fit <- match.arg(nuisance_fit)
  if (!inherits(design, "feature_sets_design")) {
    stop("feature_decode_model: design must be created by feature_sets_design().",
         call. = FALSE)
  }
  whole <- function(x, minimum) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x == as.integer(x) &&
      x >= minimum
  }
  if (!whole(purge, 0L)) stop("purge must be one non-negative integer.", call. = FALSE)
  if (!whole(n_perm, 0L)) stop("n_perm must be one non-negative integer.", call. = FALSE)
  if (!whole(min_shift, 1L)) stop("min_shift must be one positive integer.", call. = FALSE)
  if (!whole(min_voxels, 1L)) stop("min_voxels must be one positive integer.", call. = FALSE)
  if (!whole(seed, -Inf)) stop("seed must be one finite integer.", call. = FALSE)
  if (!is.logical(partial) || length(partial) != 1L || is.na(partial)) {
    stop("partial must be TRUE or FALSE.", call. = FALSE)
  }
  targets <- design$X_train
  n <- nrow(targets$X)
  if (nobs(dataset) != n) {
    stop("feature_decode_model: design rows must match dataset observations.",
         call. = FALSE)
  }
  nuisance <- .brt_validate_nuisance(nuisance, n, "feature_decode_model: nuisance",
                                     fit = nuisance_fit)
  blocks <- design$block_var_train
  if (isTRUE(design$time_series) && is.null(blocks) && purge < 1L) {
    stop("feature_decode_model: time_series = TRUE requires block_var_train or ",
         "a positive temporal purge; random row-wise CV is unsafe.", call. = FALSE)
  }
  outer_folds <- .brm_resolve_outer_folds(
    outer_crossval, n, blocks, as.integer(purge), as.integer(seed)
  )
  # One band (the pattern), so only the overall penalty is tuned.
  candidates <- .banded_ridge_candidates(
    "pattern", alphas, method = "fixed", theta = c(pattern = 1)
  )
  tuning <- .brm_resolve_tuning(
    tune_crossval, outer_folds, n, blocks, as.integer(purge), as.integer(seed),
    n_available = length(candidates$alpha)
  )
  set <- as.character(targets$set)
  groups <- .fts_groups(set, targets$families)
  spec <- create_model_spec(
    "feature_decode_model",
    dataset = dataset, design = design,
    outer_folds = outer_folds, inner_folds = tuning$inner_folds,
    inner_v = tuning$inner_v, candidates = candidates,
    purge = as.integer(purge), nuisance = nuisance,
    partial = partial,
    # a group that spans every column has nothing left to partial out
    partial_groups = if (partial) {
      groups[lengths(groups) < length(set)]
    } else list(),
    n_perm = as.integer(n_perm), perm_strategy = perm_strategy,
    perm_block = perm_block, min_shift = as.integer(min_shift),
    min_voxels = as.integer(min_voxels), seed = as.integer(seed),
    compute_performance = TRUE, return_predictions = FALSE,
    ...
  )
  spec$metric_names <- .fdm_metric_names(spec)
  clash <- unique(spec$metric_names[duplicated(spec$metric_names)])
  if (length(clash)) {
    # e.g. a set named "manner_partial" next to the partial score of "manner"
    stop("feature_decode_model: set, family and contrast names produce ",
         "duplicate metric names (", paste(clash, collapse = ", "),
         "); rename them.", call. = FALSE)
  }
  spec
}

# Held-out predictions of `targets` from `pattern`, with the scores'
# observed counterpart (targets, or their nuisance-residualized version).
.fdm_predict <- function(model, pattern, targets, covariates = NULL) {
  wide <- ncol(pattern) > nrow(pattern)
  nested <- .banded_ridge_nested_cv(
    pattern, targets, model$outer_folds, model$candidates,
    feature_groups = rep("pattern", ncol(pattern)),
    inner_folds = model$inner_folds, inner_v = model$inner_v,
    blocks = model$design$block_var_train, purge = model$purge,
    metric = "mse", seed = model$seed,
    solver = if (wide) "dual_kernel" else "svd_primal",
    solver_cache = .banded_ridge_solver_cache(),
    nuisance = model$nuisance, covariates = covariates
  )
  list(predicted = nested$predictions,
       observed = if (is.null(nested$observed)) targets else nested$observed)
}

.fdm_metric_names <- function(model) {
  targets <- model$design$X_train
  labels <- c("all", names(.fts_groups(as.character(targets$set), targets$families)),
              colnames(targets$contrasts))
  main <- as.vector(rbind(paste0(labels, "_r2"), paste0(labels, "_cor")))
  partial <- as.vector(rbind(paste0(names(model$partial_groups), "_partial_r2"),
                             paste0(names(model$partial_groups), "_partial_cor")))
  if (!length(model$partial_groups)) partial <- character()
  out <- c(main, partial)
  if (model$n_perm > 0L) out <- c(out, paste0(out, "_p"))
  out
}

#' @rdname output_schema
#' @method output_schema feature_decode_model
#' @export
output_schema.feature_decode_model <- function(model) {
  stats::setNames(as.list(rep("scalar", length(model$metric_names))),
                  model$metric_names)
}

#' @rdname fit_roi
#' @method fit_roi feature_decode_model
#' @export
fit_roi.feature_decode_model <- function(model, roi_data, context, ...) {
  indices <- roi_data$indices
  failed <- function(message) {
    roi_result(metrics = NULL, indices = indices, id = context$id,
               error = TRUE, error_message = paste0("feature_decode_model: ", message))
  }
  pattern <- roi_data$train_data
  fs <- model$design$X_train
  targets <- fs$X
  if (is.null(colnames(targets))) {
    colnames(targets) <- paste0("feature_", seq_len(ncol(targets)))
  }
  if (nrow(pattern) != nrow(targets)) {
    return(failed("pattern rows do not match the design."))
  }
  pattern <- pattern[, colSums(!is.finite(pattern)) == 0L, drop = FALSE]
  pattern <- pattern[, apply(pattern, 2L, stats::sd) > 0, drop = FALSE]
  if (ncol(pattern) < model$min_voxels) {
    return(failed("too few usable voxels."))
  }
  colnames(pattern) <- paste0("v", seq_len(ncol(pattern)))
  set <- as.character(fs$set)

  result <- tryCatch({
    .brt_with_seed(model$seed, {
      fit <- .fdm_predict(model, pattern, targets)
      scores <- .feature_target_scores(fit$observed, fit$predicted, set,
                                       fs$families, fs$contrasts)
      null <- if (model$n_perm > 0L) {
        .feature_target_null(fit$observed, fit$predicted, set, fs$families,
                             fs$contrasts, n_perm = model$n_perm,
                             strategy = model$perm_strategy,
                             block_size = model$perm_block,
                             min_shift = model$min_shift)
      } else NULL
      for (label in names(model$partial_groups)) {
        columns <- model$partial_groups[[label]]
        part <- .fdm_predict(
          model, pattern, targets[, columns, drop = FALSE],
          covariates = targets[, -columns, drop = FALSE]
        )
        partial_names <- paste0(label, c("_partial_r2", "_partial_cor"))
        own <- rep("x", length(columns))
        scores[partial_names] <-
          .feature_target_scores(part$observed, part$predicted, own)[c("all_r2", "all_cor")]
        if (!is.null(null)) {
          part_null <- .feature_target_null(
            part$observed, part$predicted, own, n_perm = model$n_perm,
            strategy = model$perm_strategy, block_size = model$perm_block,
            min_shift = model$min_shift
          )[, c("all_r2", "all_cor"), drop = FALSE]
          colnames(part_null) <- partial_names
          null <- cbind(null, part_null)
        }
      }
      if (!is.null(null)) {
        p <- .fts_p_values(scores, null[, names(scores), drop = FALSE])
        scores <- c(scores, stats::setNames(p, paste0(names(scores), "_p")))
      }
      scores
    })
  }, error = function(e) e)
  if (inherits(result, "error")) return(failed(conditionMessage(result)))
  roi_result(metrics = result[model$metric_names], indices = indices,
             id = context$id)
}
