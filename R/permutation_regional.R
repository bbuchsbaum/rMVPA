#' Permutation-Based Inference for Regional MVPA
#'
#' Computes permutation p-values for every region of a regional MVPA analysis,
#' with the same machinery as \code{\link{run_permutation_searchlight}}: the
#' same \code{\link{permutation_control}} object, the same relabelling
#' (\code{\link{permute_labels}} with \code{shuffle} and seed
#' \code{seed + i} for permutation \code{i}), the same per-model adaptation to
#' a permuted design, and the same p-value convention.
#'
#' Unlike the searchlight, every region is refitted under every permutation, so
#' each region gets its own exact null distribution: nothing is pooled across
#' regions, and \code{null_method}, \code{adjust_by}, \code{n_bins} and
#' \code{subsample} do not apply. For region \eqn{r} and metric \eqn{m},
#' \deqn{p_r = \frac{1 + \#\{b : T^{(b)}_r \ge T_r\}}{1 + B_r}}
#' where \eqn{B_r} is the number of permutations that produced a value for
#' region \eqn{r}. Larger metric values are taken as more extreme, as in the
#' searchlight. The smallest attainable p-value is \eqn{1/(B+1)}.
#'
#' Two corrections across regions are reported, separately for each metric:
#' \describe{
#'   \item{\code{p_adj}}{Benjamini-Hochberg FDR when
#'     \code{perm_ctrl$correction = "fdr"} (the default), otherwise \code{p}.}
#'   \item{\code{p_fwer}}{Family-wise error by the single-step minP method
#'     (Westfall and Young): each permutation's smallest per-region null
#'     p-value forms the null of the minimum. Because every region is
#'     scored under the same relabelling in each permutation, this respects the
#'     dependence between regions; and because it works on p-values, regions of
#'     different sizes and metric scales are comparable.}
#' }
#'
#' Efficiency: for models with the exact aggregate engine (correlation and
#' naive Bayes classifiers with fixed folds, as for the searchlight) the data
#' and folds are prepared once and each permutation only rescores the labels.
#' Otherwise region data are extracted once per batch of regions and reused
#' for every permutation; only the model is refitted. The permuted designs are
#' drawn once, so every batch and every region sees the same permutations.
#' Regions within a permutation run in parallel through the active
#' \pkg{future} plan, as in \code{\link{run_regional}}.
#'
#' @param model_spec A model specification accepted by \code{\link{run_regional}}.
#' @param region_mask Region mask (a \code{NeuroVol} or \code{NeuroSurface}
#'   with integer region ids), as for \code{\link{run_regional}}.
#' @param observed Optional \code{regional_mvpa_result} for the unpermuted
#'   data; computed with \code{run_regional(model_spec, region_mask, ...)} when
#'   \code{NULL}.
#' @param perm_ctrl A \code{\link{permutation_control}} object. Used fields:
#'   \code{n_perm}, \code{shuffle}, \code{seed}, \code{correction},
#'   \code{rsa_null}.
#' @param metric Optional character vector of metric names to test. By default
#'   every predictor of an \code{rsa_model}, otherwise every numeric column of
#'   the observed performance table.
#' @param batch_size Number of regions whose data are held in memory at once
#'   (all permutations are run on one batch before the next is extracted).
#'   \code{NULL} (default) holds all regions.
#' @param keep_null Keep the permutations x regions x metrics array of null
#'   values in the result.
#' @param ... Passed to \code{run_regional()} for the observed pass.
#'
#' @return An object of class \code{"regional_permutation_result"}: a list with
#'   \code{table} (one row per region and metric: \code{roinum}, \code{metric},
#'   \code{observed}, \code{p}, \code{p_adj}, \code{p_fwer}, \code{n_perm}),
#'   \code{vol_results} (per metric, maps of \code{p}, \code{p_adj} and
#'   \code{p_fwer}; untested voxels are \code{NA}), \code{observed},
#'   \code{perm_ctrl}, \code{metrics}, \code{p_floor}, \code{null_hypothesis},
#'   \code{engine} (\code{"aggregate_fast"} or \code{"per_region"}), and
#'   \code{null} when \code{keep_null = TRUE}.
#'
#' @examples
#' \dontrun{
#' ds <- gen_sample_dataset(c(6, 6, 6), nobs = 48, blocks = 4)
#' regions <- neuroim2::NeuroVol(sample(1:3, 216, replace = TRUE),
#'                               neuroim2::space(ds$dataset$mask))
#' cv <- blocked_cross_validation(ds$design$block_var)
#' ms <- mvpa_model(load_model("corclass"), ds$dataset, ds$design,
#'                  "classification", crossval = cv)
#' res <- run_permutation_regional(ms, regions,
#'                                 perm_ctrl = permutation_control(n_perm = 200, seed = 1))
#' res$table
#' }
#' @seealso \code{\link{run_permutation_searchlight}}, \code{\link{permutation_control}}
#' @export
run_permutation_regional <- function(model_spec,
                                     region_mask,
                                     observed   = NULL,
                                     perm_ctrl  = permutation_control(),
                                     metric     = NULL,
                                     batch_size = NULL,
                                     keep_null  = FALSE,
                                     ...) {
  assertthat::assert_that(inherits(perm_ctrl, "permutation_control"),
                          msg = "perm_ctrl must be a permutation_control object")
  rsa_null <- .resolve_rsa_permutation_null(model_spec, perm_ctrl)
  # Same refusal-before-fitting contract as the searchlight.
  invisible(.permuted_model_spec(model_spec, model_spec$design))

  if (is.null(observed)) {
    futile.logger::flog.info("Running observed regional analysis ...")
    observed <- run_regional(model_spec, region_mask, ...)
  }
  assertthat::assert_that(inherits(observed, "regional_mvpa_result"),
                          msg = "observed must be a regional_mvpa_result")

  metrics <- .resolve_regional_permutation_metrics(model_spec, observed, metric)
  obs_tab <- as.data.frame(observed$performance_table)

  prepped <- prep_regional(model_spec, region_mask)
  rois    <- prepped$region_set
  n_roi   <- length(rois)
  B       <- perm_ctrl$n_perm

  # Draw every permuted design once: each batch, and each region within it,
  # sees the same B relabellings (needed for the family-wise correction).
  futile.logger::flog.info("Drawing %d permutations (%s) ...", B, perm_ctrl$shuffle)
  perm_specs <- lapply(seq_len(B), function(i) {
    perm_seed <- if (!is.null(perm_ctrl$seed)) perm_ctrl$seed + i else NULL
    .permuted_model_spec(model_spec, permute_labels(model_spec$design,
                                                    method = perm_ctrl$shuffle,
                                                    seed = perm_seed))
  })

  null <- array(NA_real_, c(B, n_roi, length(metrics)),
                dimnames = list(NULL, as.character(rois), metrics))

  # Exact prepared engine (correlation / naive Bayes classifiers with fixed
  # folds): region class means are sums over region membership, so the data
  # and folds are prepared once and each permutation only rescores labels.
  fast <- .regional_permutation_engine(model_spec, prepped, metrics, list(...))
  if (!is.null(fast)) {
    for (i in seq_len(B)) {
      vals <- fast$score(perm_specs[[i]])
      if (is.null(vals)) { fast <- NULL; null[] <- NA_real_; break }
      null[i, rownames(vals), ] <- vals
    }
  }
  batch_size <- if (is.null(batch_size)) n_roi else max(1L, as.integer(batch_size))
  batches <- split(seq_len(n_roi), ceiling(seq_len(n_roi) / batch_size))

  for (bi in if (is.null(fast)) seq_along(batches) else integer()) {
    pos <- batches[[bi]]
    frame <- .prepare_batch_frame(
      mod_spec = model_spec, dset = model_spec$dataset,
      vox_batch = prepped$vox_iter[pos], batch_positions = pos, batch_rnums = rois[pos],
      analysis_type = "regional", use_shard_backend = FALSE, min_voxels_required = 2L
    )$frame
    if (!nrow(frame)) next
    futile.logger::flog.info("Regions batch %d / %d: %d regions x %d permutations ...",
                             bi, length(batches), nrow(frame), B)
    for (i in seq_len(B)) {
      res <- tryCatch(
        run_future(perm_specs[[i]], frame, processor = NULL, verbose = FALSE,
                   analysis_type = "regional"),
        error = function(e) {
          futile.logger::flog.warn("Permutation %d failed: %s", i, conditionMessage(e))
          NULL
        })
      if (is.null(res) || !nrow(res)) next
      vals <- .extract_perf_matrix(res, metrics, strict = TRUE)
      if (!nrow(vals)) next
      null[i, rownames(vals), ] <- vals
    }
  }

  # Score each metric: per-region exact p, BH-FDR and minP FWER across regions.
  rows <- lapply(metrics, function(m) {
    obs <- obs_tab[[m]][match(rois, obs_tab$roinum)]
    nm  <- null[, , m, drop = FALSE]
    dim(nm) <- c(B, n_roi)
    sc  <- .score_regional_null(obs, nm)
    p_adj <- if (identical(perm_ctrl$correction, "fdr")) stats::p.adjust(sc$p, method = "BH") else sc$p
    data.frame(roinum = rois, metric = m, observed = obs, p = sc$p, p_adj = p_adj,
               p_fwer = sc$p_fwer, n_perm = sc$n, stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL

  vols <- lapply(stats::setNames(metrics, metrics), function(m) {
    t <- tab[tab$metric == m, ]
    keep <- !is.na(t$p)
    lapply(c(p = "p", p_adj = "p_adj", p_fwer = "p_fwer"), function(col) {
      .map_region_values(region_mask, t$roinum[keep], t[[col]][keep], fill = NA_real_)
    })
  })

  null_hypothesis <- if (identical(rsa_null, "joint")) {
    "Joint no association with any design predictor, including nuisance predictors."
  } else if (identical(rsa_null, "individual")) {
    "No marginal association with each predictor under item exchangeability."
  } else {
    sprintf("Exchangeable labels (%s shuffle): no relation between labels and each region's patterns.",
            perm_ctrl$shuffle)
  }

  structure(list(
    table = tab, vol_results = vols, observed = observed, perm_ctrl = perm_ctrl,
    metrics = metrics, p_floor = 1 / (B + 1), null_hypothesis = null_hypothesis,
    engine = if (!is.null(fast)) fast$engine else "per_region",
    null = if (isTRUE(keep_null)) null else NULL
  ), class = "regional_permutation_result")
}

#' Exact per-region permutation p-values and the minP family-wise correction
#'
#' @param obs Observed values, one per region (NA for regions without one).
#' @param null A B x R matrix of null values (NA where a permutation produced
#'   no value for a region).
#' @return list(p, p_fwer, n) per region.
#' @keywords internal
#' @noRd
.score_regional_null <- function(obs, null) {
  B <- nrow(null); R <- ncol(null)
  n <- colSums(!is.na(null))
  ge <- vapply(seq_len(R), function(r) sum(null[, r] >= obs[r], na.rm = TRUE), numeric(1))
  p <- (1 + ge) / (1 + n)
  p[is.na(obs) | n == 0L] <- NA_real_

  # Null p-value of each permutation within its region's null (counting the
  # observed value with it, the standard resampling convention), then the
  # minimum across regions per permutation.
  null_p <- matrix(NA_real_, B, R)
  for (r in seq_len(R)) {
    v <- null[, r]
    ok <- !is.na(v)
    if (!any(ok)) next
    pool <- c(v[ok], obs[r])
    pool <- pool[!is.na(pool)]
    s <- sort(pool)
    # count of pool values >= v: length - (number strictly below v)
    below <- findInterval(v[ok], s, left.open = TRUE)
    null_p[ok, r] <- (length(s) - below) / (n[r] + 1)
  }
  min_null <- suppressWarnings(apply(null_p, 1, min, na.rm = TRUE))
  min_null <- min_null[is.finite(min_null)]
  p_fwer <- vapply(p, function(pr) if (is.na(pr)) NA_real_ else (1 + sum(min_null <= pr)) / (1 + length(min_null)), numeric(1))
  list(p = p, p_fwer = p_fwer, n = as.integer(n))
}

#' Which metrics a regional permutation run tests
#' @keywords internal
#' @noRd
.resolve_regional_permutation_metrics <- function(model_spec, observed, metric) {
  tab <- as.data.frame(observed$performance_table)
  available <- setdiff(names(tab)[vapply(tab, is.numeric, logical(1))], "roinum")
  metrics <- if (!is.null(metric)) {
    if (!is.character(metric) || !length(metric) || anyNA(metric)) {
      stop("`metric` must be NULL or a character vector of metric names.", call. = FALSE)
    }
    unique(metric)
  } else if (inherits(model_spec, "rsa_model")) {
    model_spec$design$model_predictors %||% names(model_spec$design$model_mat)
  } else {
    available
  }
  missing <- setdiff(metrics, available)
  if (length(missing)) {
    stop(sprintf("metric '%s' not found in the observed regional result; available: %s.",
                 paste(missing, collapse = "', '"), paste(available, collapse = ", ")), call. = FALSE)
  }
  if (!length(metrics)) stop("The observed regional result has no numeric metric to test.", call. = FALSE)
  metrics
}

#' @export
print.regional_permutation_result <- function(x, ...) {
  cat("Regional permutation result\n")
  cat(sprintf("  %d regions x %d metric(s), %d permutations (%s); p floor %.4g\n",
              length(unique(x$table$roinum)), length(x$metrics), x$perm_ctrl$n_perm,
              x$perm_ctrl$shuffle, x$p_floor))
  cat("  null:", x$null_hypothesis, "\n")
  for (m in x$metrics) {
    t <- x$table[x$table$metric == m, ]
    cat(sprintf("  %s: %d region(s) with p_adj < .05 (%s), %d with p_fwer < .05\n", m,
                sum(t$p_adj < 0.05, na.rm = TRUE), x$perm_ctrl$correction, sum(t$p_fwer < 0.05, na.rm = TRUE)))
  }
  invisible(x)
}

#' The exact aggregate engine over regions, when it applies
#'
#' Mirrors .permutation_engine() for the searchlight: same eligibility (the
#' engine the model resolves to, fixed-fold cross-validation), with the
#' region membership matrix in place of the spheres. Returns NULL when the
#' general per-region path must run instead.
#' @keywords internal
#' @noRd
.regional_permutation_engine <- function(model_spec, prepped, metrics, dots) {
  if (!identical(dots$engine %||% "auto", "auto")) return(NULL)
  engine <- tryCatch(.resolve_searchlight_engine(model_spec, "standard", "auto"), error = function(e) "legacy")
  if (!identical(engine, "aggregate_fast")) return(NULL)
  if (!(class(model_spec$crossval)[1L] %in% c("blocked_cross_validation", "custom_cross_validation"))) return(NULL)
  nb <- .aggregate_region_neighbourhoods(model_spec$dataset, prepped$vox_iter, prepped$region_set)
  prep <- tryCatch(.aggregate_prepare(model_spec, radius = NA_real_, neighbourhoods = nb),
                   rmvpa_engine_ineligible = function(e) NULL)
  if (is.null(prep)) return(NULL)
  futile.logger::flog.info("Regional permutations use the exact 'aggregate_fast' engine (prepared once).")
  list(engine = engine, score = function(perm_spec) {
    perf <- tryCatch(.aggregate_score(prep, y_train(perm_spec))$perf, rmvpa_engine_ineligible = function(e) NULL)
    if (is.null(perf) || !all(metrics %in% colnames(perf))) return(NULL)
    out <- perf[, metrics, drop = FALSE]
    rownames(out) <- as.character(prepped$region_set)
    out
  })
}
