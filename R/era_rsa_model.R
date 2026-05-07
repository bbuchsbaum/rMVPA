#' ERA-RSA: Encoding-Retrieval Similarity and ER Geometry
#'
#' Combines first-order encoding-retrieval similarity (ERA) with second-order
#' RSA between encoding and retrieval representational geometries. Works with
#' \code{run_regional()} and \code{run_searchlight()} using the standard rMVPA
#' iterators.
#'
#' Key outputs per ROI/searchlight sphere include:
#' - First-order ERA: top-1 accuracy, diagonal mean, diagonal-minus-off.
#' - Second-order geometry: correlation between encoding and retrieval RDMs.
#' - Optional confound-aware metrics and diagnostics (block/lag/run).
#'
#' \code{era_rsa_model()} itself returns a model specification. The scalar
#' outputs below are emitted when that specification is evaluated with
#' \code{\link{run_regional}()} or \code{\link{run_searchlight}()}: regional
#' analyses place them in \code{result$performance_table}; searchlight analyses
#' expose them as metric maps listed in \code{result$metrics}.
#'
#' @section Metrics:
#' For each ROI / searchlight center, \code{era_rsa_model} emits a set of
#' scalar metrics that are turned into spatial maps by \code{run_regional()}
#' and \code{run_searchlight()}:
#' \describe{
#'   \item{n_items}{
#'     Number of unique item keys contributing to this ROI/sphere
#'     (i.e., length of the common encoding-retrieval item set).
#'   }
#'   \item{era_top1_acc}{
#'     Top-1 encoding\eqn{\rightarrow}retrieval accuracy at the item level:
#'     fraction of retrieval trials whose most similar encoding pattern
#'     (over items) has the same \code{key_var}.
#'   }
#'   \item{era_diag_mean}{
#'     Mean encoding-retrieval similarity for matching items
#'     (mean of the diagonal of the encoding\eqn{\times}retrieval similarity
#'     matrix).
#'   }
#'   \item{era_diag_minus_off}{
#'     Diagonal-minus-off-diagonal ERA contrast: \code{era_diag_mean} minus
#'     the mean similarity to all non-matching items, capturing how much
#'     same-item pairs stand out from other pairs.
#'   }
#'   \item{geom_cor}{
#'     Correlation between encoding and retrieval representational geometries:
#'     correlation (Pearson or Spearman, per \code{rsa_simfun}) between the
#'     vectorised lower triangles of the encoding and retrieval RDMs.
#'   }
#'   \item{geom_cor_partial}{
#'     General nuisance-partial ER geometry correlation. This residualizes both
#'     encoding and retrieval geometry vectors against the confounds selected
#'     by \code{partial_against}, then correlates the residuals. Run confounds
#'     can be supplied either as \code{confound_rdms$run_enc}/\code{run_ret} or
#'     derived from \code{item_run_enc}/\code{item_run_ret}.
#'   }
#'   \item{era_diag_minus_off_same_block}{
#'     Block-limited ERA contrast when \code{item_block} is supplied:
#'     diagonal ERA minus the mean similarity to other items in the same
#'     block (e.g., run/condition), averaged over items.
#'   }
#'   \item{era_diag_minus_off_diff_block}{
#'     Cross-block ERA contrast when \code{item_block} is supplied:
#'     diagonal ERA minus the mean similarity to items in different blocks.
#'   }
#'   \item{era_lag_cor}{
#'     Lag-ERA correlation when \code{item_lag} is supplied: Spearman
#'     correlation between diagonal ERA values and the per-item lag
#'     (e.g., retrieval minus encoding onset), using complete cases only.
#'   }
#'   \item{geom_cor_run_partial}{
#'     Run-partial ER geometry correlation, when run-level confounds are
#'     supplied via \code{confound_rdms$run_enc}/\code{run_ret} or derivable
#'     from \code{item_run_enc}/\code{item_run_ret}. Computed as the
#'     correlation between encoding and retrieval RDMs after regressing out
#'     those run RDMs.
#'   }
#'   \item{geom_cor_xrun}{
#'     Cross-run-only ER geometry correlation, when \code{item_run_enc} and
#'     \code{item_run_ret} are supplied: correlation between encoding and
#'     retrieval RDMs restricted to item pairs that differ in both encoding
#'     and retrieval run.
#'   }
#'   \item{beta_*}{
#'     When \code{confound_rdms} are provided, additional \code{beta_<name>}
#'     terms give the regression coefficients from a linear model predicting
#'     retrieval geometry from encoding geometry and the confound RDMs
#'     (one coefficient per nuisance/model RDM).
#'   }
#'   \item{sp_*}{
#'     If \code{run_lm_semipartial()} is available, \code{sp_<name>} terms
#'     provide semi-partial R\eqn{^2}-like diagnostics for each confound RDM,
#'     quantifying unique variance explained in retrieval geometry.
#'   }
#' }
#'
#' @param dataset An mvpa_dataset with train_data (encoding) and test_data (retrieval).
#' @param design  An mvpa_design describing trial structure with train/test designs.
#' @param key_var Column name or formula giving the item key that links encoding
#'   and retrieval trials (e.g., ~ ImageID).
#' @param phase_var Column name or formula giving phase labels (must include
#'   encoding and retrieval levels if using a single-phase dataset; for the
#'   default two-dataset usage, this is still parsed for consistency but not
#'   required for operations).
#' @param encoding_level Level of phase_var to treat as encoding (default: first level).
#' @param retrieval_level Level of phase_var to treat as retrieval (default: second level).
#' @param distfun A distfun used to compute within-phase RDMs (e.g., cordist("pearson")).
#' @param rsa_simfun Character: similarity for comparing RDMs, one of "pearson" or "spearman".
#' @param confound_rdms Optional named list of KxK matrices or "dist" objects
#'   encoding item-level nuisance/model RDMs (e.g., block, run, time). Rows
#'   and columns should correspond to item keys (levels of \code{key_var}).
#'   When \code{run_enc} and \code{run_ret} entries are present they are used
#'   to compute \code{geom_cor_run_partial}, the ER geometry correlation after
#'   regressing out these run confounds. The same list also supplies candidate
#'   nuisance RDMs for \code{geom_cor_partial}, selected by
#'   \code{partial_against}.
#' @param partial_against Character vector selecting nuisance groups or exact
#'   \code{confound_rdms} names used for the general \code{geom_cor_partial}
#'   metric. Recognized groups are \code{"run"}, \code{"time"},
#'   \code{"block"}, \code{"category"}, \code{"global"}, and \code{"all"}.
#'   Exact confound names such as \code{"time_enc"} are also allowed. The
#'   default \code{"run"} preserves the legacy run-partial interpretation while
#'   exposing the result under the more general \code{geom_cor_partial} name.
#'   If \code{global_nuisance} is enabled and \code{partial_against} is not
#'   supplied explicitly, the effective default is \code{c("run", "global")}.
#' @param include_diag Logical retained for API compatibility. Diagonal ERA
#'   metrics are always retained, and the off-diagonal mean used by
#'   \code{era_diag_minus_off} always excludes matching-item diagonal entries.
#' @param item_block Optional factor of per-item blocks, aligned to item keys.
#'   Typically derived by aggregating a trial-level block/run variable in
#'   \code{design$train_design} to the item level (e.g., modal block per item).
#'   Used to compute block-limited ERA contrasts
#'   (\code{era_diag_minus_off_same_block} / \code{era_diag_minus_off_diff_block}).
#' @param item_lag Optional numeric per-item retrieval-minus-encoding lag,
#'   aligned to item keys. Often derived from trial onsets (e.g., mean
#'   retrieval onset minus mean encoding onset per item). Used to compute
#'   \code{era_lag_cor}, the correlation between ERA diagonal and lag.
#' @param item_run_enc Optional factor of per-item encoding runs, aligned to
#'   item keys (e.g., modal run for encoding trials of each item). Combined
#'   with \code{item_run_ret} to compute \code{geom_cor_xrun}, ER geometry
#'   restricted to item pairs differing in both encoding and retrieval run.
#' @param item_run_ret Optional factor of per-item retrieval runs, aligned to
#'   item keys. See \code{item_run_enc} for how it is used.
#' @param global_nuisance Logical or pre-supplied list controlling whole-mask
#'   global similarity nuisance. \code{FALSE} (default) disables it. \code{TRUE}
#'   computes K x K item-level RDMs (\code{global_enc}, \code{global_ret}) over
#'   the full \code{dataset$mask} once at construction time and adds them to
#'   \code{confound_rdms}. They are then picked up by \code{geom_cor_partial}
#'   when \code{partial_against} includes \code{"global"} (or \code{"all"}).
#'   A pre-computed list with elements \code{D_enc}/\code{enc} and
#'   \code{D_ret}/\code{ret} can be supplied directly for non-standard dataset
#'   backends. Caveat: each ROI/sphere is part of the global mask, so for large
#'   regional ROIs covering most of the mask the residualization partially
#'   removes the local signal.
#' @param require_run_metadata Logical; if \code{TRUE}, missing item-level run
#'   or block metadata becomes an error rather than a warning. Use this when a
#'   downstream analysis depends on \code{era_diag_minus_off_same_block},
#'   \code{era_diag_minus_off_diff_block}, or \code{geom_cor_xrun} — the
#'   constructor will refuse to silently produce schemas where those metrics
#'   are guaranteed to be \code{NA}. Default \code{FALSE} (warn only).
#' @param ... Additional fields stored on the model spec.
#'
#' @section Trial-level vs. item-level metadata:
#' \code{block_var} on \code{\link{mvpa_design}()} is \emph{trial-level}
#' metadata (one entry per row of the design table) and is used by
#' cross-validation, not by the ERA item-level metrics. The block/run-aware
#' metrics (\code{era_diag_minus_off_same_block},
#' \code{era_diag_minus_off_diff_block}, \code{geom_cor_xrun}) are computed
#' from \emph{item-level} vectors indexed by levels of \code{key_var}:
#' \code{item_block}, \code{item_run_enc}, and \code{item_run_ret}. These must
#' be supplied directly here; passing \code{block_var = ~run} to
#' \code{mvpa_design()} alone will not enable them, and the model will warn
#' that those metrics will be \code{NA}. See \code{\link{era_rsa_design}()}
#' for a helper that builds these vectors from the design table.
#'
#' For external train/test phases (encoding vs. retrieval), run labels often
#' collide across phases (both phases may have runs \code{1, 2, 3} that
#' correspond to different scans). When \code{item_run_enc} and
#' \code{item_run_ret} share atomic values that are not phase-scoped, supply
#' phase-prefixed labels such as \code{enc_1} / \code{ret_1} so cross-phase
#' equality tests are not spuriously satisfied.
#'
#' @return A model spec of class \code{"era_rsa_model"} compatible with
#'   \code{\link{run_regional}()} and \code{\link{run_searchlight}()}. When fit,
#'   the spec emits the scalar metrics documented in \strong{Metrics}.
#' @examples
#' \dontrun{
#'   # See vignette for complete ERA-RSA workflow
#' }
#' @export
era_rsa_model <- function(dataset,
                          design,
                          key_var,
                          phase_var,
                          encoding_level  = NULL,
                          retrieval_level = NULL,
                          distfun         = cordist(method = "pearson"),
                          rsa_simfun      = c("pearson", "spearman"),
                          confound_rdms   = NULL,
                          partial_against = "run",
                          include_diag    = TRUE,
                          item_block      = NULL,
                          item_lag        = NULL,
                          item_run_enc    = NULL,
                          item_run_ret    = NULL,
                          global_nuisance = FALSE,
                          require_run_metadata = FALSE,
                          ...) {

  rsa_simfun <- match.arg(rsa_simfun)

  stopifnot(inherits(dataset, "mvpa_dataset"))
  stopifnot(inherits(design,  "mvpa_design"))

  .era_check_item_metadata(
    where        = "era_rsa_model",
    item_block   = item_block,
    item_run_enc = item_run_enc,
    item_run_ret = item_run_ret,
    strict       = isTRUE(require_run_metadata),
    metric_for_block = c("era_diag_minus_off_same_block",
                         "era_diag_minus_off_diff_block"),
    metric_for_run   = "geom_cor_xrun"
  )

  # Normalize distfun
  if (is.character(distfun)) {
    distfun <- create_dist(distfun)
  }
  stopifnot(inherits(distfun, "distfun"))

  # Parse design variables from train_design for stable levels/order
  phase_vec <- parse_variable(phase_var, design$train_design)
  key_vec   <- parse_variable(key_var,   design$train_design)
  phase_fac <- factor(phase_vec)
  key_fac   <- factor(key_vec)
  phase_lev <- levels(phase_fac)

  if (is.null(encoding_level))  encoding_level  <- phase_lev[1L]
  if (is.null(retrieval_level)) retrieval_level <- phase_lev[2L]

  if (!is.null(global_nuisance) && !isFALSE(global_nuisance) && missing(partial_against)) {
    partial_against <- c("run", "global")
  }

  # Optional whole-mask global similarity nuisance: compute K x K RDMs once
  # over the full dataset mask and merge them into confound_rdms so they can
  # be picked up by `partial_against = "global"` (and the all/run/etc groups).
  global_rdms <- .era_resolve_global_nuisance(
    global_nuisance, dataset, design, key_var, distfun
  )
  if (!is.null(global_rdms)) {
    if (is.null(confound_rdms)) confound_rdms <- list()
    if (!is.null(global_rdms$D_enc) && is.null(confound_rdms$global_enc)) {
      confound_rdms$global_enc <- global_rdms$D_enc
    }
    if (!is.null(global_rdms$D_ret) && is.null(confound_rdms$global_ret)) {
      confound_rdms$global_ret <- global_rdms$D_ret
    }
  }

  # Normalize confound RDMs into matrices with names if provided
  if (!is.null(confound_rdms)) {
    stopifnot(is.list(confound_rdms), !is.null(names(confound_rdms)))
    items_all <- levels(key_fac)
    confound_rdms <- lapply(confound_rdms, function(M) {
      if (inherits(M, "dist")) {
        lab <- attr(M, "Labels")
        if (is.null(lab)) lab <- items_all[seq_len(attr(M, "Size"))]
        M <- as.matrix(M)
        rownames(M) <- colnames(M) <- lab
      } else {
        M <- as.matrix(M)
        if (is.null(rownames(M)) || is.null(colnames(M))) {
          rn <- items_all[seq_len(min(length(items_all), nrow(M)))]
          rownames(M) <- colnames(M) <- rn
        }
      }
      M
    })
  }

  create_model_spec(
    "era_rsa_model",
    dataset  = dataset,
    design   = design,
    # store parsed vectors to maintain levels/order
    key      = key_fac,
    phase    = phase_fac,
    key_var  = key_var,
    phase_var= substitute(phase_var),
    encoding_level  = encoding_level,
    retrieval_level = retrieval_level,
    distfun        = distfun,
    rsa_simfun     = rsa_simfun,
    confound_rdms  = confound_rdms,
    partial_against = partial_against,
    include_diag   = include_diag,
    item_block     = item_block,
    item_lag       = item_lag,
    item_run_enc   = item_run_enc,
    item_run_ret   = item_run_ret,
    compute_performance = TRUE,
    return_predictions  = FALSE,
    ...
  )
}


#' Output schema for era_rsa_model
#'
#' Declares all scalar metrics emitted by \code{fit_roi.era_rsa_model} so that
#' \code{combine_schema_standard} can pre-allocate output maps with consistent
#' length across ROIs.
#'
#' @param model An era_rsa_model object.
#' @return A named character vector mapping each emitted metric name to
#'   \code{"scalar"}. Base metrics are \code{n_items}, \code{era_top1_acc},
#'   \code{era_diag_mean}, \code{era_diag_minus_off}, \code{geom_cor},
#'   \code{era_diag_minus_off_same_block},
#'   \code{era_diag_minus_off_diff_block}, \code{era_lag_cor},
#'   \code{geom_cor_partial}, \code{geom_cor_run_partial}, and
#'   \code{geom_cor_xrun}. If
#'   \code{confound_rdms} is supplied, the schema also includes
#'   \code{beta_enc_geom}, one \code{beta_<name>} per confound RDM,
#'   \code{sp_enc_geom}, and one \code{sp_<name>} per confound RDM.
#' @keywords internal
#' @export
output_schema.era_rsa_model <- function(model) {
  base_names <- c("n_items", "era_top1_acc", "era_diag_mean", "era_diag_minus_off",
                  "geom_cor", "era_diag_minus_off_same_block",
                  "era_diag_minus_off_diff_block", "era_lag_cor",
                  "geom_cor_partial", "geom_cor_run_partial", "geom_cor_xrun")

  if (!is.null(model$confound_rdms)) {
    conf_names <- names(model$confound_rdms)
    beta_names <- paste0("beta_", c("enc_geom", conf_names))
    sp_names   <- paste0("sp_",   c("enc_geom", conf_names))
    base_names <- c(base_names, beta_names, sp_names)
  }

  setNames(rep("scalar", length(base_names)), base_names)
}


#' fit_roi method for era_rsa_model
#'
#' Computes first-order ERA metrics and second-order encoding-retrieval geometry
#' for a single ROI/searchlight sphere.
#'
#' @param model An era_rsa_model object.
#' @param roi_data List with train_data, test_data, indices, train_roi, test_roi.
#' @param context List with design, cv_spec, id, center_global_id.
#' @param ... Additional arguments (unused).
#' @return An \code{roi_result} object.
#' @keywords internal
#' @export
fit_roi.era_rsa_model <- function(model, roi_data, context, ...) {
  Xenc <- roi_data$train_data
  ind  <- roi_data$indices
  id   <- context$id

  # Require external test set for retrieval
  if (is.null(roi_data$test_data)) {
    return(roi_result(
      metrics = NULL, indices = ind, id = id,
      error = TRUE,
      error_message = "era_rsa_model requires dataset$test_data + design$y_test (external retrieval set)"
    ))
  }

  Xret <- roi_data$test_data
  des  <- model$design

  # Guard: non-degenerate ROI
  if (ncol(Xenc) < 1L || ncol(Xret) < 1L || nrow(Xenc) < 2L || nrow(Xret) < 2L) {
    return(roi_result(
      metrics = NULL, indices = ind, id = id,
      error = TRUE,
      error_message = "ERA-RSA: ROI too small (need >=2 trials and >=1 voxel per phase)"
    ))
  }

  # Extract keys for encoding/retrieval from design
  key_tr_full <- parse_variable(model$key_var, des$train_design)
  key_te_full <- parse_variable(model$key_var, des$test_design)
  key_tr_full <- factor(key_tr_full, levels = levels(model$key))
  key_te_full <- factor(key_te_full, levels = levels(model$key))

  # Build item-level prototypes via group means
  E_full <- group_means(Xenc, margin = 1, group = key_tr_full)
  R_full <- group_means(Xret, margin = 1, group = key_te_full)

  keys_enc    <- rownames(E_full)
  keys_ret    <- rownames(R_full)
  common_keys <- sort(intersect(keys_enc, keys_ret))
  K           <- length(common_keys)

  if (K < 2L) {
    return(roi_result(
      metrics = NULL, indices = ind, id = id,
      error = TRUE,
      error_message = sprintf("ERA-RSA: need >=2 items with both encoding and retrieval (found %d)", K)
    ))
  }

  E <- E_full[common_keys, , drop = FALSE]
  R <- R_full[common_keys, , drop = FALSE]

  # Drop voxels constant across both phases
  sd_E     <- apply(E, 2, stats::sd, na.rm = TRUE)
  sd_R     <- apply(R, 2, stats::sd, na.rm = TRUE)
  keep_vox <- (sd_E > 0) | (sd_R > 0)
  if (!any(keep_vox)) {
    return(roi_result(
      metrics = NULL, indices = ind, id = id,
      error = TRUE,
      error_message = "ERA-RSA: all voxels zero-variance"
    ))
  }
  E <- E[, keep_vox, drop = FALSE]
  R <- R[, keep_vox, drop = FALSE]

  # First-order ERA: similarity and diagnostics
  S             <- suppressWarnings(stats::cor(t(E), t(R), use = "pairwise.complete.obs"))
  diag_sim      <- diag(S)
  era_diag_mean <- mean(diag_sim, na.rm = TRUE)
  off           <- S; diag(off) <- NA_real_
  era_off_mean  <- mean(off, na.rm = TRUE)
  era_diag_minus_off <- era_diag_mean - era_off_mean
  max_enc_idx   <- apply(S, 2, function(col) if (all(is.na(col))) NA_integer_ else which.max(col))
  era_top1_acc  <- mean((!is.na(max_enc_idx)) & (max_enc_idx == seq_len(K)))

  # Second-order ER geometry
  D_enc    <- pairwise_dist(model$distfun, E)
  D_ret    <- pairwise_dist(model$distfun, R)
  dE       <- as.numeric(D_enc[lower.tri(D_enc)])
  dR       <- as.numeric(D_ret[lower.tri(D_ret)])
  geom_cor <- suppressWarnings(stats::cor(dE, dR, method = model$rsa_simfun, use = "complete.obs"))

  # Optional geometry regression with confounds
  beta_vec <- c(); sp_vec <- c()
  if (!is.null(model$confound_rdms) && K >= 3L) {
    mm <- list(enc_geom = dE)
    for (nm in names(model$confound_rdms)) {
      M <- as.matrix(model$confound_rdms[[nm]])
      if (!is.null(rownames(M)) && all(common_keys %in% rownames(M))) {
        M <- M[common_keys, common_keys, drop = FALSE]
      } else {
        M <- M[seq_len(K), seq_len(K), drop = FALSE]
      }
      mm[[nm]] <- as.numeric(M[lower.tri(M)])
    }
    df  <- as.data.frame(mm)
    fit <- try(stats::lm(dR ~ ., data = df), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      cf <- stats::coef(summary(fit))
      if (nrow(cf) > 1L) {
        keep     <- rownames(cf) != "(Intercept)"
        beta_vec <- cf[keep, "Estimate"]
        names(beta_vec) <- paste0("beta_", rownames(cf)[keep])
      }
      if (exists("run_lm_semipartial", mode = "function")) {
        tmp <- list(design = list(model_mat = mm))
        sr  <- try(run_lm_semipartial(dvec = dR, obj = tmp), silent = TRUE)
        if (!inherits(sr, "try-error") && is.numeric(sr)) {
          names(sr) <- paste0("sp_", names(mm))
          sp_vec <- sr
        }
      }
    }
  }

  # Block-limited ERA if item_block available
  era_diag_minus_off_same_block <- NA_real_
  era_diag_minus_off_diff_block <- NA_real_
  if (!is.null(model$item_block)) {
    ib <- model$item_block
    if (!is.null(names(ib))) ib <- ib[match(common_keys, names(ib))] else ib <- ib[seq_along(common_keys)]
    same_vals <- c(); diff_vals <- c()
    for (i in seq_len(K)) {
      same_idx <- setdiff(which(ib == ib[i]), i)
      diff_idx <- which(ib != ib[i])
      if (length(same_idx)) same_vals <- c(same_vals, S[i, same_idx])
      if (length(diff_idx)) diff_vals <- c(diff_vals, S[i, diff_idx])
    }
    if (length(same_vals)) era_diag_minus_off_same_block <- mean(diag_sim, na.rm = TRUE) - mean(same_vals, na.rm = TRUE)
    if (length(diff_vals)) era_diag_minus_off_diff_block <- mean(diag_sim, na.rm = TRUE) - mean(diff_vals, na.rm = TRUE)
  }

  # Lag-specific ERA if item_lag available
  era_lag_cor <- NA_real_
  if (!is.null(model$item_lag)) {
    lag <- model$item_lag
    if (!is.null(names(lag))) lag <- lag[match(common_keys, names(lag))] else lag <- lag[seq_along(common_keys)]
    if (length(lag) == length(diag_sim) && any(!is.na(lag))) {
      era_lag_cor <- suppressWarnings(stats::cor(diag_sim, lag, method = "spearman", use = "complete.obs"))
    }
  }

  # Partial ER geometry and cross-run-only geometry if nuisance/run info available
  all_confounds <- .era_geometry_confound_rdms(
    confound_rdms = model$confound_rdms,
    item_run_enc = model$item_run_enc,
    item_run_ret = model$item_run_ret,
    keys = common_keys
  )
  geom_cor_partial <- .era_partial_geometry_cor(
    dE = dE,
    dR = dR,
    confound_rdms = all_confounds,
    keys = common_keys,
    partial_against = model$partial_against %||% "run",
    method = model$rsa_simfun
  )
  geom_cor_run_partial <- NA_real_
  geom_cor_xrun        <- NA_real_
  geom_cor_run_partial <- .era_partial_geometry_cor(
    dE = dE,
    dR = dR,
    confound_rdms = all_confounds,
    keys = common_keys,
    partial_against = "run",
    method = model$rsa_simfun
  )
  if (!is.null(model$item_run_enc) && !is.null(model$item_run_ret)) {
    ire <- model$item_run_enc
    irr <- model$item_run_ret
    if (!is.null(names(ire))) ire <- ire[match(common_keys, names(ire))] else ire <- ire[seq_along(common_keys)]
    if (!is.null(names(irr))) irr <- irr[match(common_keys, names(irr))] else irr <- irr[seq_along(common_keys)]
    same_enc <- outer(ire, ire, "==")
    same_ret <- outer(irr, irr, "==")
    mask     <- lower.tri(D_enc) & !(same_enc | same_ret)
    if (any(mask)) {
      geom_cor_xrun <- suppressWarnings(stats::cor(D_enc[mask], D_ret[mask], method = model$rsa_simfun, use = "complete.obs"))
    }
  }

  # Assemble base perf vector
  perf <- c(
    n_items                       = K,
    era_top1_acc                  = era_top1_acc,
    era_diag_mean                 = era_diag_mean,
    era_diag_minus_off            = era_diag_minus_off,
    geom_cor                      = geom_cor,
    era_diag_minus_off_same_block = era_diag_minus_off_same_block,
    era_diag_minus_off_diff_block = era_diag_minus_off_diff_block,
    era_lag_cor                   = era_lag_cor,
    geom_cor_partial              = geom_cor_partial,
    geom_cor_run_partial          = geom_cor_run_partial,
    geom_cor_xrun                 = geom_cor_xrun
  )

  # Always include confound-derived metrics if schema declares them
  if (!is.null(model$confound_rdms)) {
    conf_names          <- names(model$confound_rdms)
    expected_beta_names <- paste0("beta_", c("enc_geom", conf_names))
    expected_sp_names   <- paste0("sp_",   c("enc_geom", conf_names))

    full_beta <- setNames(rep(NA_real_, length(expected_beta_names)), expected_beta_names)
    full_sp   <- setNames(rep(NA_real_, length(expected_sp_names)),   expected_sp_names)

    if (length(beta_vec)) {
      matched <- intersect(names(beta_vec), names(full_beta))
      full_beta[matched] <- beta_vec[matched]
    }
    if (length(sp_vec)) {
      matched <- intersect(names(sp_vec), names(full_sp))
      full_sp[matched] <- sp_vec[matched]
    }

    perf <- c(perf, full_beta, full_sp)
  }

roi_result(metrics = perf, indices = ind, id = id)
}


#' @noRd
.era_geometry_confound_rdms <- function(confound_rdms = NULL,
                                        item_run_enc = NULL,
                                        item_run_ret = NULL,
                                        keys) {
  out <- confound_rdms %||% list()
  if (length(out) && is.null(names(out))) {
    names(out) <- paste0("confound_", seq_along(out))
  }

  run_rdms <- .era_run_confound_rdms(item_run_enc, item_run_ret, keys)
  for (nm in names(run_rdms)) {
    if (is.null(out[[nm]])) {
      out[[nm]] <- run_rdms[[nm]]
    }
  }
  out
}

#' @noRd
.era_run_confound_rdms <- function(item_run_enc = NULL, item_run_ret = NULL, keys) {
  if (is.null(item_run_enc) || is.null(item_run_ret)) {
    return(list())
  }

  align <- function(x) {
    if (!is.null(names(x))) {
      x[match(keys, names(x))]
    } else {
      x[seq_along(keys)]
    }
  }
  ire <- align(item_run_enc)
  irr <- align(item_run_ret)

  Renc <- outer(ire, ire, FUN = function(a, b) as.numeric(a == b))
  Rret <- outer(irr, irr, FUN = function(a, b) as.numeric(a == b))
  rownames(Renc) <- colnames(Renc) <- keys
  rownames(Rret) <- colnames(Rret) <- keys
  list(run_enc = Renc, run_ret = Rret)
}

#' @noRd
.era_partial_geometry_cor <- function(dE,
                                      dR,
                                      confound_rdms,
                                      keys,
                                      partial_against = "run",
                                      method = c("pearson", "spearman")) {
  method <- match.arg(method)
  if (is.null(confound_rdms) || !length(confound_rdms) || length(partial_against) == 0L) {
    return(NA_real_)
  }

  nms <- names(confound_rdms)
  selected <- .era_select_confound_names(nms, partial_against)
  if (!length(selected)) {
    return(NA_real_)
  }

  conf <- lapply(confound_rdms[selected], .era_vectorize_geometry_confound, keys = keys)
  lens_ok <- vapply(conf, length, integer(1L)) == length(dE)
  conf <- conf[lens_ok]
  if (!length(conf)) {
    return(NA_real_)
  }

  dat <- data.frame(.dE = as.numeric(dE), .dR = as.numeric(dR), as.data.frame(conf, optional = TRUE))
  names(dat) <- c(".dE", ".dR", make.names(names(conf), unique = TRUE))
  keep <- stats::complete.cases(dat)
  keep <- keep & apply(dat, 1L, function(row) all(is.finite(row)))
  dat <- dat[keep, , drop = FALSE]
  if (nrow(dat) < 3L || stats::sd(dat$.dE) == 0 || stats::sd(dat$.dR) == 0) {
    return(NA_real_)
  }

  conf_names <- setdiff(names(dat), c(".dE", ".dR"))
  nonconstant <- vapply(dat[conf_names], function(x) stats::sd(x) > 0, logical(1L))
  conf_names <- conf_names[nonconstant]
  if (!length(conf_names)) {
    return(NA_real_)
  }

  rhs <- paste(conf_names, collapse = " + ")
  fit_e <- try(stats::lm(stats::as.formula(paste(".dE ~", rhs)), data = dat), silent = TRUE)
  fit_r <- try(stats::lm(stats::as.formula(paste(".dR ~", rhs)), data = dat), silent = TRUE)
  if (inherits(fit_e, "try-error") || inherits(fit_r, "try-error")) {
    return(NA_real_)
  }

  e_res <- stats::resid(fit_e)
  r_res <- stats::resid(fit_r)
  suppressWarnings(stats::cor(e_res, r_res, method = method, use = "complete.obs"))
}

#' @noRd
.era_select_confound_names <- function(nms, partial_against) {
  if (is.null(nms) || !length(nms)) {
    return(character())
  }
  partial_against <- unique(as.character(partial_against))
  if ("all" %in% partial_against) {
    return(nms)
  }

  selected <- partial_against[partial_against %in% nms]
  group_patterns <- list(
    run = "run",
    time = "time|temporal|lag",
    block = "block",
    category = "category|cat",
    global = "^global"
  )
  for (grp in intersect(names(group_patterns), partial_against)) {
    selected <- c(selected, nms[grepl(group_patterns[[grp]], nms, ignore.case = TRUE)])
  }
  unique(selected)
}

#' @noRd
.era_vectorize_geometry_confound <- function(x, keys) {
  K <- length(keys)
  lower_n <- K * (K - 1L) / 2L
  if (inherits(x, "dist") || is.matrix(x) || is.data.frame(x)) {
    M <- as.matrix(x)
    if (!is.null(rownames(M)) && !is.null(colnames(M)) &&
        all(keys %in% rownames(M)) && all(keys %in% colnames(M))) {
      M <- M[keys, keys, drop = FALSE]
    } else {
      M <- M[seq_len(K), seq_len(K), drop = FALSE]
    }
    return(as.numeric(M[lower.tri(M)]))
  }

  v <- as.numeric(x)
  if (length(v) == K * K) {
    M <- matrix(v, nrow = K)
    as.numeric(M[lower.tri(M)])
  } else if (length(v) == lower_n) {
    v
  } else {
    v
  }
}


#' Validate item-level metadata for ERA models
#'
#' Used by \code{era_rsa_model()} and \code{era_partition_model()} to warn (or
#' error, when \code{strict=TRUE}) about missing item-level vectors that would
#' otherwise produce silently-NA metrics or unused nuisance regressors. Also
#' detects an easy-to-miss namespace collision where encoding and retrieval
#' run labels share atomic values without being phase-scoped.
#'
#' @param where Character; calling function name used in messages.
#' @param item_block,item_run_enc,item_run_ret Optional item-level vectors.
#' @param strict Logical; if \code{TRUE}, missing metadata becomes an error.
#' @param metric_for_block Character vector of metric names that depend on
#'   \code{item_block}.
#' @param metric_for_run Character vector of metric names that depend on
#'   \code{item_run_enc} / \code{item_run_ret}.
#' @return Invisibly \code{NULL}; called for side effects.
#' @keywords internal
.era_check_item_metadata <- function(where,
                                     item_block   = NULL,
                                     item_run_enc = NULL,
                                     item_run_ret = NULL,
                                     strict       = FALSE,
                                     metric_for_block = character(),
                                     metric_for_run   = character()) {
  emit <- if (isTRUE(strict)) {
    function(msg) stop(sprintf("%s: %s", where, msg), call. = FALSE)
  } else {
    function(msg) warning(sprintf("%s: %s", where, msg), call. = FALSE)
  }

  if (is.null(item_block) && length(metric_for_block)) {
    emit(sprintf(
      "`item_block` is NULL; %s will be NA. Supply a per-item block vector named by levels of `key_var`, e.g. via `era_rsa_design(..., block_var = ~ run)$item_block`.",
      paste(sprintf("`%s`", metric_for_block), collapse = " and ")
    ))
  }

  run_missing <- is.null(item_run_enc) || is.null(item_run_ret)
  if (run_missing && length(metric_for_run)) {
    emit(sprintf(
      "%s requires both `item_run_enc` and `item_run_ret` (per-item run vectors named by `key_var` levels); %s will be NA.",
      paste(sprintf("`%s`", metric_for_run), collapse = " / "),
      paste(sprintf("`%s`", metric_for_run), collapse = " and ")
    ))
  }

  if (!is.null(item_run_enc) && !is.null(item_run_ret)) {
    enc_vals <- as.character(unique(stats::na.omit(item_run_enc)))
    ret_vals <- as.character(unique(stats::na.omit(item_run_ret)))
    overlap  <- intersect(enc_vals, ret_vals)
    if (length(overlap)) {
      warning(sprintf(
        "%s: `item_run_enc` and `item_run_ret` share label(s) %s. If encoding and retrieval are different scans, use phase-scoped labels (e.g. `enc_1` and `ret_1`) so cross-phase equality is not spuriously satisfied.",
        where, paste(sprintf("'%s'", utils::head(overlap, 6L)), collapse = ", ")
      ), call. = FALSE)
    }
  }

  invisible(NULL)
}
