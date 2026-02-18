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
#'     supplied via \code{confound_rdms$run_enc} and \code{confound_rdms$run_ret}.
#'     Computed as the correlation between encoding and retrieval RDMs
#'     after regressing out those run RDMs.
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
#'   regressing out these run confounds.
#' @param include_diag Logical; if TRUE (default) ERA off-diagonal mean excludes diagonal
#'   by setting it to NA first; diagonal metrics are always retained.
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
#' @param ... Additional fields stored on the model spec.
#'
#' @return A model spec of class "era_rsa_model" compatible with run_regional()/run_searchlight().
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
                          include_diag    = TRUE,
                          item_block      = NULL,
                          item_lag        = NULL,
                          item_run_enc    = NULL,
                          item_run_ret    = NULL,
                          ...) {

  rsa_simfun <- match.arg(rsa_simfun)

  stopifnot(inherits(dataset, "mvpa_dataset"))
  stopifnot(inherits(design,  "mvpa_design"))

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
    key_var  = substitute(key_var),
    phase_var= substitute(phase_var),
    encoding_level  = encoding_level,
    retrieval_level = retrieval_level,
    distfun        = distfun,
    rsa_simfun     = rsa_simfun,
    confound_rdms  = confound_rdms,
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
#' @return A named character vector mapping metric names to \code{"scalar"}.
#' @keywords internal
#' @export
output_schema.era_rsa_model <- function(model) {
  base_names <- c("n_items", "era_top1_acc", "era_diag_mean", "era_diag_minus_off",
                  "geom_cor", "era_diag_minus_off_same_block",
                  "era_diag_minus_off_diff_block", "era_lag_cor",
                  "geom_cor_run_partial", "geom_cor_xrun")

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

  # Run-partial ER geometry and cross-run-only geometry if run info available
  geom_cor_run_partial <- NA_real_
  geom_cor_xrun        <- NA_real_
  if (!is.null(model$confound_rdms$run_enc) && !is.null(model$confound_rdms$run_ret)) {
    Renc    <- as.numeric(as.matrix(model$confound_rdms$run_enc)[common_keys, common_keys][lower.tri(D_enc)])
    Rret    <- as.numeric(as.matrix(model$confound_rdms$run_ret)[common_keys, common_keys][lower.tri(D_ret)])
    conf_df <- data.frame(enc_run = Renc, ret_run = Rret)
    dE_res  <- stats::resid(stats::lm(dE ~ ., data = conf_df))
    dR_res  <- stats::resid(stats::lm(dR ~ ., data = conf_df))
    geom_cor_run_partial <- suppressWarnings(stats::cor(dE_res, dR_res, method = model$rsa_simfun, use = "complete.obs"))
  }
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


