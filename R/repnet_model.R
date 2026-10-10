# Representational Connectivity (ReNA-RC): repnet_model

#' Representational connectivity model (ReNA-RC)
#'
#' For each ROI or searchlight, computes representational connectivity between
#' the ROI RDM and a seed RDM, optionally controlling for confound RDMs (block,
#' lag, behavior, etc.).
#'
#' @param dataset An \code{mvpa_dataset}.
#' @param design An \code{mvpa_design}.
#' @param repnet_des Output of \code{\link{repnet_design}}.
#' @param distfun Distance function for ROI RDM (e.g. \code{cordist(method = "pearson")}).
#' @param simfun Character: similarity metric between ROI and seed RDMs
#'   (\code{"pearson"} or \code{"spearman"}).
#' @param ... Extra fields stored on the model spec.
#'
#' @return A model spec of class \code{"repnet_model"} compatible with
#'   \code{run_regional()} and \code{run_searchlight()}.
#' @examples
#' ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 16, blocks = 4, nlevels = 2)
#' # Seed RDM over trials (e.g. from a seed region or a behavioural model)
#' key_ids <- ds$design$train_design$.rownum
#' seed_rdm <- as.matrix(dist(key_ids))
#' rownames(seed_rdm) <- colnames(seed_rdm) <- as.character(key_ids)
#' rn_des <- repnet_design(ds$design, key_var = ~ .rownum, seed_rdm = seed_rdm)
#' model <- repnet_model(ds$dataset, ds$design, rn_des)
#' @export
repnet_model <- function(dataset,
                         design,
                         repnet_des,
                         distfun = cordist(method = "pearson"),
                         simfun = c("pearson", "spearman"),
                         ...) {

  assertthat::assert_that(inherits(dataset, "mvpa_dataset"))
  assertthat::assert_that(inherits(design,  "mvpa_design"))

  simfun <- match.arg(simfun)

  if (is.character(distfun)) {
    distfun <- create_dist(distfun)
  }
  assertthat::assert_that(
    inherits(distfun, "distfun"),
    msg = "distfun must be a 'distfun' or a shortcut understood by create_dist()."
  )

  create_model_spec(
    "repnet_model",
    dataset       = dataset,
    design        = design,
    key           = repnet_des$key,
    seed_rdm      = repnet_des$seed_rdm,
    confound_rdms = repnet_des$confound_rdms,
    distfun       = distfun,
    simfun        = simfun,
    # We compute per-ROI performance directly in process_roi
    compute_performance = TRUE,
    return_predictions  = FALSE,
    .fixed_design = .repnet_fixed_design(repnet_des$key, repnet_des$seed_rdm,
                                         repnet_des$confound_rdms),
    ...
  )
}


# Fixed-design helpers for repnet_model and repmed_model. The seed/confound
# design does not depend on the ROI, so its QR decomposition is built once on
# the spec. fit_roi() reuses it only when the current design matrix is
# identical, so edits to the spec after construction are still honoured.

#' @keywords internal
#' @noRd
.rep_X <- function(P) {
  cbind(`(Intercept)` = 1, P)
}

#' @keywords internal
#' @noRd
.rep_cached_decomp <- function(X, cache) {
  if (!is.null(cache) && !is.null(cache$X) && identical(cache$X, X)) cache$decomp else NULL
}

# OLS of y on X (X carries the intercept column) through a QR decomposition.
# Returns coefficients, standard errors, RSS, residual df and TSS, or NULL when
# the caller must use lm(): non-finite entries, rank deficiency (lm() reports
# NA for aliased terms) or a degenerate fit.
#' @keywords internal
#' @noRd
.rep_ols_fit <- function(y, X, decomp = NULL) {
  y <- as.numeric(y)
  if (!all(is.finite(X)) || !all(is.finite(y))) return(NULL)
  p <- ncol(X)
  n <- nrow(X)
  if (n <= p) return(NULL)
  if (is.null(decomp)) decomp <- qr(X)
  if (decomp$rank < p || !identical(as.integer(decomp$pivot), seq_len(p))) return(NULL)

  coef <- qr.coef(decomp, y)
  rss <- sum(qr.resid(decomp, y)^2)
  df <- n - p
  tss <- sum((y - mean(y))^2)
  sigma2 <- rss / df
  if (!is.finite(tss) || tss <= 0 || !is.finite(sigma2) || sigma2 <= 0) return(NULL)

  # diag((X'X)^-1) = rowSums(R^-1 ^ 2) for X = Q R
  Rinv <- backsolve(qr.R(decomp), diag(p))
  se <- sqrt(rowSums(Rinv^2) * sigma2)
  names(coef) <- names(se) <- colnames(X)
  list(coef = coef, se = se, rss = rss, df = df, tss = tss)
}

# Coefficients and standard errors for one regression. Falls back to lm() with
# the same data layout (and NA handling) the reference code used.
#' @keywords internal
#' @noRd
.rep_path_fit <- function(y, P, cache = NULL) {
  X <- .rep_X(P)
  fit <- .rep_ols_fit(y, X, .rep_cached_decomp(X, cache))
  if (!is.null(fit)) {
    return(list(coef = fit$coef, se = fit$se))
  }
  dat <- data.frame(.y = as.numeric(y), P, check.names = FALSE)
  fit_lm <- try(stats::lm(.y ~ ., data = dat), silent = TRUE)
  if (inherits(fit_lm, "try-error")) return(NULL)
  cf <- suppressWarnings(stats::coef(summary(fit_lm)))
  list(coef = cf[, "Estimate"], se = cf[, "Std. Error"])
}

# Seed and confound predictors over the lower triangle of the common items.
#' @keywords internal
#' @noRd
.repnet_predictors <- function(seed_sub, confound_rdms, common_items) {
  mm <- list(seed = as.numeric(seed_sub[lower.tri(seed_sub)]))
  for (nm in names(confound_rdms)) {
    M <- confound_rdms[[nm]]
    # Require labeled RDMs that include all common_items to avoid silent misalignment
    if (is.null(rownames(M)) || is.null(colnames(M)) ||
        !all(common_items %in% rownames(M)) || !all(common_items %in% colnames(M))) {
      return(list(error = sprintf("repnet_model: confound RDM '%s' lacks required item labels for current ROI.", nm)))
    }
    M_sub <- M[common_items, common_items, drop = FALSE]
    mm[[nm]] <- as.numeric(M_sub[lower.tri(M_sub)])
  }
  list(mm = mm)
}

# Fixed seed + confound design, built at construction. NULL when it cannot be
# built (fit_roi() then reports the problem as before).
#' @keywords internal
#' @noRd
.repnet_fixed_design <- function(key, seed_rdm, confound_rdms) {
  tryCatch({
    if (is.null(key) || is.null(seed_rdm)) return(NULL)
    items_key <- rownames(group_means(matrix(0, nrow = length(key), ncol = 1L),
                                      margin = 1, group = key))
    common_items <- sort(intersect(items_key, rownames(seed_rdm)))
    if (length(common_items) < 3L) return(NULL)
    seed_sub <- seed_rdm[common_items, common_items, drop = FALSE]
    pred <- .repnet_predictors(seed_sub, confound_rdms, common_items)
    if (!is.null(pred$error)) return(NULL)
    X <- .rep_X(as.matrix(as.data.frame(pred$mm, stringsAsFactors = FALSE)))
    list(X = X, decomp = if (all(is.finite(X))) qr(X) else NULL)
  }, error = function(e) NULL)
}


#' fit_roi method for repnet_model
#'
#' Computes representational connectivity metrics for a single ROI/searchlight.
#'
#' @rdname fit_roi
#' @method fit_roi repnet_model
#' @export
fit_roi.repnet_model <- function(model, roi_data, context, ...) {

  X   <- roi_data$train_data
  ind <- roi_data$indices

  n_obs <- nrow(X)
  n_vox <- ncol(X)

  if (n_obs < 2L || n_vox < 1L) {
    return(roi_result(
      metrics       = NULL,
      indices       = ind,
      id            = context$id,
      error         = TRUE,
      error_message = "repnet_model: insufficient observations or voxels."
    ))
  }

  key <- model$key
  if (length(key) != n_obs) {
    return(roi_result(
      metrics       = NULL,
      indices       = ind,
      id            = context$id,
      error         = TRUE,
      error_message = "repnet_model: key length mismatch with ROI data."
    ))
  }

  # Item-level prototypes
  E_full <- group_means(X, margin = 1, group = key)
  items_roi <- rownames(E_full)

  # Align with seed RDM
  seed_mat <- model$seed_rdm
  items_seed <- rownames(seed_mat)
  common_items <- intersect(items_roi, items_seed)
  K <- length(common_items)
  .log_debug(
    sprintf("repnet_model: items(roi)=%d, items(seed)=%d, common=%d",
            length(items_roi), length(items_seed), K)
  )

  if (K < 3L) {
    return(roi_result(
      metrics       = NULL,
      indices       = ind,
      id            = context$id,
      error         = TRUE,
      error_message = "repnet_model: need at least 3 common items to define RDM connectivity."
    ))
  }

  common_items <- sort(common_items)
  E <- E_full[common_items, , drop = FALSE]
  seed_sub <- seed_mat[common_items, common_items, drop = FALSE]

  # Drop constant voxels
  sd_E <- apply(E, 2L, stats::sd, na.rm = TRUE)
  keep_vox <- sd_E > 0
  if (!any(keep_vox)) {
    return(roi_result(
      metrics       = NULL,
      indices       = ind,
      id            = context$id,
      error         = TRUE,
      error_message = "repnet_model: all voxels have zero variance in this ROI."
    ))
  }
  E <- E[, keep_vox, drop = FALSE]

  # ROI RDM and seed vector
  D_roi <- pairwise_dist(model$distfun, E)
  d_roi  <- as.numeric(D_roi[lower.tri(D_roi)])
  d_seed <- as.numeric(seed_sub[lower.tri(seed_sub)])

  # Raw connectivity (simple RSA-style similarity)
  conn_raw <- suppressWarnings(
    stats::cor(d_roi, d_seed, method = model$simfun, use = "complete.obs")
  )
  conn_na <- is.na(conn_raw)

  # Build regression predictors (seed + confounds)
  pred <- .repnet_predictors(seed_sub, model$confound_rdms, common_items)
  if (!is.null(pred$error)) {
    return(roi_result(
      metrics       = NULL,
      indices       = ind,
      id            = context$id,
      error         = TRUE,
      error_message = pred$error
    ))
  }
  mm <- pred$mm

  beta_vec <- numeric(0)
  sp_vec   <- numeric(0)

  # Multiple regression on the seed + confound design. beta_<term> are the
  # regression coefficients; sp_<term> the semipartial correlations
  # t * sqrt(RSS / (TSS * df)). Both come from one QR of the fixed design.
  df <- as.data.frame(mm, stringsAsFactors = FALSE)
  X_fix <- .rep_X(as.matrix(df))
  fit <- .rep_ols_fit(d_roi, X_fix, .rep_cached_decomp(X_fix, model$.fixed_design))
  if (!is.null(fit)) {
    term <- names(fit$coef) != "(Intercept)"
    beta_vec <- fit$coef[term]
    sp_vec <- (fit$coef[term] / fit$se[term]) * sqrt(fit$rss / (fit$tss * fit$df))
    names(sp_vec) <- paste0("sp_", names(beta_vec))
  } else {
    # Reference path (non-finite rows, rank deficiency): lm() drops incomplete rows
    fit_lm <- try(stats::lm(d_roi ~ ., data = df), silent = TRUE)
    if (!inherits(fit_lm, "try-error")) {
      cf <- stats::coef(summary(fit_lm))
      if (nrow(cf) > 1L) {
        keep <- rownames(cf) != "(Intercept)"
        beta_vec <- cf[keep, "Estimate"]
        names(beta_vec) <- rownames(cf)[keep]
      }
      sp <- try(run_lm_semipartial(dvec = d_roi,
                                   obj = list(design = list(model_mat = df))),
                silent = TRUE)
      if (!inherits(sp, "try-error") && is.numeric(sp)) {
        names(sp) <- paste0("sp_", names(sp))
        sp_vec <- sp
      }
    }
  }

  perf <- c(
    n_items  = K,
    n_pairs  = length(d_roi),
    conn_raw = conn_raw
  )
  if (length(beta_vec)) {
    names(beta_vec) <- paste0("beta_", sub("^beta_", "", names(beta_vec)))
    perf <- c(perf, beta_vec)
  }
  if (length(sp_vec)) {
    perf <- c(perf, sp_vec)
  }

  roi_result(
    metrics = perf,
    indices = ind,
    id      = context$id
  )
}


#' @rdname output_schema
#' @method output_schema repnet_model
#' @export
output_schema.repnet_model <- function(model) {
  NULL
}
