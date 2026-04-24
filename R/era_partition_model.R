#' ERA Variance-Partition Model
#'
#' Compares first-order encoding-retrieval transfer with second-order geometry
#' preservation using matched variance-partition models. The model builds
#' item-level prototypes in a source state (`dataset$train_data`) and a target
#' state (`dataset$test_data`), then estimates how much variance is uniquely
#' explained by same-item transfer and by source-state geometry after optional
#' nuisance pair models.
#'
#' @param dataset An `mvpa_dataset` with `train_data` and `test_data`.
#' @param design An `mvpa_design` with train/test design tables.
#' @param key_var Column name or formula giving the item key shared across
#'   source and target states.
#' @param distfun Distance function used for within-state RDMs.
#' @param rsa_simfun Correlation method for the raw geometry summary.
#' @param first_order_nuisance Optional named list of `K x K` matrices or
#'   length `K^2` vectors for cross-state similarity nuisance regressors.
#'   Matrices are interpreted as target rows by source columns.
#' @param second_order_nuisance Optional named list of `K x K` matrices or
#'   lower-triangle vectors for geometry nuisance regressors.
#' @param item_block_enc,item_block_ret Optional item-level block labels for
#'   source and target states. Named vectors are matched to item keys.
#' @param item_time_enc,item_time_ret Optional item-level time/order values for
#'   source and target states. Named vectors are matched to item keys.
#' @param item_category Optional item-level category labels used to add a
#'   same-category nuisance model to both first- and second-order regressions.
#' @param compute_xdec_performance Logical; compute trial-level naive
#'   cross-decoding performance using the same prototype scorer as
#'   \code{\link{naive_xdec_model}}.
#' @param xdec_link_by Optional column name used to define source/target labels
#'   for the trial-level cross-decoding metrics. If \code{NULL}, \code{key_var}
#'   is used.
#' @param include_procrustes Logical; compute leave-one-item-out orthogonal
#'   Procrustes cross-decoding metrics.
#' @param procrustes_center Logical; center source and target prototypes using
#'   only alignment-training items before fitting Procrustes maps.
#' @param min_procrustes_train_items Minimum number of paired items allowed for
#'   each leave-one-item-out Procrustes alignment.
#' @param return_matrices Logical; store prototype/similarity matrices in each
#'   ROI result for diagnostics.
#' @param return_xdec_predictions Logical; store the trial-level
#'   \code{classification_result} produced by the direct cross-decoder in each
#'   ROI result.
#' @param ... Additional fields stored on the model spec.
#'
#' @return A model spec of class `era_partition_model` for `run_regional()` or
#'   `run_searchlight()`.
#' @export
era_partition_model <- function(dataset,
                                design,
                                key_var,
                                distfun = cordist(method = "pearson"),
                                rsa_simfun = c("pearson", "spearman"),
                                first_order_nuisance = NULL,
                                second_order_nuisance = NULL,
                                item_block_enc = NULL,
                                item_block_ret = NULL,
                                item_time_enc = NULL,
                                item_time_ret = NULL,
                                item_category = NULL,
                                compute_xdec_performance = TRUE,
                                xdec_link_by = NULL,
                                include_procrustes = TRUE,
                                procrustes_center = TRUE,
                                min_procrustes_train_items = 3L,
                                return_matrices = FALSE,
                                return_xdec_predictions = FALSE,
                                ...) {

  rsa_simfun <- match.arg(rsa_simfun)

  if (!inherits(dataset, "mvpa_dataset")) {
    stop("era_partition_model: `dataset` must inherit from 'mvpa_dataset'.", call. = FALSE)
  }
  if (!inherits(design, "mvpa_design")) {
    stop("era_partition_model: `design` must inherit from 'mvpa_design'.", call. = FALSE)
  }
  if (is.null(dataset$test_data) || is.null(design$test_design) || is.null(design$y_test)) {
    stop("era_partition_model requires an external test set (dataset$test_data + design$y_test).", call. = FALSE)
  }

  if (is.character(distfun)) {
    distfun <- create_dist(distfun)
  }
  if (!inherits(distfun, "distfun")) {
    stop("era_partition_model: `distfun` must be a distfun object or constructor name.", call. = FALSE)
  }

  key_vec <- factor(parse_variable(key_var, design$train_design))
  xdec_link_by_eff <- .era_partition_effective_xdec_link(
    design = design,
    link_by = xdec_link_by
  )
  xdec_key_train <- .era_partition_resolve_xdec_key(
    key_var = substitute(key_var),
    design_table = design$train_design,
    link_by = xdec_link_by_eff
  )
  xdec_levels <- levels(factor(xdec_key_train))
  xdec_performance <- if (length(xdec_levels) > 2L) {
    get_multiclass_perf(design$split_groups, class_metrics = FALSE)
  } else {
    get_binary_perf(design$split_groups)
  }

  create_model_spec(
    "era_partition_model",
    dataset = dataset,
    design = design,
    key = key_vec,
    key_var = substitute(key_var),
    distfun = distfun,
    rsa_simfun = rsa_simfun,
    first_order_nuisance = first_order_nuisance,
    second_order_nuisance = second_order_nuisance,
    item_block_enc = item_block_enc,
    item_block_ret = item_block_ret,
    item_time_enc = item_time_enc,
    item_time_ret = item_time_ret,
    item_category = item_category,
    compute_xdec_performance = isTRUE(compute_xdec_performance),
    xdec_link_by = xdec_link_by_eff,
    xdec_link_by_requested = xdec_link_by,
    xdec_performance = xdec_performance,
    include_procrustes = isTRUE(include_procrustes),
    procrustes_center = isTRUE(procrustes_center),
    min_procrustes_train_items = as.integer(min_procrustes_train_items),
    return_matrices = isTRUE(return_matrices),
    return_xdec_predictions = isTRUE(return_xdec_predictions),
    compute_performance = TRUE,
    return_predictions = FALSE,
    ...
  )
}

#' @export
print.era_partition_model <- function(x, ...) {
  cat("ERA Partition Model\n")
  cat("  key_var:                 ", deparse(x$key_var), "\n")
  cat("  direct xdec:             ", if (isTRUE(x$compute_xdec_performance)) "yes" else "no", "\n")
  cat("  xdec_link_by:            ", x$xdec_link_by %||% "(key_var)", "\n")
  cat("  second-order similarity: ", x$rsa_simfun, "\n")
  cat("  Procrustes decoder:      ", if (isTRUE(x$include_procrustes)) "yes" else "no", "\n")
  cat("  min Procrustes items:    ", x$min_procrustes_train_items, "\n")
  invisible(x)
}

#' @rdname fit_roi
#' @method fit_roi era_partition_model
#' @export
fit_roi.era_partition_model <- function(model, roi_data, context, ...) {
  if (!has_test_set(model)) {
    return(roi_result(
      metrics = NULL,
      indices = roi_data$indices,
      id = context$id,
      error = TRUE,
      error_message = "era_partition_model requires external test set"
    ))
  }

  Xenc <- as.matrix(roi_data$train_data)
  Xret <- as.matrix(roi_data$test_data)
  id <- context$id
  ind <- roi_data$indices

  if (is.null(Xret) || nrow(Xenc) < 2L || nrow(Xret) < 2L || ncol(Xenc) < 1L) {
    return(roi_result(
      metrics = NULL,
      indices = ind,
      id = id,
      error = TRUE,
      error_message = "era_partition_model: insufficient data"
    ))
  }

  des <- model$design
  key_enc <- factor(parse_variable(model$key_var, des$train_design), levels = levels(model$key))
  key_ret <- factor(parse_variable(model$key_var, des$test_design), levels = levels(model$key))

  E_full <- group_means(Xenc, margin = 1, group = key_enc)
  R_full <- group_means(Xret, margin = 1, group = key_ret)

  common_keys <- sort(intersect(rownames(E_full), rownames(R_full)))
  K <- length(common_keys)
  if (K < 3L) {
    return(roi_result(
      metrics = NULL,
      indices = ind,
      id = id,
      error = TRUE,
      error_message = sprintf("era_partition_model: need >=3 common items (found %d)", K)
    ))
  }

  E <- E_full[common_keys, , drop = FALSE]
  R <- R_full[common_keys, , drop = FALSE]

  sd_E <- apply(E, 2, stats::sd, na.rm = TRUE)
  sd_R <- apply(R, 2, stats::sd, na.rm = TRUE)
  keep_vox <- (sd_E > 0) | (sd_R > 0)
  if (!any(keep_vox)) {
    return(roi_result(
      metrics = NULL,
      indices = ind,
      id = id,
      error = TRUE,
      error_message = "era_partition_model: all voxels zero-variance"
    ))
  }
  E <- E[, keep_vox, drop = FALSE]
  R <- R[, keep_vox, drop = FALSE]

  # First-order cross-state transfer: target rows by source columns.
  S_cross <- .nx_row_cor(R, E)
  rownames(S_cross) <- common_keys
  colnames(S_cross) <- common_keys
  diag_sim <- diag(S_cross)
  off <- S_cross
  diag(off) <- NA_real_

  max_source_idx <- apply(S_cross, 1, function(row) if (all(is.na(row))) NA_integer_ else which.max(row))
  naive_top1_acc <- mean((!is.na(max_source_idx)) & (max_source_idx == seq_len(K)))
  naive_diag_mean <- mean(diag_sim, na.rm = TRUE)
  naive_diag_minus_off <- naive_diag_mean - mean(off, na.rm = TRUE)

  first_nuis <- .era_partition_first_nuisance(model, common_keys)
  first_part <- .era_partition_delta_r2(
    y = as.numeric(S_cross),
    signal = as.numeric(diag(K)),
    nuisance = first_nuis
  )

  xdec_perf <- .era_partition_empty_xdec_metrics(model)
  xdec_result <- NULL
  if (isTRUE(model$compute_xdec_performance)) {
    xdec_fit <- try(.era_partition_direct_xdec(model, Xenc, Xret), silent = TRUE)
    if (!inherits(xdec_fit, "try-error") && !is.null(xdec_fit)) {
      xdec_result <- xdec_fit$result
      xdec_vals <- try(model$xdec_performance(xdec_result), silent = TRUE)
      if (!inherits(xdec_vals, "try-error") && is.numeric(xdec_vals)) {
        names(xdec_vals) <- paste0("xdec_", names(xdec_vals))
        matched <- intersect(names(xdec_perf), names(xdec_vals))
        xdec_perf[matched] <- xdec_vals[matched]
      }
    }
  }

  # Second-order geometry preservation.
  D_enc <- pairwise_dist(model$distfun, E)
  D_ret <- pairwise_dist(model$distfun, R)
  dE <- as.numeric(D_enc[lower.tri(D_enc)])
  dR <- as.numeric(D_ret[lower.tri(D_ret)])
  geom_cor <- suppressWarnings(stats::cor(dE, dR, method = model$rsa_simfun, use = "complete.obs"))

  second_nuis <- .era_partition_second_nuisance(model, common_keys)
  second_part <- .era_partition_delta_r2(
    y = dR,
    signal = dE,
    nuisance = second_nuis
  )

  pro <- if (isTRUE(model$include_procrustes)) {
    .era_partition_procrustes_scores(
      E = E,
      R = R,
      center = isTRUE(model$procrustes_center),
      min_train_items = model$min_procrustes_train_items
    )
  } else {
    list(scores = matrix(NA_real_, K, K), train_n = rep(NA_real_, K))
  }

  pro_scores <- pro$scores
  rownames(pro_scores) <- common_keys
  colnames(pro_scores) <- common_keys
  pro_diag <- diag(pro_scores)
  pro_off <- pro_scores
  diag(pro_off) <- NA_real_
  pro_idx <- apply(pro_scores, 1, function(row) if (all(is.na(row))) NA_integer_ else which.max(row))
  pro_valid <- !is.na(pro_idx)
  procrustes_top1_acc <- if (any(pro_valid)) {
    mean(pro_idx[pro_valid] == seq_len(K)[pro_valid])
  } else {
    NA_real_
  }
  pro_diag_mean <- mean(pro_diag, na.rm = TRUE)
  pro_diag_minus_off <- pro_diag_mean - mean(pro_off, na.rm = TRUE)
  if (is.nan(pro_diag_mean)) pro_diag_mean <- NA_real_
  if (is.nan(pro_diag_minus_off)) pro_diag_minus_off <- NA_real_
  pro_train_n_mean <- mean(pro$train_n, na.rm = TRUE)
  if (is.nan(pro_train_n_mean)) pro_train_n_mean <- NA_real_

  perf <- c(
    n_items = K,
    first_order_delta_r2 = first_part$delta_r2,
    first_order_partial_r2 = first_part$partial_r2,
    first_order_full_r2 = first_part$full_r2,
    first_order_nuisance_r2 = first_part$nuisance_r2,
    first_order_beta = first_part$beta,
    first_order_n_obs = first_part$n_obs,
    second_order_delta_r2 = second_part$delta_r2,
    second_order_partial_r2 = second_part$partial_r2,
    second_order_full_r2 = second_part$full_r2,
    second_order_nuisance_r2 = second_part$nuisance_r2,
    second_order_beta = second_part$beta,
    second_order_n_obs = second_part$n_obs,
    naive_top1_acc = naive_top1_acc,
    naive_diag_mean = naive_diag_mean,
    naive_diag_minus_off = naive_diag_minus_off,
    xdec_perf,
    geom_cor = geom_cor,
    procrustes_top1_acc = procrustes_top1_acc,
    procrustes_diag_mean = pro_diag_mean,
    procrustes_diag_minus_off = pro_diag_minus_off,
    procrustes_train_n_mean = pro_train_n_mean,
    nuisance_first_order_n = length(first_nuis),
    nuisance_second_order_n = length(second_nuis)
  )

  diag_result <- list()
  if (isTRUE(model$return_matrices)) {
    diag_result <- c(diag_result, list(
      keys = common_keys,
      source_prototypes = E,
      target_prototypes = R,
      cross_similarity = S_cross,
      source_rdm = D_enc,
      target_rdm = D_ret,
      procrustes_similarity = pro_scores
    ))
  }
  if (isTRUE(model$return_xdec_predictions)) {
    diag_result$xdec_result <- xdec_result
  }
  if (!length(diag_result)) {
    diag_result <- NULL
  }

  roi_result(metrics = perf, indices = ind, id = id, result = diag_result)
}

#' @rdname output_schema
#' @method output_schema era_partition_model
#' @export
output_schema.era_partition_model <- function(model) {
  nms <- c(
    "n_items",
    "first_order_delta_r2",
    "first_order_partial_r2",
    "first_order_full_r2",
    "first_order_nuisance_r2",
    "first_order_beta",
    "first_order_n_obs",
    "second_order_delta_r2",
    "second_order_partial_r2",
    "second_order_full_r2",
    "second_order_nuisance_r2",
    "second_order_beta",
    "second_order_n_obs",
    "naive_top1_acc",
    "naive_diag_mean",
    "naive_diag_minus_off",
    .era_partition_xdec_metric_names(model),
    "geom_cor",
    "procrustes_top1_acc",
    "procrustes_diag_mean",
    "procrustes_diag_minus_off",
    "procrustes_train_n_mean",
    "nuisance_first_order_n",
    "nuisance_second_order_n"
  )
  setNames(rep("scalar", length(nms)), nms)
}

#' @keywords internal
.era_partition_effective_xdec_link <- function(design, link_by = NULL) {
  if (!is.null(link_by) && length(link_by) == 1L &&
      link_by %in% colnames(design$train_design) &&
      link_by %in% colnames(design$test_design)) {
    return(link_by)
  }
  NULL
}

#' @keywords internal
.era_partition_resolve_xdec_key <- function(key_var, design_table, link_by = NULL) {
  if (!is.null(link_by) && length(link_by) == 1L && link_by %in% colnames(design_table)) {
    design_table[[link_by]]
  } else {
    parse_variable(key_var, design_table)
  }
}

#' @keywords internal
.era_partition_xdec_keys <- function(model) {
  des <- model$design
  key_tr <- factor(.era_partition_resolve_xdec_key(
    key_var = model$key_var,
    design_table = des$train_design,
    link_by = model$xdec_link_by
  ))
  key_te <- factor(.era_partition_resolve_xdec_key(
    key_var = model$key_var,
    design_table = des$test_design,
    link_by = model$xdec_link_by
  ), levels = levels(key_tr))

  list(key_tr = key_tr, key_te = key_te)
}

#' @keywords internal
.era_partition_direct_xdec <- function(model, Xenc, Xret) {
  keys <- .era_partition_xdec_keys(model)
  if (length(levels(keys$key_tr)) < 2L) {
    return(NULL)
  }

  core <- .naive_xdec_fit_core(
    Xtr = Xenc,
    Xte = Xret,
    kernel = NULL,
    key_tr = keys$key_tr,
    key_te = keys$key_te
  )

  list(
    result = classification_result(
      core$obs,
      core$pred,
      core$probs,
      testind = seq_len(nrow(Xret)),
      test_design = model$design$test_design,
      predictor = NULL
    )
  )
}

#' @keywords internal
.era_partition_xdec_metric_names <- function(model) {
  if (!isTRUE(model$compute_xdec_performance)) {
    return(character())
  }

  base <- c("Accuracy", "AUC")
  split_groups <- model$design$split_groups
  if (!is.null(split_groups) && length(split_groups) > 0L) {
    split_names <- unlist(lapply(names(split_groups), function(tag) {
      paste0(base, "_", tag)
    }), use.names = FALSE)
    base <- c(base, split_names)
  }

  paste0("xdec_", base)
}

#' @keywords internal
.era_partition_empty_xdec_metrics <- function(model) {
  nms <- .era_partition_xdec_metric_names(model)
  setNames(rep(NA_real_, length(nms)), nms)
}

#' @keywords internal
.era_partition_delta_r2 <- function(y, signal, nuisance = NULL) {
  y <- as.numeric(y)
  signal <- as.numeric(signal)
  if (length(y) != length(signal)) {
    stop(".era_partition_delta_r2: `y` and `signal` must have matching length.", call. = FALSE)
  }

  dat <- data.frame(.y = y, .signal = signal)
  nuisance <- .era_partition_nuisance_dataframe(nuisance, length(y))
  if (length(nuisance)) {
    dat <- cbind(dat, nuisance)
  }

  keep <- stats::complete.cases(dat)
  keep <- keep & apply(dat, 1, function(row) all(is.finite(row)))
  dat <- dat[keep, , drop = FALSE]

  out_na <- list(
    delta_r2 = NA_real_,
    partial_r2 = NA_real_,
    full_r2 = NA_real_,
    nuisance_r2 = NA_real_,
    beta = NA_real_,
    n_obs = nrow(dat)
  )

  if (nrow(dat) < 3L || stats::sd(dat$.y) == 0 || stats::sd(dat$.signal) == 0) {
    return(out_na)
  }

  nuisance_names <- setdiff(names(dat), c(".y", ".signal"))
  full_terms <- c(".signal", nuisance_names)
  full_formula <- stats::as.formula(paste(".y ~", paste(full_terms, collapse = " + ")))
  nuisance_formula <- if (length(nuisance_names)) {
    stats::as.formula(paste(".y ~", paste(nuisance_names, collapse = " + ")))
  } else {
    stats::as.formula(".y ~ 1")
  }

  fit_full <- try(stats::lm(full_formula, data = dat), silent = TRUE)
  fit_nuis <- try(stats::lm(nuisance_formula, data = dat), silent = TRUE)
  if (inherits(fit_full, "try-error") || inherits(fit_nuis, "try-error")) {
    return(out_na)
  }

  r2_full <- .era_partition_r2(fit_full, dat$.y)
  r2_nuis <- .era_partition_r2(fit_nuis, dat$.y)
  delta <- r2_full - r2_nuis
  denom <- 1 - r2_nuis
  partial <- if (is.finite(denom) && denom > 0) delta / denom else NA_real_
  beta <- unname(stats::coef(fit_full)[".signal"])
  if (is.na(beta)) beta <- NA_real_

  list(
    delta_r2 = delta,
    partial_r2 = partial,
    full_r2 = r2_full,
    nuisance_r2 = r2_nuis,
    beta = beta,
    n_obs = nrow(dat)
  )
}

#' @keywords internal
.era_partition_r2 <- function(fit, y) {
  rss <- sum(stats::resid(fit)^2)
  tss <- sum((y - mean(y))^2)
  if (!is.finite(tss) || tss <= 0) {
    return(NA_real_)
  }
  1 - rss / tss
}

#' @keywords internal
.era_partition_nuisance_dataframe <- function(nuisance, n) {
  if (is.null(nuisance) || length(nuisance) == 0L) {
    return(data.frame())
  }
  if (!is.list(nuisance)) {
    stop("Nuisance regressors must be supplied as a named list.", call. = FALSE)
  }
  nms <- names(nuisance)
  if (is.null(nms) || any(!nzchar(nms))) {
    nms <- paste0("nuisance_", seq_along(nuisance))
  }
  nms <- make.names(nms, unique = TRUE)
  out <- lapply(nuisance, function(x) {
    x <- as.numeric(x)
    if (length(x) != n) {
      stop("Nuisance regressor length does not match the analysis vector.", call. = FALSE)
    }
    x
  })
  names(out) <- nms
  as.data.frame(out, optional = TRUE)
}

#' @keywords internal
.era_partition_first_nuisance <- function(model, keys) {
  K <- length(keys)
  out <- list()

  block_enc <- .era_partition_align_item_vector(model$item_block_enc, keys)
  block_ret <- .era_partition_align_item_vector(model$item_block_ret, keys)
  if (!is.null(block_enc) && !is.null(block_ret)) {
    out$same_block_cross <- as.numeric(outer(block_ret, block_enc, "=="))
  }

  time_enc <- .era_partition_align_item_vector(model$item_time_enc, keys)
  time_ret <- .era_partition_align_item_vector(model$item_time_ret, keys)
  if (!is.null(time_enc) && !is.null(time_ret)) {
    out$enc_time <- matrix(rep(as.numeric(time_enc), each = K), nrow = K)
    out$ret_time <- matrix(rep(as.numeric(time_ret), times = K), nrow = K)
    out$abs_lag <- abs(outer(as.numeric(time_ret), as.numeric(time_enc), "-"))
  }

  category <- .era_partition_align_item_vector(model$item_category, keys)
  if (!is.null(category)) {
    out$same_category <- as.numeric(outer(category, category, "=="))
  }

  user <- .era_partition_cross_user_nuisance(model$first_order_nuisance, keys)
  c(out, user)
}

#' @keywords internal
.era_partition_second_nuisance <- function(model, keys) {
  out <- list()

  block_enc <- .era_partition_align_item_vector(model$item_block_enc, keys)
  block_ret <- .era_partition_align_item_vector(model$item_block_ret, keys)
  if (!is.null(block_enc)) {
    M <- outer(block_enc, block_enc, "==")
    out$same_block_enc <- as.numeric(M[lower.tri(M)])
  }
  if (!is.null(block_ret)) {
    M <- outer(block_ret, block_ret, "==")
    out$same_block_ret <- as.numeric(M[lower.tri(M)])
  }

  time_enc <- .era_partition_align_item_vector(model$item_time_enc, keys)
  time_ret <- .era_partition_align_item_vector(model$item_time_ret, keys)
  if (!is.null(time_enc)) {
    M <- abs(outer(as.numeric(time_enc), as.numeric(time_enc), "-"))
    out$temporal_distance_enc <- as.numeric(M[lower.tri(M)])
  }
  if (!is.null(time_ret)) {
    M <- abs(outer(as.numeric(time_ret), as.numeric(time_ret), "-"))
    out$temporal_distance_ret <- as.numeric(M[lower.tri(M)])
  }

  category <- .era_partition_align_item_vector(model$item_category, keys)
  if (!is.null(category)) {
    M <- outer(category, category, "==")
    out$same_category <- as.numeric(M[lower.tri(M)])
  }

  user <- .era_partition_geometry_user_nuisance(model$second_order_nuisance, keys)
  c(out, user)
}

#' @keywords internal
.era_partition_align_item_vector <- function(x, keys) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.null(names(x))) {
    x <- x[match(keys, names(x))]
  } else {
    x <- x[seq_along(keys)]
  }
  x
}

#' @keywords internal
.era_partition_cross_user_nuisance <- function(nuisance, keys) {
  if (is.null(nuisance)) return(list())
  K <- length(keys)
  lapply(nuisance, function(x) {
    if (is.matrix(x) || is.data.frame(x) || inherits(x, "dist")) {
      M <- as.matrix(x)
      if (!is.null(rownames(M)) && !is.null(colnames(M)) &&
          all(keys %in% rownames(M)) && all(keys %in% colnames(M))) {
        M <- M[keys, keys, drop = FALSE]
      } else {
        M <- M[seq_len(K), seq_len(K), drop = FALSE]
      }
      as.numeric(M)
    } else {
      as.numeric(x)
    }
  })
}

#' @keywords internal
.era_partition_geometry_user_nuisance <- function(nuisance, keys) {
  if (is.null(nuisance)) return(list())
  K <- length(keys)
  lower_n <- K * (K - 1L) / 2L
  lapply(nuisance, function(x) {
    if (is.matrix(x) || is.data.frame(x) || inherits(x, "dist")) {
      M <- as.matrix(x)
      if (!is.null(rownames(M)) && !is.null(colnames(M)) &&
          all(keys %in% rownames(M)) && all(keys %in% colnames(M))) {
        M <- M[keys, keys, drop = FALSE]
      } else {
        M <- M[seq_len(K), seq_len(K), drop = FALSE]
      }
      as.numeric(M[lower.tri(M)])
    } else {
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
  })
}

#' @keywords internal
.era_partition_procrustes_fit <- function(E, R, train_idx, center = TRUE) {
  Etr <- as.matrix(E[train_idx, , drop = FALSE])
  Rtr <- as.matrix(R[train_idx, , drop = FALSE])
  mu_E <- if (isTRUE(center)) colMeans(Etr) else rep(0, ncol(Etr))
  mu_R <- if (isTRUE(center)) colMeans(Rtr) else rep(0, ncol(Rtr))
  Ec <- sweep(Etr, 2, mu_E, "-")
  Rc <- sweep(Rtr, 2, mu_R, "-")
  s <- svd(crossprod(Ec, Rc))
  Tmat <- s$u %*% t(s$v)
  list(T = Tmat, mu_E = mu_E, mu_R = mu_R, train_idx = train_idx)
}

#' @keywords internal
.era_partition_apply_procrustes <- function(E, fit) {
  sweep(sweep(as.matrix(E), 2, fit$mu_E, "-") %*% fit$T, 2, fit$mu_R, "+")
}

#' @keywords internal
.era_partition_procrustes_scores <- function(E, R, center = TRUE, min_train_items = 3L) {
  E <- as.matrix(E)
  R <- as.matrix(R)
  K <- nrow(E)
  scores <- matrix(NA_real_, K, K)
  train_n <- rep(NA_real_, K)
  min_train_items <- as.integer(min_train_items)
  if (is.na(min_train_items) || min_train_items < 1L) {
    min_train_items <- 1L
  }

  for (i in seq_len(K)) {
    train_idx <- setdiff(seq_len(K), i)
    train_n[i] <- length(train_idx)
    if (length(train_idx) < min_train_items) {
      next
    }
    fit <- try(.era_partition_procrustes_fit(E, R, train_idx, center = center), silent = TRUE)
    if (inherits(fit, "try-error")) {
      next
    }
    E_aligned <- .era_partition_apply_procrustes(E, fit)
    scores[i, ] <- .nx_row_cor(R[i, , drop = FALSE], E_aligned)
  }
  list(scores = scores, train_n = train_n)
}
