#' Subspace Alignment cross-decoder
#'
#' Fast unsupervised domain-adaptation baseline following Fernando et al. (ICCV 2013).
#' Learns separate PCA subspaces for source (train) and target (test), aligns them
#' with a closed-form map \eqn{M = X_S^T X_T}, projects both domains, and
#' classifies target trials via correlation to source class prototypes in the
#' aligned space. Requires an external test set but no target labels for fitting.
#'
#' @param dataset mvpa_dataset with `train_data` (source) and `test_data` (target).
#' @param design mvpa_design with `y_train` (source labels) and `y_test` (for evaluation).
#' @param d Integer subspace dimension; capped automatically by samples/features.
#' @param center,scale Logical flags for per-domain z-normalization prior to PCA.
#' @param return_predictions logical; keep per-ROI predictions (default TRUE).
#' @param ... Additional arguments stored on the model spec.
#'
#' @return A model spec of class `subspace_alignment_model` for use with
#'   `run_regional()` / `run_searchlight()`.
#'
#' @examples
#' \dontrun{
#'   ds <- gen_sample_dataset(c(5,5,5), 20, external_test=TRUE)
#'   model <- subspace_alignment_model(ds$dataset, ds$design, d=10)
#' }
#' @export
subspace_alignment_model <- function(dataset,
                                     design,
                                     d = 20L,
                                     center = TRUE,
                                     scale = TRUE,
                                     return_predictions = TRUE,
                                     ...) {

  if (is.null(design$y_test) || is.null(dataset$test_data)) {
    stop("subspace_alignment_model requires an external test set (dataset$test_data + design$y_test)")
  }

  if (is.numeric(design$y_train)) {
    stop("subspace_alignment_model currently supports classification (factor y_train) only")
  }

  perf_fun <- if (length(levels(design$y_train)) > 2) {
    get_multiclass_perf(design$split_groups, class_metrics = FALSE)
  } else {
    get_binary_perf(design$split_groups)
  }

  create_model_spec(
    "subspace_alignment_model",
    dataset = dataset,
    design = design,
    d = as.integer(d),
    center = center,
    scale = scale,
    performance = perf_fun,
    compute_performance = TRUE,
    return_predictions = return_predictions,
    ...
  )
}

#' @export
compute_performance.subspace_alignment_model <- function(obj, result) {
  obj$performance(result)
}


#' @rdname fit_roi
#' @method fit_roi subspace_alignment_model
#' @export
fit_roi.subspace_alignment_model <- function(model, roi_data, context, ...) {
  if (!has_test_set(model)) {
    return(roi_result(
      metrics = NULL,
      indices = roi_data$indices,
      id = context$id,
      error = TRUE,
      error_message = "subspace_alignment_model requires external test set"
    ))
  }

  Xtr <- roi_data$train_data
  Xte <- roi_data$test_data
  des <- model$design
  ytr <- factor(des$y_train)
  yte <- factor(des$y_test)

  if (is.null(Xte) || nrow(Xtr) < 2L || ncol(Xtr) < 1L || nrow(Xte) < 1L) {
    return(roi_result(
      metrics = NULL,
      indices = roi_data$indices,
      id = context$id,
      error = TRUE,
      error_message = "subspace_alignment_model: insufficient samples or features"
    ))
  }

  # Per-domain z-normalization
  Xtr <- scale(Xtr, center = model$center, scale = model$scale)
  Xte <- scale(Xte, center = model$center, scale = model$scale)

  # Cap d
  d_req <- if (is.null(model$d)) 20L else as.integer(model$d)
  d_cap <- min(d_req, ncol(Xtr), ncol(Xte), nrow(Xtr) - 1L, nrow(Xte) - 1L)
  if (is.na(d_cap) || d_cap < 1L) {
    return(roi_result(
      metrics = NULL,
      indices = roi_data$indices,
      id = context$id,
      error = TRUE,
      error_message = "subspace_alignment_model: subspace dimension < 1"
    ))
  }

  # PCA
  pca_safe <- function(mat, ncomp) {
    tryCatch({
      irlba::prcomp_irlba(mat, n = ncomp, center = FALSE, scale. = FALSE)
    }, error = function(e) {
      stats::prcomp(mat, center = FALSE, scale. = FALSE, rank. = ncomp)
    })
  }

  pcs <- pca_safe(Xtr, d_cap)
  pct <- pca_safe(Xte, d_cap)

  Xs <- pcs$rotation[, seq_len(d_cap), drop = FALSE]
  Xt <- pct$rotation[, seq_len(d_cap), drop = FALSE]

  # Alignment
  M <- crossprod(Xs, Xt)
  Xa <- Xs %*% M

  S_aligned <- Xtr %*% Xa
  T_proj    <- Xte %*% Xt

  # Class prototypes
  proto <- .nx_rowsum_mean(S_aligned, ytr)
  levs <- rownames(proto)
  obs  <- factor(as.character(yte), levels = levs)

  scores <- cor(t(T_proj), t(proto))
  if (is.null(colnames(scores))) colnames(scores) <- levs
  scores[is.na(scores)] <- 0
  probs <- .nx_softmax(scores)
  colnames(probs) <- levs
  pred <- factor(levs[max.col(scores, ties.method = "first")], levels = levs)

  align_err <- sqrt(sum((Xa - Xt)^2))

  cres <- classification_result(obs, pred, probs,
                                testind = seq_len(nrow(Xte)),
                                test_design = des$test_design,
                                predictor = list(d_used = d_cap,
                                                 alignment_frob = align_err))

  perf <- compute_performance(model, cres)

  roi_result(
    metrics = c(perf, d_used = d_cap, alignment_frob = align_err),
    indices = roi_data$indices,
    id = context$id,
    result = cres
  )
}

#' @rdname output_schema
#' @method output_schema subspace_alignment_model
#' @export
output_schema.subspace_alignment_model <- function(model) {
  NULL
}

