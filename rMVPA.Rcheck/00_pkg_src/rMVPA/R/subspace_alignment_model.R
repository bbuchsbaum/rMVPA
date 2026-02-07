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


#' Per-ROI Subspace Alignment processing
#'
#' Computes PCA subspaces on source/target, aligns them, projects data, and
#' classifies target trials by correlation to source class prototypes in the
#' aligned space.
#'
#' @keywords internal
#' @export
process_roi.subspace_alignment_model <- function(mod_spec, roi, rnum, ...) {
  if (!has_test_set(mod_spec)) {
    return(tibble::tibble(
      result = list(NULL),
      indices = list(neuroim2::indices(roi$train_roi)),
      performance = list(NULL),
      id = rnum,
      error = TRUE,
      error_message = "subspace_alignment_model requires external test set"
    ))
  }

  Xtr <- as.matrix(neuroim2::values(roi$train_roi))  # n_train x p
  Xte <- as.matrix(neuroim2::values(roi$test_roi))   # n_test  x p
  des <- mod_spec$design
  ytr <- factor(des$y_train)
  yte <- factor(des$y_test)

  if (nrow(Xtr) < 2L || ncol(Xtr) < 1L || nrow(Xte) < 1L) {
    return(tibble::tibble(
      result = list(NULL),
      indices = list(neuroim2::indices(roi$train_roi)),
      performance = list(NULL),
      id = rnum,
      error = TRUE,
      error_message = "subspace_alignment_model: insufficient samples or features"
    ))
  }

  # Per-domain z-normalization (separately for source/target)
  Xtr <- scale(Xtr, center = mod_spec$center, scale = mod_spec$scale)
  Xte <- scale(Xte, center = mod_spec$center, scale = mod_spec$scale)

  # Cap d to avoid ill-posed PCA
  d_req <- if (is.null(mod_spec$d)) 20L else as.integer(mod_spec$d)
  d_cap <- min(d_req,
               ncol(Xtr), ncol(Xte),
               nrow(Xtr) - 1L, nrow(Xte) - 1L)
  if (is.na(d_cap) || d_cap < 1L) {
    return(tibble::tibble(
      result = list(NULL),
      indices = list(neuroim2::indices(roi$train_roi)),
      performance = list(NULL),
      id = rnum,
      error = TRUE,
      error_message = "subspace_alignment_model: subspace dimension < 1"
    ))
  }

  # PCA helpers with fallback to stats::prcomp
  pca_safe <- function(mat, ncomp) {
    tryCatch({
      irlba::prcomp_irlba(mat, n = ncomp, center = FALSE, scale. = FALSE)
    }, error = function(e) {
      stats::prcomp(mat, center = FALSE, scale. = FALSE, rank. = ncomp)
    })
  }

  pcs <- pca_safe(Xtr, d_cap)
  pct <- pca_safe(Xte, d_cap)

  Xs <- pcs$rotation[, seq_len(d_cap), drop = FALSE]  # p x d
  Xt <- pct$rotation[, seq_len(d_cap), drop = FALSE]  # p x d

  # Alignment: closed-form M = Xs^T Xt
  M <- crossprod(Xs, Xt)               # d x d
  Xa <- Xs %*% M                       # p x d (aligned source basis)

  S_aligned <- Xtr %*% Xa              # n_train x d
  T_proj    <- Xte %*% Xt              # n_test  x d

  # Class prototypes in aligned space
  proto <- .nx_rowsum_mean(S_aligned, ytr)  # classes x d
  levs <- rownames(proto)
  obs  <- factor(as.character(yte), levels = levs)

  # Correlation to prototypes; softmax for probabilities
  scores <- cor(t(T_proj), t(proto))
  if (is.null(colnames(scores))) colnames(scores) <- levs
  # Handle possible NA correlations (e.g., zero-variance rows)
  scores[is.na(scores)] <- 0
  probs <- .nx_softmax(scores)
  colnames(probs) <- levs
  pred <- factor(levs[max.col(scores, ties.method = "first")], levels = levs)

  # Alignment diagnostics
  align_err <- sqrt(sum((Xa - Xt)^2))

  cres <- classification_result(obs, pred, probs,
                                testind = seq_len(nrow(Xte)),
                                test_design = des$test_design,
                                predictor = list(d_used = d_cap,
                                                 alignment_frob = align_err))

  perf <- compute_performance(mod_spec, cres)

  tibble::tibble(
    result = list(cres),
    indices = list(neuroim2::indices(roi$train_roi)),
    performance = list(c(perf, d_used = d_cap, alignment_frob = align_err)),
    id = rnum,
    error = FALSE,
    error_message = "~"
  )
}
