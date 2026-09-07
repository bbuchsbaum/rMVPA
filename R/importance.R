#' @rdname extract_weights
#' @export
extract_weights.sda <- function(object, ...) {
  # sda$beta is K x P (classes x features); transpose to P x K
  if (is.null(object$beta)) {
    stop("extract_weights.sda: fitted sda object has no 'beta' component. ",
         "Did you pass a raw sda fit?")
  }
  t(as.matrix(object$beta))
}

#' @rdname extract_weights
#' @export
extract_weights.glmnet <- function(object, ...) {
  s <- object$opt_lambda
  if (is.null(s)) {
    stop("extract_weights.glmnet: fitted glmnet object has no 'opt_lambda'. ",
         "Expected a model fitted via glmnet_opt.")
  }

  cf <- coef(object, s = s)

  if (is.list(cf)) {
    # Multinomial: cf is a list of sparse column vectors, one per class
    # Each element is (P+1) x 1 sparse matrix (intercept + features)
    W <- do.call(cbind, lapply(cf, function(cc) {
      as.numeric(cc[-1, , drop = TRUE])
    }))
    colnames(W) <- names(cf)
  } else {
    # Binomial: cf is a single (P+1) x 1 sparse vector
    W <- matrix(as.numeric(cf[-1, , drop = TRUE]), ncol = 1)
  }

  W
}

#' @rdname extract_weights
#' @export
extract_weights.spacenet_fit <- function(object, ...) {
  beta <- object$beta
  if (is.matrix(beta)) {
    beta
  } else {
    matrix(as.numeric(beta), ncol = 1)
  }
}

#' @rdname extract_weights
#' @export
extract_weights.model_fit <- function(object, ...) {
  fit <- object$fit
  if (is.null(fit)) {
    stop("extract_weights.model_fit: object has no underlying $fit component.")
  }
  extract_weights(fit, ...)
}

#' @rdname extract_weights
#' @export
extract_weights.default <- function(object, ...) {
  stop("extract_weights: no method for class '", paste(class(object), collapse = "/"),
       "'. Implement an extract_weights method for this model type.")
}

#' Haufe Feature Importance (Activation Patterns)
#'
#' Computes per-feature importance using the Haufe et al. (2014) transformation
#' from decoding weights to encoding (activation) patterns:
#' \code{A = Sigma_x \%*\% W \%*\% solve(t(W) \%*\% Sigma_x \%*\% W)}.
#'
#' @details
#' Two equivalent computation paths are available. Supplying \code{Sigma_x}
#' uses the explicit \eqn{P \times P} covariance. Supplying the training
#' observations \code{X} instead uses the matrix-free identity
#' \deqn{Z = X_c W}
#' \deqn{A = X_c^\top Z (Z^\top Z)^{+}}
#' where \eqn{X_c} is the column-centred data. The covariance normalisation
#' cancels exactly, so the two paths agree to numerical precision, but the
#' matrix-free path never forms a \eqn{P \times P} matrix: it needs only
#' \eqn{n \times D}, \eqn{P \times D}, and \eqn{D \times D} quantities. Use it
#' for whole-brain feature counts, where a dense covariance is infeasible.
#' Because \eqn{Z} is centred, \eqn{X_c^\top Z = X^\top Z}, so \code{X} itself
#' is never copied; \code{block_size} additionally limits the number of feature
#' columns multiplied at once.
#'
#' @param W A P x D weight matrix (features x discriminant directions).
#' @param Sigma_x A P x P covariance matrix of the training features. Ignored
#'   when \code{X} is supplied.
#' @param summary_fun A function applied to the rows of A to produce a scalar
#'   importance per feature. Defaults to the L2 norm across discriminants.
#' @param X Optional n x P matrix of training observations. When supplied, the
#'   matrix-free path is used and \code{Sigma_x} is not needed.
#' @param center Logical; centre the columns of \code{X} before computing the
#'   patterns (default \code{TRUE}, matching \code{cov()}).
#' @param block_size Optional integer; when \code{X} is supplied, compute the
#'   activation patterns in column blocks of at most this many features.
#'
#' @return A list with components:
#'   \describe{
#'     \item{A}{The P x D activation pattern matrix.}
#'     \item{importance}{A numeric vector of length P with per-feature importance.}
#'   }
#'
#' @references
#' Haufe, S., Meinecke, F., Goergen, K., Daehne, S., Haynes, J.D.,
#' Blankertz, B., & Biessmann, F. (2014). On the interpretation of weight
#' vectors of linear models in multivariate neuroimaging. NeuroImage, 87, 96-110.
#'
#' @examples
#' \donttest{
#'   W <- matrix(rnorm(10*2), 10, 2)
#'   X <- matrix(rnorm(50*10), 50, 10)
#'   haufe_importance(W, cov(X))
#'   # identical, without forming cov(X):
#'   haufe_importance(W, X = X)
#' }
#' @export
haufe_importance <- function(W, Sigma_x = NULL,
                              summary_fun = function(A) sqrt(rowSums(A^2)),
                              X = NULL, center = TRUE, block_size = NULL) {
  W <- as.matrix(W)

  if (!is.null(X)) {
    A <- .haufe_patterns_from_data(W, X, center = center, block_size = block_size)
  } else {
    if (is.null(Sigma_x)) {
      stop("haufe_importance: supply either 'Sigma_x' or 'X'.", call. = FALSE)
    }
    Sigma_x <- as.matrix(Sigma_x)

    stopifnot(nrow(W) == nrow(Sigma_x),
              nrow(Sigma_x) == ncol(Sigma_x))

    # W'*Sigma_x*W  (D x D)
    WtSW <- crossprod(W, Sigma_x %*% W)

    # Use ginv for robustness to singular matrices
    WtSW_inv <- tryCatch(
      solve(WtSW),
      error = function(e) MASS::ginv(WtSW)
    )

    # A = Sigma_x * W * (W'*Sigma_x*W)^{-1}
    A <- Sigma_x %*% W %*% WtSW_inv
  }

  importance <- summary_fun(A)

  list(A = A, importance = as.numeric(importance))
}

#' @keywords internal
#' @noRd
.haufe_patterns_from_data <- function(W, X, center = TRUE, block_size = NULL) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.numeric(X)) stop("haufe_importance: 'X' must be a numeric matrix.", call. = FALSE)
  n <- nrow(X)
  P <- ncol(X)
  if (nrow(W) != P) {
    stop(sprintf("haufe_importance: nrow(W) = %d must equal ncol(X) = %d.", nrow(W), P),
         call. = FALSE)
  }
  if (n < 2L) stop("haufe_importance: 'X' needs at least two observations.", call. = FALSE)

  # Z = X_c W computed without materialising X_c: X W - 1 (mu' W)
  Z <- X %*% W
  if (center) {
    mu <- colMeans(X)
    Z <- sweep(Z, 2L, as.numeric(crossprod(mu, W)), "-")
  }
  # Scale by (n - 1) so that the D x D matrix matches W' cov(X) W exactly and the
  # ginv fallback behaves identically to the covariance path.
  denom <- n - 1
  ZtZ <- crossprod(Z) / denom
  ZtZ_inv <- tryCatch(
    solve(ZtZ),
    error = function(e) MASS::ginv(ZtZ)
  )

  # X_c' Z = X' Z because Z is column-centred (when center = TRUE); when not
  # centring, X' Z is the uncentred cross-product, matching crossprod(X)/denom.
  D <- ncol(W)
  if (is.null(block_size) || block_size >= P) {
    XtZ <- crossprod(X, Z) / denom
    A <- XtZ %*% ZtZ_inv
  } else {
    block_size <- max(1L, as.integer(block_size))
    A <- matrix(0, P, D)
    starts <- seq.int(1L, P, by = block_size)
    for (s in starts) {
      e <- min(P, s + block_size - 1L)
      A[s:e, ] <- (crossprod(X[, s:e, drop = FALSE], Z) / denom) %*% ZtZ_inv
    }
  }
  dimnames(A) <- list(colnames(X), colnames(W))
  A
}


# ---------- model_importance methods ----------

#' @rdname model_importance
#' @param summary_fun Optional function to summarize the activation pattern matrix rows. Default NULL uses L2 norm.
#' @export
model_importance.sda <- function(object, X_train, summary_fun = NULL, ...) {
  W <- extract_weights(object)
  haufe_args <- list(W = W, X = X_train)
  if (!is.null(summary_fun)) {
    haufe_args$summary_fun <- summary_fun
  }
  res <- do.call(haufe_importance, haufe_args)
  res$importance
}

#' @rdname model_importance
#' @param summary_fun Optional function to summarize the activation pattern matrix rows. Default NULL uses L2 norm.
#' @export
model_importance.glmnet <- function(object, X_train, summary_fun = NULL, ...) {
  W <- extract_weights(object)
  haufe_args <- list(W = W, X = X_train)
  if (!is.null(summary_fun)) {
    haufe_args$summary_fun <- summary_fun
  }
  res <- do.call(haufe_importance, haufe_args)
  res$importance
}

#' @rdname model_importance
#' @param summary_fun Optional function to summarize the activation pattern matrix rows. Default NULL uses L2 norm.
#' @export
model_importance.spacenet_fit <- function(object, X_train, summary_fun = NULL, ...) {
  W <- extract_weights(object)
  haufe_args <- list(W = W, X = X_train)
  if (!is.null(summary_fun)) {
    haufe_args$summary_fun <- summary_fun
  }
  res <- do.call(haufe_importance, haufe_args)
  res$importance
}

#' @rdname model_importance
#' @details
#' \code{model_importance.randomForest} returns Gini-based importance
#' (MeanDecreaseGini) by default, or permutation-based importance
#' (MeanDecreaseAccuracy) when the forest was trained with
#' \code{importance = TRUE}.
#' Both are \strong{backward} (decoding) measures and may assign high
#' importance to suppressor variables that do not carry neural signal.
#' For neuroscience interpretation, prefer \code{\link{haufe_importance}}
#' with a linear model or use \code{\link{region_importance}} for a
#' model-agnostic alternative bounded by out-of-sample performance.
#' @export
model_importance.randomForest <- function(object, X_train, ...) {
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    return(NULL)
  }
  imp <- randomForest::importance(object)
  if ("MeanDecreaseAccuracy" %in% colnames(imp)) {
    as.numeric(imp[, "MeanDecreaseAccuracy"])
  } else if ("MeanDecreaseGini" %in% colnames(imp)) {
    as.numeric(imp[, "MeanDecreaseGini"])
  } else if ("IncNodePurity" %in% colnames(imp)) {
    as.numeric(imp[, "IncNodePurity"])
  } else {
    as.numeric(imp[, 1])
  }
}

#' @rdname model_importance
#' @export
model_importance.default <- function(object, X_train, ...) {
  NULL
}


# ---------- cv_score_global ----------

#' Compute a scalar CV performance score for a feature matrix
#'
#' Runs the shared \code{cv_run_global} loop on the supplied feature matrix
#' and returns a single scalar performance value. Used internally by
#' \code{region_importance}.
#'
#' @param mspec An \code{mvpa_model} specification (may have dataset stripped).
#' @param X A T x P feature matrix (the subset of columns for this iteration).
#' @param y Response vector.
#' @param crossval A cross-validation specification.
#' @param feature_ids Integer vector of feature IDs for \code{train_model}.
#' @param metric Character name of the metric to extract, or NULL for the first.
#' @return A single numeric performance value.
#' @keywords internal
#' @noRd
cv_score_global <- function(mspec, X, y, crossval, feature_ids, metric = NULL) {
  cv_out <- cv_run_global(mspec, X, y, crossval, feature_ids,
                          return_fits = FALSE, extract_weights_fn = FALSE)

  design <- mspec$design
  cres <- tryCatch(
    wrap_result(cv_out$result_table, design),
    error = function(e) {
      futile.logger::flog.warn("cv_score_global: wrap_result failed: %s", e$message)
      NULL
    }
  )

  if (is.null(cres)) return(NA_real_)

  perf <- tryCatch(
    compute_performance(mspec, cres),
    error = function(e) {
      futile.logger::flog.warn("cv_score_global: compute_performance failed: %s", e$message)
      NULL
    }
  )

  if (is.null(perf) || length(perf) == 0) return(NA_real_)

  if (!is.null(metric)) {
    if (metric %in% names(perf)) perf[[metric]] else perf[[1]]
  } else {
    perf[[1]]
  }
}


# ---------- region_importance ----------

#' @rdname region_importance
#' @param model_spec An \code{mvpa_model} specification.
#' @param n_iter Number of random subset iterations (default 200).
#' @param subset_fraction Fraction of features sampled per iteration (default 0.5).
#' @param metric Character name of performance metric to use (e.g. "Accuracy").
#'   NULL uses the first metric returned by \code{compute_performance}.
#' @param aggregation How to aggregate multi-basis feature importance (default "mean").
#' @param ... Additional arguments (currently unused).
#' @return A \code{region_importance_result} object.
#' @export
region_importance.mvpa_model <- function(model_spec, n_iter = 200,
                                          subset_fraction = 0.5,
                                          metric = NULL, aggregation = c("mean", "sum", "maxabs"), ...) {
  aggregation <- match.arg(aggregation)
  dataset <- model_spec$dataset
  design  <- model_spec$design

  X <- get_feature_matrix(dataset)
  stopifnot(is.matrix(X))

  P <- ncol(X)
  N <- nrow(X)
  y <- y_train(design)
  stopifnot(length(y) == N)

  feature_ids <- feature_ids_for_dataset(dataset, P)
  subset_size <- max(1L, floor(P * subset_fraction))

  # Validate subset_fraction
  if (subset_fraction <= 0 || subset_fraction >= 1) {
    stop("subset_fraction must be between 0 and 1 (exclusive).")
  }

  # Create lightweight spec for parallel workers (excludes large dataset)
  mspec_lite <- as_worker_spec(model_spec)

  crossval <- model_spec$crossval

  # Parallel iterations
  iter_results <- furrr::future_map(seq_len(n_iter), function(i) {
    selected <- sort(sample.int(P, subset_size))
    X_sub <- X[, selected, drop = FALSE]
    fids_sub <- feature_ids[selected]

    perf_i <- cv_score_global(mspec_lite, X_sub, y, crossval, fids_sub, metric)

    list(perf = perf_i, included = selected)
  }, .options = furrr::furrr_options(seed = TRUE))

  # Collect performance and inclusion info
  perfs    <- vapply(iter_results, function(r) r$perf, numeric(1))
  included <- lapply(iter_results, function(r) r$included)

  # Build iteration log
  iteration_log <- tibble::tibble(
    iter        = seq_len(n_iter),
    performance = perfs,
    n_features  = vapply(included, length, integer(1))
  )

  # Pre-compute membership matrix (logical, n_iter x P)
  membership <- matrix(FALSE, nrow = n_iter, ncol = P)
  for (i in seq_len(n_iter)) {
    membership[i, included[[i]]] <- TRUE
  }

  # Aggregate per-feature
  importance  <- numeric(P)
  p_values    <- numeric(P)
  ci_lower    <- numeric(P)
  ci_upper    <- numeric(P)
  n_in_vec    <- integer(P)
  n_out_vec   <- integer(P)
  mean_in_vec <- numeric(P)
  mean_out_vec <- numeric(P)

  for (k in seq_len(P)) {
    in_mask  <- membership[, k]
    out_mask <- !in_mask

    in_perfs  <- perfs[in_mask]
    out_perfs <- perfs[out_mask]

    n_in_k  <- sum(in_mask)
    n_out_k <- sum(out_mask)
    n_in_vec[k]  <- n_in_k
    n_out_vec[k] <- n_out_k

    if (n_in_k == 0 || n_out_k == 0) {
      importance[k]  <- NA_real_
      p_values[k]    <- NA_real_
      ci_lower[k]    <- NA_real_
      ci_upper[k]    <- NA_real_
      mean_in_vec[k] <- if (n_in_k > 0) mean(in_perfs) else NA_real_
      mean_out_vec[k] <- if (n_out_k > 0) mean(out_perfs) else NA_real_
      next
    }

    mean_in  <- mean(in_perfs)
    mean_out <- mean(out_perfs)
    mean_in_vec[k]  <- mean_in
    mean_out_vec[k] <- mean_out

    importance[k] <- mean_in - mean_out

    p_values[k] <- tryCatch(
      suppressWarnings(stats::wilcox.test(in_perfs, out_perfs)$p.value),
      error = function(e) NA_real_
    )

    # Approximate CI: quantiles of per-ROI importance scores within region k
    # (not a formal bootstrap CI for mean_in - mean_out)
    diffs <- in_perfs - mean_out
    ci_lower[k] <- stats::quantile(diffs, 0.025, na.rm = TRUE)
    ci_upper[k] <- stats::quantile(diffs, 0.975, na.rm = TRUE)
  }

  # Warn if any features were never included or never excluded
  n_never_in  <- sum(n_in_vec == 0)
  n_never_out <- sum(n_out_vec == 0)
  if (n_never_in > 0) {
    warning(sprintf("%d feature(s) were never sampled; their importance is NA.", n_never_in))
  }
  if (n_never_out > 0) {
    warning(sprintf("%d feature(s) were always sampled (never excluded); their importance is NA.", n_never_out))
  }

  # Stats table
  stats_table <- tibble::tibble(
    feature_id = feature_ids,
    importance = importance,
    p_value    = p_values,
    ci_lower   = ci_lower,
    ci_upper   = ci_upper,
    n_in       = n_in_vec,
    n_out      = n_out_vec,
    mean_in    = mean_in_vec,
    mean_out   = mean_out_vec
  )

  # Build spatial maps
  importance_map <- tryCatch(
    build_output_map(dataset, importance, feature_ids, aggregation = aggregation),
    error = function(e) {
      futile.logger::flog.warn("region_importance: build_output_map failed for importance: %s", e$message)
      NULL
    }
  )

  neglog10p <- -log10(pmax(p_values, .Machine$double.xmin))
  p_value_map <- tryCatch(
    build_output_map(dataset, neglog10p, feature_ids, aggregation = aggregation),
    error = function(e) {
      futile.logger::flog.warn("region_importance: build_output_map failed for p-values: %s", e$message)
      NULL
    }
  )

  region_importance_result(
    importance            = importance,
    importance_map        = importance_map,
    p_values              = p_values,
    p_value_map           = p_value_map,
    stats_table           = stats_table,
    iteration_log         = iteration_log,
    model_spec            = model_spec,
    importance_vector_raw = importance
  )
}


# ---------- region_importance_result constructor ----------

#' Construct a Region Importance Result
#'
#' @param importance Numeric vector of per-feature importance scores.
#' @param importance_map Spatial object of importance values.
#' @param p_values Numeric vector of Wilcoxon p-values.
#' @param p_value_map Spatial object of -log10(p) values.
#' @param stats_table Tibble with per-feature statistics.
#' @param iteration_log Tibble with per-iteration performance.
#' @param model_spec The input model specification.
#' @param importance_vector_raw Raw importance vector before spatial mapping (default NULL).
#'
#' @return An S3 object of class \code{region_importance_result}.
#' @examples
#' \dontrun{
#'   # Typically created by region_importance(), not directly
#'   result <- region_importance_result(
#'     importance = c(0.1, 0.2),
#'     importance_map = NULL,
#'     p_values = c(0.05, 0.01),
#'     p_value_map = NULL,
#'     stats_table = tibble::tibble(feature_id = 1:2),
#'     iteration_log = tibble::tibble(iter = 1:10),
#'     model_spec = list()
#'   )
#' }
#' @export
region_importance_result <- function(importance, importance_map,
                                      p_values, p_value_map,
                                      stats_table, iteration_log,
                                      model_spec, importance_vector_raw = NULL) {
  structure(
    list(
      importance            = importance,
      importance_map        = importance_map,
      p_values              = p_values,
      p_value_map           = p_value_map,
      stats_table           = stats_table,
      iteration_log         = iteration_log,
      model_spec            = model_spec,
      importance_vector_raw = importance_vector_raw
    ),
    class = "region_importance_result"
  )
}


#' @export
#' @method print region_importance_result
print.region_importance_result <- function(x, ...) {
  cat("\n  Region Importance Result\n\n")

  P <- length(x$importance)
  n_valid <- sum(!is.na(x$importance))
  cat("  Features:", P, "(", n_valid, "with valid importance)\n")
  cat("  Iterations:", nrow(x$iteration_log), "\n")

  if (n_valid > 0) {
    imp <- x$importance[!is.na(x$importance)]
    cat("  Importance range: [", round(min(imp), 4), ",",
        round(max(imp), 4), "]\n")

    # Top 5 features by importance
    top_idx <- order(x$importance, decreasing = TRUE, na.last = TRUE)
    n_show <- min(5, n_valid)
    cat("  Top", n_show, "features:\n")
    for (i in seq_len(n_show)) {
      k <- top_idx[i]
      cat("    ", x$stats_table$feature_id[k], ": importance =",
          round(x$importance[k], 4),
          ", p =", format.pval(x$p_values[k], digits = 3), "\n")
    }
  }

  if (!is.null(x$importance_map)) {
    cat("  Importance map:", class(x$importance_map)[1], "\n")
  }
  if (!is.null(x$p_value_map)) {
    cat("  P-value map:", class(x$p_value_map)[1], "\n")
  }

  # Mean performance across iterations
  valid_perf <- x$iteration_log$performance[!is.na(x$iteration_log$performance)]
  if (length(valid_perf) > 0) {
    cat("  Mean iteration performance:", round(mean(valid_perf), 4),
        "(SD:", round(sd(valid_perf), 4), ")\n")
  }

  cat("\n")
}
