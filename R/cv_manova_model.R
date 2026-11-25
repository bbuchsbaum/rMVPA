#' Cross-validated MANOVA pattern distinctness
#'
#' Estimate a signed pattern distinctness (D) for each specified effect across
#' independent runs, using the leave-one-run-out estimator of Allefeld and
#' Haynes (2014), doi:10.1016/j.neuroimage.2013.11.043. Run the model with
#' [run_regional()] or [run_searchlight()]. Classical MANOVA is available
#' separately through [manova_model()].
#'
#' @param dataset An [mvpa_dataset()] with rows in design order. For time
#'   series, supply data and regressors already filtered and temporally
#'   whitened together. The estimator assumes independent residual observations
#'   within runs after preprocessing and independence between runs.
#' @param design A [manova_design()] with a one-sided formula. A common
#'   model matrix defines the coefficient coordinates. Coefficients, including
#'   all nuisance effects, are estimated separately within each run. Unordered
#'   factors default to sum coding unless they carry explicit contrasts.
#' @param runs Vector identifying the independent acquisition run of each row.
#' @param effects Formula term names to score (default: all non-intercept
#'   terms), or a named list of coefficient contrasts. Each list entry is a
#'   vector or matrix with one column per design coefficient and one contrast
#'   per row. Named coefficients are aligned by name. A multi-row contrast
#'   tests their joint row space. All tested contrasts must be estimable in
#'   every run. Other formula terms remain in the GLM as nuisance variables.
#' @param precision Spatial noise model. `"full"` uses the inverse residual
#'   scatter pooled across training runs and requires more than p + 1 training
#'   residual degrees of freedom for p voxels. `"diagonal"` sums univariate
#'   distinctness over voxels, with the corresponding univariate inverse-
#'   variance correction and a requirement of more than two training residual
#'   degrees of freedom. This changes the spatial metric.
#' @param residual_df Optional numeric vector of residual degrees of freedom
#'   per run after temporal filtering, in order of first occurrence of each
#'   run (or named by run). Defaults to run size minus design rank.
#' @return A model specification. Regional and searchlight results contain one
#'   scalar per named effect. The score is D summed over voxels and contrast
#'   dimensions. A negative observation is valid. Scores are neither p-values
#'   nor classification accuracies. Divide by voxel count for per-voxel units.
#' @examples
#' d <- gen_sample_dataset(c(3, 3, 3), nobs = 48, nlevels = 2, blocks = 4)
#' trials <- expand.grid(place = factor(c("front", "back")),
#'                       vowel = factor(c("a", "i", "u")), rep = 1:2, run = 1:4)
#' model <- cv_manova_model(d$dataset, manova_design(~ place + vowel, trials),
#'                          runs = trials$run, precision = "diagonal")
#' @md
#' @export
cv_manova_model <- function(dataset, design, runs, effects = NULL,
                            precision = c("full", "diagonal"), residual_df = NULL) {
  precision <- match.arg(precision)
  if (!inherits(design, "manova_design") || length(design$formula) != 2L) {
    stop("Supply a manova_design with a one-sided formula.")
  }
  frame <- stats::model.frame(design$formula, as.data.frame(design$data),
                              na.action = stats::na.fail)
  frame[] <- lapply(frame, function(x) if (is.character(x)) factor(x) else x)
  use_sum <- vapply(frame, function(x) {
    is.factor(x) && !is.ordered(x) && is.null(attr(x, "contrasts"))
  }, logical(1))
  coding <- lapply(frame[use_sum], function(x) stats::contr.sum(levels(x)))
  X <- stats::model.matrix(attr(frame, "terms"), frame,
                           contrasts.arg = if (length(coding)) coding else NULL)
  n <- nrow(X)
  if (n != nobs(dataset) || !all(is.finite(X))) {
    stop("Design rows must match dataset observations and be finite.")
  }
  if (!is.atomic(runs) || !is.null(dim(runs)) || length(runs) != n ||
      anyNA(runs) || length(unique(runs)) < 2L) {
    stop("Supply one run ID per observation and at least two independent runs.")
  }
  design$block_var <- factor(runs, levels = unique(runs))
  cv <- blocked_cross_validation(design$block_var)
  groups <- split(seq_len(n), cv$block_var)
  bases <- .cv_manova_effects(X, attr(frame, "terms"), effects)
  plans <- lapply(groups, function(rows) .cv_manova_run(X, rows, bases))
  df <- vapply(plans, function(z) length(z$rows) - ncol(z$u), integer(1))
  if (!is.null(residual_df)) {
    if (!is.null(names(residual_df))) {
      if (anyDuplicated(names(residual_df)) || !setequal(names(residual_df), names(df))) {
        stop("Named residual_df must identify each run exactly once.")
      }
      residual_df <- residual_df[names(df)]
    }
    if (!is.numeric(residual_df) || length(residual_df) != length(df) ||
        any(!is.finite(residual_df)) || any(residual_df > df) || any(residual_df < 0)) {
      stop("residual_df must be between zero and the available degrees of freedom in each run.")
    }
    df[] <- residual_df
  }
  create_model_spec("cv_manova_model", dataset, design, crossval = cv,
                    run_designs = plans, residual_df = df, precision = precision,
                    metric_names = names(bases), compute_performance = TRUE)
}

.cv_manova_effects <- function(X, terms, effects) {
  labels <- attr(terms, "term.labels")
  if (is.null(effects)) effects <- labels
  if (is.character(effects)) {
    if (!all(effects %in% labels)) stop("Unknown effect in design formula.")
    assign <- attr(X, "assign")
    bases <- lapply(match(effects, labels), function(i) {
      diag(ncol(X))[, assign == i, drop = FALSE]
    })
    names(bases) <- effects
  } else {
    if (!is.list(effects)) stop("effects must be term names or a named list of contrasts.")
    bases <- lapply(effects, function(C) {
      if (is.numeric(C) && is.null(dim(C))) C <- rbind(C)
      if (!is.matrix(C) || !is.numeric(C) || !nrow(C) || ncol(C) != ncol(X) ||
          any(!is.finite(C))) stop("Contrasts must be finite with one column per coefficient.")
      if (!is.null(colnames(C))) {
        if (!setequal(colnames(C), colnames(X)) || anyDuplicated(colnames(C))) {
          stop("Contrast column names must match the design coefficients.")
        }
        C <- C[, colnames(X), drop = FALSE]
      }
      q <- qr(t(C))
      if (!q$rank) stop("Contrast has rank zero.")
      qr.Q(q)[, seq_len(q$rank), drop = FALSE]
    })
  }
  nm <- names(bases)
  if (!length(bases) || is.null(nm) || anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) {
    stop("Specify at least one effect, with unique nonempty names.")
  }
  bases
}

.cv_manova_run <- function(X, rows, bases) {
  x <- X[rows, , drop = FALSE]
  dec <- svd(x)
  keep <- dec$d > max(dim(x)) * .Machine$double.eps * max(dec$d)
  u <- dec$u[, keep, drop = FALSE]
  v <- dec$v[, keep, drop = FALSE]
  inverse <- tcrossprod(sweep(v, 2L, dec$d[keep], "/"), u)
  effects <- lapply(names(bases), function(name) {
    b <- bases[[name]]
    if (max(abs(b - v %*% crossprod(v, b))) > 1e-7) {
      stop("Effect '", name, "' is not estimable in the run starting at row ", rows[1L], ".")
    }
    list(weights = crossprod(b, inverse), gram = crossprod(x %*% b))
  })
  list(rows = rows, u = u, effects = effects)
}

.cv_manova_scores <- function(model, Y) {
  Y <- as.matrix(Y)
  if (nrow(Y) != length(model$crossval$block_var) || !ncol(Y) || any(!is.finite(Y))) {
    stop("ROI data must be finite with one row per observation and at least one voxel.")
  }
  diagonal <- model$precision == "diagonal"
  p <- ncol(Y)
  fold_df <- sum(model$residual_df) - model$residual_df
  correction <- fold_df - (if (diagonal) 1L else p) - 1L
  if (any(correction <= 0)) stop("Insufficient training residual degrees of freedom for the requested precision.")
  fits <- lapply(model$run_designs, function(run) {
    y <- Y[run$rows, , drop = FALSE]
    e <- y - run$u %*% crossprod(run$u, y)
    # Exact fits can leave projection roundoff instead of a zero residual.
    floor <- (nrow(y) * .Machine$double.eps)^2 * colSums(y^2)
    e[, colSums(e^2) <= floor] <- 0
    list(beta = lapply(run$effects, function(h) h$weights %*% y),
         scatter = if (diagonal) colSums(e^2) else crossprod(e))
  })
  n <- vapply(model$run_designs, function(run) length(run$rows), integer(1))
  scores <- numeric(length(model$metric_names))
  for (l in seq_along(fits)) {
    train <- fits[-l]
    scatter <- Reduce("+", lapply(train, "[[", "scatter"))
    if (diagonal) {
      if (any(scatter <= 0)) stop("Zero training residual variance.")
    } else {
      chol_scatter <- tryCatch(chol(scatter), error = function(e) NULL)
      if (is.null(chol_scatter) || rcond(chol_scatter) < sqrt(.Machine$double.eps)) {
        stop("Singular training residual covariance.")
      }
    }
    for (i in seq_along(scores)) {
      train_beta <- Reduce("+", lapply(train, function(f) f$beta[[i]]))
      test_beta <- model$run_designs[[l]]$effects[[i]]$gram %*% fits[[l]]$beta[[i]]
      if (diagonal) {
        product <- sum(colSums(train_beta * test_beta) / scatter)
      } else {
        # For E = R'R, whiten both contrast patterns by R^{-T}.
        a <- backsolve(chol_scatter, t(train_beta), transpose = TRUE)
        b <- backsolve(chol_scatter, t(test_beta), transpose = TRUE)
        product <- sum(a * b)
      }
      scores[i] <- scores[i] + correction[l] * product / sum(n[-l])
    }
  }
  stats::setNames(scores / length(fits), model$metric_names)
}

#' @export
output_schema.cv_manova_model <- function(model) {
  stats::setNames(as.list(rep("scalar", length(model$metric_names))), model$metric_names)
}

#' @export
fit_roi.cv_manova_model <- function(model, roi_data, context, ...) {
  tryCatch({
    scores <- .cv_manova_scores(model, roi_data$train_data)
    roi_result(metrics = scores, indices = roi_data$indices, id = context$id)
  }, error = function(e) {
    roi_result(metrics = NULL, indices = roi_data$indices, id = context$id,
               error = TRUE, error_message = conditionMessage(e))
  })
}

# MANOVA models a multivariate response and has no classification labels.
#' @export
y_train.manova_design <- function(obj) NULL
