# Pure banded-ridge reference core -------------------------------------------

# The helpers in this file deliberately have no dependency on rMVPA dataset,
# ROI, searchlight, or spatial classes.  They define the numerical contract
# against which optimized and integrated implementations are compared.
#
# For training rows, let Xs be X centered and divided columnwise by the sample
# standard deviation, and let Yc be Y centered but not scaled.  For feature
# groups g, the estimator minimizes
#
#   ||Yc - Xs B||_F^2 + sum_g lambda_g ||B_g||_F^2.
#
# The normal equations are therefore (crossprod(Xs) + diag(lambda)) B =
# crossprod(Xs, Yc), with no hidden division by n.  The alternative
# parameterization is lambda_g = alpha / theta_g, where named theta values lie
# on the non-negative simplex.  theta_g = 0 is represented by excluding that
# block from the solve (the exact lambda_g -> Inf limit); its coefficients are
# identically zero.  Centering and scaling are learned only on training rows.

.brc_response_matrix <- function(Y, n_expected = NULL, context = "Y") {
  if (is.numeric(Y) && is.null(dim(Y))) {
    Y <- matrix(Y, ncol = 1L)
  }

  if (!is.matrix(Y) || !is.numeric(Y)) {
    stop(context, " must be a numeric matrix or numeric vector.", call. = FALSE)
  }
  if (nrow(Y) < 1L || ncol(Y) < 1L) {
    stop(context, " must have at least one row and one response column.", call. = FALSE)
  }
  if (!is.null(n_expected) && nrow(Y) != n_expected) {
    stop(context, " must have the same number of rows as X.", call. = FALSE)
  }
  if (any(!is.finite(Y))) {
    stop(context, " must contain only finite values.", call. = FALSE)
  }

  response_names <- colnames(Y)
  if (is.null(response_names)) {
    response_names <- paste0("response_", seq_len(ncol(Y)))
  } else {
    if (anyNA(response_names) || any(!nzchar(response_names))) {
      stop(context, " column names must be non-empty and non-missing when supplied.", call. = FALSE)
    }
    response_names <- make.unique(response_names)
  }
  colnames(Y) <- response_names
  Y
}

.brc_numeric_matrix <- function(X, context, n_expected = NULL) {
  if (!is.matrix(X) || !is.numeric(X)) {
    stop(context, " must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(X) < 1L || ncol(X) < 1L) {
    stop(context, " must have at least one row and one column.", call. = FALSE)
  }
  if (!is.null(n_expected) && nrow(X) != n_expected) {
    stop(context, " must have the same number of rows as the other feature bands.", call. = FALSE)
  }
  if (any(!is.finite(X))) {
    stop(context, " must contain only finite values.", call. = FALSE)
  }
  X
}

.brc_validate_group_names <- function(x, context = "groups") {
  if (anyNA(x) || any(!nzchar(x))) {
    stop(context, " must contain non-empty, non-missing names.", call. = FALSE)
  }
  x
}

.brc_feature_names <- function(X, group) {
  raw <- colnames(X)
  if (is.null(raw)) {
    raw <- stats::ave(seq_along(group), group, FUN = seq_along)
  } else {
    raw <- as.character(raw)
    bad <- is.na(raw) | !nzchar(raw)
    if (any(bad)) {
      raw[bad] <- stats::ave(seq_along(group), group, FUN = seq_along)[bad]
    }
  }
  make.unique(paste0(group, "::", raw))
}

.brc_grouped_matrix <- function(X,
                                groups = NULL,
                                expected = NULL,
                                context = "X") {
  expected_layout <- if (is.null(expected)) NULL else expected$layout

  if (is.data.frame(X)) {
    stop(context, " must be a numeric matrix or a non-empty named list of numeric matrices.",
         call. = FALSE)
  }

  if (is.matrix(X)) {
    if (!is.null(expected_layout) && expected_layout != "matrix") {
      stop(context, " must use the same matrix/list layout as the fitted data.", call. = FALSE)
    }
    X <- .brc_numeric_matrix(X, context)

    if (is.null(groups)) {
      group <- if (is.null(expected)) rep("all", ncol(X)) else expected$group
    } else {
      if (!(is.character(groups) || is.factor(groups)) || length(groups) != ncol(X)) {
        stop("groups must be a character/factor vector with one value per X column.", call. = FALSE)
      }
      group <- .brc_validate_group_names(as.character(groups))
    }

    if (length(group) != ncol(X)) {
      stop(context, " has a different number of columns than the fitted data.", call. = FALSE)
    }
    group_names <- unique(group)
    list_sizes <- as.integer(tabulate(match(group, group_names), nbins = length(group_names)))
    names(list_sizes) <- group_names

    if (!is.null(expected)) {
      if (!identical(group, expected$group)) {
        stop(context, " group membership/order does not match the fitted data.", call. = FALSE)
      }
      if (!is.null(expected$input_colnames) &&
          !identical(colnames(X), expected$input_colnames)) {
        stop(context, " column names/order do not match the fitted data.", call. = FALSE)
      }
      feature_names <- expected$feature_names
    } else {
      feature_names <- .brc_feature_names(X, group)
    }

    return(list(
      X = X,
      group = group,
      group_names = group_names,
      feature_names = feature_names,
      contract = list(
        layout = "matrix",
        group = group,
        group_names = group_names,
        list_sizes = list_sizes,
        input_colnames = colnames(X),
        list_colnames = NULL,
        feature_names = feature_names
      )
    ))
  }

  if (!is.list(X) || length(X) < 1L) {
    stop(context, " must be a numeric matrix or a non-empty named list of numeric matrices.", call. = FALSE)
  }
  if (!is.null(expected_layout) && expected_layout != "list") {
    stop(context, " must use the same matrix/list layout as the fitted data.", call. = FALSE)
  }
  if (!is.null(groups)) {
    stop("groups must be NULL when X is a named list; list names define the bands.", call. = FALSE)
  }

  band_names <- names(X)
  if (is.null(band_names)) {
    stop(context, " list must be named with one unique name per feature band.", call. = FALSE)
  }
  band_names <- .brc_validate_group_names(as.character(band_names), paste0(context, " list names"))
  if (anyDuplicated(band_names)) {
    stop(context, " list band names must be unique.", call. = FALSE)
  }

  if (!is.null(expected)) {
    if (!setequal(band_names, expected$group_names)) {
      stop(context, " list band names do not match the fitted data.", call. = FALSE)
    }
    X <- X[expected$group_names]
    band_names <- names(X)
  }

  mats <- vector("list", length(X))
  n_rows <- NULL
  input_colnames <- vector("list", length(X))
  names(input_colnames) <- band_names
  for (ii in seq_along(X)) {
    nm <- band_names[[ii]]
    mats[[ii]] <- .brc_numeric_matrix(
      X[[ii]],
      paste0(context, "[['", nm, "']]"),
      n_expected = n_rows
    )
    if (is.null(n_rows)) n_rows <- nrow(mats[[ii]])
    input_colnames[[ii]] <- colnames(mats[[ii]])
  }

  list_sizes <- vapply(mats, ncol, integer(1))
  names(list_sizes) <- band_names
  if (!is.null(expected) && !identical(list_sizes, expected$list_sizes)) {
    stop(context, " list band column counts do not match the fitted data.", call. = FALSE)
  }
  if (!is.null(expected)) {
    for (nm in band_names) {
      expected_names <- expected$list_colnames[[nm]]
      if (!is.null(expected_names) && !identical(input_colnames[[nm]], expected_names)) {
        stop(context, " column names/order for band '", nm,
             "' do not match the fitted data.", call. = FALSE)
      }
    }
  }

  Xmat <- do.call(cbind, mats)
  group <- rep(band_names, times = list_sizes)
  feature_names <- if (is.null(expected)) {
    .brc_feature_names(Xmat, group)
  } else {
    expected$feature_names
  }

  list(
    X = Xmat,
    group = group,
    group_names = band_names,
    feature_names = feature_names,
    contract = list(
      layout = "list",
      group = group,
      group_names = band_names,
      list_sizes = list_sizes,
      input_colnames = NULL,
      list_colnames = input_colnames,
      feature_names = feature_names
    )
  )
}

.brc_named_group_values <- function(x, group_names, name) {
  if (!is.numeric(x) || is.null(names(x)) || length(x) != length(group_names)) {
    stop(name, " must be a named numeric vector with one value per feature band.", call. = FALSE)
  }
  if (anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x)) ||
      !setequal(names(x), group_names)) {
    stop(name, " names must match the feature-band names exactly.", call. = FALSE)
  }
  x <- as.numeric(x[group_names])
  names(x) <- group_names
  x
}

.brc_penalty_spec <- function(group,
                              group_names,
                              lambdas = NULL,
                              alpha = NULL,
                              theta = NULL,
                              simplex_tol = sqrt(.Machine$double.eps)) {
  lambda_mode <- !is.null(lambdas)
  alpha_theta_mode <- !is.null(alpha) || !is.null(theta)
  if (lambda_mode == alpha_theta_mode) {
    stop("Supply exactly one penalty form: `lambdas`, or both `alpha` and `theta`.", call. = FALSE)
  }

  if (lambda_mode) {
    lambdas <- .brc_named_group_values(lambdas, group_names, "lambdas")
    if (any(!is.finite(lambdas)) || any(lambdas < 0)) {
      stop("lambdas must contain finite non-negative values.", call. = FALSE)
    }
    included_by_group <- setNames(rep(TRUE, length(group_names)), group_names)
    theta_out <- NULL
    alpha_out <- NULL
  } else {
    if (is.null(alpha) || is.null(theta)) {
      stop("Both `alpha` and `theta` are required for the alpha/theta penalty form.", call. = FALSE)
    }
    if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha < 0) {
      stop("alpha must be a finite non-negative scalar.", call. = FALSE)
    }
    theta <- .brc_named_group_values(theta, group_names, "theta")
    if (any(!is.finite(theta)) || any(theta < 0)) {
      stop("theta must contain finite non-negative values.", call. = FALSE)
    }
    theta_sum <- sum(theta)
    if (!is.finite(theta_sum) || abs(theta_sum - 1) > simplex_tol) {
      stop("theta must sum to one within the simplex tolerance.", call. = FALSE)
    }
    theta <- theta / theta_sum
    included_by_group <- theta > 0
    lambdas <- rep(Inf, length(group_names))
    names(lambdas) <- group_names
    lambdas[included_by_group] <- alpha / theta[included_by_group]
    if (any(!is.finite(lambdas[included_by_group]))) {
      stop("alpha/theta produced a non-finite penalty; positive theta values are too small.", call. = FALSE)
    }
    theta_out <- theta
    alpha_out <- unname(alpha)
  }

  included <- unname(included_by_group[group])
  penalty <- unname(lambdas[group])
  list(
    lambdas = lambdas,
    alpha = alpha_out,
    theta = theta_out,
    included_by_group = included_by_group,
    included = included,
    penalty = penalty
  )
}

#' Column-wise centering and scaling without sweep()
#'
#' Elementwise `X[i, j] - center[j]` and `/ scale[j]` written with recycled
#' vectors. This is the same IEEE arithmetic as `sweep()` (bit-identical
#' results, dimnames preserved) without sweep's transposed intermediate arrays,
#' and runs in roughly half the time on n x p matrices. The statistics are
#' unnamed first: `rep()` would otherwise replicate feature names n times.
#'
#' @keywords internal
#' @noRd
.brc_center_columns <- function(X, center) {
  X - rep(unname(center), each = nrow(X))
}

.brc_standardize_columns <- function(X, center, scale) {
  n <- nrow(X)
  (X - rep(unname(center), each = n)) / rep(unname(scale), each = n)
}

.brc_offset_columns <- function(X, offset) {
  X + rep(unname(offset), each = nrow(X))
}

.brc_column_sds <- function(X) {
  n <- nrow(X)
  centered <- .brc_center_columns(X, colMeans(X))
  sqrt(colSums(centered * centered) / (n - 1L))
}

#' Training-split column centering and scaling used by every banded-ridge fit
#'
#' Constant columns keep a unit scale so they contribute zero after centering.
#' All solver paths and the tuning preprocessing receipt call this helper so the
#' receipt is exactly what the fits use.
#'
#' @keywords internal
#' @noRd
.brc_column_standardization <- function(X) {
  center <- colMeans(X)
  scale <- .brc_column_sds(X)
  constant <- !is.finite(scale) | scale == 0
  scale[constant] <- 1
  list(center = center, scale = scale, constant = constant)
}

.brc_symmetric_solve <- function(A, B, tolerance = NULL) {
  A <- (A + t(A)) / 2
  chol_A <- tryCatch(chol(A), error = function(e) NULL)
  if (!is.null(chol_A)) {
    solution <- backsolve(chol_A, forwardsolve(t(chol_A), B))
    return(list(
      solution = solution,
      method = "direct_cholesky",
      rank = ncol(A),
      tolerance = 0,
      eigenvalues = NULL
    ))
  }

  ee <- eigen(A, symmetric = TRUE)
  scale <- max(abs(ee$values), 1)
  if (is.null(tolerance)) {
    tolerance <- max(dim(A)) * .Machine$double.eps * scale
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance < 0) {
    stop("solver tolerance must be a finite non-negative scalar.", call. = FALSE)
  }
  if (any(ee$values < -tolerance)) {
    stop("penalized normal equations are not positive semidefinite within tolerance.", call. = FALSE)
  }
  keep <- ee$values > tolerance
  if (any(keep)) {
    projected <- crossprod(ee$vectors[, keep, drop = FALSE], B)
    projected <- projected / ee$values[keep]
    solution <- ee$vectors[, keep, drop = FALSE] %*% projected
  } else {
    solution <- matrix(0, nrow(A), ncol(B))
  }
  list(
    solution = solution,
    method = "direct_symmetric_pseudoinverse",
    rank = sum(keep),
    tolerance = tolerance,
    eigenvalues = ee$values
  )
}

#' Fit the direct reference banded-ridge estimator
#'
#' This intentionally slow, transparent backend implements the objective
#' described at the top of this file. It is the numerical oracle for optimized
#' backends, not the scalable execution path.
#'
#' @keywords internal
#' @noRd
.banded_ridge_fit_direct <- function(X,
                                     Y,
                                     groups = NULL,
                                     lambdas = NULL,
                                     alpha = NULL,
                                     theta = NULL,
                                     solver_tolerance = NULL) {
  gx <- .brc_grouped_matrix(X, groups = groups, context = "X")
  Y <- .brc_response_matrix(Y, n_expected = nrow(gx$X))
  if (nrow(gx$X) < 2L) {
    stop("X and Y must contain at least two training rows.", call. = FALSE)
  }

  penalty <- .brc_penalty_spec(
    group = gx$group,
    group_names = gx$group_names,
    lambdas = lambdas,
    alpha = alpha,
    theta = theta
  )

  standardization <- .brc_column_standardization(gx$X)
  x_center <- standardization$center
  x_scale <- standardization$scale
  constant_x <- standardization$constant
  Xs <- .brc_standardize_columns(gx$X, x_center, x_scale)

  y_center <- colMeans(Y)
  Ys <- .brc_center_columns(Y, y_center)

  included <- penalty$included
  Xi <- Xs[, included, drop = FALSE]
  pi <- penalty$penalty[included]
  A <- crossprod(Xi) + diag(pi, nrow = length(pi), ncol = length(pi))
  B <- crossprod(Xi, Ys)
  solved <- .brc_symmetric_solve(A, B, tolerance = solver_tolerance)

  beta_standardized <- matrix(0, ncol(gx$X), ncol(Y))
  beta_standardized[included, ] <- solved$solution
  dimnames(beta_standardized) <- list(gx$feature_names, colnames(Y))

  coefficients <- beta_standardized / x_scale
  dimnames(coefficients) <- dimnames(beta_standardized)
  intercept <- y_center - drop(crossprod(x_center, coefficients))
  names(intercept) <- colnames(Y)

  fit <- list(
    coefficients = coefficients,
    standardized_coefficients = beta_standardized,
    intercept = intercept,
    x_center = setNames(x_center, gx$feature_names),
    x_scale = setNames(x_scale, gx$feature_names),
    y_center = setNames(y_center, colnames(Y)),
    constant_x = setNames(constant_x, gx$feature_names),
    group = setNames(gx$group, gx$feature_names),
    group_names = gx$group_names,
    feature_names = gx$feature_names,
    response_names = colnames(Y),
    penalty_by_group = penalty$lambdas,
    penalty_by_feature = setNames(penalty$penalty, gx$feature_names),
    alpha = penalty$alpha,
    theta = penalty$theta,
    excluded_groups = names(penalty$included_by_group)[!penalty$included_by_group],
    layout_contract = gx$contract,
    solver = list(
      backend = "direct",
      linear_solver = solved$method,
      rank = solved$rank,
      tolerance = solved$tolerance,
      n = nrow(gx$X),
      p = ncol(gx$X),
      v = ncol(Y),
      included_p = sum(included)
    )
  )
  class(fit) <- c("banded_ridge_direct_fit", "list")
  fit
}

#' Predict from the direct reference banded-ridge estimator
#'
#' @keywords internal
#' @noRd
.banded_ridge_predict_direct <- function(object, newx, groups = NULL) {
  if (!inherits(object, "banded_ridge_direct_fit")) {
    stop("object must be a banded_ridge_direct_fit.", call. = FALSE)
  }
  gx <- .brc_grouped_matrix(
    newx,
    groups = groups,
    expected = object$layout_contract,
    context = "newx"
  )
  out <- gx$X %*% object$coefficients
  out <- .brc_offset_columns(out, object$intercept)
  colnames(out) <- object$response_names
  out
}

.brc_split_indices <- function(x, n, name) {
  if (!is.numeric(x) || length(x) < 1L || any(!is.finite(x)) ||
      any(x != as.integer(x))) {
    stop(name, " must be a non-empty vector of finite integer row indices.", call. = FALSE)
  }
  x <- as.integer(x)
  if (any(x < 1L | x > n)) {
    stop(name, " contains out-of-range row indices.", call. = FALSE)
  }
  if (anyDuplicated(x)) {
    stop(name, " must not contain duplicated row indices.", call. = FALSE)
  }
  x
}

.brc_subset_rows <- function(X, rows) {
  if (is.matrix(X)) {
    return(X[rows, , drop = FALSE])
  }
  lapply(X, function(x) x[rows, , drop = FALSE])
}

#' Fit on explicit training rows and predict explicit test rows
#'
#' @keywords internal
#' @noRd
.banded_ridge_fit_predict_direct <- function(X,
                                             Y,
                                             train,
                                             test,
                                             groups = NULL,
                                             lambdas = NULL,
                                             alpha = NULL,
                                             theta = NULL,
                                             solver_tolerance = NULL) {
  n <- if (is.matrix(X)) nrow(X) else if (is.list(X) && length(X)) nrow(X[[1L]]) else 0L
  train <- .brc_split_indices(train, n, "train")
  test <- .brc_split_indices(test, n, "test")
  if (length(intersect(train, test))) {
    stop("train and test row indices must be disjoint.", call. = FALSE)
  }
  if (length(train) < 2L) {
    stop("train must contain at least two rows.", call. = FALSE)
  }
  Y <- .brc_response_matrix(Y, n_expected = n)
  fit <- .banded_ridge_fit_direct(
    X = .brc_subset_rows(X, train),
    Y = Y[train, , drop = FALSE],
    groups = groups,
    lambdas = lambdas,
    alpha = alpha,
    theta = theta,
    solver_tolerance = solver_tolerance
  )
  predictions <- .banded_ridge_predict_direct(
    fit,
    .brc_subset_rows(X, test),
    groups = groups
  )
  rownames(predictions) <- as.character(test)
  list(fit = fit, predictions = predictions, train = train, test = test)
}

#' Score finite multi-response predictions
#'
#' @keywords internal
#' @noRd
.banded_ridge_score <- function(observed, predicted) {
  observed <- .brc_response_matrix(observed, context = "observed")
  predicted <- .brc_response_matrix(predicted, context = "predicted")
  if (!identical(dim(observed), dim(predicted))) {
    stop("observed and predicted must have identical dimensions.", call. = FALSE)
  }
  if (!identical(colnames(observed), colnames(predicted))) {
    stop("observed and predicted response names/order must match.", call. = FALSE)
  }

  residual <- observed - predicted
  mse <- colMeans(residual * residual)
  observed_centered <- .brc_center_columns(observed, colMeans(observed))
  predicted_centered <- .brc_center_columns(predicted, colMeans(predicted))
  ss_tot <- colSums(observed_centered * observed_centered)
  pred_ss <- colSums(predicted_centered * predicted_centered)
  denom <- sqrt(ss_tot * pred_ss)
  correlation <- rep(NA_real_, ncol(observed))
  ok_cor <- nrow(observed) >= 2L & is.finite(denom) & denom > 0
  correlation[ok_cor] <- colSums(
    observed_centered[, ok_cor, drop = FALSE] *
      predicted_centered[, ok_cor, drop = FALSE]
  ) / denom[ok_cor]
  r2 <- rep(NA_real_, ncol(observed))
  ok_r2 <- is.finite(ss_tot) & ss_tot > 0
  r2[ok_r2] <- 1 - colSums(residual[, ok_r2, drop = FALSE]^2) / ss_tot[ok_r2]

  data.frame(
    response = colnames(observed),
    n_obs = nrow(observed),
    mse = unname(mse),
    correlation = correlation,
    r2 = r2,
    stringsAsFactors = FALSE
  )
}
