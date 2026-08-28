# Certified scalable solvers for the banded-ridge objective -----------------

.brs_validate_alphas <- function(alphas) {
  if (!is.numeric(alphas) || !length(alphas) || any(!is.finite(alphas)) ||
      any(alphas < 0)) {
    stop("alphas must contain finite non-negative values.", call. = FALSE)
  }
  as.numeric(alphas)
}

.brs_batch_size <- function(x, total, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0 ||
      (!is.infinite(x) && x != as.integer(x))) {
    stop(name, " must be a positive integer or Inf.", call. = FALSE)
  }
  if (is.infinite(x)) total else min(as.integer(x), total)
}

.brs_batches <- function(total, batch_size) {
  split(seq_len(total), ceiling(seq_len(total) / batch_size))
}

.brs_prepare <- function(X, Y, groups, theta) {
  gx <- .brc_grouped_matrix(X, groups = groups, context = "X")
  Y <- .brc_response_matrix(Y, n_expected = nrow(gx$X))
  if (nrow(gx$X) < 2L) {
    stop("X and Y must contain at least two training rows.", call. = FALSE)
  }
  theta <- .brc_named_group_values(theta, gx$group_names, "theta")
  if (any(!is.finite(theta)) || any(theta < 0) ||
      abs(sum(theta) - 1) > sqrt(.Machine$double.eps)) {
    stop("theta must be finite, non-negative, and sum to one.", call. = FALSE)
  }
  theta <- theta / sum(theta)
  included_group <- theta > 0
  included <- unname(included_group[gx$group])

  x_center <- colMeans(gx$X)
  x_scale <- .brc_column_sds(gx$X)
  constant_x <- !is.finite(x_scale) | x_scale == 0
  x_scale[constant_x] <- 1
  Xs <- sweep(sweep(gx$X, 2L, x_center, "-"), 2L, x_scale, "/")
  y_center <- colMeans(Y)
  Yc <- sweep(Y, 2L, y_center, "-")
  theta_feature <- unname(theta[gx$group])

  list(
    Xs = Xs,
    Yc = Yc,
    Y = Y,
    n = nrow(Xs),
    p = ncol(Xs),
    v = ncol(Y),
    theta = theta,
    theta_feature = theta_feature,
    included_group = included_group,
    included = included,
    x_center = x_center,
    x_scale = x_scale,
    y_center = y_center,
    constant_x = constant_x,
    group = gx$group,
    group_names = gx$group_names,
    feature_names = gx$feature_names,
    response_names = colnames(Y),
    layout_contract = gx$contract
  )
}

#' Create an exact-verification decomposition cache
#'
#' A caller-supplied key is never trusted by itself. Reuse requires byte-exact
#' identity of standardized training X, group membership, and theta.
#'
#' @keywords internal
#' @noRd
.banded_ridge_solver_cache <- function() {
  cache <- new.env(parent = emptyenv())
  cache$entries <- list()
  cache$decomposition_count <- 0L
  cache$hit_count <- 0L
  class(cache) <- c("banded_ridge_solver_cache", "environment")
  cache
}

.brs_cache_signature <- function(prepared) {
  list(Xs = prepared$Xs, group = prepared$group,
       theta = prepared$theta, included = prepared$included)
}

.brs_decomposition <- function(prepared,
                               backend,
                               cache,
                               cache_key,
                               solver_options,
                               builder) {
  signature <- .brs_cache_signature(prepared)
  signature$solver_options <- solver_options
  if (is.null(cache)) {
    return(list(value = builder(), cache_hit = FALSE,
                decomposition_count = 1L, cache_key = NULL))
  }
  if (!inherits(cache, "banded_ridge_solver_cache")) {
    stop("cache must be created by .banded_ridge_solver_cache().", call. = FALSE)
  }
  if (!is.character(cache_key) || length(cache_key) != 1L ||
      is.na(cache_key) || !nzchar(cache_key)) {
    stop("cache_key must be one non-empty string when cache is supplied.", call. = FALSE)
  }
  entry_key <- paste(backend, cache_key, sep = "::")
  existing <- cache$entries[[entry_key]]
  if (!is.null(existing)) {
    if (!identical(existing$signature, signature)) {
      stop("cache_key collision: training data, groups, theta, or solver options changed.",
           call. = FALSE)
    }
    cache$hit_count <- cache$hit_count + 1L
    return(list(value = existing$value, cache_hit = TRUE,
                decomposition_count = 0L, cache_key = cache_key))
  }
  value <- builder()
  cache$entries[[entry_key]] <- list(signature = signature, value = value)
  cache$decomposition_count <- cache$decomposition_count + 1L
  list(value = value, cache_hit = FALSE,
       decomposition_count = 1L, cache_key = cache_key)
}

.brs_spectral_tolerance <- function(values, dimension, tolerance = NULL) {
  scale <- max(abs(values), 1)
  if (is.null(tolerance)) {
    return(dimension * .Machine$double.eps * scale)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance < 0) {
    stop("solver_tolerance must be a finite non-negative scalar.", call. = FALSE)
  }
  tolerance
}

.brs_make_fit <- function(prepared,
                          alpha,
                          beta_standardized,
                          backend,
                          solver,
                          coefficients = TRUE,
                          dual_weights = NULL,
                          training_xs = NULL) {
  if (is.null(beta_standardized)) {
    raw_coefficients <- NULL
    intercept <- NULL
  } else {
    dimnames(beta_standardized) <- list(prepared$feature_names,
                                        prepared$response_names)
    raw_coefficients <- sweep(beta_standardized, 1L, prepared$x_scale, "/")
    dimnames(raw_coefficients) <- dimnames(beta_standardized)
    intercept <- prepared$y_center - drop(crossprod(prepared$x_center,
                                                     raw_coefficients))
    names(intercept) <- prepared$response_names
  }
  fit <- list(
    coefficients = if (coefficients) raw_coefficients else NULL,
    standardized_coefficients = if (coefficients) beta_standardized else NULL,
    intercept = intercept,
    x_center = setNames(prepared$x_center, prepared$feature_names),
    x_scale = setNames(prepared$x_scale, prepared$feature_names),
    y_center = setNames(prepared$y_center, prepared$response_names),
    constant_x = setNames(prepared$constant_x, prepared$feature_names),
    group = setNames(prepared$group, prepared$feature_names),
    group_names = prepared$group_names,
    feature_names = prepared$feature_names,
    response_names = prepared$response_names,
    alpha = alpha,
    theta = prepared$theta,
    penalty_by_group = {
      z <- rep(Inf, length(prepared$theta)); names(z) <- names(prepared$theta)
      z[prepared$theta > 0] <- alpha / prepared$theta[prepared$theta > 0]
      z
    },
    excluded_groups = names(prepared$theta)[prepared$theta == 0],
    layout_contract = prepared$layout_contract,
    dual_weights = dual_weights,
    training_xs = training_xs,
    solver = c(list(backend = backend, jitter = 0,
                    estimator_change = FALSE), solver)
  )
  class(fit) <- c("banded_ridge_optimized_fit", "list")
  fit
}

.brs_svd_decomposition <- function(prepared, solver_tolerance) {
  included <- prepared$included
  weight <- sqrt(prepared$theta_feature[included])
  Z <- sweep(prepared$Xs[, included, drop = FALSE], 2L, weight, "*")
  decomposition <- svd(Z, nu = min(dim(Z)), nv = min(dim(Z)))
  tolerance <- .brs_spectral_tolerance(
    decomposition$d, max(dim(Z)), solver_tolerance
  )
  list(U = decomposition$u, d = decomposition$d, V = decomposition$v,
       weight = weight, included = included, tolerance = tolerance,
       rank = sum(decomposition$d > tolerance),
       dimensions = c(U = length(decomposition$u),
                      d = length(decomposition$d),
                      V = length(decomposition$v), Z = length(Z)))
}

.brs_svd_fit_one <- function(prepared,
                             decomposition,
                             alpha,
                             response_batch_size,
                             recover_dual) {
  beta <- matrix(0, prepared$p, prepared$v)
  dual <- if (recover_dual) matrix(0, prepared$n, prepared$v) else NULL
  batches <- .brs_batches(prepared$v, response_batch_size)
  d <- decomposition$d
  keep <- d > decomposition$tolerance
  for (batch in batches) {
    Yb <- prepared$Yc[, batch, drop = FALSE]
    UY <- crossprod(decomposition$U, Yb)
    multiplier <- if (alpha == 0) {
      out <- numeric(length(d)); out[keep] <- 1 / d[keep]; out
    } else {
      d / (d * d + alpha)
    }
    gamma <- decomposition$V %*% sweep(UY, 1L, multiplier, "*")
    beta[decomposition$included, batch] <-
      sweep(gamma, 1L, decomposition$weight, "*")

    if (recover_dual) {
      if (alpha == 0) {
        dual_multiplier <- numeric(length(d))
        dual_multiplier[keep] <- 1 / (d[keep] * d[keep])
        dual[, batch] <- decomposition$U %*%
          sweep(UY, 1L, dual_multiplier, "*")
      } else {
        projected <- decomposition$U %*% UY
        dual[, batch] <- decomposition$U %*%
          sweep(UY, 1L, 1 / (d * d + alpha), "*") +
          (Yb - projected) / alpha
      }
    }
  }
  .brs_make_fit(
    prepared, alpha, beta, "svd_primal",
    solver = list(
      rank = decomposition$rank,
      tolerance = decomposition$tolerance,
      spectral_fallback = if (alpha == 0) "minimum_norm_pseudoinverse" else "none",
      peak_intermediate_dimensions = c(
        U = length(decomposition$U), V = length(decomposition$V),
        response_projection = length(d) * response_batch_size,
        coefficient_batch = prepared$p * response_batch_size
      )
    ),
    coefficients = TRUE, dual_weights = dual
  )
}

#' Fit an alpha path from one economy-SVD decomposition
#'
#' @keywords internal
#' @noRd
.banded_ridge_fit_svd_path <- function(X,
                                       Y,
                                       alphas,
                                       theta,
                                       groups = NULL,
                                       alpha_batch_size = Inf,
                                       response_batch_size = Inf,
                                       recover_dual = FALSE,
                                       solver_tolerance = NULL,
                                       cache = NULL,
                                       cache_key = NULL) {
  alphas <- .brs_validate_alphas(alphas)
  prepared <- .brs_prepare(X, Y, groups, theta)
  alpha_batch <- .brs_batch_size(alpha_batch_size, length(alphas),
                                 "alpha_batch_size")
  response_batch <- .brs_batch_size(response_batch_size, prepared$v,
                                    "response_batch_size")
  decomposition <- .brs_decomposition(
    prepared, "svd_primal", cache, cache_key,
    list(solver_tolerance = solver_tolerance),
    function() .brs_svd_decomposition(prepared, solver_tolerance)
  )
  fits <- vector("list", length(alphas))
  for (alpha_indices in .brs_batches(length(alphas), alpha_batch)) {
    for (ii in alpha_indices) {
      fits[[ii]] <- .brs_svd_fit_one(
        prepared, decomposition$value, alphas[[ii]], response_batch,
        isTRUE(recover_dual)
      )
      fits[[ii]]$solver$cache_hit <- decomposition$cache_hit
      fits[[ii]]$solver$cache_key <- decomposition$cache_key
    }
  }
  names(fits) <- paste0("alpha_", sprintf("%.17g", alphas))
  out <- list(
    fits = fits,
    alphas = alphas,
    theta = prepared$theta,
    backend = "svd_primal",
    decomposition_count = decomposition$decomposition_count,
    cache_hit = decomposition$cache_hit,
    cache_key = decomposition$cache_key,
    alpha_batch_size = alpha_batch_size,
    response_batch_size = response_batch_size,
    decomposition_dimensions = decomposition$value$dimensions
  )
  class(out) <- c("banded_ridge_solver_path", "list")
  out
}

.brs_dual_decomposition <- function(prepared, solver_tolerance) {
  band_kernels <- lapply(prepared$group_names, function(group_name) {
    idx <- prepared$group == group_name
    Xg <- prepared$Xs[, idx, drop = FALSE]
    tcrossprod(Xg)
  })
  names(band_kernels) <- prepared$group_names
  K <- matrix(0, prepared$n, prepared$n)
  for (group_name in prepared$group_names) {
    K <- K + prepared$theta[[group_name]] * band_kernels[[group_name]]
  }
  K <- (K + t(K)) / 2
  spectral <- eigen(K, symmetric = TRUE)
  tolerance <- .brs_spectral_tolerance(spectral$values, prepared$n,
                                       solver_tolerance)
  if (any(spectral$values < -tolerance)) {
    stop("weighted training kernel is not positive semidefinite within tolerance.",
         call. = FALSE)
  }
  values <- pmax(spectral$values, 0)
  list(
    vectors = spectral$vectors,
    values = values,
    tolerance = tolerance,
    rank = sum(values > tolerance),
    band_kernels = band_kernels,
    dimensions = c(
      weighted_kernel = length(K),
      band_kernels = sum(vapply(band_kernels, length, integer(1))),
      eigenvectors = length(spectral$vectors),
      eigenvalues = length(values)
    )
  )
}

.brs_dual_fit_one <- function(prepared,
                              decomposition,
                              alpha,
                              response_batch_size,
                              recover_primal) {
  dual <- matrix(0, prepared$n, prepared$v)
  beta <- if (recover_primal) matrix(0, prepared$p, prepared$v) else NULL
  values <- decomposition$values
  keep <- values > decomposition$tolerance
  for (batch in .brs_batches(prepared$v, response_batch_size)) {
    Yb <- prepared$Yc[, batch, drop = FALSE]
    projected <- crossprod(decomposition$vectors, Yb)
    multiplier <- if (alpha == 0) {
      out <- numeric(length(values)); out[keep] <- 1 / values[keep]; out
    } else {
      1 / (values + alpha)
    }
    weights <- decomposition$vectors %*% sweep(projected, 1L, multiplier, "*")
    dual[, batch] <- weights
    if (recover_primal) {
      beta[, batch] <- sweep(
        crossprod(prepared$Xs, weights), 1L, prepared$theta_feature, "*"
      )
    }
  }
  .brs_make_fit(
    prepared, alpha, beta, "dual_kernel",
    solver = list(
      rank = decomposition$rank,
      tolerance = decomposition$tolerance,
      spectral_fallback = if (alpha == 0) "minimum_norm_pseudoinverse" else "none",
      band_kernel_dimensions = vapply(decomposition$band_kernels, dim, integer(2)),
      peak_intermediate_dimensions = c(
        weighted_kernel = prepared$n * prepared$n,
        band_kernel = prepared$n * prepared$n,
        response_projection = prepared$n * response_batch_size,
        dual_batch = prepared$n * response_batch_size,
        primal_batch = if (recover_primal) prepared$p * response_batch_size else 0
      )
    ),
    coefficients = isTRUE(recover_primal),
    dual_weights = dual,
    training_xs = prepared$Xs
  )
}

#' Fit an alpha path from one weighted-kernel eigendecomposition
#'
#' @keywords internal
#' @noRd
.banded_ridge_fit_dual_path <- function(X,
                                        Y,
                                        alphas,
                                        theta,
                                        groups = NULL,
                                        alpha_batch_size = Inf,
                                        response_batch_size = Inf,
                                        recover_primal = TRUE,
                                        solver_tolerance = NULL,
                                        cache = NULL,
                                        cache_key = NULL) {
  alphas <- .brs_validate_alphas(alphas)
  prepared <- .brs_prepare(X, Y, groups, theta)
  alpha_batch <- .brs_batch_size(alpha_batch_size, length(alphas),
                                 "alpha_batch_size")
  response_batch <- .brs_batch_size(response_batch_size, prepared$v,
                                    "response_batch_size")
  decomposition <- .brs_decomposition(
    prepared, "dual_kernel", cache, cache_key,
    list(solver_tolerance = solver_tolerance),
    function() .brs_dual_decomposition(prepared, solver_tolerance)
  )
  fits <- vector("list", length(alphas))
  for (alpha_indices in .brs_batches(length(alphas), alpha_batch)) {
    for (ii in alpha_indices) {
      fits[[ii]] <- .brs_dual_fit_one(
        prepared, decomposition$value, alphas[[ii]], response_batch,
        isTRUE(recover_primal)
      )
      fits[[ii]]$solver$cache_hit <- decomposition$cache_hit
      fits[[ii]]$solver$cache_key <- decomposition$cache_key
    }
  }
  names(fits) <- paste0("alpha_", sprintf("%.17g", alphas))
  out <- list(
    fits = fits,
    alphas = alphas,
    theta = prepared$theta,
    backend = "dual_kernel",
    decomposition_count = decomposition$decomposition_count,
    cache_hit = decomposition$cache_hit,
    cache_key = decomposition$cache_key,
    alpha_batch_size = alpha_batch_size,
    response_batch_size = response_batch_size,
    decomposition_dimensions = decomposition$value$dimensions
  )
  class(out) <- c("banded_ridge_solver_path", "list")
  out
}

#' Predict from a certified SVD-primal or dual-kernel fit
#'
#' @keywords internal
#' @noRd
.banded_ridge_predict_optimized <- function(object, newx, groups = NULL) {
  if (!inherits(object, "banded_ridge_optimized_fit")) {
    stop("object must be a banded_ridge_optimized_fit.", call. = FALSE)
  }
  gx <- .brc_grouped_matrix(
    newx, groups = groups, expected = object$layout_contract, context = "newx"
  )
  if (object$solver$backend == "svd_primal") {
    out <- gx$X %*% object$coefficients
    out <- sweep(out, 2L, object$intercept, "+")
  } else if (object$solver$backend == "dual_kernel") {
    Xs_new <- sweep(sweep(gx$X, 2L, object$x_center, "-"),
                    2L, object$x_scale, "/")
    cross_kernel <- matrix(0, nrow(Xs_new), nrow(object$training_xs))
    for (group_name in object$group_names) {
      idx <- unname(object$group) == group_name
      cross_kernel <- cross_kernel + object$theta[[group_name]] *
        (Xs_new[, idx, drop = FALSE] %*%
           t(object$training_xs[, idx, drop = FALSE]))
    }
    out <- cross_kernel %*% object$dual_weights
    out <- sweep(out, 2L, object$y_center, "+")
  } else {
    stop("Unknown optimized solver backend.", call. = FALSE)
  }
  colnames(out) <- object$response_names
  out
}

.brs_memory_estimates <- function(n, p, v, n_bands,
                                  alpha_batch_size, response_batch_size) {
  eight <- 8
  r <- min(n, p)
  c(
    direct = eight * (p * p + p * response_batch_size),
    svd_primal = eight * (n * r + p * r + r +
                            r * response_batch_size +
                            p * response_batch_size * alpha_batch_size),
    dual_kernel = eight * (n * n * (n_bands + 2) +
                             n * response_batch_size * alpha_batch_size +
                             p * response_batch_size)
  )
}

#' Plan a solver without exceeding its estimated intermediate storage contract
#'
#' @keywords internal
#' @noRd
.banded_ridge_solver_plan <- function(n,
                                      p,
                                      v,
                                      n_bands,
                                      n_alphas = 1L,
                                      alpha_batch_size = 1L,
                                      response_batch_size = v,
                                      memory_limit_mb = 1024,
                                      solver = c("auto", "direct", "svd_primal",
                                                 "dual_kernel")) {
  solver <- match.arg(solver)
  values <- c(n = n, p = p, v = v, n_bands = n_bands,
              n_alphas = n_alphas, alpha_batch_size = alpha_batch_size,
              response_batch_size = response_batch_size)
  if (any(!is.numeric(values)) || any(!is.finite(values)) ||
      any(values <= 0) || any(values != as.integer(values))) {
    stop("shape, band, alpha, and batch counts must be positive integers.",
         call. = FALSE)
  }
  if (!is.numeric(memory_limit_mb) || length(memory_limit_mb) != 1L ||
      !is.finite(memory_limit_mb) || memory_limit_mb <= 0) {
    stop("memory_limit_mb must be one positive finite value.", call. = FALSE)
  }
  estimates <- .brs_memory_estimates(
    as.integer(n), as.integer(p), as.integer(v), as.integer(n_bands),
    min(as.integer(alpha_batch_size), as.integer(n_alphas)),
    min(as.integer(response_batch_size), as.integer(v))
  )
  limit <- memory_limit_mb * 1024^2
  fits <- estimates <= limit
  if (solver != "auto") {
    chosen <- solver
    if (!fits[[chosen]]) {
      stop("Requested solver exceeds the configured estimated storage contract.",
           call. = FALSE)
    }
    reason <- "explicit override within memory contract"
  } else {
    if (p <= 64L && fits[["direct"]]) {
      chosen <- "direct"
      reason <- "measured small-p direct threshold p <= 64 and storage fits"
    } else if (p >= 4L * n && fits[["dual_kernel"]]) {
      chosen <- "dual_kernel"
      reason <- "conservative p >= 4n dual threshold and storage fits"
    } else if (n_alphas > 1L && fits[["svd_primal"]]) {
      chosen <- "svd_primal"
      reason <- "alpha-path SVD reuse and storage fits"
    } else if (fits[["direct"]]) {
      chosen <- "direct"
      reason <- "small-path direct reference and storage fits"
    } else {
      available <- names(estimates)[fits]
      if (!length(available)) {
        stop("No solver fits the configured estimated storage contract.",
             call. = FALSE)
      }
      chosen <- available[[which.min(estimates[available])]]
      reason <- "smallest estimated allocation among feasible solvers"
    }
  }
  list(
    solver = chosen,
    reason = reason,
    estimated_bytes = estimates,
    estimated_mib = estimates / 1024^2,
    memory_limit_mb = memory_limit_mb,
    fits_memory_contract = fits,
    heuristic_version = "br-auto-v2-measured-p64-conservative-p4n"
  )
}
