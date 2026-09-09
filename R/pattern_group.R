# Group summaries use explicit target and spatial identities. No component
# matching, interpolation of SE maps, or voxelwise covariance is inferred.

.pattern_group_basis <- function(x) {
  if (!inherits(x, "pattern_basis")) x <- pattern_basis(x)
  B <- x$matrix
  if (!is.matrix(B) || !is.numeric(B) || any(!is.finite(B)) ||
      ncol(B) < 1L || ncol(B) > nrow(B) ||
      min(svd(B, nu = 0, nv = 0)$d) <= 1e-9 * max(svd(B, nu = 0, nv = 0)$d))
    stop("reference_basis must have a finite full-rank target matrix.", call. = FALSE)
  x$target_ids <- .pattern_ids(x$target_ids, nrow(B), "basis target_ids")
  # Group identity concerns the actual raw coordinates, including a caller's
  # explicit reference rotation; never reuse a stale supplied hash.
  x$source_basis_id <- x$source_basis_id %||% x$basis_id
  order <- order(x$target_ids)
  x$basis_id <- digest::digest(list(matrix = unname(B[order, , drop = FALSE]),
                                  target_ids = x$target_ids[order], type = x$type))
  x
}

.pattern_group_transform <- function(subject, reference) {
  b <- .pattern_group_basis(subject$basis)
  if (!identical(b$type, reference$type) || !setequal(b$target_ids, reference$target_ids) ||
      ncol(b$matrix) != ncol(reference$matrix))
    stop("Subject and reference target bases are incompatible.", call. = FALSE)
  B <- b$matrix[match(reference$target_ids, b$target_ids), , drop = FALSE]
  H <- qr.solve(reference$matrix, B, tol = 1e-9)
  if (sqrt(sum((reference$matrix %*% H - B)^2)) > 1e-8 * sqrt(sum(B^2)) || min(svd(H, nu = 0, nv = 0)$d) < 1e-9 * max(svd(H, nu = 0, nv = 0)$d))
    stop("Subject and reference bases must span the same target subspace.", call. = FALSE)
  H
}

.pattern_psd <- function(V) {
  e <- eigen((V + t(V))/2, symmetric = TRUE)
  tcrossprod(sweep(e$vectors, 2L, sqrt(pmax(e$values, 0)), "*"))
}

#' Combine independent subject confirmations in shared target coordinates
#'
#' @param subject_confirmations List of \code{pattern_confirmation} objects,
#'   each with a distinct subject_id and independent confirmation observations.
#' @param reference_basis A \code{pattern_basis} or discovery fit defining the
#'   shared raw target coordinates. Must be fixed independently of confirmation
#'   outcomes. Every subject basis must span exactly this target subspace.
#' @param spatial_mapping Optional list (one element per subject) of named
#'   character vectors: names are shared feature IDs, values are that subject's
#'   source feature IDs. Mappings must be one-to-one and cover the same shared
#'   features. With NULL, all subjects must have the same feature-ID set.
#'   This is correspondence, not spatial interpolation or parcel averaging.
#' @param effects \code{"random"} estimates the equally weighted population
#'   mean across subjects; \code{"fixed"} estimates a common effect using full
#'   inverse-covariance weighting.
#' @return A \code{pattern_group_result} with mean estimates, full mean
#'   covariance, component tests, invariant omnibus tests and effect norms,
#'   moment estimates of between-subject covariance, aligned subject estimates,
#'   leave-one-subject-out descriptive expression, prediction summaries, and
#'   provenance. Save the complete object with saveRDS.
#' @details Subject coefficients and their full covariance are transformed to
#'   the reference basis before pooling. Unequal target subspaces are rejected;
#'   Procrustes approximation would change the estimand. Original feature units,
#'   target names/units, nuisance meaning, and preprocessing recipes must agree.
#'   The preprocessing identifier and nuisance column names are checked; their
#'   scientific equivalence remains the caller's responsibility. Target sampling
#'   and omitted target effects must also permit a common conditional estimand.
#'
#'   Random effects use the arithmetic subject mean and sample coefficient
#'   covariance divided by number of subjects. Between-subject covariance is
#'   the positive-semidefinite part of sample covariance minus average sampling
#'   covariance (a moment estimate, not REML). Sampling error remains in the mean
#'   covariance; it is not added a second time. Component t tests have s-1 df;
#'   the omnibus uses Hotelling's T-squared transformed to F(r, s-r). These are
#'   exact for iid Gaussian subject estimates with a common total covariance,
#'   and approximate with heterogeneous within-subject precision or nonnormal
#'   effects. Requires s > r+1. Subject count, not observation count, sets df.
#'
#'   Fixed effects use multivariate GLS and asymptotic normal/chi-squared tests
#'   treating estimated within-subject covariances as known; heterogeneity Q has
#'   r(s-1) df under the common-effect null. They do not support population
#'   generalization. Holm correction is separate for all component and all
#'   omnibus tests. Singular mean covariance yields NA omnibus inference.
#'
#'   Norms, omnibus tests, heterogeneity trace, and expression are invariant to
#'   orthogonal reference rotations. Component estimates and tests are basis
#'   dependent. Arbitrary scaling of reference axes changes norms. Expression
#'   projects a subject's coefficient vector onto the normalized mean of the
#'   other subjects; it is descriptive, not an independent group prediction.
#'   Interpolation would need cross-feature covariance, which confirmations do
#'   not store, so many-to-one spatial mappings are deliberately rejected.
#' @seealso \code{pattern_confirm}, \code{pattern_basis}
#' @export
pattern_group <- function(subject_confirmations, reference_basis,
                          spatial_mapping = NULL, effects = c("random", "fixed")) {
  effects <- match.arg(effects)
  subjects <- subject_confirmations; s <- length(subjects)
  if (!is.list(subjects) || s < 2L || !all(vapply(subjects, inherits, logical(1), "pattern_confirmation")))
    stop("Supply at least two subject confirmations.", call. = FALSE)
  reference <- .pattern_group_basis(reference_basis); r <- ncol(reference$matrix)
  if (effects == "random" && s <= r + 1L)
    stop("Random-effects inference needs more subjects than target dimensions plus one.", call. = FALSE)
  ids <- vapply(subjects, function(x) .pattern_string(x$provenance$subject_id, "subject_id"), character(1))
  if (anyDuplicated(ids)) stop("Each subject_id must be unique.", call. = FALSE)
  observations <- unlist(lapply(subjects, function(x) x$provenance$observation_ids), use.names = FALSE)
  discoveries <- unlist(lapply(subjects, function(x) x$provenance$discovery_ids), use.names = FALSE)
  if (anyDuplicated(observations) || length(intersect(observations, discoveries)))
    stop("Subject confirmation IDs overlap each other or a discovery set.", call. = FALSE)
  recipes <- vapply(subjects, function(x) .pattern_string(x$provenance$preprocessing_id, "preprocessing_id"), character(1))
  nuisances <- lapply(subjects, function(x) colnames(x$nuisance))
  if (length(unique(recipes)) != 1L || !all(vapply(nuisances, identical, logical(1), nuisances[[1]])))
    stop("Subject preprocessing recipes and nuisance column definitions must match.", call. = FALSE)
  if (is.null(spatial_mapping)) {
    common <- subjects[[1]]$feature_ids
    if (!all(vapply(subjects, function(x) setequal(x$feature_ids, common), logical(1))))
      stop("Feature sets differ; supply explicit spatial_mapping.", call. = FALSE)
    spatial_mapping <- unname(lapply(subjects, function(x) stats::setNames(common, common)))
  }
  if (!is.list(spatial_mapping) || length(spatial_mapping) != s)
    stop("spatial_mapping must have one named character vector per subject.", call. = FALSE)
  # Named lists are matched by subject identity, never by accidental list order.
  if (!is.null(names(spatial_mapping))) {
    if (anyDuplicated(names(spatial_mapping)) || !setequal(names(spatial_mapping), ids))
      stop("Named spatial_mapping entries must match subject IDs.", call. = FALSE)
    spatial_mapping <- spatial_mapping[match(ids, names(spatial_mapping))]
  }
  common <- names(spatial_mapping[[1]])
  if (!length(common)) stop("Spatial mappings need shared feature names.", call. = FALSE)
  common <- .pattern_ids(common, length(common), "shared feature IDs"); p <- length(common)
  estimates <- array(NA_real_, c(p, r, s), dimnames = list(common, paste0("component", seq_len(r)), ids))
  positions <- transforms <- vector("list", s)
  for (j in seq_len(s)) {
    x <- subjects[[j]]; map <- spatial_mapping[[j]]
    x$feature_ids <- .pattern_ids(x$feature_ids, nrow(x$estimate), "subject feature_ids")
    if (!is.character(map) || anyDuplicated(names(map)) || !setequal(names(map), common) ||
        length(map) != p || anyNA(map) || anyDuplicated(map) || !all(map %in% x$feature_ids))
      stop("Spatial mappings must be complete one-to-one feature correspondences.", call. = FALSE)
    rows <- match(map[match(common, names(map))], x$feature_ids)
    H <- .pattern_group_transform(x, reference); transforms[[j]] <- H
    estimates[, , j] <- x$estimate[rows, , drop = FALSE] %*% t(H)
    positions[[j]] <- rows
  }
  if (any(!is.finite(estimates))) stop("Subject estimates must be finite.", call. = FALSE)
  mean_estimate <- se <- matrix(NA_real_, p, r, dimnames = dimnames(estimates)[1:2])
  mean_covariance <- between <- array(NA_real_, c(r, r, p))
  omnibus <- rep(NA_real_, p); p_omnibus <- rep(NA_real_, p)
  Qval <- Qp <- rep(NA_real_, p)
  expression <- matrix(NA_real_, p, s, dimnames = list(common, ids))
  for (v in seq_len(p)) {
    Y <- t(matrix(estimates[v, , ], r, s))
    # Transform one feature at a time: do not duplicate every subject's
    # rank-by-rank covariance array in the group workspace.
    Vs <- lapply(seq_len(s), function(j) {
      H <- transforms[[j]]
      V <- H %*% .pattern_confirmation_cov(subjects[[j]], positions[[j]][v]) %*% t(H)
      if (any(!is.finite(V)) || max(abs(V - t(V))) > 1e-10 * max(abs(V)) ||
          min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) < -1e-10 * max(abs(V)))
        stop("Subject sampling covariances must be finite symmetric positive semidefinite.", call. = FALSE)
      V
    })
    have_precision <- all(vapply(Vs, .pattern_group_pd, logical(1)))
    if (have_precision) {
      W <- lapply(Vs, solve); Vfixed <- solve(Reduce(`+`, W))
      fixed_mu <- as.vector(Vfixed %*% Reduce(`+`, lapply(seq_len(s), function(j) W[[j]] %*% Y[j, ])))
    }
    sample_cov <- matrix(stats::cov(Y), r, r)
    between[, , v] <- .pattern_psd(sample_cov - Reduce(`+`, Vs) / s)
    if (effects == "random") {
      mu <- colMeans(Y); Vmean <- sample_cov / s
    } else {
      if (!have_precision)
        stop("Fixed-effects pooling requires nonsingular subject covariances.", call. = FALSE)
      Vmean <- Vfixed; mu <- fixed_mu
    }
    mean_estimate[v, ] <- mu; mean_covariance[, , v] <- Vmean
    se[v, ] <- sqrt(pmax(diag(Vmean), 0))
    if (.pattern_group_pd(Vmean)) {
      wald <- sum(mu * solve(Vmean, mu))
      omnibus[v] <- if (effects == "random") wald * (s - r) / (r * (s - 1)) else wald
      p_omnibus[v] <- if (effects == "random") stats::pf(omnibus[v], r, s - r, lower.tail = FALSE) else
        stats::pchisq(omnibus[v], r, lower.tail = FALSE)
    }
    if (have_precision) {
      Qval[v] <- sum(vapply(seq_len(s), function(j) {
        delta <- Y[j, ] - fixed_mu
        sum(delta * (W[[j]] %*% delta))
      }, numeric(1)))
      Qp[v] <- stats::pchisq(Qval[v], r * (s - 1), lower.tail = FALSE)
    }
    for (j in seq_len(s)) {
      other <- colMeans(Y[-j, , drop = FALSE]); norm <- sqrt(sum(other^2))
      if (norm > 0) expression[v, j] <- sum(Y[j, ] * other) / norm
    }
  }
  stat <- mean_estimate / se; stat[!is.finite(stat) | se <= 0] <- NA_real_
  pval <- if (effects == "random") 2 * stats::pt(-abs(stat), s - 1L) else 2 * stats::pnorm(-abs(stat))
  structure(list(estimate = mean_estimate, se = se, statistic = stat, p = pval,
    p_holm = matrix(stats::p.adjust(pval, "holm"), p, r, dimnames = dimnames(pval)),
    covariance = mean_covariance, between_covariance = between,
    omnibus = data.frame(feature_id = common, statistic = omnibus, df1 = r,
      df2 = if (effects == "random") s - r else Inf, p = p_omnibus,
      p_holm = stats::p.adjust(p_omnibus, "holm")),
    effect_norm = sqrt(rowSums(mean_estimate^2)),
    heterogeneity = data.frame(feature_id = common,
      trace = vapply(seq_len(p), function(v) sum(diag(matrix(between[, , v], r, r))), numeric(1)),
      Q = Qval, df = r * (s - 1), p = Qp),
    subject_estimates = estimates, subject_expression = expression,
    prediction = .pattern_group_prediction(subjects, ids),
    reference_basis = reference, transforms = stats::setNames(transforms, ids),
    effects = effects, n_subjects = s, feature_ids = common,
    provenance = list(subject_ids = ids, confirmations = lapply(subjects, `[[`, "provenance"),
                      spatial_mapping = spatial_mapping, preprocessing_id = recipes[1])),
    class = "pattern_group_result")
}


.pattern_group_prediction <- function(subjects, ids) {
  tables <- lapply(seq_along(subjects), function(j) {
    tab <- subjects[[j]]$prediction
    if (is.null(tab)) return(NULL)
    tab$subject_id <- ids[j]; tab
  })
  if (any(vapply(tables, is.null, logical(1)))) return(NULL)
  rows <- do.call(rbind, tables)
  keys <- interaction(rows$response, rows$metric, drop = TRUE)
  summary <- do.call(rbind, lapply(split(rows, keys), function(x) {
    finite <- is.finite(x$value)
    data.frame(response = x$response[1], metric = x$metric[1],
      mean = if (all(finite)) mean(x$value) else NA_real_, n_subjects = sum(finite))
  }))
  rownames(summary) <- NULL
  list(subjects = rows, summary = summary)
}

#' @rdname pattern_group
#' @param x A group result to print.
#' @param ... Reserved.
#' @export
print.pattern_group_result <- function(x, ...) {
  cat("pattern_group_result:", x$n_subjects, "subjects,", nrow(x$estimate), "features,", ncol(x$estimate), "components\n")
  cat("  effects:", x$effects, "in the supplied reference target basis\n")
  cat("  subject expression is descriptive; component tests depend on the reference axes\n")
  invisible(x)
}


.pattern_group_pd <- function(V) {
  ev <- eigen((V + t(V))/2, symmetric = TRUE, only.values = TRUE)$values
  min(ev) > 1e-10 * max(ev)
}
