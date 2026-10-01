# Exact sphere-aggregation searchlight engine.
#
# For classifiers whose per-sphere computation reduces to sums of per-voxel
# quantities, every sphere can be scored at once. Per fold, the per-voxel
# statistics (class means, products with test rows) are computed once. A
# sparse centre x voxel neighbourhood matrix then sums them over every sphere
# in a few sparse products. This is the same estimator as the general
# per-sphere path, not an approximation: the same voxel screening, the same
# class means, the same Pearson correlations, softmax, zapsmall rounding,
# fold pooling and metrics, up to floating-point summation order.
#
# Supported: corclass (method = "pearson", robust = FALSE) and naive_bayes.
# Data outside the proven regime (missing or non-finite values, identical
# voxel columns, a class absent from a training fold, and for naive_bayes a
# zero within-class variance) raise an "rmvpa_engine_ineligible" condition,
# and the run falls back to the general-purpose iterator.

#' @keywords internal
#' @noRd
.is_aggregate_fast_path <- function(model_spec, method) {
  if (!inherits(model_spec, "mvpa_model")) return(FALSE)
  if (!identical(method, "standard")) return(FALSE)
  label <- model_spec$model$label
  if (!(identical(label, "corclass") || identical(label, "naive_bayes"))) return(FALSE)
  if (!isTRUE(has_crossval(model_spec)) || isTRUE(has_test_set(model_spec))) return(FALSE)

  ds <- model_spec$dataset
  if (!inherits(ds, "mvpa_image_dataset") ||
      inherits(ds, "mvpa_multibasis_image_dataset")) return(FALSE)
  if (!is.null(model_spec$feature_selector)) return(FALSE)

  y <- y_train(model_spec)
  if (!is.factor(y) || nlevels(y) < 2L) return(FALSE)

  grid <- model_spec$tune_grid
  if (!is.null(grid)) {
    if (!is.data.frame(grid) || nrow(grid) != 1L) return(FALSE)
    if (identical(label, "corclass") &&
        (!identical(as.character(grid$method), "pearson") ||
         !identical(as.logical(grid$robust), FALSE))) return(FALSE)
  }

  if (!isTRUE(model_spec$compute_performance)) return(FALSE)
  perf <- model_spec$performance
  kind <- attr(perf, "rmvpa_perf_kind", exact = TRUE)
  if (is.null(kind) || !(kind %in% c("multiclass", "binary"))) return(FALSE)
  if (!is.null(attr(perf, "rmvpa_split_list", exact = TRUE))) return(FALSE)
  if (identical(kind, "binary") != (nlevels(y) == 2L)) return(FALSE)
  TRUE
}

#' @keywords internal
#' @noRd
.aggregate_ineligible <- function(reason) {
  structure(
    class = c("rmvpa_engine_ineligible", "error", "condition"),
    list(message = paste0("aggregate_fast: ", reason), call = NULL)
  )
}

#' Centre x voxel neighbourhood matrix, columns indexed by mask position
#' @keywords internal
#' @noRd
.aggregate_neighbourhoods <- function(ds, radius) {
  space_obj <- resolve_volume_space(ds)
  dims <- spatial_dim_shape(space_obj)
  spacing <- neuroim2::spacing(space_obj)[1:3]

  mask_indices <- ds$mask_indices
  if (is.null(mask_indices)) mask_indices <- compute_mask_indices(ds$mask)
  centers <- intersect(get_center_ids(ds), mask_indices)

  n_space <- prod(dims)
  mask_active <- logical(n_space)
  mask_active[mask_indices] <- TRUE
  col_lookup <- integer(n_space)
  col_lookup[mask_indices] <- seq_along(mask_indices)

  offsets <- .searchlight_offsets(radius, spacing = spacing)
  center_coords <- neuroim2::index_to_grid(space_obj, centers)
  cols <- lapply(seq_along(centers), function(i) {
    ids <- .offset_to_indices(center_coords[i, ], offsets, dims, mask_active)
    cc <- col_lookup[ids]
    unique(cc[cc > 0L])
  })
  S <- Matrix::sparseMatrix(
    i = rep.int(seq_along(cols), lengths(cols)),
    j = unlist(cols, use.names = FALSE),
    x = 1,
    dims = c(length(centers), length(mask_indices))
  )
  list(S = S, centers = centers, mask_indices = mask_indices)
}

#' Mann-Whitney AUC for each column of `score` (ties count one half), the
#' value yardstick::roc_auc_vec() returns with event_level = "second".
#' @keywords internal
#' @noRd
.aggregate_col_auc <- function(score, positive) {
  n_pos <- sum(positive)
  n_neg <- length(positive) - n_pos
  if (n_pos == 0L || n_neg == 0L) return(rep(NA_real_, ncol(score)))
  ranks <- matrixStats::colRanks(score, ties.method = "average", preserveShape = TRUE)
  (colSums(ranks[positive, , drop = FALSE]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

#' Per-fold corclass probabilities for a block of centres
#'
#' Returns a matrix with one row per (test row, centre) pair, with the test
#' row varying fastest, and one column per class. Each centre's rows are the
#' probabilities prob_corsimFit() would produce for that sphere: a softmax of
#' the Pearson correlations to the class means, then zapsmall() over the
#' centre's (test row x class) block. Centres with fewer than two valid voxels
#' in this fold get NA rows.
#' @keywords internal
#' @noRd
.aggregate_corclass_fold <- function(x_train, y_train, x_test, tS, valid, classes) {
  n_te <- nrow(x_test)
  K <- length(classes)
  n_c <- ncol(tS)

  # Screen out invalid voxels by zeroing their rows of t(S).
  tS_f <- Matrix::Diagonal(x = as.numeric(valid)) %*% tS
  p <- Matrix::colSums(tS_f)

  M <- group_means(x_train, 1, factor(y_train, levels = classes))
  M <- M[classes, , drop = FALSE]

  agg <- function(A) as.matrix(A %*% tS_f)
  Sx <- agg(x_test)
  Sxx <- agg(x_test * x_test)
  Sm <- agg(M)
  Smm <- agg(M * M)

  p_rep <- rep(p, each = n_te)
  pm1 <- p_rep - 1
  ss_x <- Sxx - Sx * Sx / p_rep
  sd_x <- pmax(sqrt(pmax(ss_x / pm1, 0)), .Machine$double.eps)
  # Conditioning of the aggregated centred sums: (mean^2 * p) / (sum of
  # squared deviations). It bounds how much summation error is amplified.
  kappa <- matrixStats::colMaxs(matrix((Sx * Sx / p_rep) / pmax(ss_x, .Machine$double.xmin), n_te))

  r <- matrix(NA_real_, n_te * n_c, K)
  for (k in seq_len(K)) {
    Sxm <- agg(x_test * rep(M[k, ], each = n_te))
    sm_rep <- rep(Sm[k, ], each = n_te)
    ss_m <- Smm[k, ] - Sm[k, ] * Sm[k, ] / p
    kappa <- pmax(kappa, (Sm[k, ] * Sm[k, ] / p) / pmax(ss_m, .Machine$double.xmin))
    sd_m <- pmax(sqrt(pmax(rep(ss_m, each = n_te) / pm1, 0)), .Machine$double.eps)
    r[, k] <- as.vector((Sxm - Sx * sm_rep / p_rep) / (pm1 * sd_x * sd_m))
  }
  err_bound <- 16 * p * .Machine$double.eps * (1 + kappa)

  bad <- p < 2
  r[rep(bad, each = n_te), ] <- 0
  e <- exp(r - matrixStats::rowMaxs(r))
  probs <- e / rowSums(e)

  # zapsmall() per centre over its (test row x class) block. Centres whose
  # aggregated values may fall on the wrong side of a rounding boundary,
  # given their error bound, are flagged for exact recomputation.
  zs <- .aggregate_zapsmall_blocks(probs, n_te, err_bound)
  probs <- zs$probs
  probs[rep(bad, each = n_te), ] <- NA_real_
  list(probs = probs, flag = zs$flag & !bad)
}

#' Rounding digits zapsmall() uses for a block whose largest |value| is mx
#'
#' zapsmall() has changed across R versions (R 4.3 rounds to
#' digits - log10(mx), non-integer; later versions take the ceiling). The rule
#' of the running R is identified once, by comparing candidate rules with
#' base::zapsmall() on probe matrices. If neither matches, NULL means "call
#' zapsmall() per block".
#' @keywords internal
#' @noRd
.aggregate_zapsmall_rule <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    rules <- list(
      plain = function(mx, d) ifelse(mx > 0, pmax(0, d - log10(mx)), d),
      ceiling = function(mx, d) ifelse(mx > 0, pmax(0, d - ceiling(log10(mx))), d)
    )
    probe_rng <- function(seed) {
      old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
      on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
      set.seed(seed)
      lapply(c(1e-3, 0.0137, 0.126, 0.5, 0.97, 3.3, 47), function(scale) {
        matrix(stats::runif(12) * scale, 3)
      })
    }
    probes <- probe_rng(20261001L)
    d <- getOption("digits")
    for (nm in names(rules)) {
      ok <- all(vapply(probes, function(m) {
        identical(base::zapsmall(m), round(m, digits = rules[[nm]](max(abs(m)), d)))
      }, logical(1)))
      if (ok) { cached <<- rules[[nm]]; return(cached) }
    }
    cached <<- FALSE
    cached
  }
})

#' zapsmall() applied separately to each centre's block of rows
#' @keywords internal
#' @noRd
.aggregate_zapsmall_blocks <- function(probs, n_te, err_bound) {
  K <- ncol(probs)
  n_c <- nrow(probs) %/% n_te
  rule <- .aggregate_zapsmall_rule()
  if (isFALSE(rule)) {
    for (c in seq_len(n_c)) {
      rows <- (c - 1L) * n_te + seq_len(n_te)
      probs[rows, ] <- base::zapsmall(probs[rows, , drop = FALSE])
    }
    # Without a known rule the boundaries cannot be checked: recompute all.
    return(list(probs = probs, flag = rep(TRUE, n_c)))
  }
  digits <- getOption("digits")
  block_max <- matrixStats::colMaxs(matrix(matrixStats::rowMaxs(abs(probs)), n_te))
  zd <- rule(block_max, digits)
  # round() uses the nearest integer number of digits.
  D <- floor(zd + 0.5)
  D_el <- rep(rep(D, each = n_te), times = K)
  scaled <- probs * 10^D_el
  boundary_gap <- abs(scaled - floor(scaled) - 0.5) / 10^D_el
  min_gap <- matrixStats::colMins(matrix(matrixStats::rowMins(matrix(boundary_gap, nrow(probs))), n_te))
  # The digits themselves change where zd crosses k + 0.5.
  zd_gap <- abs(zd - (floor(zd) + 0.5))
  zd_eb <- err_bound / (pmax(block_max, .Machine$double.xmin) * log(10))
  flag <- (min_gap <= err_bound) | (zd_gap <= zd_eb)
  list(probs = round(probs, digits = rep(rep(zd, each = n_te), times = K)), flag = flag)
}

#' Per-voxel Gaussian naive Bayes fit, identical to MVPAModels$naive_bayes$fit
#' on any subset of columns (every statistic is per voxel).
#' @keywords internal
#' @noRd
.aggregate_nb_fit <- function(x, y, classes) {
  y <- factor(y, levels = classes)
  mus <- vars <- matrix(NA_real_, length(classes), ncol(x), dimnames = list(classes, NULL))
  for (k in seq_along(classes)) {
    samples <- x[y == classes[k], , drop = FALSE]
    nk <- nrow(samples)
    mus[k, ] <- colMeans(samples)
    vars[k, ] <- matrixStats::colVars(samples) * (nk - 1) / nk
  }
  list(mus = mus, vars = vars, log_priors = as.numeric(log(table(y) / length(y))))
}

#' Per-fold naive Bayes probabilities for a block of centres
#'
#' Same layout as .aggregate_corclass_fold(). Per-voxel log densities are the
#' same dnorm() values the per-sphere path computes. Only their summation over
#' each sphere differs in order, so each row carries an error bound on its
#' probabilities, used to flag near-ties for exact recomputation.
#' @keywords internal
#' @noRd
.aggregate_nb_fold <- function(nb_fit, x_test, tS, valid) {
  n_te <- nrow(x_test)
  K <- nrow(nb_fit$mus)
  tS_f <- Matrix::Diagonal(x = as.numeric(valid)) %*% tS
  p <- Matrix::colSums(tS_f)
  n_c <- ncol(tS)

  lp <- matrix(NA_real_, n_te * n_c, K)
  abs_sum <- numeric(n_te * n_c)
  for (k in seq_len(K)) {
    ll <- stats::dnorm(x_test, mean = rep(nb_fit$mus[k, ], each = n_te),
                       sd = rep(sqrt(nb_fit$vars[k, ]), each = n_te), log = TRUE)
    ll[!is.finite(ll)] <- -1e100
    lp[, k] <- as.vector(as.matrix(ll %*% tS_f)) + nb_fit$log_priors[k]
    abs_sum <- pmax(abs_sum, as.vector(as.matrix(abs(ll) %*% tS_f)))
  }
  bad <- p < 2
  lp[rep(bad, each = n_te), ] <- 0
  shifted <- exp(lp - matrixStats::rowMaxs(lp))
  probs <- shifted / rowSums(shifted)
  probs[rep(bad, each = n_te), ] <- NA_real_
  # Summation error of the log posteriors bounds the probability error
  # (softmax is 1-Lipschitz in each coordinate).
  err <- 16 * rep(pmax(p, 1), each = n_te) * .Machine$double.eps * (abs_sum + 1)
  list(probs = probs, flag = rep(FALSE, n_c), err = err)
}

#' Centres whose pooled probabilities have a near-tie the error bound cannot
#' resolve
#'
#' `err` bounds the absolute error of each row's log posteriors, so each
#' probability carries a relative error of at most about 2 * err. A centre is
#' flagged when an observation's top two probabilities, or two observations'
#' one-vs-rest AUC scores, are within their combined error. Exact ties count
#' too, except between saturated observations: all non-top probabilities
#' below 1e-17. Those produce the same exact scores (1 and -1/(K-1), or 1 and
#' 0 for two classes) under any perturbation within the bound.
#' @keywords internal
#' @noRd
.aggregate_near_tie_flags <- function(pooled, err, n_obs, ok) {
  K <- ncol(pooled)
  P <- pooled / rowSums(pooled)
  P[is.na(P)] <- 1 / K  # centres with too few voxels; masked out by `ok`
  rel <- 4 * err
  top <- matrixStats::rowMaxs(P)
  rest <- rowSums(P) - top
  sat <- rest < 1e-17

  # First-order softmax sensitivity: |dp_k| <= 2 * err * min(p_k, 1 - p_k),
  # plus a floor for the rounding of the score arithmetic itself.
  floor_eps <- 16 * .Machine$double.eps
  second <- apply(P, 1, function(v) sort(v, decreasing = TRUE)[2])
  top_err <- rel * (pmin(top, 1 - top) + pmin(second, 1 - second)) + floor_eps
  top_flag <- !sat & (top - second) <= top_err
  flag <- matrixStats::colAnys(matrix(top_flag, n_obs))

  for (k in seq_len(K)) {
    # The one-vs-rest score p_k - mean(p_-k) equals (K p_k - 1) / (K - 1).
    score <- if (K == 2L) P[, k] else P[, k] - rowMeans(P[, -k, drop = FALSE])
    e <- (if (K == 2L) 1 else K / (K - 1)) * rel * pmin(P[, k], 1 - P[, k]) + floor_eps
    e[sat] <- 0
    S <- matrix(score, n_obs); E <- matrix(e, n_obs); SAT <- matrix(sat, n_obs)
    for (c in which(!flag & ok)) {
      o <- order(S[, c])
      s_c <- S[o, c]; e_c <- E[o, c]; sat_c <- SAT[o, c]
      gap <- diff(s_c)
      lo <- seq_len(n_obs - 1L); hi <- lo + 1L
      near <- gap <= e_c[lo] + e_c[hi]
      robust <- gap == 0 & sat_c[lo] & sat_c[hi]
      if (any(near & !robust)) flag[c] <- TRUE
    }
  }
  flag & ok
}

#' Mask-column indices of each centre's sphere in the general path's voxel
#' order (as get_searchlight() returns them), for exact recomputation.
#' @keywords internal
#' @noRd
.aggregate_generic_order <- function(ds, radius, centres, mask_indices) {
  sl <- get_searchlight(ds, "standard", radius)
  sp <- neuroim2::space(ds$mask)
  parent <- vapply(sl, function(w) as.integer(w@parent_index), 1L)
  lapply(centres, function(cid) {
    w <- sl[[match(cid, parent)]]
    match(as.integer(neuroim2::grid_to_index(sp, w@coords)), mask_indices)
  })
}

#' @keywords internal
#' @noRd
run_searchlight_aggregate_fast <- function(model_spec, radius, verbose = FALSE,
                                           block_size = 4000L, ...) {
  ds <- model_spec$dataset
  y_all <- y_train(model_spec)
  classes <- levels(y_all)
  K <- length(classes)
  kind <- attr(model_spec$performance, "rmvpa_perf_kind", exact = TRUE)
  class_metrics <- isTRUE(attr(model_spec$performance, "rmvpa_class_metrics", exact = TRUE))

  nb <- .aggregate_neighbourhoods(ds, radius)
  if (length(nb$centers) == 0L) return(empty_searchlight_result(ds))

  x_all <- as.matrix(neuroim2::series(ds$train_data, nb$mask_indices))
  if (nrow(x_all) != length(y_all)) {
    stop("aggregate_fast: mismatch between train rows and y_train length.")
  }
  if (!all(is.finite(x_all))) {
    stop(.aggregate_ineligible("data contain missing or non-finite values"))
  }
  label <- model_spec$model$label
  is_nb <- identical(label, "naive_bayes")
  # Pearson correlation across voxels is invariant to a common shift; removing
  # the grand mean limits cancellation in the aggregated sums. The original
  # values are kept for exact recomputation of flagged centres. Naive Bayes
  # is not shift invariant and needs no shift (its per-voxel terms are exact).
  x_orig <- x_all
  if (!is_nb) x_all <- x_all - mean(x_all)

  folds <- generate_folds(model_spec$crossval, tibble::tibble(.row = seq_len(nrow(x_all))), y_all)
  fold_list <- lapply(seq_len(nrow(folds)), function(i) {
    tr <- as.integer(.extract_sample_indices(folds$train[[i]]))
    te <- as.integer(.extract_sample_indices(folds$test[[i]]))
    ytr <- factor(y_all[tr], levels = classes)
    if (any(table(ytr) == 0L)) {
      stop(.aggregate_ineligible("a class is absent from a training fold"))
    }
    valid <- nonzeroVarianceColumns2(x_all[tr, , drop = FALSE])
    if (anyDuplicated(t(x_all[tr, valid, drop = FALSE])) > 0L) {
      stop(.aggregate_ineligible("identical voxel columns in a training fold"))
    }
    nb_fit <- NULL
    if (is_nb) {
      # naive_bayes floors a zero within-class variance at a value that depends
      # on the other voxels in the sphere; such data are left to the general path.
      nb_fit <- .aggregate_nb_fit(x_all[tr, , drop = FALSE], ytr, classes)
      if (any(nb_fit$vars[, valid, drop = FALSE] <= .Machine$double.eps)) {
        stop(.aggregate_ineligible("zero within-class variance in a training fold"))
      }
    }
    list(train = tr, test = te, valid = valid, nb_fit = nb_fit)
  })
  testind <- sort(unique(unlist(lapply(fold_list, `[[`, "test"))))
  observed <- y_all[testind]
  n_obs <- length(testind)

  tS_all <- Matrix::t(nb$S)
  n_centres <- length(nb$centers)
  blocks <- split(seq_len(n_centres), ceiling(seq_len(n_centres) / block_size))
  metric_names <- c("Accuracy", "AUC",
                    if (class_metrics && identical(kind, "multiclass")) paste0("AUC_", classes))
  perf <- matrix(NA_real_, n_centres, length(metric_names), dimnames = list(NULL, metric_names))
  n_repaired <- 0L

  for (blk in blocks) {
    n_b <- length(blk)
    tS <- tS_all[, blk, drop = FALSE]
    pooled <- matrix(0, n_obs * n_b, K)
    ok <- rep(TRUE, n_b)
    flagged <- rep(FALSE, n_b)
    err <- numeric(n_obs * n_b)
    offsets <- (seq_len(n_b) - 1L) * n_obs
    for (f in fold_list) {
      fr <- if (is_nb) {
        .aggregate_nb_fold(f$nb_fit, x_all[f$test, , drop = FALSE], tS, f$valid)
      } else {
        .aggregate_corclass_fold(x_all[f$train, , drop = FALSE], y_all[f$train],
                                 x_all[f$test, , drop = FALSE], tS, f$valid, classes)
      }
      pr <- fr$probs
      ok <- ok & !is.na(pr[(seq_len(n_b) - 1L) * length(f$test) + 1L, 1])
      flagged <- flagged | fr$flag
      dest <- as.vector(outer(match(f$test, testind), offsets, "+"))
      pooled[dest, ] <- pooled[dest, , drop = FALSE] + pr
      if (!is.null(fr$err)) err[dest] <- pmax(err[dest], fr$err)
    }
    if (is_nb) {
      flagged <- flagged | .aggregate_near_tie_flags(pooled, err, n_obs, ok)
    }
    # Exact recomputation, with the per-sphere classifier code, for centres
    # whose aggregated values were too close to a rounding boundary.
    repair <- which(flagged & ok)
    generic_cols <- if (length(repair)) {
      .aggregate_generic_order(ds, radius, nb$centers[blk[repair]], nb$mask_indices)
    }
    for (ri in seq_along(repair)) {
      b <- repair[ri]
      cols <- generic_cols[[ri]]
      rows_b <- offsets[b] + seq_len(n_obs)
      pooled[rows_b, ] <- 0
      for (f in fold_list) {
        vc <- cols[f$valid[cols]]
        ytr <- factor(y_all[f$train], levels = classes)
        pr <- if (is_nb) {
          nbm <- MVPAModels$naive_bayes
          fit <- nbm$fit(x_orig[f$train, vc, drop = FALSE], ytr, NULL, NULL, classes, NULL, NULL, TRUE)
          nbm$prob(fit, x_orig[f$test, vc, drop = FALSE])
        } else {
          fit <- corsimFit(x_orig[f$train, vc, drop = FALSE], ytr, "pearson", FALSE)
          prob_corsimFit(fit, x_orig[f$test, vc, drop = FALSE])
        }
        dest <- offsets[b] + match(f$test, testind)
        pooled[dest, ] <- pooled[dest, , drop = FALSE] + pr
      }
    }
    n_repaired <- n_repaired + sum(flagged & ok)
    bad_rows <- rep(!ok, each = n_obs)
    pooled[bad_rows, ] <- 1
    # wrap_result(): renormalise pooled rows, then the exact argmax.
    pooled <- pooled / rowSums(pooled)
    pred <- matrix(max.col(pooled, ties.method = "first"), n_obs)
    acc <- colMeans(pred == as.integer(observed))

    if (identical(kind, "binary")) {
      auc <- 2 * .aggregate_col_auc(matrix(pooled[, 2], n_obs), observed == classes[2]) - 1
      vals <- cbind(Accuracy = acc, AUC = auc)
    } else {
      auc_k <- vapply(seq_len(K), function(k) {
        score <- pooled[, k] - rowMeans(pooled[, -k, drop = FALSE])
        2 * .aggregate_col_auc(matrix(score, n_obs), observed == classes[k]) - 1
      }, numeric(n_b))
      auc_k <- matrix(auc_k, nrow = n_b)
      vals <- cbind(Accuracy = acc, AUC = rowMeans(auc_k, na.rm = TRUE))
      if (class_metrics) vals <- cbind(vals, auc_k)
    }
    vals[!ok, ] <- NA_real_
    perf[blk, ] <- vals
  }

  good <- !is.na(perf[, "Accuracy"])
  if (!any(good)) return(empty_searchlight_result(ds))
  out <- wrap_out(perf[good, , drop = FALSE], ds, ids = nb$centers[good])
  attr(out, "bad_results") <- tibble::tibble()
  attr(out, "aggregate_recomputed_centres") <- n_repaired
  out
}
