# Confirmation is conditional on a frozen discovery basis. No direction,
# feature selection, or preprocessing parameter is learned from these rows.

.pattern_ids <- function(x, n, name) {
  if (!is.atomic(x) || length(x) != n || anyNA(x) ||
      any(!nzchar(as.character(x))) || anyDuplicated(x)) {
    stop(name, " must contain one unique, nonempty ID per row/feature.", call. = FALSE)
  }
  as.character(x)
}

.pattern_string <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x))
    stop(name, " must be one nonempty string.", call. = FALSE)
  x
}

#' Specify a confirmation error model
#'
#' @param error Independent Gaussian rows, a CR1 sandwich over independent
#'   blocks, or a restricted-residual block wild bootstrap with Rademacher signs.
#' @param n_resamples Number of Monte Carlo sign draws (at least 99).
#' @param seed Nonnegative integer seed. Resampling restores the caller's RNG.
#' @return A \code{pattern_confirmation_plan}.
#' @details The independent model uses residual degrees of freedom. Block
#'   sandwich t and Wald F tests use number of blocks minus one and are
#'   small-sample approximations. The sign-flip option uses the same sandwich
#'   statistics and refits residuals under each null. With estimated nuisance
#'   effects this is an approximate wild bootstrap, not an exact permutation
#'   test. Independent blocks and enough blocks for the tested dimension are
#'   required. No option establishes independence from row labels alone.
#' @export
confirmation_plan <- function(error = c("independent", "block_robust", "sign_flip"),
                              n_resamples = 999L, seed = 1L) {
  error <- match.arg(error)
  if (length(n_resamples) != 1L || !is.finite(n_resamples) ||
      n_resamples < 99 || n_resamples != floor(n_resamples) || n_resamples > .Machine$integer.max)
    stop("n_resamples must be an integer of at least 99.", call. = FALSE)
  if (length(seed) != 1L || !is.finite(seed) || seed < 0 ||
      seed != floor(seed) || seed > .Machine$integer.max)
    stop("seed must be a nonnegative integer.", call. = FALSE)
  structure(list(error = error, n_resamples = as.integer(n_resamples),
                 seed = as.integer(seed)), class = "pattern_confirmation_plan")
}

#' Describe the frozen target coordinates of a pattern fit
#'
#' @param fit A pattern fit, view, or global result with a refit. Views use the
#'   underlying fitted basis, not their independent display rotations.
#' @return A \code{pattern_basis} containing the raw-target-to-score matrix,
#'   target IDs/type, centering, and a hash of the complete frozen transform.
#' @details The intercept in confirmation absorbs target centering. The matrix
#'   still includes training scaling and whitening. Equality of C alone does
#'   not establish compatible coordinates. Target names and physical units
#'   must have the same meaning in every subject.
#' @export
pattern_basis <- function(fit) {
  fit <- .pattern_base_fit(fit)
  yt <- fit$y_transform
  B <- (yt$Wy / (yt$sd %||% rep(1, yt$q))) %*% fit$C
  rownames(B) <- yt$response_ids
  structure(list(matrix = B, target_ids = yt$response_ids, type = yt$type,
                 center = yt$mu,
                 basis_id = digest::digest(list(C = fit$C, y_transform = yt))),
            class = "pattern_basis")
}

.pattern_confirmation_data <- function(fit, dataset, design) {
  X <- if (is.matrix(dataset)) dataset else get_feature_matrix(dataset)
  X <- as.matrix(X)
  if (!is.numeric(X) || ncol(X) != fit$p_input || nrow(X) < 3L || any(!is.finite(X)))
    stop("Confirmation requires finite numeric data with all input feature columns and at least three rows.", call. = FALSE)
  values <- if (is.factor(design) || is.atomic(design) || is.matrix(design)) design else {
    tt <- model_targets(design, "train")
    if (!is.null(tt$row_weights) && any(tt$row_weights != 1))
      stop("Confirmation does not support nonuniform row weights.", call. = FALSE)
    tt$values
  }
  if (fit$y_transform$type == "continuous") {
    values <- as.matrix(values)
    if (!is.numeric(values) || any(!is.finite(values))) stop("Targets must be finite numeric values.", call. = FALSE)
    if (ncol(values) > 1L && is.null(colnames(values)))
      stop("Multivariate confirmation targets require frozen response IDs as column names.", call. = FALSE)
    if (!is.null(colnames(values)) && !identical(colnames(values), fit$y_transform$response_ids))
      stop("Target columns must match the frozen response IDs and order.", call. = FALSE)
  }
  Tm <- .pattern_targets_apply(fit$y_transform, values) %*% fit$C
  if (nrow(Tm) != nrow(X) || any(!is.finite(Tm))) stop("Targets and data must have matching finite rows.", call. = FALSE)
  list(X = X, targets = values, T = Tm)
}

.pattern_nuisance <- function(nuisance, n) {
  if (is.null(nuisance)) return(matrix(1, n, 1L, dimnames = list(NULL, "intercept")))
  N <- as.matrix(nuisance)
  if (!is.numeric(N) || nrow(N) != n || any(!is.finite(N)))
    stop("nuisance must be a finite numeric matrix with one row per observation; omit the intercept.", call. = FALSE)
  if (is.null(colnames(N))) colnames(N) <- paste0("nuisance", seq_len(ncol(N)))
  if (anyDuplicated(colnames(N)) || any(!nzchar(colnames(N))))
    stop("Nuisance column names must be unique and nonempty.", call. = FALSE)
  cbind(intercept = 1, N)
}

.pattern_lm_setup <- function(D, r, blocks = NULL) {
  n <- nrow(D); d <- ncol(D)
  q <- qr(D, tol = 1e-9)
  if (q$rank != d || n <= d) stop("Confirmation design is rank deficient or has no residual degrees of freedom.", call. = FALSE)
  # QR rather than normal-equation inversion; L maps observations to slopes.
  L <- backsolve(qr.R(q), t(qr.Q(q)))[order(q$pivot), , drop = FALSE][seq_len(r), , drop = FALSE]
  g <- if (is.null(blocks)) n else length(unique(blocks))
  if (!is.null(blocks) && g <= r + 1L)
    stop("Need more independent blocks than tested dimensions plus one.", call. = FALSE)
  list(qr = q, D = D, L = L, bread = tcrossprod(L), inverse_bread = solve(tcrossprod(L)),
       r = r, n = n, d = d,
       blocks = blocks, g = g, df = if (is.null(blocks)) n - d else g - 1L)
}

.pattern_mass_lm <- function(Y, setup) {
  s <- setup; r <- s$r; p <- ncol(Y)
  beta <- t(s$L %*% Y)
  E <- qr.resid(s$qr, Y)
  sigma2 <- colSums(E^2) / (s$n - s$d)
  covariance <- NULL
  if (is.null(s$blocks)) {
    se <- sqrt(outer(sigma2, diag(s$bread)))
    Fval <- rowSums((beta %*% s$inverse_bread) * beta) / (r * sigma2)
  } else {
    covariance <- array(0, c(r, r, p))
    correction <- s$g / (s$g - 1) * (s$n - 1) / (s$n - s$d)
    for (rows in split(seq_len(s$n), s$blocks)) {
      scores <- s$L[, rows, drop = FALSE] %*% E[rows, , drop = FALSE]
      for (a in seq_len(r)) for (b in seq_len(r))
        covariance[a, b, ] <- covariance[a, b, ] + scores[a, ] * scores[b, ]
    }
    covariance <- covariance * correction
    se <- matrix(NA_real_, p, r); Fval <- rep(NA_real_, p)
    if (r == 1L) {
      se[, 1] <- sqrt(pmax(covariance[1, 1, ], 0))
      Fval <- beta[, 1]^2 / covariance[1, 1, ]
    } else if (r == 2L) {
      a <- covariance[1, 1, ]; b <- covariance[1, 2, ]; c <- covariance[2, 2, ]
      se[, 1] <- sqrt(pmax(a, 0)); se[, 2] <- sqrt(pmax(c, 0))
      det <- a * c - b^2
      # Eigenvalue-relative criterion is orthogonally invariant.
      largest <- (a + c + sqrt((a - c)^2 + 4 * b^2)) / 2
      ok <- det > 1e-10 * largest^2
      Fval[ok] <- (c[ok] * beta[ok, 1]^2 - 2*b[ok]*beta[ok, 1]*beta[ok, 2] +
                    a[ok] * beta[ok, 2]^2) / (2 * det[ok])
    } else {
      for (v in seq_len(p)) {
        V <- matrix(covariance[, , v], r, r)
        se[v, ] <- sqrt(pmax(diag(V), 0))
        ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
        if (min(ev) > 1e-10 * max(ev)) Fval[v] <- sum(beta[v, ] * solve(V, beta[v, ])) / r
      }
    }
  }
  tval <- beta / se
  # Zero residual variation cannot support a sampling distribution.
  tval[!is.finite(tval) | se <= 0] <- NA_real_
  Fval[!is.finite(Fval) | sigma2 <= 0] <- NA_real_
  list(estimate = beta, se = se, t = tval,
       p = 2 * stats::pt(-abs(tval), s$df), F = Fval,
       p_omnibus = stats::pf(Fval, r, s$df, lower.tail = FALSE),
       covariance = covariance, bread = s$bread, sigma2 = sigma2, df = s$df)
}

.pattern_with_seed <- function(seed, fun) {
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (existed) assign(".Random.seed", old, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  fun()
}

.pattern_wild_test <- function(Y, setup, observed, inference) {
  s <- setup; r <- s$r; p <- ncol(Y); B <- inference$n_resamples
  group <- match(s$blocks, unique(s$blocks))
  signs <- .pattern_with_seed(inference$seed, function()
    matrix(sample(c(-1, 1), s$g * B, replace = TRUE), s$g, B))
  p_raw <- p_max <- matrix(NA_real_, p, r)
  p_F <- p_Fmax <- rep(NA_real_, p)
  # One restricted null for the omnibus; one for each component conditional
  # on all other scores. Each null shares draws across all feature columns.
  for (k in 0:r) {
    excluded <- if (k == 0L) seq_len(r) else k
    reduced <- qr(s$D[, -excluded, drop = FALSE], tol = 1e-9)
    E0 <- qr.resid(reduced, Y); mu0 <- Y - E0
    obs <- if (k == 0L) observed$F else abs(observed$t[, k])
    valid <- is.finite(obs)
    count <- count_max <- rep(1, p)
    for (b in seq_len(B)) {
      draw <- .pattern_mass_lm(mu0 + E0 * signs[group, b], s)
      statistic <- if (k == 0L) draw$F else abs(draw$t[, k])
      # A singular bootstrap covariance is not evidence against the null.
      statistic[!is.finite(statistic)] <- Inf
      count <- count + as.integer(statistic >= obs - 1e-12)
      count_max <- count_max + as.integer(max(statistic[valid], -Inf) >= obs - 1e-12)
    }
    raw <- count / (B + 1); adjusted <- count_max / (B + 1)
    raw[!valid] <- adjusted[!valid] <- NA_real_
    if (k == 0L) { p_F <- raw; p_Fmax <- adjusted } else {
      p_raw[, k] <- raw
      # Max over voxels within a component, Bonferroni over components.
      p_max[, k] <- pmin(1, r * adjusted)
    }
  }
  list(p = p_raw, p_max = p_max, p_omnibus = p_F, p_omnibus_max = p_Fmax)
}

.pattern_confirmation_cov <- function(object, v) {
  if (is.null(object$covariance)) object$covariance_factor * object$residual_variance[v] else
    matrix(object$covariance[, , v], ncol(object$estimate), ncol(object$estimate))
}

#' Confirm a frozen pattern basis on independent observations
#'
#' @param fit A pattern fit, view, or global result with a refit.
#' @param dataset Confirmation observations by all input features, or an MVPA
#'   dataset whose training partition contains only confirmation observations.
#' @param design Confirmation targets (factor or numeric matrix), or an MVPA
#'   design with those targets in its training partition.
#' @param block_var Row-aligned independent run/block labels. Required for
#'   block-based error models; these labels do not add nuisance intercepts.
#' @param inference A \code{confirmation_plan()}.
#' @param observation_ids Globally meaningful, unique confirmation row IDs.
#' @param discovery_ids Complete discovery row IDs in the same namespace,
#'   including all rows used for selection, tuning, or preprocessing. Required:
#'   the legacy fit's positional train/test IDs are not adequate provenance.
#' @param nuisance Numeric nuisance columns, without an intercept. Categorical
#'   nuisance variables must be coded by the caller using model.matrix.
#' @param feature_ids Unique IDs for all input columns, in their input order.
#' @param preprocessing_id A stable identifier for measurement units and the
#'   preprocessing recipe; matching strings are a caller assertion, not proof.
#' @param subject_id Optional unique participant ID, required for group analysis.
#' @return A \code{pattern_confirmation}: original-feature-unit unpenalized
#'   estimates and SEs, t and omnibus F tests, Holm-adjusted p values, complete
#'   within-feature coefficient covariance (factorized for independent errors),
#'   frozen basis and provenance. Sign-flip plans also return bootstrap p and
#'   max-statistic adjusted p values. Nonestimable sampling distributions are NA.
#' @details Regresses every input feature on frozen target scores and nuisance
#'   columns, including an intercept. Discovery feature screening is not repeated
#'   or used to select the confirmation hypothesis family. All target scores must
#'   be estimable after nuisance adjustment. Independent Gaussian homoskedastic
#'   errors give exact t/F tests; block methods are approximate. Holm correction
#'   covers all feature-component tests, and separately all omnibus tests.
#'   Component columns depend on the frozen basis; omnibus tests are invariant
#'   to nonsingular changes of score coordinates. Rank selection is not a rank
#'   hypothesis test and this API does not report a supported population rank.
#'   IDs detect overlap but the caller must establish genuine independence,
#'   including no shared preprocessing fit or correlated repeated measurements.
#' @seealso \code{pattern_component_tests}, \code{pattern_basis}
#' @export
pattern_confirm <- function(fit, dataset, design, block_var = NULL,
                            inference = confirmation_plan(), observation_ids,
                            discovery_ids, nuisance = NULL, feature_ids,
                            preprocessing_id, subject_id = NULL) {
  fit <- .pattern_base_fit(fit)
  if (!inherits(inference, "pattern_confirmation_plan")) stop("Use confirmation_plan().", call. = FALSE)
  inference <- do.call(confirmation_plan, unclass(inference))
  dat <- .pattern_confirmation_data(fit, dataset, design)
  n <- nrow(dat$X); p <- ncol(dat$X); r <- ncol(dat$T)
  ids <- .pattern_ids(observation_ids, n, "observation_ids")
  discovery_ids <- .pattern_ids(discovery_ids, length(discovery_ids), "discovery_ids")
  if (!length(discovery_ids) || length(discovery_ids) < fit$n_train)
    stop("discovery_ids must cover at least all fitted training rows.", call. = FALSE)
  if (length(intersect(ids, discovery_ids))) stop("Discovery and confirmation IDs overlap.", call. = FALSE)
  feature_ids <- .pattern_ids(feature_ids, p, "feature_ids")
  if (!is.null(colnames(dat$X)) && !identical(colnames(dat$X), feature_ids))
    stop("Data column names must match feature_ids in order.", call. = FALSE)
  preprocessing_id <- .pattern_string(preprocessing_id, "preprocessing_id")
  if (!is.null(subject_id)) subject_id <- .pattern_string(subject_id, "subject_id")
  blocks <- NULL
  if (inference$error != "independent") {
    if (!is.atomic(block_var) || length(block_var) != n || anyNA(block_var) || any(!nzchar(as.character(block_var))))
      stop("block_var must identify an independent block for every confirmation row.", call. = FALSE)
    blocks <- as.character(block_var)
  } else if (!is.null(block_var)) {
    stop("Independent inference does not use block_var; choose a block error model or omit it.", call. = FALSE)
  }
  N <- .pattern_nuisance(nuisance, n)
  setup <- .pattern_lm_setup(cbind(dat$T, N), r, blocks)
  out <- .pattern_mass_lm(dat$X, setup)
  bootstrap <- if (inference$error == "sign_flip") .pattern_wild_test(dat$X, setup, out, inference) else NULL
  dimnames(out$estimate) <- dimnames(out$se) <- dimnames(out$t) <- dimnames(out$p) <-
    list(feature_ids, paste0("component", seq_len(r)))
  if (!is.null(bootstrap)) {
    dimnames(bootstrap$p) <- dimnames(bootstrap$p_max) <- dimnames(out$p)
  }
  p_primary <- if (is.null(bootstrap)) out$p else bootstrap$p
  omnibus_primary <- if (is.null(bootstrap)) out$p_omnibus else bootstrap$p_omnibus
  basis <- pattern_basis(fit)
  structure(list(estimate = out$estimate, se = out$se, t = out$t,
                 p = p_primary, p_holm = matrix(stats::p.adjust(p_primary, "holm"), p, r,
                                               dimnames = dimnames(out$p)),
                 omnibus = data.frame(feature_id = feature_ids, F = out$F,
                                      df1 = r, df2 = out$df, p = omnibus_primary,
                                      p_holm = stats::p.adjust(omnibus_primary, "holm")),
                 bootstrap = bootstrap, covariance = out$covariance,
                 covariance_factor = out$bread, residual_variance = out$sigma2,
                 df = out$df, basis = basis, feature_ids = feature_ids,
                 data_hash = digest::digest(list(dat$X, dat$targets)),
                 fit_hash = digest::digest(fit), nuisance = N,
                 prediction = .pattern_confirmation_prediction(fit, dat),
                 inference = inference,
                 provenance = list(discovery_ids = discovery_ids, observation_ids = ids,
                   subject_id = subject_id, preprocessing_id = preprocessing_id,
                   preprocessing_hash = digest::digest(list(preprocessing_id, fit$x_transform, fit$y_transform)),
                   basis_id = basis$basis_id, block_ids = blocks,
                   exchangeability_unit = if (is.null(blocks)) "row" else "block",
                   nuisance_hash = digest::digest(N))), class = "pattern_confirmation")
}

#' Test frozen component association and incremental decoding value
#'
#' @param fit The discovery fit used for confirmation.
#' @param confirmation A \code{pattern_confirmation} from that fit.
#' @param dataset,design The exact confirmation data and targets used by
#'   \code{pattern_confirm}. Content hashes are checked before computing scores.
#' @param calibration Optional list with \code{dataset}, \code{design}, and
#'   \code{observation_ids}. Its IDs must be included in the declared discovery
#'   set. Full and leave-one-component-out ordinary least-squares decoding heads
#'   are fitted on these rows alone. Without calibration only association is tested.
#' @return A list with an association table and, when calibration is supplied,
#'   an incremental-loss table, fitted decoding heads, and row-level losses.
#' @details Association is the partial correlation of each frozen decoded score
#'   with its corresponding frozen target score, adjusting for confirmation
#'   nuisance columns (not for other components). Tests use the confirmation
#'   error model and Holm correction over components.
#'
#'   Incremental gain is reduced-head minus full-head squared prediction error,
#'   averaged over raw target columns, then within independent blocks, then
#'   equally across blocks. Categorical targets use one-hot squared loss;
#'   linear heads are not probability models. Positive gain favors retaining
#'   the component. The paired block-mean t test is approximate; sign-flip plans
#'   use a centered, studentized block wild bootstrap, also approximate. These
#'   tests concern separately fitted decoding heads, not the original Gaussian
#'   posterior decoder or a supported latent rank. Component tests depend on the
#'   frozen basis. Score degeneracy is reported for association and rejected
#'   when a decoding head cannot be estimated.
#' @export
pattern_component_tests <- function(fit, confirmation, dataset, design, calibration = NULL) {
  fit <- .pattern_base_fit(fit)
  if (!inherits(confirmation, "pattern_confirmation")) stop("Expected a pattern_confirmation.", call. = FALSE)
  dat <- .pattern_confirmation_data(fit, dataset, design)
  if (!identical(pattern_basis(fit)$basis_id, confirmation$basis$basis_id) ||
      !identical(digest::digest(list(dat$X, dat$targets)), confirmation$data_hash) ||
      !identical(digest::digest(fit), confirmation$fit_hash))
    stop("The fit, data, and targets must match the confirmation record.", call. = FALSE)
  Z <- predict(fit, dat$X, type = "scores")
  r <- ncol(Z); plan <- confirmation$inference
  blocks <- confirmation$provenance$block_ids
  N <- confirmation$nuisance
  association <- data.frame(component = seq_len(r), correlation = NA_real_,
                            t = NA_real_, df = NA_real_, p = NA_real_)
  Nqr <- qr(N)
  for (k in seq_len(r)) {
    z <- qr.resid(Nqr, Z[, k]); target <- qr.resid(Nqr, dat$T[, k])
    if (sum(z^2) <= .Machine$double.eps * max(1, sum(Z[, k]^2))) next
    s <- .pattern_lm_setup(cbind(dat$T[, k], N), 1L, blocks)
    m <- .pattern_mass_lm(Z[, k, drop = FALSE], s)
    association$correlation[k] <- sum(z * target) / sqrt(sum(z^2) * sum(target^2))
    association$t[k] <- m$t[1]; association$df[k] <- m$df
    association$p[k] <- if (plan$error == "sign_flip")
      .pattern_wild_test(Z[, k, drop = FALSE], s, m, plan)$p[1] else m$p[1]
  }
  association$p_holm <- stats::p.adjust(association$p, "holm")
  result <- list(association = association, incremental = NULL,
                 basis_id = confirmation$basis$basis_id,
                 provenance = confirmation$provenance)
  if (is.null(calibration)) return(result)
  if (!is.list(calibration) || !all(c("dataset", "design", "observation_ids") %in% names(calibration)))
    stop("calibration requires dataset, design, and observation_ids.", call. = FALSE)
  cal <- .pattern_confirmation_data(fit, calibration$dataset, calibration$design)
  cal_ids <- .pattern_ids(calibration$observation_ids, nrow(cal$X), "calibration observation_ids")
  if (length(intersect(cal_ids, confirmation$provenance$observation_ids)))
    stop("Calibration and confirmation IDs overlap.", call. = FALSE)
  if (!all(cal_ids %in% confirmation$provenance$discovery_ids))
    stop("Calibration IDs must be included in the declared discovery set.", call. = FALSE)
  raw <- function(y) {
    if (fit$y_transform$type == "continuous") return(as.matrix(y))
    out <- matrix(0, length(y), fit$y_transform$K)
    out[cbind(seq_along(y), match(as.character(y), fit$y_transform$levels))] <- 1
    out
  }
  Ycal <- raw(cal$targets); Ytest <- raw(dat$targets)
  Dcal <- cbind(1, predict(fit, cal$X, type = "scores")); Dtest <- cbind(1, Z)
  head <- function(cols) {
    q <- qr(Dcal[, cols, drop = FALSE], tol = 1e-9)
    if (q$rank < length(cols) || nrow(Dcal) <= length(cols))
      stop("Calibration decoding head is rank deficient or lacks residual degrees of freedom.", call. = FALSE)
    qr.coef(q, Ycal)
  }
  full <- head(seq_len(r + 1L))
  full_loss <- rowMeans((Ytest - Dtest %*% full)^2)
  reduced <- vector("list", r); losses <- matrix(NA_real_, nrow(Z), r)
  for (k in seq_len(r)) {
    cols <- setdiff(seq_len(r + 1L), k + 1L)
    reduced[[k]] <- head(cols)
    losses[, k] <- rowMeans((Ytest - Dtest[, cols, drop = FALSE] %*% reduced[[k]])^2)
  }
  unit <- blocks %||% as.character(seq_len(nrow(Z)))
  index <- split(seq_len(nrow(Z)), unit)
  gain <- t(matrix(vapply(index, function(i) colMeans(losses[i, , drop = FALSE] - full_loss[i]), numeric(r)),
                   nrow = r, ncol = length(index)))
  g <- nrow(gain); avg <- colMeans(gain)
  se <- apply(gain, 2L, stats::sd) / sqrt(g)
  stat <- avg / se; stat[!is.finite(stat) | se <= 0] <- NA_real_
  pval <- 2 * stats::pt(-abs(stat), g - 1L)
  if (plan$error == "sign_flip") {
    pval <- .pattern_with_seed(plan$seed, function() {
      count <- rep(1, r); centered <- sweep(gain, 2L, avg, "-")
      for (b in seq_len(plan$n_resamples)) {
        draw <- centered * sample(c(-1, 1), g, replace = TRUE)
        tb <- colMeans(draw) / (apply(draw, 2L, stats::sd) / sqrt(g))
        tb[!is.finite(tb)] <- Inf
        count <- count + as.integer(abs(tb) >= abs(stat) - 1e-12)
      }
      count / (plan$n_resamples + 1)
    })
  }
  result$incremental <- data.frame(component = seq_len(r), gain = avg, se = se,
                                   t = stat, df = g - 1L, p = pval,
                                   p_holm = stats::p.adjust(pval, "holm"))
  result$heads <- list(full = full, reduced = reduced, calibration_ids = cal_ids)
  result$losses <- list(full = full_loss, reduced = losses, unit = unit, block_gain = gain)
  result
}


.pattern_confirmation_prediction <- function(fit, dat) {
  yt <- fit$y_transform
  if (yt$type == "categorical") {
    prob <- predict(fit, dat$X, type = "prob")
    index <- match(as.character(dat$targets), yt$levels)
    truth <- matrix(0, nrow(prob), ncol(prob)); truth[cbind(seq_len(nrow(prob)), index)] <- 1
    return(data.frame(response = "class", metric = c("Accuracy", "logloss", "Brier"),
      value = c(mean(max.col(prob, ties.method = "first") == index),
                -mean(log(pmax(prob[cbind(seq_len(nrow(prob)), index)], .Machine$double.xmin))),
                mean(rowSums((prob - truth)^2)))))
  }
  prediction <- predict(fit, dat$X, type = "decode")
  data.frame(response = rep(yt$response_ids, 2), metric = rep(c("MSE", "cor"), each = yt$q),
    value = c(colMeans((prediction - dat$targets)^2),
      vapply(seq_len(yt$q), function(j) {
        if (stats::sd(prediction[, j]) == 0 || stats::sd(dat$targets[, j]) == 0) return(NA_real_)
        stats::cor(prediction[, j], dat$targets[, j])
      }, numeric(1))))
}

#' @rdname pattern_confirm
#' @param x A confirmation result to print.
#' @param ... Reserved.
#' @export
print.pattern_confirmation <- function(x, ...) {
  cat("pattern_confirmation:", nrow(x$estimate), "features x", ncol(x$estimate), "frozen components\n")
  cat("  error model:", x$inference$error, "; denominator df:", x$df, "\n")
  cat("  independent confirmation rows:", length(x$provenance$observation_ids), "(declared)\n")
  cat("  estimates in original feature units; component and omnibus tests in separate families\n")
  invisible(x)
}
