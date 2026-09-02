# Leakage-safe nested tuning for banded ridge -------------------------------

.brt_with_seed <- function(seed, code) {
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("seed must be one finite integer.", call. = FALSE)
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(code)
}

.brt_validate_rows <- function(rows, n, context) {
  .brc_split_indices(rows, n, context)
}

.brt_validate_blocks <- function(blocks, n) {
  if (is.null(blocks)) return(NULL)
  if (length(blocks) != n || anyNA(blocks)) {
    stop("blocks must have one non-missing value per observation.", call. = FALSE)
  }
  as.character(blocks)
}

#' Drop training rows within a temporal purge gap of any test row
#'
#' @keywords internal
#' @noRd
.brt_purge_train <- function(train, test, purge) {
  if (purge > 0L && length(train) && length(test)) {
    keep <- vapply(train, function(i) all(abs(i - test) > purge), logical(1))
    train <- train[keep]
  }
  train
}

#' Construct deterministic ordinary, blocked, or purged CV folds
#'
#' @keywords internal
#' @noRd
.banded_ridge_make_folds <- function(rows,
                                     n_folds,
                                     blocks = NULL,
                                     purge = 0L,
                                     seed = 1L) {
  if (!is.numeric(rows) || !length(rows) || any(!is.finite(rows)) ||
      any(rows != as.integer(rows)) || anyDuplicated(rows) || any(rows < 1L)) {
    stop("rows must be unique positive integer observation indices.", call. = FALSE)
  }
  rows <- as.integer(rows)
  n <- max(rows)
  if (!is.numeric(n_folds) || length(n_folds) != 1L || !is.finite(n_folds) ||
      n_folds != as.integer(n_folds) || n_folds < 2L) {
    stop("n_folds must be an integer of at least two.", call. = FALSE)
  }
  n_folds <- as.integer(n_folds)
  if (!is.numeric(purge) || length(purge) != 1L || !is.finite(purge) ||
      purge != as.integer(purge) || purge < 0L) {
    stop("purge must be one non-negative integer.", call. = FALSE)
  }
  purge <- as.integer(purge)

  if (!is.null(blocks)) {
    if (length(blocks) < n || anyNA(blocks[rows])) {
      stop("blocks must cover all requested rows without missing values.", call. = FALSE)
    }
    blocks <- as.character(blocks)
    units <- unique(blocks[rows])
    if (length(units) < n_folds) {
      stop("There are fewer independent blocks than requested folds.", call. = FALSE)
    }
    shuffled <- .brt_with_seed(seed, sample(units, length(units), replace = FALSE))
    unit_sizes <- vapply(shuffled, function(z) sum(blocks[rows] == z), integer(1))
    fold_units <- vector("list", n_folds)
    fold_sizes <- integer(n_folds)
    for (ii in order(unit_sizes, decreasing = TRUE, method = "radix")) {
      target <- which.min(fold_sizes)
      fold_units[[target]] <- c(fold_units[[target]], shuffled[[ii]])
      fold_sizes[[target]] <- fold_sizes[[target]] + unit_sizes[[ii]]
    }
    test_sets <- lapply(fold_units, function(z) rows[blocks[rows] %in% z])
  } else {
    if (length(rows) < n_folds) {
      stop("There are fewer rows than requested folds.", call. = FALSE)
    }
    shuffled <- .brt_with_seed(seed, sample(rows, length(rows), replace = FALSE))
    assignment <- rep(seq_len(n_folds), length.out = length(rows))
    test_sets <- lapply(seq_len(n_folds), function(k) sort(shuffled[assignment == k]))
  }

  folds <- lapply(seq_len(n_folds), function(k) {
    test <- sort(test_sets[[k]])
    train <- .brt_purge_train(setdiff(rows, test), test, purge)
    if (length(train) < 2L) {
      stop("Fold ", k, " has fewer than two training rows after blocking/purge.",
           call. = FALSE)
    }
    list(id = paste0("fold_", k), train = sort(train), test = test)
  })
  folds
}

.brt_validate_folds <- function(folds,
                                allowed_rows,
                                n,
                                context,
                                blocks = NULL,
                                purge = 0L,
                                require_coverage = TRUE) {
  if (!is.list(folds) || length(folds) < 2L) {
    stop(context, " must be a list containing at least two folds.", call. = FALSE)
  }
  allowed_rows <- sort(.brt_validate_rows(allowed_rows, n, paste0(context, " allowed_rows")))
  blocks <- .brt_validate_blocks(blocks, n)
  out <- vector("list", length(folds))
  ids <- character(length(folds))

  for (ii in seq_along(folds)) {
    fold <- folds[[ii]]
    if (!is.list(fold) || is.null(fold$train) || is.null(fold$test)) {
      stop(context, " fold ", ii, " must contain train and test indices.", call. = FALSE)
    }
    train <- .brt_validate_rows(fold$train, n, paste0(context, " fold ", ii, " train"))
    test <- .brt_validate_rows(fold$test, n, paste0(context, " fold ", ii, " test"))
    if (length(train) < 2L) {
      stop(context, " fold ", ii, " must have at least two training rows.", call. = FALSE)
    }
    if (length(intersect(train, test))) {
      stop(context, " fold ", ii, " train and test rows must be disjoint.", call. = FALSE)
    }
    if (!all(c(train, test) %in% allowed_rows)) {
      stop(context, " fold ", ii, " contains rows outside its allowed outer-training set.",
           call. = FALSE)
    }
    if (!is.null(blocks) && length(intersect(blocks[train], blocks[test]))) {
      stop(context, " fold ", ii, " splits at least one declared block.", call. = FALSE)
    }
    if (purge > 0L && any(vapply(train, function(i) any(abs(i - test) <= purge), logical(1)))) {
      stop(
        context, " fold ", ii, " has training rows within purge = ", purge,
        " of its test rows. Explicit fold lists are validated as supplied, not ",
        "purged: remove those training rows yourself, or supply a fold count ",
        "(or, for outer_crossval, a cross_validation object) so the purge is ",
        "applied automatically.",
        call. = FALSE
      )
    }
    ids[[ii]] <- if (is.null(fold$id)) paste0("fold_", ii) else as.character(fold$id)
    if (length(ids[[ii]]) != 1L || is.na(ids[[ii]]) || !nzchar(ids[[ii]])) {
      stop(context, " fold IDs must be non-empty strings.", call. = FALSE)
    }
    out[[ii]] <- list(id = ids[[ii]], train = sort(train), test = sort(test))
  }
  if (anyDuplicated(ids)) stop(context, " fold IDs must be unique.", call. = FALSE)
  if (require_coverage) {
    all_test <- unlist(lapply(out, `[[`, "test"), use.names = FALSE)
    counts <- tabulate(match(all_test, allowed_rows), nbins = length(allowed_rows))
    if (length(all_test) != length(allowed_rows) || any(counts != 1L)) {
      stop(context, " test rows must cover every allowed row exactly once.", call. = FALSE)
    }
  }
  out
}

.brt_simplex_grid <- function(n_groups, grid_points) {
  if (!is.numeric(grid_points) || length(grid_points) != 1L ||
      !is.finite(grid_points) || grid_points != as.integer(grid_points) ||
      grid_points < 2L) {
    stop("grid_points must be an integer of at least two.", call. = FALSE)
  }
  denominator <- as.integer(grid_points) - 1L
  recurse <- function(remaining, dimensions) {
    if (dimensions == 1L) return(matrix(remaining, nrow = 1L))
    do.call(rbind, lapply(0:remaining, function(x) {
      cbind(x, recurse(remaining - x, dimensions - 1L))
    }))
  }
  recurse(denominator, n_groups) / denominator
}

.brt_theta_matrix <- function(theta, group_names, context = "theta") {
  if (is.numeric(theta) && is.null(dim(theta))) {
    if (is.null(names(theta))) {
      stop(context, " vector must be named by feature band.", call. = FALSE)
    }
    if (anyNA(names(theta)) || any(!nzchar(names(theta))) ||
        anyDuplicated(names(theta)) || !setequal(names(theta), group_names)) {
      stop(context, " vector names must match feature-band names exactly.",
           call. = FALSE)
    }
    theta <- matrix(theta[group_names], nrow = 1L,
                    dimnames = list(NULL, group_names))
  }
  if (!is.matrix(theta) || !is.numeric(theta) || nrow(theta) < 1L ||
      ncol(theta) != length(group_names)) {
    stop(context, " must be a numeric matrix with one column per feature band.",
         call. = FALSE)
  }
  if (is.null(colnames(theta)) || anyDuplicated(colnames(theta)) ||
      !setequal(colnames(theta), group_names)) {
    stop(context, " column names must match feature-band names exactly.", call. = FALSE)
  }
  theta <- theta[, group_names, drop = FALSE]
  if (any(!is.finite(theta)) || any(theta < 0)) {
    stop(context, " must contain finite non-negative values.", call. = FALSE)
  }
  if (any(abs(rowSums(theta) - 1) > sqrt(.Machine$double.eps))) {
    stop(context, " rows must sum to one.", call. = FALSE)
  }
  theta <- theta / rowSums(theta)
  theta
}

.brt_candidate_ids <- function(alpha, theta, group_names) {
  vapply(seq_along(alpha), function(ii) {
    pieces <- paste0(group_names, "=", sprintf("%.17g", theta[ii, ]))
    paste(c(paste0("alpha=", sprintf("%.17g", alpha[[ii]])), pieces), collapse = "|")
  }, character(1))
}

.brt_canonical_candidates <- function(alpha, theta, group_names) {
  if (!is.numeric(alpha) || !length(alpha) || any(!is.finite(alpha)) || any(alpha < 0)) {
    stop("candidate alpha values must be finite and non-negative.", call. = FALSE)
  }
  theta <- .brt_theta_matrix(theta, group_names, "candidate theta")
  if (nrow(theta) != length(alpha)) {
    stop("candidate alpha and theta rows must have equal lengths.", call. = FALSE)
  }
  ord_args <- c(list(alpha), lapply(seq_len(ncol(theta)), function(j) theta[, j]),
                list(method = "radix"))
  ord <- do.call(order, ord_args)
  alpha <- as.numeric(alpha[ord])
  theta <- theta[ord, , drop = FALSE]
  ids <- .brt_candidate_ids(alpha, theta, group_names)
  if (anyDuplicated(ids)) stop("candidate alpha/theta combinations must be unique.", call. = FALSE)
  out <- list(alpha = alpha, theta = theta, candidate_id = ids,
              group_names = group_names)
  class(out) <- c("banded_ridge_candidates", "list")
  out
}

#' Create a canonical alpha-by-theta banded-ridge candidate manifest
#'
#' @keywords internal
#' @noRd
.banded_ridge_candidates <- function(group_names,
                                     alphas,
                                     method = c("fixed", "grid", "random"),
                                     theta = NULL,
                                     grid_points = 3L,
                                     n_theta = 20L,
                                     seed = 1L) {
  method <- match.arg(method)
  if (!is.character(group_names) || !length(group_names) || anyNA(group_names) ||
      any(!nzchar(group_names)) || anyDuplicated(group_names)) {
    stop("group_names must be unique non-empty feature-band names.", call. = FALSE)
  }
  if (!is.numeric(alphas) || !length(alphas) || any(!is.finite(alphas)) ||
      any(alphas < 0)) {
    stop("alphas must contain finite non-negative values.", call. = FALSE)
  }
  alphas <- sort(unique(as.numeric(alphas)))
  g <- length(group_names)
  theta_rows <- switch(method,
    fixed = {
      if (is.null(theta)) stop("theta is required for method = 'fixed'.", call. = FALSE)
      .brt_theta_matrix(theta, group_names)
    },
    grid = .brt_simplex_grid(g, grid_points),
    random = {
      if (!is.numeric(n_theta) || length(n_theta) != 1L || !is.finite(n_theta) ||
          n_theta != as.integer(n_theta) || n_theta < 1L) {
        stop("n_theta must be one positive integer.", call. = FALSE)
      }
      draws <- .brt_with_seed(seed, matrix(stats::rexp(as.integer(n_theta) * g),
                                           as.integer(n_theta), g))
      draws / rowSums(draws)
    }
  )
  colnames(theta_rows) <- group_names
  expanded_alpha <- rep(alphas, each = nrow(theta_rows))
  expanded_theta <- theta_rows[rep(seq_len(nrow(theta_rows)), times = length(alphas)),
                               , drop = FALSE]
  out <- .brt_canonical_candidates(expanded_alpha, expanded_theta, group_names)
  attr(out, "generation") <- list(method = method, seed = as.integer(seed),
                                   grid_points = as.integer(grid_points),
                                   n_theta = as.integer(n_theta))
  out
}

.brt_validate_candidates <- function(candidates, group_names) {
  if (!is.list(candidates) || is.null(candidates$alpha) || is.null(candidates$theta)) {
    stop("candidates must be created by .banded_ridge_candidates().", call. = FALSE)
  }
  generation <- attr(candidates, "generation")
  out <- .brt_canonical_candidates(candidates$alpha, candidates$theta, group_names)
  attr(out, "generation") <- generation
  out
}

.brt_scope <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be one of 'response', 'shared', 'roi', or 'fixed'.",
         call. = FALSE)
  }
  x <- match.arg(x, c("response", "shared", "roi", "fixed"))
  if (x == "roi") "shared" else x
}

.brt_theta_keys <- function(theta) {
  apply(theta, 1L, function(x) paste(sprintf("%.17g", x), collapse = ","))
}

.brt_best <- function(values, candidate_indices, tolerance) {
  ok <- is.finite(values)
  if (!any(ok)) return(NA_integer_)
  best <- min(values[ok])
  threshold <- tolerance * max(1, abs(best))
  candidate_indices[which(ok & values <= best + threshold)[[1L]]]
}

#' Candidate indices that survive the fixed alpha/theta scope constraints
#'
#' @keywords internal
#' @noRd
.brt_available_candidates <- function(candidates,
                                      alpha_scope,
                                      theta_scope,
                                      fixed_alpha = NULL,
                                      fixed_theta = NULL) {
  alpha_scope <- .brt_scope(alpha_scope, "alpha_scope")
  theta_scope <- .brt_scope(theta_scope, "theta_scope")
  keep <- rep(TRUE, length(candidates$alpha))
  if (alpha_scope == "fixed") {
    if (is.null(fixed_alpha)) {
      unique_alpha <- unique(candidates$alpha)
      if (length(unique_alpha) != 1L) {
        stop("fixed alpha scope requires fixed_alpha or exactly one candidate alpha.",
             call. = FALSE)
      }
      fixed_alpha <- unique_alpha
    }
    if (!is.numeric(fixed_alpha) || length(fixed_alpha) != 1L ||
        !is.finite(fixed_alpha) || fixed_alpha < 0) {
      stop("fixed_alpha must be one finite non-negative value.", call. = FALSE)
    }
    keep <- keep & candidates$alpha == fixed_alpha
  }
  theta_keys <- .brt_theta_keys(candidates$theta)
  if (theta_scope == "fixed") {
    if (is.null(fixed_theta)) {
      unique_theta <- unique(theta_keys)
      if (length(unique_theta) != 1L) {
        stop("fixed theta scope requires fixed_theta or exactly one candidate theta.",
             call. = FALSE)
      }
      fixed_key <- unique_theta
    } else {
      fixed_theta <- .brt_theta_matrix(fixed_theta, candidates$group_names, "fixed_theta")
      if (nrow(fixed_theta) != 1L) stop("fixed_theta must contain one simplex point.", call. = FALSE)
      fixed_key <- .brt_theta_keys(fixed_theta)
    }
    keep <- keep & theta_keys == fixed_key
  }
  available <- which(keep)
  if (!length(available)) stop("No candidate matches the fixed alpha/theta constraints.", call. = FALSE)
  available
}

.brt_select_candidates <- function(scores,
                                    candidates,
                                    alpha_scope,
                                    theta_scope,
                                    fixed_alpha = NULL,
                                    fixed_theta = NULL,
                                    tie_tolerance = 1e-12) {
  alpha_scope <- .brt_scope(alpha_scope, "alpha_scope")
  theta_scope <- .brt_scope(theta_scope, "theta_scope")
  if (!is.numeric(tie_tolerance) || length(tie_tolerance) != 1L ||
      !is.finite(tie_tolerance) || tie_tolerance < 0) {
    stop("tie_tolerance must be a finite non-negative scalar.", call. = FALSE)
  }
  available <- .brt_available_candidates(
    candidates, alpha_scope, theta_scope,
    fixed_alpha = fixed_alpha, fixed_theta = fixed_theta
  )
  theta_keys <- .brt_theta_keys(candidates$theta)
  v <- ncol(scores)

  if (alpha_scope != "shared" && theta_scope != "shared") {
    selected <- vapply(seq_len(v), function(j) {
      .brt_best(scores[available, j], available, tie_tolerance)
    }, integer(1))
  } else if (alpha_scope == "shared" && theta_scope == "shared") {
    objective <- rowMeans(scores[available, , drop = FALSE], na.rm = FALSE)
    winner <- .brt_best(objective, available, tie_tolerance)
    selected <- rep(winner, v)
  } else if (alpha_scope == "shared") {
    alpha_values <- sort(unique(candidates$alpha[available]))
    choices <- matrix(NA_integer_, length(alpha_values), v)
    objectives <- rep(Inf, length(alpha_values))
    for (aa in seq_along(alpha_values)) {
      idx <- available[candidates$alpha[available] == alpha_values[[aa]]]
      choices[aa, ] <- vapply(seq_len(v), function(j) {
        .brt_best(scores[idx, j], idx, tie_tolerance)
      }, integer(1))
      if (!anyNA(choices[aa, ])) {
        objectives[[aa]] <- mean(scores[cbind(choices[aa, ], seq_len(v))])
      }
    }
    best_alpha <- .brt_best(objectives, seq_along(alpha_values), tie_tolerance)
    selected <- if (is.na(best_alpha)) rep(NA_integer_, v) else choices[best_alpha, ]
  } else {
    theta_values <- unique(theta_keys[available])
    theta_values <- sort(theta_values)
    choices <- matrix(NA_integer_, length(theta_values), v)
    objectives <- rep(Inf, length(theta_values))
    for (tt in seq_along(theta_values)) {
      idx <- available[theta_keys[available] == theta_values[[tt]]]
      choices[tt, ] <- vapply(seq_len(v), function(j) {
        .brt_best(scores[idx, j], idx, tie_tolerance)
      }, integer(1))
      if (!anyNA(choices[tt, ])) {
        objectives[[tt]] <- mean(scores[cbind(choices[tt, ], seq_len(v))])
      }
    }
    best_theta <- .brt_best(objectives, seq_along(theta_values), tie_tolerance)
    selected <- if (is.na(best_theta)) rep(NA_integer_, v) else choices[best_theta, ]
  }
  if (anyNA(selected)) {
    stop("At least one response has no valid candidate under the requested scopes.",
         call. = FALSE)
  }
  selected
}

.brt_metric <- function(observed, predicted, metric) {
  score <- .banded_ridge_score(observed, predicted)
  unname(score[[metric]])
}

.brt_candidate_args <- function(candidates, index) {
  list(alpha = candidates$alpha[[index]],
       theta = setNames(candidates$theta[index, ], candidates$group_names))
}

.brt_fit_candidate <- function(X,
                               Y,
                               train,
                               test,
                               feature_groups,
                               args,
                               solver = "direct",
                               solver_cache = NULL,
                               cache_key = NULL) {
  tryCatch({
    if (solver == "direct") {
      result <- .banded_ridge_fit_predict_direct(
        X, Y, train = train, test = test, groups = feature_groups,
        alpha = args$alpha, theta = args$theta
      )
    } else {
      fit_function <- if (solver == "svd_primal") {
        .banded_ridge_fit_svd_path
      } else {
        .banded_ridge_fit_dual_path
      }
      path <- fit_function(
        .brc_subset_rows(X, train), Y[train, , drop = FALSE],
        alphas = args$alpha, theta = args$theta, groups = feature_groups,
        cache = solver_cache, cache_key = cache_key
      )
      fit <- path$fits[[1L]]
      predictions <- .banded_ridge_predict_optimized(
        fit, .brc_subset_rows(X, test), groups = feature_groups
      )
      rownames(predictions) <- as.character(test)
      result <- list(fit = fit, predictions = predictions,
                     train = train, test = test)
    }
    list(ok = TRUE, result = result, error = NULL)
  }, error = function(e) list(ok = FALSE, result = NULL, error = conditionMessage(e)))
}

.brt_inner_folds <- function(inner_folds,
                             outer_index,
                             outer_train,
                             n,
                             inner_v,
                             blocks,
                             purge,
                             seed) {
  if (is.null(inner_folds)) {
    raw <- .banded_ridge_make_folds(
      outer_train, n_folds = inner_v, blocks = blocks, purge = purge,
      seed = seed + outer_index - 1L
    )
  } else {
    if (!is.list(inner_folds) || length(inner_folds) < outer_index ||
        !is.list(inner_folds[[outer_index]])) {
      stop("inner_folds must contain one fold list per outer fold.", call. = FALSE)
    }
    raw <- inner_folds[[outer_index]]
  }
  .brt_validate_folds(
    raw, allowed_rows = outer_train, n = n,
    context = paste0("inner_folds[[", outer_index, "]]"), blocks = blocks,
    purge = purge, require_coverage = TRUE
  )
}

.brt_preprocessing_receipt <- function(X, Y, fold, feature_groups, group_names) {
  neutral <- setNames(rep(1, length(group_names)), group_names)
  fit <- .banded_ridge_fit_direct(
    .brc_subset_rows(X, fold$train), Y[fold$train, , drop = FALSE],
    groups = feature_groups, lambdas = neutral
  )
  list(fold_id = fold$id, train = fold$train, test = fold$test,
       x_center = fit$x_center, x_scale = fit$x_scale,
       y_center = fit$y_center, constant_x = fit$constant_x)
}

#' Leakage-safe nested per-response tuning for the banded-ridge core
#'
#' Inner candidate scores and choices use only the corresponding outer-training
#' rows. Final reported performance is computed only from predictions assembled
#' once for every outer-test row.
#'
#' @keywords internal
#' @noRd
.banded_ridge_nested_cv <- function(X,
                                    Y,
                                    outer_folds,
                                    candidates,
                                    feature_groups = NULL,
                                    inner_folds = NULL,
                                    inner_v = 5L,
                                    blocks = NULL,
                                    purge = 0L,
                                    metric = c("mse", "correlation", "r2"),
                                    alpha_scope = "response",
                                    theta_scope = "response",
                                    fixed_alpha = NULL,
                                    fixed_theta = NULL,
                                    response_batch_size = Inf,
                                    seed = 1L,
                                    tie_tolerance = 1e-12,
                                    solver = c("direct", "svd_primal", "dual_kernel"),
                                    solver_cache = NULL,
                                    cache_prefix = "nested") {
  metric <- match.arg(metric)
  solver <- match.arg(solver)
  gx <- .brc_grouped_matrix(X, groups = feature_groups, context = "X")
  Y <- .brc_response_matrix(Y, n_expected = nrow(gx$X))
  n <- nrow(Y)
  v <- ncol(Y)
  blocks <- .brt_validate_blocks(blocks, n)
  candidates <- .brt_validate_candidates(candidates, gx$group_names)
  tune <- !(length(inner_v) == 1L && (is.logical(inner_v) || is.integer(inner_v)) &&
              is.na(inner_v))
  if (tune && (!is.numeric(inner_v) || length(inner_v) != 1L ||
               !is.finite(inner_v) || inner_v != as.integer(inner_v) ||
               inner_v < 2L)) {
    stop("inner_v must be an integer of at least two, or NA to skip inner ",
         "tuning when exactly one candidate is available.", call. = FALSE)
  }
  if (!is.numeric(response_batch_size) || length(response_batch_size) != 1L ||
      is.na(response_batch_size) || response_batch_size <= 0 ||
      (!is.infinite(response_batch_size) && response_batch_size != as.integer(response_batch_size))) {
    stop("response_batch_size must be a positive integer or Inf.", call. = FALSE)
  }
  batch_size <- if (is.infinite(response_batch_size)) v else as.integer(response_batch_size)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed != as.integer(seed)) {
    stop("seed must be one finite integer.", call. = FALSE)
  }
  if (!is.character(cache_prefix) || length(cache_prefix) != 1L ||
      is.na(cache_prefix) || !nzchar(cache_prefix)) {
    stop("cache_prefix must be one non-empty string.", call. = FALSE)
  }
  if (solver == "direct" && !is.null(solver_cache)) {
    stop("solver_cache is only used by svd_primal or dual_kernel.", call. = FALSE)
  }
  if (!is.numeric(purge) || length(purge) != 1L || !is.finite(purge) ||
      purge != as.integer(purge) || purge < 0L) {
    stop("purge must be one non-negative integer.", call. = FALSE)
  }
  alpha_scope_canonical <- .brt_scope(alpha_scope, "alpha_scope")
  theta_scope_canonical <- .brt_scope(theta_scope, "theta_scope")
  available <- .brt_available_candidates(
    candidates, alpha_scope_canonical, theta_scope_canonical,
    fixed_alpha = fixed_alpha, fixed_theta = fixed_theta
  )
  if (!tune && length(available) != 1L) {
    stop("inner_v = NA skips inner tuning, which requires exactly one ",
         "candidate under the fixed alpha/theta scopes; ", length(available),
         " are available.", call. = FALSE)
  }
  outer_folds <- .brt_validate_folds(
    outer_folds, allowed_rows = seq_len(n), n = n, context = "outer_folds",
    blocks = blocks, purge = as.integer(purge), require_coverage = TRUE
  )

  oof <- matrix(NA_real_, n, v, dimnames = list(rownames(Y), colnames(Y)))
  if (is.null(rownames(oof))) rownames(oof) <- as.character(seq_len(n))
  outer_results <- vector("list", length(outer_folds))
  selection_rows <- list()
  error_rows <- list()
  max_dimensions <- c(inner_oof_rows = 0L, responses = v, candidates = length(candidates$alpha),
                      inner_folds = 0L, outer_test_rows = 0L)

  for (oo in seq_along(outer_folds)) {
    outer <- outer_folds[[oo]]
    inners <- if (tune) {
      .brt_inner_folds(
        inner_folds, oo, outer$train, n, as.integer(inner_v), blocks,
        as.integer(purge), as.integer(seed)
      )
    } else list()
    c_count <- length(candidates$alpha)
    fold_scores <- array(
      NA_real_, dim = c(length(inners), c_count, v),
      dimnames = list(vapply(inners, `[[`, character(1), "id"),
                      candidates$candidate_id, colnames(Y))
    )
    candidate_scores <- matrix(
      NA_real_, c_count, v,
      dimnames = list(candidates$candidate_id, colnames(Y))
    )
    preprocessing <- lapply(inners, function(fold) {
      .brt_preprocessing_receipt(X, Y, fold, feature_groups, gx$group_names)
    })

    for (cc in if (tune) seq_len(c_count) else integer(0)) {
      args <- .brt_candidate_args(candidates, cc)
      inner_oof <- matrix(NA_real_, length(outer$train), v,
                          dimnames = list(as.character(outer$train), colnames(Y)))
      row_lookup <- setNames(seq_along(outer$train), outer$train)
      candidate_valid <- TRUE
      for (ii in seq_along(inners)) {
        inner <- inners[[ii]]
        evaluated <- .brt_fit_candidate(
          X, Y, inner$train, inner$test, feature_groups, args,
          solver = solver, solver_cache = solver_cache,
          cache_key = paste(cache_prefix, "outer", oo, "inner", ii,
                            "theta", .brt_theta_keys(matrix(
                              args$theta, nrow = 1L
                            )), sep = "::")
        )
        if (!evaluated$ok) {
          candidate_valid <- FALSE
          error_rows[[length(error_rows) + 1L]] <- data.frame(
            outer_fold = outer$id, inner_fold = inner$id,
            candidate_id = candidates$candidate_id[[cc]],
            error = evaluated$error, stringsAsFactors = FALSE
          )
          next
        }
        pred <- evaluated$result$predictions
        fold_scores[ii, cc, ] <- .brt_metric(
          Y[inner$test, , drop = FALSE], pred, metric
        )
        inner_oof[unname(row_lookup[as.character(inner$test)]), ] <- pred
      }
      if (candidate_valid && all(is.finite(inner_oof))) {
        candidate_scores[cc, ] <- .brt_metric(
          Y[outer$train, , drop = FALSE], inner_oof, metric
        )
      }
    }

    selected <- if (tune) {
      objective <- if (metric == "mse") candidate_scores else -candidate_scores
      .brt_select_candidates(
        objective, candidates, alpha_scope_canonical, theta_scope_canonical,
        fixed_alpha = fixed_alpha, fixed_theta = fixed_theta,
        tie_tolerance = tie_tolerance
      )
    } else rep(available, v)

    outer_pred <- matrix(NA_real_, length(outer$test), v,
                         dimnames = list(as.character(outer$test), colnames(Y)))
    outer_fits <- list()
    for (cc in sort(unique(selected))) {
      response_indices <- which(selected == cc)
      batches <- split(response_indices,
                       ceiling(seq_along(response_indices) / batch_size))
      args <- .brt_candidate_args(candidates, cc)
      for (bb in seq_along(batches)) {
        response_index <- batches[[bb]]
        fitted <- .brt_fit_candidate(
          X, Y[, response_index, drop = FALSE], outer$train, outer$test,
          feature_groups, args, solver = solver, solver_cache = solver_cache,
          cache_key = paste(cache_prefix, "outer", oo, "refit", "theta",
                            .brt_theta_keys(matrix(args$theta, nrow = 1L)),
                            sep = "::")
        )
        if (!fitted$ok) {
          stop("Selected candidate failed during outer refit: ", fitted$error,
               call. = FALSE)
        }
        outer_pred[, response_index] <- fitted$result$predictions
        outer_fits[[length(outer_fits) + 1L]] <- list(
          candidate_index = cc,
          candidate_id = candidates$candidate_id[[cc]],
          response_indices = response_index,
          responses = colnames(Y)[response_index],
          fit = fitted$result$fit
        )
      }
    }
    if (any(!is.finite(outer_pred))) {
      stop("Outer predictions were not completely and finitely assembled.", call. = FALSE)
    }
    oof[outer$test, ] <- outer_pred

    selected_table <- data.frame(
      outer_fold = outer$id,
      response = colnames(Y),
      response_index = seq_len(v),
      candidate_index = selected,
      candidate_id = candidates$candidate_id[selected],
      alpha = candidates$alpha[selected],
      inner_score = candidate_scores[cbind(selected, seq_len(v))],
      stringsAsFactors = FALSE
    )
    for (gg in seq_along(gx$group_names)) {
      selected_table[[paste0("theta_", gx$group_names[[gg]])]] <-
        candidates$theta[cbind(selected, gg)]
    }
    selection_rows[[oo]] <- selected_table
    outer_results[[oo]] <- list(
      fold_id = outer$id,
      train = outer$train,
      test = outer$test,
      inner_folds = inners,
      preprocessing = preprocessing,
      inner_fold_scores = fold_scores,
      candidate_scores = candidate_scores,
      selected = selected_table,
      outer_fits = outer_fits,
      predictions = outer_pred
    )
    max_dimensions <- pmax(max_dimensions, c(
      inner_oof_rows = length(outer$train), responses = v,
      candidates = c_count, inner_folds = length(inners),
      outer_test_rows = length(outer$test)
    ))
  }
  if (any(!is.finite(oof))) {
    stop("Outer folds did not produce one finite prediction per response and row.",
         call. = FALSE)
  }
  candidate_errors <- if (length(error_rows)) do.call(rbind, error_rows) else {
    data.frame(outer_fold = character(), inner_fold = character(),
               candidate_id = character(), error = character(),
               stringsAsFactors = FALSE)
  }
  fold_performance <- do.call(rbind, lapply(seq_along(outer_folds), function(oo) {
    rows <- outer_folds[[oo]]$test
    score <- .banded_ridge_score(Y[rows, , drop = FALSE], oof[rows, , drop = FALSE])
    data.frame(outer_fold = outer_folds[[oo]]$id, score,
               stringsAsFactors = FALSE)
  }))
  out <- list(
    predictions = oof,
    performance = .banded_ridge_score(Y, oof),
    fold_performance = fold_performance,
    selections = do.call(rbind, selection_rows),
    outer_results = outer_results,
    candidates = candidates,
    candidate_errors = candidate_errors,
    response_names = colnames(Y),
    group_names = gx$group_names,
    metric = metric,
    inner_tuning = if (tune) "nested" else "none",
    alpha_scope = alpha_scope_canonical,
    theta_scope = theta_scope_canonical,
    fixed_alpha = fixed_alpha,
    fixed_theta = fixed_theta,
    seed = as.integer(seed),
    purge = as.integer(purge),
    response_batch_size = response_batch_size,
    solver = solver,
    tie_breaking = paste0(
      "minimize ", if (metric == "mse") "MSE" else paste0("negative ", metric),
      "; within relative tolerance ", sprintf("%.17g", tie_tolerance),
      ", choose canonical numeric alpha/theta order"
    ),
    peak_intermediate_dimensions = max_dimensions,
    provenance = list(
      outer_fold_ids = vapply(outer_folds, `[[`, character(1), "id"),
      inner_fold_ids = lapply(outer_results, function(x) {
        vapply(x$inner_folds, `[[`, character(1), "id")
      }),
      candidate_generation = attr(candidates, "generation")
    )
  )
  class(out) <- c("banded_ridge_nested_cv", "list")
  out
}
