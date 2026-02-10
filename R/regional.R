#' Coerce a ClusteredNeuroVol to DenseNeuroVol
#'
#' Several regional-analysis helpers rely on \code{as.vector()} and logical
#' indexing (\code{x[x > 0]}), neither of which are defined for
#' \code{ClusteredNeuroVol}.
#' This helper converts to dense when needed, and is a no-op otherwise.
#'
#' @param x A \code{NeuroVol} (or subclass).
#' @return A \code{NeuroVol} that supports \code{as.vector()} and \code{[}.
#' @keywords internal
#' @noRd
.ensure_dense_vol <- function(x) {
  if (inherits(x, "ClusteredNeuroVol")) {
    neuroim2::as.dense(x)
  } else {
    x
  }
}

#' @keywords internal
get_unique_regions.NeuroVol <- function(region_mask, ...) {
  sort(unique(region_mask[region_mask > 0]))
}

#' @keywords internal
get_unique_regions.NeuroSurface <- function(region_mask, ...) {
  sort(unique(region_mask@data[region_mask@data > 0]))
}

#' Combine Regional Results
#'
#' This function combines regional results from a list into a single data frame.
#'
#' @param results A list of regional results.
#' @return A data frame with combined regional results.
#' @details
#' The function is used to combine regional results from a list into a single data frame.
#' It handles both factor and non-factor observed values and creates a combined data frame with the corresponding columns.
#' @keywords internal
#' @noRd
combine_regional_results = function(results) {
  roinum=NULL
  .rownum=NULL

  if (!is.data.frame(results) || !"result" %in% names(results) || nrow(results) == 0) {
    return(tibble::tibble())
  }
  # Check if the observed values are factors (for categorical data)
  if (!is.null(results$result[[1]]) && is.factor(results$result[[1]]$observed)) {
    results %>% dplyr::rowwise() %>% dplyr::do( {
      
      # Use test indices if provided in the classification_result to preserve
      # alignment with the original test design (important when some rows are dropped)
      row_index <- if (!is.null(.$result$testind)) .$result$testind else seq_along(.$result$observed)

      # Create a tibble containing observed, predicted, and additional information
      tib1 <- tibble::tibble(
        .rownum=row_index,
        roinum=rep(.$id, length(row_index)),
        observed=.$result$observed,
        pobserved=sapply(seq_along(.$result$observed), function(i) .$result$probs[i, .$result$observed[i]]),
        predicted=.$result$predicted,
        correct=as.character(.$result$observed) == as.character(.$result$predicted)
      )
      
      # Create a tibble with the probabilities for each class
      tib2 <- tibble::as_tibble(.$result$probs, .name_repair=.name_repair)
      names(tib2) <- paste0("prob_", names(tib2))
      
      # Combine tib1 and tib2
      cbind(tib1, tib2)
    })
  } else if (!is.null(results$result[[1]])) {
    # For non-factor observed values (for continuous data)
    results %>% dplyr::rowwise() %>% dplyr::do({
      row_index <- if (!is.null(.$result$testind)) .$result$testind else seq_along(.$result$observed)
      tibble::tibble(
        .rownum=row_index,
        roinum=rep(.$id, length(row_index)),
        observed=.$result$observed,
        predicted=.$result$predicted)
    })
  } else {
    tibble::tibble()
  }
}

#' Combine prediction tables
#'
#' Combines multiple prediction tables (e.g., from different models or regions) into a single table.
#' Supports weighted combination and collapsing regions.
#'
#' @param predtabs A list of prediction tables (data frames) to be combined.
#' @param wts A vector of weights, with the same length as \code{predtabs}. Default is equal weights.
#' @param collapse_regions A logical value; if TRUE, regions are collapsed in the final prediction table.
#'
#' @return A combined prediction table (data frame).
#' @import dplyr
#' @importFrom purrr map_df
#' @details
#' For classification, this function pools class-probability columns (`prob_*`)
#' and returns pooled predicted class and correctness. For regression (numeric
#' `observed`), it pools `predicted` values and returns residual diagnostics.
#' @examples
#' # Create example prediction tables
#' observed = factor(sample(letters[1:2], 10, replace = TRUE))
#' predtab1 <- data.frame(.rownum = 1:10,
#'                        roinum = rep(1, 10),
#'                        observed = observed,
#'                        prob_A = runif(10),
#'                        prob_B = runif(10))
#' predtab2 <- data.frame(.rownum = 1:10,
#'                        roinum = rep(2, 10),
#'                        observed = observed,
#'                        prob_A = runif(10),
#'                        prob_B = runif(10))
#'
#' # Combine the tables
#' combined_table <- combine_prediction_tables(list(predtab1, predtab2))
#' @export
combine_prediction_tables <- function(predtabs, wts=rep(1,length(predtabs)), collapse_regions=FALSE) {
  assert_that(length(predtabs) == length(wts))
  assert_that(sum(wts) > 0)
  assert_that(all(purrr::map_lgl(predtabs, is.data.frame)))

  wts <- wts/sum(wts)

  .weight <- NULL
  .rownum <- NULL
  roinum <-  NULL
  observed <- NULL
  predicted <- NULL
  residual <- NULL
  abs_error <- NULL
  sq_error <- NULL

  ## applies constant weight to each table and concatenates
  ptab <- map(seq_along(predtabs), function(i) predtabs[[i]] %>% mutate(.tableid=i, .weight=wts[i])) %>%
    map_df(bind_rows) %>% as_tibble(.name_repair=.name_repair)

  if (is.character(predtabs[[1]]$observed) || is.factor(predtabs[[1]]$observed)) {
    probs <- if (collapse_regions) {
      ptab %>% dplyr::group_by(.rownum, observed) %>% summarise_at(
        vars(starts_with("prob_")),
        list(~stats::weighted.mean(., w = .weight))
      ) %>% dplyr::ungroup()
    } else {
      ## groups over rownames, condition, and roinum, then compute weighted means of probabilities
      ptab %>% dplyr::group_by(.rownum, observed, roinum) %>% summarise_at(
        vars(starts_with("prob_")),
        list(~stats::weighted.mean(., w = .weight))
      ) %>% dplyr::ungroup()
    }

    p <- probs %>% dplyr::select(dplyr::starts_with("prob_"))
    pmat <- as.matrix(p)
    pobserved <- pmat[cbind(seq_len(nrow(probs)), as.integer(probs$observed))]
    mc <- max.col(pmat)
    preds <- levels(probs$observed)[mc]

    prediction_table <- tibble(
      .rownum = probs$.rownum,
      roinum = if (collapse_regions) 1 else probs$roinum,
      observed = probs$observed,
      pobserved = pobserved,
      predicted = preds,
      correct = preds == as.character(probs$observed)
    ) %>% bind_cols(p)
  } else if (is.numeric(predtabs[[1]]$observed)) {
    pred_avg <- if (collapse_regions) {
      ptab %>%
        dplyr::group_by(.rownum, observed) %>%
        dplyr::summarise(
          predicted = stats::weighted.mean(predicted, w = .weight),
          .groups = "drop"
        ) %>%
        dplyr::mutate(roinum = 1L)
    } else {
      ptab %>%
        dplyr::group_by(.rownum, observed, roinum) %>%
        dplyr::summarise(
          predicted = stats::weighted.mean(predicted, w = .weight),
          .groups = "drop"
        )
    }

    prediction_table <- pred_avg %>%
      dplyr::mutate(
        residual = observed - predicted,
        abs_error = abs(residual),
        sq_error = residual^2
      ) %>%
      dplyr::select(.rownum, roinum, observed, predicted, residual, abs_error, sq_error)
  } else {
    stop("combine_prediction_tables: `observed` must be factor/character or numeric.")
  }

  prediction_table
}

#' @keywords internal
#' @noRd
.classification_perf_from_table <- function(prediction_table) {
  prob_cols <- grep("^prob_", names(prediction_table), value = TRUE)
  if (length(prob_cols) == 0L) {
    stop("classification performance requires `prob_` columns in prediction_table.")
  }

  observed <- prediction_table$observed
  if (!is.factor(observed)) {
    observed <- factor(observed)
  }
  lvls <- levels(observed)
  probs <- as.matrix(prediction_table[, prob_cols, drop = FALSE])

  current_prob_levels <- sub("^prob_", "", prob_cols)
  if (!all(lvls %in% current_prob_levels)) {
    stop("classification performance could not match probability columns to observed levels.")
  }
  probs <- probs[, match(lvls, current_prob_levels), drop = FALSE]
  colnames(probs) <- lvls

  predicted <- prediction_table$predicted
  predicted <- factor(predicted, levels = lvls)
  cres <- classification_result(
    observed = observed,
    predicted = predicted,
    probs = probs,
    testind = prediction_table$.rownum
  )
  performance(cres)
}

#' @keywords internal
#' @noRd
.regression_perf_from_table <- function(prediction_table) {
  rres <- regression_result(
    observed = as.numeric(prediction_table$observed),
    predicted = as.numeric(prediction_table$predicted),
    testind = prediction_table$.rownum
  )
  performance(rres)
}

#' @keywords internal
#' @noRd
.compute_pooled_performance <- function(prediction_table) {
  if (is.factor(prediction_table$observed) || is.character(prediction_table$observed)) {
    .classification_perf_from_table(prediction_table)
  } else if (is.numeric(prediction_table$observed)) {
    .regression_perf_from_table(prediction_table)
  } else {
    stop(".compute_pooled_performance: unsupported observed type.")
  }
}

#' @keywords internal
#' @noRd
.stacking_frame <- function(prediction_table) {
  .rownum <- NULL
  observed <- NULL
  roinum <- NULL
  .class <- NULL
  .feature <- NULL
  .value <- NULL
  predicted <- NULL

  obs <- prediction_table %>%
    dplyr::select(.rownum, observed) %>%
    dplyr::distinct(.rownum, .keep_all = TRUE)

  if (is.factor(obs$observed) || is.character(obs$observed)) {
    prob_cols <- grep("^prob_", names(prediction_table), value = TRUE)
    if (length(prob_cols) == 0L) {
      stop("stacking for classification requires `prob_` columns.")
    }

    xwide <- prediction_table %>%
      dplyr::select(.rownum, roinum, dplyr::all_of(prob_cols)) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(prob_cols),
        names_to = ".class",
        values_to = ".value"
      ) %>%
      dplyr::mutate(.feature = paste0("roi", roinum, "__", .class)) %>%
      dplyr::select(.rownum, .feature, .value) %>%
      tidyr::pivot_wider(names_from = .feature, values_from = .value)

    dat <- dplyr::left_join(obs, xwide, by = ".rownum") %>%
      dplyr::arrange(.rownum)

    y <- dat$observed
    if (!is.factor(y)) {
      y <- factor(y)
    }
    x <- as.matrix(dat[, setdiff(names(dat), c(".rownum", "observed")), drop = FALSE])
    list(.rownum = dat$.rownum, y = y, x = x)
  } else if (is.numeric(obs$observed)) {
    xwide <- prediction_table %>%
      dplyr::select(.rownum, roinum, predicted) %>%
      dplyr::mutate(.feature = paste0("roi", roinum)) %>%
      dplyr::select(.rownum, .feature, predicted) %>%
      tidyr::pivot_wider(names_from = .feature, values_from = predicted)

    dat <- dplyr::left_join(obs, xwide, by = ".rownum") %>%
      dplyr::arrange(.rownum)

    x <- as.matrix(dat[, setdiff(names(dat), c(".rownum", "observed")), drop = FALSE])
    list(.rownum = dat$.rownum, y = as.numeric(dat$observed), x = x)
  } else {
    stop(".stacking_frame: unsupported observed type.")
  }
}

#' @keywords internal
#' @noRd
.normalize_fold_list <- function(fold_list, rownum) {
  out <- lapply(fold_list, as.integer)
  all_idx <- sort(unique(unlist(out, use.names = FALSE)))
  n <- length(rownum)

  if (all(all_idx %in% seq_len(n))) {
    return(out)
  }

  if (all(all_idx %in% rownum)) {
    return(lapply(out, function(idx) {
      mapped <- match(idx, rownum)
      mapped[!is.na(mapped)]
    }))
  }

  stop("stack_folds list must be indices in 1..nrow(unique(.rownum)) or actual `.rownum` ids.")
}

#' @keywords internal
#' @noRd
.derive_stack_folds <- function(model_spec, rownum, y, stack_folds = NULL, stack_seed = NULL) {
  n <- length(rownum)
  if (n < 2L) {
    return(list(seq_len(n)))
  }

  if (is.list(stack_folds)) {
    folds <- .normalize_fold_list(stack_folds, rownum)
    return(Filter(function(idx) length(idx) > 0L, folds))
  }

  if (is.numeric(stack_folds) && length(stack_folds) == n) {
    fold_ids <- as.integer(stack_folds)
    fold_ids[is.na(fold_ids)] <- max(fold_ids, na.rm = TRUE) + seq_len(sum(is.na(fold_ids)))
    folds <- split(seq_len(n), fold_ids)
    return(Filter(function(idx) length(idx) > 0L, folds))
  }

  block_var <- tryCatch(model_spec$design$block_var, error = function(...) NULL)
  if (!is.null(block_var) && length(block_var) >= max(rownum, na.rm = TRUE)) {
    bsub <- block_var[rownum]
    if (length(unique(bsub[!is.na(bsub)])) >= 2L) {
      folds <- split(seq_len(n), as.integer(as.factor(bsub)))
      folds <- Filter(function(idx) length(idx) > 0L, folds)
      if (length(folds) >= 2L) {
        return(folds)
      }
    }
  }

  k <- if (length(stack_folds) == 1L && !is.na(stack_folds)) as.integer(stack_folds) else min(5L, n)
  k <- max(2L, min(k, n))
  create_mvpa_folds(y, k = k, list = TRUE, seed = stack_seed)
}

#' @keywords internal
#' @noRd
.impute_train_test <- function(x_train, x_test) {
  if (ncol(x_train) == 0L) {
    return(list(train = x_train, test = x_test))
  }
  col_means <- colMeans(x_train, na.rm = TRUE)
  col_means[!is.finite(col_means)] <- 0

  train_na <- which(is.na(x_train), arr.ind = TRUE)
  if (nrow(train_na) > 0L) {
    x_train[train_na] <- col_means[train_na[, 2]]
  }

  test_na <- which(is.na(x_test), arr.ind = TRUE)
  if (nrow(test_na) > 0L) {
    x_test[test_na] <- col_means[test_na[, 2]]
  }

  list(train = x_train, test = x_test)
}

#' @keywords internal
#' @noRd
.crossfit_stack_predictions <- function(prediction_table, model_spec,
                                        stack_folds = NULL, stack_seed = NULL,
                                        stack_lambda = 1e-3) {
  sf <- .stacking_frame(prediction_table)
  y <- sf$y
  x <- sf$x
  rownum <- sf$.rownum
  n <- NROW(x)
  if (n == 0L) {
    stop(".crossfit_stack_predictions: no rows available for stacking.")
  }

  folds <- .derive_stack_folds(model_spec, rownum, y, stack_folds, stack_seed)
  folds <- Filter(function(idx) length(idx) > 0L, folds)
  if (length(folds) < 2L) {
    stop(".crossfit_stack_predictions: need at least 2 folds for cross-fitted stacking.")
  }

  if (is.factor(y)) {
    lvls <- levels(y)
    k <- length(lvls)
    prob_pred <- matrix(NA_real_, nrow = n, ncol = k, dimnames = list(NULL, lvls))
    prior <- as.numeric(table(y)[lvls]) / length(y)

    for (te in folds) {
      tr <- setdiff(seq_len(n), te)
      if (length(tr) < 2L || ncol(x) == 0L) {
        prob_pred[te, ] <- matrix(prior, nrow = length(te), ncol = k, byrow = TRUE)
        next
      }

      xtr <- x[tr, , drop = FALSE]
      xte <- x[te, , drop = FALSE]
      imp <- .impute_train_test(xtr, xte)
      xtr <- imp$train
      xte <- imp$test

      ytr <- y[tr]
      present <- lvls[lvls %in% unique(as.character(ytr))]
      if (length(present) < 2L) {
        local_prior <- rep(0, k)
        local_prior[match(present, lvls)] <- 1
        prob_pred[te, ] <- matrix(local_prior, nrow = length(te), ncol = k, byrow = TRUE)
        next
      }

      if (k == 2L) {
        ybin <- as.integer(ytr == lvls[2])
        fit <- glmnet::glmnet(
          x = xtr,
          y = ybin,
          family = "binomial",
          alpha = 0,
          lambda = stack_lambda,
          standardize = TRUE,
          intercept = TRUE
        )
        p2 <- as.numeric(stats::predict(fit, newx = xte, s = stack_lambda, type = "response"))
        p2 <- pmin(pmax(p2, 0), 1)
        prob_pred[te, 1] <- 1 - p2
        prob_pred[te, 2] <- p2
      } else {
        fit <- glmnet::glmnet(
          x = xtr,
          y = ytr,
          family = "multinomial",
          alpha = 0,
          lambda = stack_lambda,
          standardize = TRUE,
          intercept = TRUE
        )
        p <- stats::predict(fit, newx = xte, s = stack_lambda, type = "response")
        p <- if (length(dim(p)) == 3L) p[, , 1, drop = FALSE] else p
        if (length(dim(p)) == 3L) {
          p <- p[, , 1]
        }
        p <- as.matrix(p)
        if (is.null(colnames(p))) {
          colnames(p) <- lvls[seq_len(ncol(p))]
        }
        pfull <- matrix(0, nrow = nrow(p), ncol = k, dimnames = list(NULL, lvls))
        keep <- intersect(colnames(p), lvls)
        pfull[, keep] <- p[, keep, drop = FALSE]
        rs <- rowSums(pfull)
        rs[rs <= 0 | !is.finite(rs)] <- 1
        prob_pred[te, ] <- pfull / rs
      }
    }

    still_na <- which(!is.finite(prob_pred), arr.ind = TRUE)
    if (nrow(still_na) > 0L) {
      prob_pred[still_na] <- prior[still_na[, 2]]
    }
    rs <- rowSums(prob_pred)
    rs[rs <= 0 | !is.finite(rs)] <- 1
    prob_pred <- prob_pred / rs

    pred_class <- lvls[max.col(prob_pred, ties.method = "first")]
    pobs <- prob_pred[cbind(seq_len(n), as.integer(y))]
    prob_tbl <- tibble::as_tibble(prob_pred) %>%
      dplyr::rename_with(~paste0("prob_", .x))

    tibble(
      .rownum = rownum,
      roinum = 1L,
      observed = y,
      pobserved = pobs,
      predicted = pred_class,
      correct = pred_class == as.character(y)
    ) %>% bind_cols(prob_tbl)
  } else {
    pred <- rep(NA_real_, n)
    global_mean <- mean(y, na.rm = TRUE)

    for (te in folds) {
      tr <- setdiff(seq_len(n), te)
      if (length(tr) < 2L || ncol(x) == 0L) {
        pred[te] <- mean(y[tr], na.rm = TRUE)
        next
      }

      xtr <- x[tr, , drop = FALSE]
      xte <- x[te, , drop = FALSE]
      imp <- .impute_train_test(xtr, xte)
      xtr <- imp$train
      xte <- imp$test

      fit <- glmnet::glmnet(
        x = xtr,
        y = y[tr],
        family = "gaussian",
        alpha = 0,
        lambda = stack_lambda,
        standardize = TRUE,
        intercept = TRUE
      )
      pred[te] <- as.numeric(stats::predict(fit, newx = xte, s = stack_lambda, type = "response"))
    }

    pred[!is.finite(pred)] <- global_mean
    residual <- y - pred
    tibble(
      .rownum = rownum,
      roinum = 1L,
      observed = y,
      predicted = pred,
      residual = residual,
      abs_error = abs(residual),
      sq_error = residual^2
    )
  }
}

#' @keywords internal
#' @noRd
.pool_regional_predictions <- function(prediction_table, model_spec,
                                       method = c("none", "mean", "stack"),
                                       stack_folds = NULL, stack_seed = NULL,
                                       stack_lambda = 1e-3, pooled_weights = NULL) {
  method <- match.arg(method)
  if (method == "none" || is.null(prediction_table) || nrow(prediction_table) == 0L) {
    return(list(prediction_table = NULL, performance = NULL))
  }

  predtabs <- split(prediction_table, prediction_table$roinum)
  weights <- if (is.null(pooled_weights)) rep(1, length(predtabs)) else pooled_weights
  if (!is.null(pooled_weights) && length(pooled_weights) != length(predtabs)) {
    stop("`pooled_weights` must have one weight per ROI in `prediction_table`.")
  }

  pooled_tab <- if (method == "mean") {
    combine_prediction_tables(predtabs, wts = weights, collapse_regions = TRUE)
  } else {
    .crossfit_stack_predictions(
      prediction_table = prediction_table,
      model_spec = model_spec,
      stack_folds = stack_folds,
      stack_seed = stack_seed,
      stack_lambda = stack_lambda
    )
  }

  list(
    prediction_table = pooled_tab,
    performance = .compute_pooled_performance(pooled_tab)
  )
}


#'
#' Merge multiple regional MVPA results into a single result.
#'
#' @param obj A \code{regional_mvpa_result} object.
#' @param ... Additional \code{regional_mvpa_result} objects to be merged.
#'
#' @return A merged \code{regional_mvpa_result} object.
#' @rdname merge_results-methods
#' @method merge_results regional_mvpa_result
#' @export
merge_results.regional_mvpa_result <- function(obj, ...) {
  rlist <- list(obj, ...)
  combine_prediction_tables(rlist)
}


#' Create a \code{regional_mvpa_result} instance
#'
#' Constructs a regional MVPA result object that stores the results of MVPA analysis in a specific region.
#'
#' @param model_spec A model specification object.
#' @param performance_table A data frame with performance measures.
#' @param prediction_table A data frame with prediction results.
#' @param vol_results A list of voxel-level results.
#' @param fits Optional model fits.
#' @param pooled_prediction_table Optional pooled prediction table (e.g., ROI-averaged or stacked).
#' @param pooled_performance Optional named numeric vector of pooled performance metrics.
#'
#' @return A \code{regional_mvpa_result} object.
#' @examples
#' # Create example inputs
#' model_spec <- list(dataset = "Example dataset")
#' performance_table <- data.frame(accuracy = c(0.8, 0.85))
#' prediction_table <- data.frame(observed = factor(rep(letters[1:2], 5)),
#'                                 predicted = factor(rep(letters[1:2], 5)))
#' vol_results <- list(vol1 = "Example vol_result 1", vol2 = "Example vol_result 2")
#' fits <- list(fit1 = "Example fit 1", fit2 = "Example fit 2")
#'
#' # Construct a regional_mvpa_result
#' regional_result <- regional_mvpa_result(model_spec, performance_table,
#'                                         prediction_table, vol_results, fits = fits)
#' @export
regional_mvpa_result <- function(model_spec, performance_table, prediction_table, vol_results, fits=fits,
                                 pooled_prediction_table = NULL, pooled_performance = NULL) {
  ret <- list(model_spec=model_spec, 
              performance_table=performance_table,
              prediction_table=prediction_table,
              vol_results=vol_results,
              fits=fits,
              pooled_prediction_table = pooled_prediction_table,
              pooled_performance = pooled_performance)
  
  class(ret) <- c("regional_mvpa_result", "list")
  ret
  
}

#' Prepare regional data for MVPA analysis
#'
#' This function processes the input data and prepares the regions for MVPA analysis by extracting
#' voxel indices for each region of interest (ROI) specified in the region_mask.
#'
#' @param model_spec A model specification object.
#' @param region_mask A mask representing different regions in the brain image.
#'
#' @return A list containing information about the regions for further processing:
#'   * allrois: A vector of unique ROI labels.
#'   * region_vec: A vector representation of the region_mask.
#'   * region_set: A sorted vector of unique ROI labels in the region_mask.
#'   * vox_iter: A list containing voxel indices for each ROI.
#'   * lens: A vector containing the number of voxels in each ROI.
#'   * keep: A logical vector indicating if an ROI should be kept for analysis (those with more than one voxel).
#'
#' @examples
#' # Create example data
#' sample_data <- gen_sample_dataset(c(5, 5, 5), nobs = 100, blocks = 4)
#'
#' # Create a simple region mask with 3 ROIs
#' mask_vol <- sample_data$dataset$mask
#' region_mask <- neuroim2::NeuroVol(
#'   sample(1:3, size = sum(mask_vol > 0), replace = TRUE),
#'   space = neuroim2::space(mask_vol),
#'   indices = which(mask_vol > 0)
#' )
#'
#' # Create a basic model spec
#' model_spec <- list(dataset = sample_data$dataset)
#'
#' # Prepare regional data
#' regional_data <- prep_regional(model_spec, region_mask)
#' @export
prep_regional <- function(model_spec, region_mask) {
  region_mask <- .ensure_dense_vol(region_mask)
  allrois <- get_unique_regions(region_mask)
  ##allrois <- sort(unique(region_mask[region_mask>0]))
  region_vec <- as.vector(region_mask)
  region_set <- sort(as.integer(unique(region_vec[region_vec > 0])))
  if (length(region_set) < 1) {
    stop("run_regional: invalid ROI mask, number of ROIs = 0")
  }

  # Intersect region_mask with dataset$mask: only voxels present in BOTH
  # are included.  Voxels that are in a region but outside the data mask
  # (e.g. no valid fMRI signal) are silently dropped here.
  mask_vec <- as.vector(model_spec$dataset$mask)
  roi_lens  <- sapply(region_set, function(rnum) sum(region_vec == rnum))
  vox_iter  <- lapply(region_set, function(rnum) which(region_vec == rnum & mask_vec > 0))
  lens      <- sapply(vox_iter, length)

  # Report voxels lost to data-mask exclusion
  lost <- roi_lens - lens
  if (any(lost > 0)) {
    affected <- which(lost > 0)
    detail   <- sprintf("ROI %d: %d/%d voxels outside data mask",
                        region_set[affected], lost[affected], roi_lens[affected])
    futile.logger::flog.info(
      "prep_regional: %d of %d ROIs lost voxels to the data mask:\n  %s",
      length(affected), length(region_set), paste(detail, collapse = "\n  "))
  }

  keep <- lens > 1

  if (all(!keep)) {
    futile.logger::flog.error("run_regional: no ROIs have more than one voxel after intersecting with the data mask.")
    stop("run_regional: no ROIs have more than one voxel.")
  }

  if (any(lens < 2)) {
    futile.logger::flog.warn(
      "prep_regional: removing %d ROI(s) with < 2 voxels after mask intersection: %s",
      sum(!keep), paste(region_set[!keep], collapse = ", "))
    vox_iter   <- vox_iter[keep]
    region_set <- region_set[keep]
  }

  list(allrois=allrois, region_vec=region_vec, region_set=region_set,
       vox_iter=vox_iter, lens=lens, keep=keep)
  
}

#' Compile performance results and generate volumetric results
#'
#' This function compiles the performance results from the regional MVPA analysis
#' and generates volumetric results based on the region_mask.
#'
#' @param results A list containing the results from regional MVPA analysis.
#' @param region_mask A mask representing different regions in the brain image.
#'
#' @return A list containing:
#'   * vols: A list of volumetric results for each performance metric.
#'   * perf_mat: A data frame containing the compiled performance results with an added 'roinum' column.
#'
#' @importFrom neuroim2 map_values
#' @keywords internal
#' @noRd
comp_perf <- function(results, region_mask) {
  region_mask <- .ensure_dense_vol(region_mask)
  roinum <- NULL

  # If results is not a data frame/tibble or lacks expected columns, return empty
  if (!is.data.frame(results) || !all(c("performance", "id") %in% names(results))) {
    return(list(vols = list(), perf_mat = tibble::tibble(roinum = integer(0))))
  }
  ## compile performance results
  perf_mat <- tryCatch({
    do.call(rbind, results$performance)
  }, error = function(e) {
    message("Warning: Error creating performance matrix: ", e$message)
    return(NULL)
  })
  
  # Ensure we keep original names, make unique if duplicates exist
  perf_mat <- as_tibble(perf_mat, .name_repair = "unique")
  
  # Check if perf_mat is NULL or has 0 columns
  if (is.null(perf_mat) || !is.data.frame(perf_mat) || ncol(perf_mat) == 0) {
    message("Warning: Performance matrix is empty or invalid. Returning empty results.")
    id_vec <- if ("id" %in% names(results)) unlist(results$id) else integer(0)
    return(list(vols = list(), perf_mat = tibble::tibble(roinum = id_vec)))
  }
  
  ## generate volumetric results
  ## TODO fails when region_mask is an logical vol
  id_vec <- if ("id" %in% names(results)) as.integer(unlist(results$id)) else integer(0)
  vols <- lapply(seq_len(ncol(perf_mat)), function(i) {
    if (length(id_vec) == 0L) {
      return(NULL)
    }
    map_values(region_mask, cbind(id_vec, perf_mat[[i]]))
  })
  names(vols) <- names(perf_mat)
  
  perfmat <- tibble::as_tibble(perf_mat,.name_repair=.name_repair) %>% dplyr::mutate(roinum = unlist(results$id)) %>% dplyr::select(roinum, dplyr::everything())
  list(vols=vols, perf_mat=perfmat)
}

#' @rdname run_regional-methods
#' @param compute_performance Logical indicating whether to compute performance metrics (defaults to \code{model_spec$compute_performance}).
#' @param return_predictions Logical indicating whether to combine a full prediction table (defaults to \code{model_spec$return_predictions}).
#' @param return_fits Logical indicating whether to return the fitted models (defaults to \code{model_spec$return_fits}).
#' @param pool_predictions Character scalar controlling pooled outputs:
#'   \code{"none"} (default), \code{"mean"} (weighted mean pooling across ROIs),
#'   or \code{"stack"} (cheap cross-fitted stacking on OOF ROI predictions).
#' @param pooled_weights Optional numeric vector of ROI weights used when
#'   \code{pool_predictions = "mean"}. Must have one weight per ROI.
#' @param stack_folds Optional fold specification for stacking. Can be:
#'   an integer number of folds, a fold-id vector, or a list of fold indices.
#' @param stack_seed Optional seed used when auto-generating stacking folds.
#' @param stack_lambda Ridge penalty used by \pkg{glmnet} in stacking.
#' @details This function serves as the base implementation for regional analyses, orchestrating data preparation, iteration over regions, performance computation, and result aggregation. Specific `run_regional` methods for different model classes may call this function or provide specialized behavior.
#' @export
run_regional_base <- function(model_spec,
                              region_mask,
                              coalesce_design_vars = FALSE,
                              processor = NULL,
                              verbose = FALSE,
                              compute_performance = model_spec$compute_performance,
                              return_predictions = model_spec$return_predictions,
                              return_fits = model_spec$return_fits,
                              pool_predictions = c("none", "mean", "stack"),
                              pooled_weights = NULL,
                              stack_folds = NULL,
                              stack_seed = NULL,
                              stack_lambda = 1e-3,
                              ...) {

  pool_predictions <- match.arg(pool_predictions)
 
  # 1) Prepare regions
  prepped <- prep_regional(model_spec, region_mask)
  
  # 2) Iterate over regions
  results <- mvpa_iterate(
    model_spec,
    prepped$vox_iter,
    ids = prepped$region_set,
    processor = processor,
    verbose = verbose,
    analysis_type = "regional",
    ...
  )
  
  # 3) Performance computation
  perf <- if (isTRUE(compute_performance)) {
    comp_perf(results, region_mask)
  } else {
    list(vols = list(), perf_mat = tibble::tibble())
  }
  
  # 4) Predictions
  prediction_table <- NULL
  if (isTRUE(return_predictions)) {
    prediction_table <- combine_regional_results(results)
    
    if (coalesce_design_vars &&
        !is.null(prediction_table) &&
        nrow(prediction_table) > 0L &&
        ".rownum" %in% names(prediction_table)) {
      prediction_table <- coalesce_join(
        prediction_table,
        test_design(model_spec$design),
        by = ".rownum"
      )
    }
  }

  pooled_prediction_table <- NULL
  pooled_performance <- NULL
  if (isTRUE(return_predictions) &&
      !is.null(prediction_table) &&
      nrow(prediction_table) > 0L &&
      pool_predictions != "none") {
    pooled <- .pool_regional_predictions(
      prediction_table = prediction_table,
      model_spec = model_spec,
      method = pool_predictions,
      stack_folds = stack_folds,
      stack_seed = stack_seed,
      stack_lambda = stack_lambda,
      pooled_weights = pooled_weights
    )
    pooled_prediction_table <- pooled$prediction_table
    pooled_performance <- pooled$performance
  }
  
  # 5) Fits
  fits <- NULL
  if (isTRUE(return_fits)) {
    fits <- lapply(results$result, "[[", "predictor")
  }
  
  # 6) Construct and return final result
  regional_mvpa_result(
    model_spec        = model_spec,
    performance_table = perf$perf_mat,
    prediction_table  = prediction_table,
    vol_results       = perf$vols,
    fits             = fits,
    pooled_prediction_table = pooled_prediction_table,
    pooled_performance = pooled_performance
  )
}


#' Default Method for run_regional
#'
#' @rdname run_regional-methods
#' @details This is the fallback method called when no specialized `run_regional` method is found for the class of `model_spec`. It typically calls `run_regional_base`.
#' @export
run_regional.default <- function(model_spec, region_mask, ...) {
  run_regional_base(model_spec, region_mask, ...)
}


#' Regional MVPA for `mvpa_model` Objects
#'
#' @rdname run_regional-methods
#' @details This method provides the standard regional analysis pipeline for objects of class `mvpa_model` by calling `run_regional_base`.
#' @export
run_regional.mvpa_model <- function(model_spec, region_mask,
                                    coalesce_design_vars = FALSE,
                                    processor = NULL,
                                    verbose = FALSE,
                                    ...) {
  
  run_regional_base(
    model_spec,
    region_mask,
    coalesce_design_vars = coalesce_design_vars,
    processor = processor,
    verbose = verbose,
    ...
  )
}


#' Regional MVPA for `rsa_model` Objects
#'
#' @rdname run_regional-methods
#' @param return_fits Whether to return each region's fitted model (default \code{FALSE}).
#' @param compute_performance \code{logical} indicating whether to compute performance metrics (default \code{TRUE}).
#' @details For `rsa_model` objects, `return_predictions` defaults to `FALSE` as standard RSA typically doesn't produce a prediction table in the same way as classification/regression models.
#' @export
run_regional.rsa_model <- function(model_spec, region_mask,
                                   return_fits = FALSE,
                                   compute_performance = TRUE,
                                   coalesce_design_vars = FALSE,
                                   ...) {
  
  run_regional_base(
    model_spec,
    region_mask,
    coalesce_design_vars  = coalesce_design_vars,
    compute_performance   = compute_performance,
    return_fits           = return_fits,
    return_predictions    = FALSE,  # Override default for RSA
    ...
  )
}

#' Regional MVPA for `vector_rsa_model` Objects
#'
#' @rdname run_regional-methods
#' @param return_fits Logical indicating whether to return the fitted models (default \code{FALSE}).
#' @param compute_performance Logical indicating whether to compute performance metrics (default \code{TRUE}).
#' @details For `vector_rsa_model` objects, `return_predictions` defaults to `FALSE` in `run_regional_base`.
#' If `model_spec$return_predictions` is TRUE, this method will assemble an `observation_scores_table`.
#' @importFrom dplyr bind_rows rename mutate row_number left_join
#' @importFrom tidyr unnest
#' @export
run_regional.vector_rsa_model <- function(model_spec, region_mask,
                                         return_fits = FALSE,
                                         compute_performance = TRUE,
                                         coalesce_design_vars = FALSE, # Usually FALSE for RSA
                                         processor = NULL,
                                         verbose = FALSE,
                                         ...) {
  
  # 1) Prepare regions (using base helper)
  prepped <- prep_regional(model_spec, region_mask)
  
  # 2) Iterate over regions using mvpa_iterate
  # The result from merge_results.vector_rsa_model will contain:
  # - performance: list column with the summary performance matrix
  # - result: list column containing list(rsa_scores=scores_vector) or NULL
  iteration_results <- mvpa_iterate(
    model_spec,
    prepped$vox_iter,
    ids = prepped$region_set,
    processor = processor, # Use default processor unless specified
    verbose = verbose,
    analysis_type = "regional",
    ...
  )
  
  # 3) Performance computation (using base helper)
  # This extracts the 'performance' column from iteration_results
  perf <- if (isTRUE(compute_performance)) {
    comp_perf(iteration_results, region_mask)
  } else {
    list(vols = list(), perf_mat = tibble::tibble())
  }
  
  # 4) Assemble observation scores (if requested)
  prediction_table <- NULL
  if (isTRUE(model_spec$return_predictions) && "result" %in% names(iteration_results)) {
    # Filter out NULL results (where return_predictions was FALSE or errors occurred)
    valid_results <- iteration_results[!sapply(iteration_results$result, is.null), ]
    
    if (nrow(valid_results) > 0) {
      # Create a tibble: roinum | rsa_scores_list
      scores_data <- tibble::tibble(
          roinum = valid_results$id, 
          scores_list = lapply(valid_results$result, function(res) res$rsa_scores)
      )
      
      # Unnest to get a long table: roinum | observation_index | rsa_score
      prediction_table <- scores_data %>%
           mutate(observation_index = map(scores_list, seq_along)) %>% # Add observation index within ROI
           tidyr::unnest(cols = c(scores_list, observation_index)) %>% 
           dplyr::rename(rsa_score = scores_list) # Rename the scores column
           
       # Optionally merge design variables (might need adjustment based on score indices)
       if (coalesce_design_vars) {
            # We need a way to map observation_index back to the original design .rownum
            # This assumes scores are in the same order as the original y_train 
            # (which `second_order_similarity` preserves)
            # Need the original design dataframe 
            orig_design <- model_spec$design$design_table # Assuming it's stored here? Check mvpa_design
            if (!is.null(orig_design)) {
                # Add .rownum based on the original sequence
                # This relies on the assumption that the number of scores matches nrow(orig_design)
                num_obs_in_design <- nrow(orig_design)
                prediction_table <- prediction_table %>%
                   # Need to handle potential mismatch if scores length != num_obs_in_design
                   # For now, assume they match and add .rownum directly
                   dplyr::mutate(.rownum = observation_index) %>%
                   # Perform the join
                   coalesce_join(orig_design, by = ".rownum")
            } else {
                 warning("coalesce_design_vars=TRUE but original design table not found in model_spec$design$design_table")
            }
       }
           
    } else {
         warning("return_predictions=TRUE, but no observation scores were returned from processing.")
    }
  }
  
  # 5) Fits (using base logic - check if applicable for vector_rsa)
  # train_model returns scores, not a fit object, so fits will likely be NULL
  fits <- NULL
  if (isTRUE(return_fits)) {
      # The `result` column now holds scores, not fits. This needs reconsideration.
      # fits <- lapply(iteration_results$result, "[[<some_fit_element>") # This won't work
      warning("`return_fits=TRUE` requested for vector_rsa_model, but this model type does not currently return standard fit objects.")
  }
  
  # 6) Construct and return final result (using base constructor)
  regional_mvpa_result(
    model_spec        = model_spec,
    performance_table = perf$perf_mat,
    prediction_table  = prediction_table, # Add the assembled scores table
    vol_results       = perf$vols,
    fits             = fits
  )
}

#' Regional MVPA for `feature_rsa_model` Objects
#'
#' @rdname run_regional-methods
#' @details For `feature_rsa_model` objects, `return_predictions` defaults to `FALSE`
#'   (set at model creation via `create_model_spec`). This method delegates to
#'   `run_regional_base`.
#' @export
run_regional.feature_rsa_model <- function(model_spec, region_mask,
                                           coalesce_design_vars = FALSE,
                                           processor = NULL,
                                           verbose = FALSE,
                                           ...) {
  run_regional_base(
    model_spec,
    region_mask,
    coalesce_design_vars = coalesce_design_vars,
    processor = processor,
    verbose = verbose,
    ...
  )
}
