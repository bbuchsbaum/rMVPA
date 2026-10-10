#' @keywords internal
gen_id <- function(n) {
  width <- max(2, nchar(n))  # Minimum width of 2 for consistent formatting
  sprintf(paste0("%0", width, "d"), seq_len(n))
}

#' @keywords internal
#' @noRd
.fast_filter_roi_enabled <- function() {
  .fast_path_enabled("fast_filter_roi")
}

#' @keywords internal
#' @noRd
.filter_roi_stats_slow <- function(trdat) {
  nas <- apply(trdat, 2, function(v) any(is.na(v)))
  sdnonzero <- apply(trdat, 2, sd, na.rm = TRUE) > 0
  list(nas = as.logical(nas), sdnonzero = as.logical(sdnonzero))
}

#' @keywords internal
#' @noRd
.filter_roi_stats_fast <- function(trdat) {
  p <- ncol(trdat)
  if (is.null(p) || p == 0L) {
    return(list(nas = logical(0), sdnonzero = logical(0)))
  }

  if (!is.numeric(trdat) && !is.integer(trdat)) {
    return(.filter_roi_stats_slow(trdat))
  }

  nas <- colSums(is.na(trdat)) > 0L
  sdnonzero <- rep(FALSE, p)
  valid_cols <- which(!nas)

  if (length(valid_cols) > 0L && nrow(trdat) >= 2L) {
    # Per-column min/max in one vectorised pass. Columns here have no NA, so
    # the span equals the former row-by-row pmin/pmax recursion exactly
    # (including Inf and overflow cases, which are not finite and so drop).
    x <- trdat[, valid_cols, drop = FALSE]
    rng <- matrixStats::colRanges(x)
    span <- rng[, 2L] - rng[, 1L]
    sdnonzero[valid_cols] <- is.finite(span) & (span > 0)
  }

  list(nas = nas, sdnonzero = sdnonzero)
}

#' Is a preserved centre column unusable?
#'
#' A centre voxel is never preserved when it is all-NA in the training data or
#' has an infinite value in the test data. A partly-NA training column, or NA
#' in the test data, is preserved as before (retrieval NA is handled downstream).
#' `cols` are column positions of the train matrix (one voxel, or the basis
#' columns of one multibasis voxel group).
#' @keywords internal
#' @noRd
.filter_roi_centre_unusable <- function(trdat, tedat, cols) {
  if (length(cols) == 0L) return(FALSE)
  all_na_train <- colSums(is.na(trdat[, cols, drop = FALSE])) == nrow(trdat)
  if (any(all_na_train)) return(TRUE)
  if (!is.null(tedat)) {
    return(any(is.infinite(tedat[, cols, drop = FALSE])))
  }
  FALSE
}

#' @keywords internal
#' @noRd
.filter_roi_stats <- function(trdat) {
  if (.fast_filter_roi_enabled()) {
    .filter_roi_stats_fast(trdat)
  } else {
    .filter_roi_stats_slow(trdat)
  }
}

#' @keywords internal
.get_samples <- function(obj, voxlist) {
  ret <- lapply(voxlist, function(vox) {
    # Expand searchlight windows to full voxel indices when possible (duck-typed)
    vox_expanded <- tryCatch({
      coords <- tryCatch(vox@coords, error = function(...) NULL)
      sp     <- tryCatch(vox@space,  error = function(...) NULL)
      if (!is.null(coords) && !is.null(sp)) {
        neuroim2::grid_to_index(sp, coords)
      } else {
        vox
      }
    }, error = function(...) vox)

    sam <- data_sample(obj, vox_expanded)
  })
  
  n <- length(ret)
  df <- tibble::tibble(sample = ret)
  df[[".id"]] <- gen_id(n)
  df
}

#' @export
get_samples.mvpa_dataset <- function(obj, vox_list) {
  .get_samples(obj, vox_list)
}

#' @export
get_samples.mvpa_surface_dataset <- function(obj, vox_list) {
  .get_samples(obj, vox_list)
}


#' @rdname data_sample
#' @export
data_sample.mvpa_dataset <- function(obj, vox, ...) {
  structure(
    list(
      #data = obj,
      data=NULL,
      vox=vox
    ),
    class = "data_sample"
  )
}

#' @rdname data_sample
#' @export
data_sample.mvpa_multibasis_image_dataset <- function(obj, vox, ...) {
  structure(
    list(
      data = NULL,
      vox = vox
    ),
    class = c("multibasis_data_sample", "data_sample")
  )
}


#' @export
print.data_sample <- function(x, ...) {
  if (is.matrix(x$vox)) {
    cat("data sample with : ", nrow(x$vox), "features")
  } else {
    cat("data sample with : ", length(x$vox), "features")
  }
}

#' @keywords internal
#' @export
filter_roi.default <- function(roi, preserve = NULL, ...) {
  stop("Unsupported ROI type")
}

#' @keywords internal
#' @export
#' @importFrom neuroim2 ROIVec space coords values
filter_roi.ROIVec <- function(roi, preserve = NULL, min_voxels = 2, ...) {
  # Extract the train data values
  trdat <- values(roi$train_roi)
  stats <- .filter_roi_stats(trdat)
  nas <- stats$nas
  sdnonzero <- stats$sdnonzero

  # Test-set columns reach predict() unchanged. A column with an infinite test
  # value is dropped (it is never meaningful). NA in the test data is NOT a
  # drop criterion here: retrieval missing-data semantics are handled downstream
  # (fit_roi's pairwise-complete similarity fallback; see
  # test_era_rsa_searchlight_fast.R). Clean data has no such columns.
  tedat <- if (is.null(roi$test_roi)) NULL else values(roi$test_roi)
  test_bad <- if (is.null(tedat)) {
    rep(FALSE, ncol(trdat))
  } else {
    if (ncol(tedat) != ncol(trdat)) {
      stop("filter_roi: train and test ROI have different numbers of columns")
    }
    colSums(is.infinite(tedat)) > 0L
  }

  # Read basis_count: NULL for standard datasets, integer for multibasis

  basis_count <- roi$basis_count

  # --- Multibasis branch: filter by voxel group ---
  if (!is.null(basis_count) && basis_count > 1L) {
    if (ncol(trdat) %% basis_count != 0L) {
      stop(sprintf("filter_roi: ncol(train_roi) = %d is not divisible by basis_count = %d",
                   ncol(trdat), basis_count))
    }
    V_phys <- ncol(trdat) %/% basis_count
    # Column layout is basis-major: [b1_v1..b1_vV, b2_v1..b2_vV, ...]
    group_idx <- make_group_idx(V_phys, basis_count)

    # Per-column checks (same criteria as standard path, plus test-set finiteness)
    col_bad <- nas | !sdnonzero | test_bad

    # A physical voxel group is bad if ANY of its k columns fails
    group_bad <- as.logical(tapply(col_bad, group_idx, any))
    group_keep <- !group_bad

    n_removed <- sum(!group_keep)
    if (n_removed > 0) {
      # Count reasons at column level, report at group level
      group_na <- as.logical(tapply(nas, group_idx, any))
      group_zv <- as.logical(tapply(!sdnonzero, group_idx, any))
      n_na  <- sum(group_na & !group_keep)
      n_zv  <- sum(group_zv & !group_na & !group_keep)
      .log_debug(
        "filter_roi: removing %d/%d physical voxels (NA: %d, zero-variance: %d); %d remain. [basis_count=%d]",
        n_removed, V_phys, n_na, n_zv, sum(group_keep), basis_count)
    }

    # Preserve specified voxel if requested (e.g., center voxel in searchlight).
    # An unusable centre (all-NA in train, or infinite in test) is not preserved.
    if (!is.null(preserve)) {
      gi <- neuroim2::indices(roi$train_roi)
      # First V_phys indices correspond to basis-1 block
      first_block_ids <- gi[seq_len(V_phys)]
      preserve_group <- match(preserve, first_block_ids)
      if (!is.na(preserve_group) && !group_keep[preserve_group]) {
        group_cols <- which(group_idx == preserve_group)
        if (.filter_roi_centre_unusable(trdat, tedat, group_cols)) {
          .log_debug(
            "Not preserving voxel group %s: all-NA in train or infinite in test data",
            preserve)
        } else {
          .log_debug(
            "Preserving voxel group %s that would have been filtered (NA or zero variance)",
            preserve)
          group_keep[preserve_group] <- TRUE
        }
      }
    }

    # min_voxels check on physical voxel count
    if (sum(group_keep) < min_voxels) {
      group_na <- as.logical(tapply(nas, group_idx, any))
      group_zv <- as.logical(tapply(!sdnonzero, group_idx, any))
      reasons <- c()
      if (any(group_na)) reasons <- c(reasons, sprintf("NA voxel groups: %d", sum(group_na)))
      if (any(group_zv)) reasons <- c(reasons, sprintf("zero-variance voxel groups: %d", sum(group_zv)))
      .log_debug(
        "filter_roi.ROIVec: %d valid physical voxels (< %d). Reasons: %s",
        sum(group_keep), min_voxels,
        paste(reasons, collapse = "; ")
      )
      stop(sprintf("filter_roi: roi has %d valid physical voxels (need at least %d)",
                    sum(group_keep), min_voxels))
    }

    # Expand group-level keep to column-level keep
    keep <- as.logical(group_keep[group_idx])

  } else {
    # --- Standard (non-multibasis) path: unchanged ---

    # Determine columns to keep
    keep <- !nas & sdnonzero & !test_bad

    n_removed <- sum(!keep)
    if (n_removed > 0) {
      n_na  <- sum(nas)
      n_zv  <- sum(!sdnonzero & !nas)  # zero-var but not already counted as NA
      .log_debug(
        "filter_roi: removing %d/%d voxels (NA: %d, zero-variance: %d, infinite test: %d); %d remain.",
        n_removed, length(keep), n_na, n_zv, sum(test_bad & !nas & sdnonzero), sum(keep))
    }

    # Preserve specified voxel if requested (e.g., center voxel in searchlight).
    # An unusable centre (all-NA in train, or infinite in test) is not preserved.
    if (!is.null(preserve)) {
      gi <- neuroim2::indices(roi$train_roi)
      kp <- match(preserve, gi)
      if (!is.na(kp) && !keep[kp]) {
        if (.filter_roi_centre_unusable(trdat, tedat, kp)) {
          .log_debug("Not preserving voxel %s: all-NA in train or infinite in test data", preserve)
        } else {
          .log_debug("Preserving voxel %s that would have been filtered (NA or zero variance)", preserve)
          keep[kp] <- TRUE
        }
      }
    }

    # Need at least 2 valid columns for multivariate analysis
    # (even with preserve, searchlight edges need multiple voxels)
    if (sum(keep) < min_voxels) {
      reasons <- c()
      if (any(nas)) reasons <- c(reasons, sprintf("NA voxels: %d", sum(nas)))
      if (any(!sdnonzero)) reasons <- c(reasons, sprintf("zero-variance voxels: %d", sum(!sdnonzero)))
      .log_debug(
        "filter_roi.ROIVec: %d valid columns (< %d). Reasons: %s",
        sum(keep), min_voxels,
        paste(reasons, collapse = "; ")
      )
      stop(sprintf("filter_roi: roi has %d valid columns (need at least %d)", sum(keep), min_voxels))
    }
  }

  # --- Shared: build filtered ROI and return ---
  # Nothing filtered: hand back the input unchanged rather than rebuilding
  # identical ROIVec objects. Without column_positions the retained columns
  # are all of them (filter_roi.list treats a missing field as identity).
  if (all(keep)) {
    return(roi)
  }

  if (is.null(roi$test_roi)) {
    troi <- ROIVec(space(roi$train_roi),
                   coords(roi$train_roi)[keep,,drop=FALSE],
                   data=trdat[,keep,drop=FALSE])
    result <- list(train_roi=troi, test_roi=NULL)
  } else {
    # Filter train ROI data
    troi <- ROIVec(space(roi$train_roi),
                   coords(roi$train_roi)[keep,,drop=FALSE],
                   data=trdat[,keep,drop=FALSE])

    # Filter test ROI data
    teroi <- ROIVec(space(roi$test_roi),
                    coords(roi$test_roi)[keep,,drop=FALSE],
                    data=tedat[,keep,drop=FALSE])

    result <- list(train_roi=troi, test_roi=teroi)
  }

  # Carry actual retained columns; centroid coordinates need not identify clusters.
  result$column_positions <- which(keep)
  if (!is.null(roi$feature_positions)) result$feature_positions <- roi$feature_positions[keep]
  # Propagate basis_count to downstream consumers
  result$basis_count <- roi$basis_count
  result
}

#' @keywords internal
#' @export
filter_roi.list <- function(roi, preserve = NULL, min_voxels = 2, ...) {
  # Expect a list with train_roi/test_roi as produced by as_roi()
  if (!is.list(roi) || is.null(roi$train_roi)) {
    stop("filter_roi.list: expected list with train_roi/test_roi")
  }
  filtered <- filter_roi(roi$train_roi, preserve = preserve, min_voxels = min_voxels, ...)
  result <- list(train_roi = filtered$train_roi,
       test_roi  = if (!is.null(roi$test_roi)) {
         # Apply the same voxel filter to the test ROI if present
         # Use the logical mask from filtered train ROI indices
         keep_idx <- neuroim2::indices(filtered$train_roi)
         orig_idx <- neuroim2::indices(roi$train_roi)
         kp <- filtered$column_positions %||% match(keep_idx, orig_idx)
         # Guard against mismatch
         if (any(is.na(kp))) {
           warning("filter_roi.list: index mismatch between filtered train and original train ROI; ",
                   "returning unfiltered test ROI.", call. = FALSE)
           roi$test_roi
         } else {
           neuroim2::ROIVec(neuroim2::space(roi$test_roi),
                            neuroim2::coords(roi$test_roi)[kp,,drop=FALSE],
                            data = neuroim2::values(roi$test_roi)[, kp, drop = FALSE])
         }
       } else NULL)
  result$basis_count <- roi$basis_count
  if (!is.null(roi$feature_positions)) {
    # A NULL column_positions means filter_roi returned its input unchanged (identity)
    pos <- filtered$column_positions %||% seq_along(roi$feature_positions)
    result$feature_positions <- roi$feature_positions[pos]
  }
  result
}

#' @keywords internal
#' @export
filter_roi.ROISurfaceVector <- function(roi, preserve = NULL, min_voxels = 2, ...) {
  # Extract the train data values
  trdat <- roi$train_roi@data
  stats <- .filter_roi_stats(trdat)
  nas <- stats$nas
  sdnonzero <- stats$sdnonzero
  
  # Determine columns to keep
  keep <- !nas & sdnonzero
  
  # Preserve specified voxel if requested (e.g., center voxel in searchlight)
  if (!is.null(preserve)) {
    gi <- roi$train_roi@indices
    kp <- match(preserve, gi)
    if (!is.na(kp) && !keep[kp]) {
      .log_debug("Preserving voxel %s that would have been filtered (NA or zero variance)", preserve)
      keep[kp] <- TRUE
    }
  }

  # Need at least 2 valid columns for multivariate analysis
  # (even with preserve, searchlight edges need multiple voxels)
  if (sum(keep) < min_voxels) {
    stop(sprintf("filter_roi: roi has %d valid columns (need at least %d)", sum(keep), min_voxels))
  }
  
  # If there's no test ROI data, return filtered train ROI data only
  if (is.null(roi$test_roi)) {
    troi <- neurosurf::ROISurfaceVector(
      geometry = roi$train_roi@geometry,
      indices = roi$train_roi@indices[keep],
      data = trdat[, keep, drop = FALSE]
    )
    list(train_roi=troi, test_roi=NULL)
  } else {
    # Filter train ROI data
    troi <- neurosurf::ROISurfaceVector(
      geometry = roi$train_roi@geometry,
      indices = roi$train_roi@indices[keep],
      data = trdat[, keep, drop = FALSE]
    )
    
    # Filter test ROI data
    tedat <- roi$test_roi@data
    teroi <- neurosurf::ROISurfaceVector(
      geometry = roi$test_roi@geometry,
      indices = roi$test_roi@indices[keep],
      data = tedat[, keep, drop = FALSE]
    )
    
    list(train_roi=troi, test_roi=teroi)
  }
}


#' @keywords internal
#' @noRd
compute_mask_indices <- function(mask) {
  if (inherits(mask, c("NeuroVol", "NeuroVec"))) {
    vals <- neuroim2::values(mask)
    if (is.matrix(vals)) {
      vals <- vals[, 1, drop = TRUE]
    }
    which(vals != 0)
  } else {
    which(mask != 0)
  }
}

#' Logical membership lookup for mask voxel indices
#'
#' Built once per run (mvpa_iterate stores it as \code{dataset$mask_lookup}) so
#' that \code{as_roi} can test voxel membership by indexing instead of hashing
#' the whole mask for every ROI. Position \code{i} is TRUE when voxel \code{i}
#' is in \code{mask_indices}. The lookup length is \code{max(mask_indices)}:
#' any index beyond it cannot be in the mask.
#' @keywords internal
#' @noRd
.mask_membership_lookup <- function(mask_indices) {
  n <- if (length(mask_indices)) as.integer(max(mask_indices)) else 0L
  lk <- logical(n)
  lk[as.integer(mask_indices)] <- TRUE
  lk
}

#' Per-element membership of voxel indices in a mask lookup
#'
#' Equivalent to \code{match(vox, mask_indices, 0L) > 0L}. Values that are NA,
#' below 1, beyond the lookup, or non-integer are never members.
#' @keywords internal
#' @noRd
.mask_membership <- function(vox, lookup) {
  v <- as.vector(vox)
  k <- logical(length(v))
  n <- length(lookup)
  ok <- which(!is.na(v) & v >= 1 & v <= n & v == floor(v))
  if (length(ok)) {
    k[ok] <- lookup[as.integer(v[ok])]
  }
  k
}

#' Restrict ROI voxel indices to the data mask
#'
#' Replicates the former \code{setdiff}/\code{intersect} logic exactly: when every
#' voxel is inside the mask the input is returned unchanged (duplicates and order
#' kept); otherwise the result is \code{intersect(vox, mask_indices)}, i.e.
#' unique in-mask voxels in first-occurrence order, with the type of the mask
#' index vector. Uses \code{data$mask_lookup} when present, otherwise builds the
#' lookup from \code{data$mask_indices} (or the mask) for this call.
#' @return list with \code{vox} (restricted indices) and \code{n_outside} (number
#'   of distinct voxels outside the mask, as \code{length(setdiff(vox, mask_indices))}).
#' @keywords internal
#' @noRd
.restrict_vox_to_mask <- function(vox, data) {
  ref <- data$mask_indices
  lookup <- data$mask_lookup
  if (!is.logical(lookup)) {
    if (is.null(ref)) ref <- compute_mask_indices(data$mask)
    lookup <- .mask_membership_lookup(ref)
  }

  k <- .mask_membership(vox, lookup)
  if (all(k)) {
    return(list(vox = vox, n_outside = 0L))
  }

  v <- as.vector(vox)
  first <- !duplicated(v)
  n_outside <- sum(!k & first)
  if (is.null(ref)) ref <- integer()
  # intersect(x, y) in base R is c(u[!duplicated(u) & in y], y[numeric()])
  list(vox = c(v[k & first], ref[0L]), n_outside = as.integer(n_outside))
}

#' @keywords internal
#' @noRd
#' @importFrom neuroim2 series_roi
#' @method as_roi data_sample
#' @export
as_roi.data_sample <- function(obj, data, ...) {

  # Defensive validation: ensure all voxel indices fall within the data mask.
  # For regional analysis this is normally a no-op (prep_regional already

  # intersects region_mask with dataset$mask).  For searchlight analysis
  # edge spheres may extend beyond the brain mask, so this catch is needed.
  if (!is.null(data$mask)) {
    restricted <- .restrict_vox_to_mask(obj$vox, data)

    if (restricted$n_outside > 0L) {
      .log_debug(
        "as_roi.data_sample: %d voxel(s) outside data mask (expected for edge searchlights; should not occur for regional).",
        restricted$n_outside)
      obj$vox <- restricted$vox

      # Don't fail for <2 voxels -- let filter_roi handle the minimum check
      if (length(obj$vox) == 0) {
        .log_debug("as_roi.data_sample: 0 voxels remain after mask filtering.")
        train_roi <- structure(list(message = "No voxels remaining after mask filtering"),
                              class = "try-error")
        return(list(train_roi = train_roi, test_roi = NULL))
      }
    }
  }

  train_roi <- try(series_roi(data$train_data, obj$vox), silent = TRUE)

  test_roi <- if (has_test_set(data)) {
    try(series_roi(data$test_data, obj$vox), silent = TRUE)
  }

  if (is.null(test_roi)) {
    list(train_roi = train_roi, test_roi = NULL)
  } else {
    list(train_roi = train_roi, test_roi = test_roi)
  }
}

#' @keywords internal
#' @noRd
#' @importFrom neuroim2 ROIVec coords values space series_roi
combine_multibasis_roi <- function(vec_list, vox) {
  roi_list <- lapply(vec_list, function(v) try(neuroim2::series_roi(v, vox), silent = TRUE))
  err_idx <- which(vapply(roi_list, function(x) inherits(x, "try-error"), logical(1)))
  if (length(err_idx) > 0) {
    return(roi_list[[err_idx[1]]])
  }

  ref <- roi_list[[1]]
  coords_ref <- neuroim2::coords(ref)
  space_ref <- neuroim2::space(ref)

  all_data <- do.call(cbind, lapply(roi_list, neuroim2::values))
  all_coords <- do.call(rbind, replicate(length(roi_list), coords_ref, simplify = FALSE))
  neuroim2::ROIVec(space_ref, all_coords, data = all_data)
}

#' @keywords internal
#' @noRd
#' @method as_roi multibasis_data_sample
#' @export
as_roi.multibasis_data_sample <- function(obj, data, ...) {
  if (!is.null(data$mask)) {
    restricted <- .restrict_vox_to_mask(obj$vox, data)
    if (restricted$n_outside > 0L) {
      .log_debug(
        "as_roi.multibasis_data_sample: %d voxel(s) outside data mask.",
        restricted$n_outside
      )
      obj$vox <- restricted$vox
      if (length(obj$vox) == 0) {
        train_roi <- structure(
          list(message = "No voxels remaining after mask filtering"),
          class = "try-error"
        )
        return(list(train_roi = train_roi, test_roi = NULL))
      }
    }
  }

  train_roi <- combine_multibasis_roi(data$train_data, obj$vox)
  test_roi <- if (has_test_set(data)) combine_multibasis_roi(data$test_data, obj$vox) else NULL

  list(train_roi = train_roi, test_roi = test_roi, basis_count = length(data$train_data))
}

#' @keywords internal
#' @noRd
#' @method as_roi numeric
#' @export
as_roi.numeric <- function(obj, data, ...) {
  ds <- data_sample(data, obj)
  as_roi(ds, data, ...)
}

#' @keywords internal
#' @noRd
#' @method as_roi integer
#' @export
as_roi.integer <- function(obj, data, ...) {
  ds <- data_sample(data, obj)
  as_roi(ds, data, ...)
}

#' @keywords internal
#' @noRd
#' @importFrom neuroim2 space series series_roi
as.data.frame.data_sample <- function(x, row.names = NULL, optional = FALSE, ..., data) {
  train_mat <- neuroim2::series(data$train_data, x$vox)
  
  test_mat <- if (has_test_set(data)) {
    neuroim2::series(data$test_data, x$vox)
  }
  
  cds <- if (is.vector(x$vox)) {
    cds <- neuroim2::index_to_grid(space(data$mask), x$vox)
  } else {
    x$vox
  }
  
  if (!is.null(test_mat)) {
    .type <- rep(c("train", "test"), c(nrow(train_mat), nrow(test_mat)))
    ret <- as.data.frame(rbind(train_mat, test_mat))
    ret$.type <- .type
    attr(ret, "vox") <- cds
    ret
  } else {
    .type <- rep(c("train"), nrow(train_mat))
    ret <- as.data.frame(train_mat)
    ret$.type <- .type
    attr(ret, "vox") <- cds
    ret
  }
  
}
  
 
