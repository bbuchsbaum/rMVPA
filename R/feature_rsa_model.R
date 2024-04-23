
feature_rsa_design <- function(S, labels, k=0) {
  assertthat::assert_that(is.matrix(S))
  assertthat::assert_that(nrow(S) == length(labels))
  assertthat::assert_that(isSymmetric(S))
  
  S <- (S + t(S))/2
  
  if (k == 0) {
    eres <- eigen(S)
    k <- max(which(eres$values > 1))
    k <- max(k, 2)
    F <- eres$vectors[,1:k, drop=FALSE]
  } else {
    assertthat::assert_that(k > 0 && k <= nrow(S))
    eres <- eigen(S)
    F <- eres$vectors[, 1:k, drop=FALSE]
  }
  
  ret <- list(
    S=S,
    F=F,
    labels=labels
  )
  
  class(ret) <- "feature_rsa_design"
  ret
}



feature_rsa_model <- function(dataset,
                       design,
                       method=c("scca", "pls", "pca"),
                       crossval=NULL) {
  
  assert_that(inherits(dataset, "mvpa_dataset"))
  assert_that(inherits(design, "feature_rsa_design"))
  
  if (is.null(crossval) && !is.null(design$block_var)) {
    crossval <- blocked_cross_validation(design$block_var)
  }
  
  assertthat::assert_that(!is.null(crossval))
  
  ret <- list(
              dataset=dataset,
              design=design,
              crossval=crossval)
  
  class(ret) <- "feature_rsa_model"
  ret
}


#' Train an RSA Model
#'
#' This function trains an RSA (representational similarity analysis) model using the specified method and distance calculation.
#'
#' @param obj An object of class \code{featurersa_model}.
#' @param train_dat The training data.
#' @param indices The indices of the training data.
#' @param ... Additional arguments passed to the training method.
#' @return The trained model.
train_model.feature_rsa_model <- function(obj, train_dat, indices, ...) {
  X <- as.matrix(train_dat)
  whitening::scca(X, obj$design$F)
  
}



#' @importFrom neuroim2 indices values
do_feature_rsa <- function(roi, mod_spec, rnum) {
  xtrain <- tibble::as_tibble(neuroim2::values(roi$train_roi), .name_repair=.name_repair)
  ind <- indices(roi$train_roi)
  ret <- try(train_model(mod_spec, xtrain, ind))
  if (inherits(ret, "try-error")) {
    tibble::tibble(result=list(NULL), indices=list(ind), performance=list(ret), id=rnum, error=TRUE, error_message=attr(ret, "condition")$message)
  } else {
    tibble::tibble(result=list(NULL), indices=list(ind), performance=list(ret), id=rnum, error=FALSE, error_message="~")
  }
}

  

#' feature_rsa_iterate
#'
#' Runs feature-based representational similarity analysis (RSA) for each voxel set in a list.
#'
#' @param mod_spec An object of class \code{rsa_model} specifying the RSA model.
#' @param vox_list A \code{list} of voxel indices or coordinates for each voxel set.
#' @param ids A \code{vector} of IDs for each voxel set (defaults to 1:length(vox_list)).
#' @param permute Logical, whether to permute the labels (defaults to FALSE).
#' @importFrom dplyr do rowwise
#' @inheritParams mvpa_iterate
feature_rsa_iterate <- function(mod_spec, vox_list, ids=1:length(vox_list),  permute=FALSE) {
  assert_that(length(ids) == length(vox_list), msg=paste("length(ids) = ", length(ids), "::", "length(vox_list) =", length(vox_list)))
  sframe <- get_samples(mod_spec$dataset, vox_list)
  
  ## iterate over searchlights using parallel futures
  sf <- sframe %>% dplyr::mutate(rnum=ids) 
  fut_feature_rsa(mod_spec,sf, regtype, distmethod)
}


#' @keywords internal
fut_feature_rsa <- function(mod_spec, sf) {
  #mod_spec$dataset <- NULL
  gc()
  
  #browser()
  sf %>% furrr::future_pmap(function(sample, rnum, .id) {
    ## extract_roi?
    do_feature_rsa(as_roi(sample, mod_spec$dataset), mod_spec, rnum)
  }, .options = furrr::furrr_options(seed = T)) %>% dplyr::bind_rows()
  
  
}

  