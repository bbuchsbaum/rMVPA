

run_rfit <- function(dvec, obj) {
  form <- paste("dvec", "~", paste(names(obj$design$model_mat), collapse = " + "))
  obj$design$model_mat$dvec <- dvec
  
  res <- if (!is.null(obj$design$include)) {
    rfit(form, data=obj$design$model_mat, subset=obj$design$include)
  } else {
    rfit(form, data=obj$design$model_mat)
  }
  
  coef(res)[-1]
}

run_lm <- function(dvec, obj) {
  form <- paste("dvec", "~", paste(names(obj$design$model_mat), collapse = " + "))
  vnames <- names(obj$design$model_mat)
  obj$design$model_mat$dvec <- dvec
  
  res <- if (!is.null(obj$design$include)) {
    rfit(form, data=obj$design$model_mat, subset=obj$design$include)
  } else {
    rfit(form, data=obj$design$model_mat)
  }
  
  res <- coef(summary(res))[-1,3]
  names(res) <- vnames
  res
}

run_cor <- function(dvec, obj, method) {
  res <- if (!is.null(obj$design$include)) {
    sapply(obj$design$model_mat, function(x) cor(dvec[obj$design$include], x[obj$design$include]))
  } else {
    sapply(obj$design$model_mat, function(x) cor(dvec, x, method=method))
  }
  
  names(res) <- names(obj$design$model_mat)
  res
}

train_model.rsa_model <- function(obj, train_dat, indices, wts=NULL, method=c("lm", "rfit", "pearson", "spearman"), distmethod=c("pearson", "spearman")) {
  method <- match.arg(method)
  distmethod <- match.arg(distmethod)
  
  dtrain <- 1 - cor(t(train_dat), method=distmethod)
  dvec <- dtrain[lower.tri(dtrain)]
  
  switch(method,
         rfit=run_rfit(dvec, obj),
         lm=run_lm(dvec,obj),
         pearson=run_cor(dvec,obj,"pearson"),
         spearman=run_cor(dvec,obj,"spearman"))
  
  
 
}



do_rsa <- function(roi, mod_spec, rnum, method, distmethod) {
  xtrain <- tibble::as_tibble(values(roi$train_roi))
  ind <- indices(roi$train_roi)
  ret <- train_model(mod_spec, xtrain, ind, method, distmethod)
  tibble::tibble(result=list(NULL), indices=list(ind), performance=list(ret), id=rnum, error=FALSE, error_message="~")
}


#' rsa_iterate
#' 
#' Run RSA analysis for each of a list of voxels sets
#' 
#' @param mod_spec a class of type \code{mvpa_model}
#' @param vox_list a \code{list} of voxel indices/coordinates
#' @param ids a \code{vector} of ids for each voxel set
#' @param regtype the analysis method, one of: \code{lm}, \code{rfit}, \code{pearson}, \code{spearman}
#' @param distmethod the method used to computer distances between oservations, one of: \code{pearson}, \code{spearman}
#' @importFrom dplyr do rowwise
#' @export
rsa_iterate <- function(mod_spec, vox_list, ids=1:length(vox_iter), regtype=c("rfit", "lm", "pearson", "spearman"), distmethod=c("spearman", "pearson")) {
  distmethod <- match.arg(distmethod)
  regtype <- match.arg(regtype)
  
  assert_that(length(ids) == length(vox_list))
  sframe <- get_samples(mod_spec$dataset, vox_list)
  
  ret <- sframe %>% dplyr::mutate(rnum=ids) %>% 
    dplyr::rowwise() %>% 
    dplyr::do(do_rsa(as_roi(.$sample), mod_spec, .$rnum, regtype, distmethod))
}
