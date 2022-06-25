

combine_regional_results = function(results) {
  roinum=NULL
  .rownum=NULL
  if (is.factor(results$result[[1]]$observed)) {
    results %>% dplyr::rowwise() %>% dplyr::do( {

      tib1 <- tibble::tibble(
          .rownum=seq_along(.$result$observed),
          roinum=rep(.$id, length(.$result$observed)),
          observed=.$result$observed,
          pobserved=sapply(seq_along(.$result$observed), function(i) .$result$probs[i, .$result$observed[i]]),
          predicted=.$result$predicted,
          correct=as.character(.$result$observed) == as.character(.$result$predicted)
      )
      #tib2 <- tibble::as_tibble(.$result$probs, .name_repair="universal")
      tib2 <- tibble::as_tibble(.$result$probs, .name_repair=.name_repair)
      names(tib2) <- paste0("prob_", names(tib2))
      cbind(tib1, tib2)
    })
  } else {
   
    results %>% dplyr::rowwise() %>% dplyr::do(
      tibble::tibble(
        .rownum=seq_along(.$result$observed),
        roinum=rep(.$id, length(.$result$observed)),
        observed=.$result$observed,
        predicted=.$result$predicted)
    )
  }
}



#' @import dplyr
#' @importFrom purrr map_df
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
  
  if (is.character(predtabs[[1]]$observed) || is.factor(predtabs[[1]]$observed)) {
    ptab <- map(1:length(predtabs), function(i) predtabs[[i]] %>% mutate(.tableid=i, .weight=wts[i])) %>% 
      map_df(bind_rows) %>% as_tibble(.name_repair=.name_repair)
    
    probs <- if (collapse_regions) {
      ptab %>% dplyr::group_by(.rownum,observed) %>% summarise_at(vars(starts_with("prob_")), 
                                                             funs(stats::weighted.mean(., w=.weight)))
    } else {
      ptab %>% dplyr::group_by(.rownum,observed,roinum) %>% summarise_at(vars(starts_with("prob_")), 
                                                                    funs(stats::weighted.mean(., w=.weight)))
    }
    
    p <- probs %>% ungroup() %>% dplyr::select(dplyr::starts_with("prob_"))
    pmat <- as.matrix(p)
      
    pobserved <- pmat[cbind(seq(1,nrow(probs)),as.integer(probs$observed))]
    mc <- max.col(pmat)
    preds <- levels(probs$observed)[mc]
    
    
    prediction_table <- tibble(
      .rownum=probs$.rownum,
      roinum=if (collapse_regions) 1 else probs$roinum,
      observed=probs$observed,
      pobserved=pobserved,
      predicted=preds,
      correct = predicted == probs$observed,
    ) %>% bind_cols(p)
  } else if (is.numeric(predtabs[[1]]$observed)) {
    stop("combining continuous predictions not implemented")
  }
}

merge_results.regional_mvpa_result <- function(x, ...) {
  rlist <- list(x,...)
  combine_prediction_tables(rlist)
}

regional_mvpa_result <- function(model_spec, performance_table, prediction_table, vol_results, fits=fits) {
  ret <- list(model_spec=model_spec, 
              performance_table=performance_table,
              prediction_table=prediction_table,
              vol_results=vol_results,
              fits=fits)
  
  class(ret) <- c("regional_mvpa_result", "list")
  ret
  
}


prep_regional <- function(model_spec, region_mask) {
  allrois <- sort(unique(region_mask[region_mask>0]))
  region_vec <- as.vector(region_mask)
  region_set <- sort(as.integer(unique(region_vec[region_vec > 0])))
  if (length(region_set) < 1) {
    stop("run_regional: invalid ROI mask, number of ROIs = 0")
  }
  
  vox_iter <- lapply(region_set, function(rnum) which(region_vec == rnum & model_spec$dataset$mask > 0))
  lens <- sapply(vox_iter, length)
  keep <- lens > 1
  
  if (all(!keep)) {
    futile.logger::flog.error("run_regional: no ROIs have more than one voxel.")
    stop()
  }
  
  if (any(lens < 2)) {
    futile.logger::flog.warn(paste("some ROIs have less than two voxels, removing them from list: ", paste(region_set[lens < 2], collapse=",")))
    vox_iter <- vox_iter[keep]
    region_set <- region_set[keep]
  }
  
  list(allrois=allrois, region_vec=region_vec, region_set=region_set,
       vox_iter=vox_iter, lens=lens, keep=keep)
  
}

#' @importFrom neuroim2 map_values
#' @keywords internal
comp_perf <- function(results, region_mask) {
  roinum <- NULL
  ## compile performance results
  perf_mat <- do.call(rbind, results$performance)
  
  ## generate volumetric results
  ## TODO fails when region_mask is an logical vol
  vols <- lapply(1:ncol(perf_mat), function(i) map_values(region_mask, 
                                                          cbind(as.integer(results$id), 
                                                                perf_mat[,i])))
  names(vols) <- colnames(perf_mat)
  
  perfmat <- tibble::as_tibble(perf_mat,.name_repair=.name_repair) %>% dplyr::mutate(roinum = unlist(results$id)) %>% dplyr::select(roinum, dplyr::everything())
  list(vols=vols, perf_mat=perfmat)
}

#' @param return_fits whether to return model fit for every ROI (default is \code{FALSE} to save memory)
#' @param return_predictions whether to return full prediction table with per trial probabilities (can be a large table, set \code{FALSE} to limit memory use)
#' @param compute_performance \code{logical} indicating whether to compute performance measures (e.g. Accuracy, AUC) 
#' @param coalesce_design_vars concatenate additional design variables with output stored in `prediction_table`
#' @import itertools 
#' @import foreach
#' @import doParallel
#' @import parallel
#' @return a \code{list} of type \code{regional_mvpa_result} containing a named list of \code{NeuroVol} objects, where each element contains performance metric and is 
#' labelled according to the metric used (e.g. Accuracy, AUC)
#' @export
#' @rdname run_regional
run_regional.mvpa_model <- function(model_spec, region_mask, return_fits=FALSE, 
                                    return_predictions = TRUE, compute_performance=TRUE, 
                                    coalesce_design_vars=FALSE, ...) {  
  ## to get rid of package check warnings
  roinum=NULL
  ###
  
  prepped <- prep_regional(model_spec, region_mask)
  flog.info("model is: %s", model_spec$model$label)
  
  ## run mvpa for each region
  results <- mvpa_iterate(model_spec, prepped$vox_iter, ids=prepped$region_set,
                          compute_performance=compute_performance, 
                          return_fits = return_fits)

  perf <- if (compute_performance) comp_perf(results, region_mask) else list(vols=list(), perf_mat=tibble())
   
  
  ## compile full prediction table
  prediction_table <- if (return_predictions) {
    combine_regional_results(results) 
  }
  
  if (coalesce_design_vars && return_predictions) {
    prediction_table <- coalesce_join(prediction_table, test_design(model_spec$design), 
                                      by=".rownum")
  }
  
  fits <- if (return_fits) {
    lapply(results$result, "[[", "predictor")
  }
  
  regional_mvpa_result(model_spec=model_spec, performance_table=perf$perf_mat, 
                       prediction_table=prediction_table, vol_results=perf$vols, fits=fits)
}


#' @export
#' @rdname run_regional
#' @param distmethod the distance cmputing method ("pearson" or "spearman")
#' @param regtype the regression method ("pearson", "spearman", "lm", or "rfit")
run_regional.rsa_model <- function(model_spec, region_mask, return_fits=FALSE, 
                                   compute_performance=TRUE,regtype=c("pearson", "spearman", "lm", "rfit"), 
                                   distmethod=c("pearson", "spearman"), coalesce_design_vars=FALSE, ...) {  
  regtype <- match.arg(regtype)
  distmethod <- match.arg(distmethod)
  
  prepped <- prep_regional(model_spec, region_mask)
  
  results <- rsa_iterate(model_spec, prepped$vox_iter, ids=prepped$region_set, regtype=regtype, distmethod=distmethod)
  perf <- if (compute_performance) comp_perf(results, region_mask) else list(vols=list(), perf_mat=tibble())
  
  regional_mvpa_result(model_spec=model_spec, performance_table=perf$perf_mat, 
                       prediction_table=NULL, vol_results=perf$vols, fits=NULL)
  
}
  


