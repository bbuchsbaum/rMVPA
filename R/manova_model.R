#' Create a MANOVA Design
#'
#' This function creates a MANOVA design object containing a formula expression and a named list of data.
#'
#' @param formula A formula expression specifying the MANOVA regression model.
#' @param data A named list containing the dissimilarity matrices and any other auxiliary variables.
#' @return A MANOVA design object with class attributes "manova_design" and "list".
#' @details
#' The function takes a formula expression and a named list of data as input, and returns a MANOVA design object.
#' The object is a list that contains the formula expression and the named list of data with class attributes "manova_design" and "list".
#' This object can be further used for MANOVA analysis or other related multivariate statistical methods.
#' @importFrom assertthat assert_that
#' @importFrom purrr is_formula
#' @importFrom ffmanova ffmanova
#' @examples
#' # Create a MANOVA design
#' formula <- y ~ x1 + x2
#' data_list <- list(y = dissimilarity_matrix_y, x1 = dissimilarity_matrix_x1, x2 = dissimilarity_matrix_x2)
#' manova_design_obj <- manova_design(formula, data_list)
#' @export
manova_design <- function(formula, data) {
  assert_that(purrr::is_formula(formula))
  
  des <- list(
    formula=formula,
    data=data
  )
  class(des) <- c("manova_design", "list")
  des
}


#' Create a MANOVA Model
#'
#' This function creates a MANOVA model object containing an `mvpa_dataset` instance and a `manova_design` instance.
#'
#' @param dataset An \code{mvpa_dataset} instance.
#' @param design A \code{manova_design} instance.
#' @return A MANOVA model object with class attributes "manova_model" and "list".
#' @details
#' The function takes an `mvpa_dataset` instance and a `manova_design` instance as input, and returns a MANOVA model object.
#' The object is a list that contains the dataset and the design with class attributes "manova_model" and "list".
#' This object can be used for further multivariate statistical analysis using the MANOVA method.
#' @importFrom assertthat assert_that
#' @importFrom purrr is_formula
#' @examples
#' # Create a MANOVA model
#' dataset <- create_mvpa_dataset(data_matrix, labels, subject_ids)
#' formula <- y ~ x1 + x2
#' data_list <- list(y = dissimilarity_matrix_y, x1 = dissimilarity_matrix_x1, x2 = dissimilarity_matrix_x2)
#' design <- manova_design(formula, data_list)
#' manova_model_obj <- manova_model(dataset, design)
#' @export
manova_model <- function(dataset,
                      design) {
  
  assert_that(inherits(dataset, "mvpa_dataset"))
  assert_that(inherits(design, "manova_design"))
  
  create_model_spec("manova_model", dataset, design)
  
  
}


#' Train a MANOVA Model
#'
#' This function trains a multivariate analysis of variance (MANOVA) model using the specified design.
#'
#' @param obj An object of class \code{manova_model}.
#' @param train_dat The training data.
#' @param indices The indices of the training data.
#' @param ... Additional arguments passed to the training method.
#' @return A named numeric vector of -log(p-values) for each predictor in the MANOVA model.
#' @importFrom stats as.formula
train_model.manova_model <- function(obj, train_dat, indices, ...) {
  dframe <- obj$design$data
  dframe$response <- as.matrix(train_dat)
  form <- stats::as.formula(paste("response", paste(as.character(obj$design$formula), collapse='')))
  
  fres=ffmanova(form, data=dframe)
  pvals=fres$pValues
  names(pvals) <- sanitize(names(pvals))   
  lpvals <- -log(pvals)
  lpvals
}

# process_roi.manova_model <- function( mod_spec, roi, rnum,...) {
#   xtrain <- tibble::as_tibble(neuroim2::values(roi$train_roi), .name_repair=.name_repair)
#   ind <- indices(roi$train_roi)
#   ret <- try(train_model(mod_spec, xtrain, ind))
#   if (inherits(ret, "try-error")) {
#     flog.warn("error fitting model %s : %s", rnum, attr(ret, "condition")$message)
#     ## error encountered, store error messages
#     emessage <- if (is.null(attr(ret, "condition")$message)) "" else attr(ret, "condition")$message
#     tibble::tibble(result=list(NULL), indices=list(ind), performance=list(NULL), 
#                    id=rnum, error=TRUE, error_message=emessage)
#   } else {
#     tibble::tibble(result=list(NULL), indices=list(ind), performance=list(ret), 
#                    id=rnum, error=FALSE, error_message="~")
#   }
#   
# }


  


#' @importFrom neuroim2 indices values
#' @keywords internal
# do_manova <- function(roi, mod_spec, rnum) {
#   xtrain <- tibble::as_tibble(neuroim2::values(roi$train_roi), .name_repair=.name_repair)
#   ind <- indices(roi$train_roi)
#   ret <- try(train_model(mod_spec, xtrain, ind))
#   if (inherits(ret, "try-error")) {
#     flog.warn("error fitting model %s : %s", rnum, attr(ret, "condition")$message)
#     ## error encountered, store error messages
#     emessage <- if (is.null(attr(ret, "condition")$message)) "" else attr(ret, "condition")$message
#     tibble::tibble(result=list(NULL), indices=list(ind), performance=list(NULL), 
#                    id=rnum, error=TRUE, error_message=emessage)
#   } else {
#     tibble::tibble(result=list(NULL), indices=list(ind), performance=list(ret), 
#                    id=rnum, error=FALSE, error_message="~")
#   }
# }

# MANOVA Iteration for Voxel Sets
#
# This function runs a MANOVA analysis for each of a list of voxel sets.
#
# @param mod_spec A \code{mvpa_model} object representing the model specification.
# @param vox_list A \code{list} of voxel indices or coordinates.
# @param ids A \code{vector} of IDs for each voxel set.
# @param batch_size An \code{integer} specifying the batch size for processing (default is 10% of the total number of IDs).
#
# @return A \code{data.frame} containing the MANOVA results for each voxel set.
#
# @importFrom dplyr do rowwise
# @references Langsrud, O. (2000). Fifty-fifty MANOVA: multivariate analysis of variance for collinear responses. Proceedings of The Industrial Statistics in Action, 2000(2), 250-264.
# @export
# manova_iterate <- function(mod_spec, vox_list, ids=1:length(vox_list), batch_size=as.integer(.1*length(ids))) {
#   assert_that(length(ids) == length(vox_list), msg=paste("length(ids) = ", length(ids), "::", "length(vox_list) =", length(vox_list)))
#   futile.logger::flog.info("manova_iterate: extracting voxel ids")
#   
#   batch_size <- max(1, batch_size)
#   nbatches <- as.integer(length(ids)/batch_size)
#   batch_group <- sort(rep(1:nbatches, length.out=length(ids)))
#   batch_ids <- split(1:length(ids), batch_group)
#   rnums <- split(ids, batch_group)
#   
#   dset <- mod_spec$dataset
#   ##mod_spec$dataset <- NULL
#   
#   tot <- length(ids)
#   
#   result <- purrr::map(1:length(batch_ids), function(i) {
#     futile.logger::flog.info("manova_iterate: compute manovas ...")
#     sf <- get_samples(mod_spec$dataset, vox_list[batch_ids[[i]]]) %>% mutate(.id=batch_ids[[i]], rnum=rnums[[i]])
#     if (nrow(coords(sf$sample[[1]]$vox)) > 1) {
#       sf <- sf %>% rowwise() %>% mutate(roi=list(extract_roi(sample,dset))) %>% select(-sample)
#       fut_manova(mod_spec, sf)
#     }
#   }) %>% bind_rows()
#   
#   result
#   
# }



# MANOVA Iteration for Voxel Sets
#
# This function runs a MANOVA analysis for each of a list of voxel sets.
#
# @param mod_spec A \code{mvpa_model} object representing the model specification.
# @param vox_list A \code{list} of voxel indices or coordinates.
# @param ids A \code{vector} of IDs for each voxel set.
# @param batch_size An \code{integer} specifying the batch size for processing (default is 10% of the total number of IDs).
#
# @return A \code{data.frame} containing the MANOVA results for each voxel set.
#
# @importFrom dplyr do rowwise
# @references Langsrud, O. (2000). Fifty-fifty MANOVA: multivariate analysis of variance for collinear responses. Proceedings of The Industrial Statistics in Action, 2000(2), 250-264.
# @export
# manova_iterate <- function(mod_spec, vox_list, ids=1:length(vox_list),   batch_size=as.integer(.1*length(ids))) {
#   assert_that(length(ids) == length(vox_list), msg=paste("length(ids) = ", length(ids), "::", "length(vox_list) =", length(vox_list)))
#   futile.logger::flog.info("manova_iterate: extracting voxel ids")
#   
#   batch_size <- max(1, batch_size)
#   nbatches <- as.integer(length(ids)/batch_size)
#   batch_group <- sort(rep(1:nbatches, length.out=length(ids)))
#   batch_ids <- split(1:length(ids), batch_group)
#   rnums <- split(ids, batch_group)
#   
#   dset <- mod_spec$dataset
#   ##mod_spec$dataset <- NULL
#   
#   tot <- length(ids)
#   
#   result <- purrr::map(1:length(batch_ids), function(i) {
#     futile.logger::flog.info("manova_iterate: compute manovas ...")
#     sf <- get_samples(mod_spec$dataset, vox_list[batch_ids[[i]]]) %>% mutate(.id=batch_ids[[i]], rnum=rnums[[i]])
#     if (nrow(coords(sf$sample[[1]]$vox)) > 1) {
#       sf <- sf %>% rowwise() %>% mutate(roi=list(extract_roi(sample,dset))) %>% select(-sample)
#       fut_manova(mod_spec, sf)
#     }
#   }) %>% bind_rows()
#   
#   result
#   
# }

#'@keywords internal
# fut_manova <- function(mod_spec, sf) {
#   mod_spec$dataset <- NULL
#   gc()
#   sf %>% furrr::future_pmap(function(.id, rnum, roi) {
#     do_manova(roi, mod_spec, rnum)
#   },.options = furrr::furrr_options(seed = T)) %>% purrr::discard(is.null) %>% dplyr::bind_rows()
# }

#' @export
#' @method print manova_model
print.manova_model <- function(x, ...) {
  # Ensure crayon is available
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty printing. Please install it.")
  }
  
  # Define color scheme
  header_style <- crayon::bold$cyan
  section_style <- crayon::yellow
  info_style <- crayon::white
  number_style <- crayon::green
  formula_style <- crayon::italic$blue
  var_style <- crayon::magenta
  
  # Print header
  cat("\n", header_style("█▀▀ MANOVA Model ▀▀█"), "\n\n")
  
  # Formula section
  cat(section_style("├─ Model Specification"), "\n")
  cat(info_style("│  └─ Formula: "), formula_style(deparse(x$design$formula)), "\n")
  
  # Dataset information
  cat(section_style("├─ Dataset"), "\n")
  dims <- dim(x$dataset$train_data)
  dim_str <- paste0(paste(dims[-length(dims)], collapse=" × "), 
                   " × ", number_style(dims[length(dims)]), " observations")
  cat(info_style("│  ├─ Dimensions: "), dim_str, "\n")
  cat(info_style("│  └─ Type: "), class(x$dataset$train_data)[1], "\n")
  
  # Variables section
  cat(section_style("├─ Variables"), "\n")
  predictors <- all.vars(x$design$formula[[3]])  # Get predictor names from RHS of formula
  response <- all.vars(x$design$formula[[2]])    # Get response name from LHS of formula
  cat(info_style("│  ├─ Response: "), var_style(response), "\n")
  cat(info_style("│  └─ Predictors: "), var_style(paste(predictors, collapse=", ")), "\n")
  
  # Data structure
  cat(section_style("└─ Data Structure"), "\n")
  
  # Check if there's a test set
  has_test <- !is.null(x$dataset$test_data)
  cat(info_style("   ├─ Test Set: "), 
      if(has_test) crayon::green("Present") else crayon::red("None"), "\n")
  
  # Check if there's a mask
  if (!is.null(x$dataset$mask)) {
    mask_sum <- sum(x$dataset$mask > 0)
    cat(info_style("   └─ Active Voxels/Vertices: "), 
        number_style(format(mask_sum, big.mark=",")), "\n")
  } else {
    cat(info_style("   └─ Mask: "), crayon::red("None"), "\n")
  }
  
  cat("\n")
}

#' @export
#' @method print manova_design
print.manova_design <- function(x, ...) {
  # Ensure crayon is available
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty printing. Please install it.")
  }
  
  # Define color scheme
  header_style <- crayon::bold$cyan
  section_style <- crayon::yellow
  info_style <- crayon::white
  formula_style <- crayon::italic$blue
  var_style <- crayon::magenta
  
  # Print header
  cat("\n", header_style("█▀▀ MANOVA Design ▀▀█"), "\n\n")
  
  # Formula section
  cat(section_style("├─ Formula"), "\n")
  cat(info_style("│  └─ "), formula_style(deparse(x$formula)), "\n")
  
  # Data section
  cat(section_style("└─ Variables"), "\n")
  var_names <- names(x$data)
  cat(info_style("   ├─ Total Variables: "), crayon::green(length(var_names)), "\n")
  cat(info_style("   └─ Names: "), var_style(paste(var_names, collapse=", ")), "\n")
  
  cat("\n")
}


