

sanitize <- function(name) {
  name <- gsub(":", ".", name)
  name <- gsub(" ", "", name)
  name <- gsub("[\\(\\)]", ".", name, perl=TRUE)
  name <- gsub(",", "_", name)
  name <- gsub("\\.$", "", name)
  name
}

#' rsa_design
#' 
#' @param formula a formula expression specifying the dissimilarity-based regression function
#' @param data a named list containing the dissimilarity matrices and any other auxiliary variables
#' @param block_var an optional \code{formula}, \code{character} name or \code{integer} vector designating the block structure.
#' @param split_by an optional \code{formula} indicating grouping structure for evaluating test performance.
#' @param keep_intra_run a \code{logical} indicating whether to exclude within-run comparisons
#' @importFrom assertthat assert_that
rsa_design <- function(formula, data, block_var=NULL, split_by=NULL, keep_intra_run=FALSE) {
  assert_that(purrr::is_formula(formula))
  
  ## check that all variables are either matrices, "dist", or vectors
  nr <- sapply(data, function(x) {
    if (is.matrix(x)) {
      nrow(x)
    } else if (inherits(x, "dist")) {
      attr(x, "Size")
    } else if (is.vector(x)) {
      length(x)
    } else {
      stop(paste("illegal variable type", class(x)))
    }
  })
  
  assert_that(all(nr == nr[1]), msg="all elements in 'data' must have the same number of rows")
  
  check_split <- function(split_var) {
    minSplits <- min(table(split_var))
    if (minSplits < 3) {
      stop(paste("error: splitting condition results in fewer than 3 observations in at least one set"))
    }
  }
  
  split_groups <- if (!is.null(split_by)) {
    split_var <- parse_variable(split_by, data)
    split(seq_along(split_var), split_var)
  } 
  
  #time_var <- if (!is.null(time_var)) {
  #  parse_variable(time_var, data)
  #} 
  
  block_var <- if (!is.null(block_var)) {
    parse_variable(block_var, data)
  }
  
  include <- if (!is.null(block_var) && !keep_intra_run) {
    as.vector(dist(block_var)) != 0
  }
  
  
  des <- list(
    formula=formula,
    data=data,
    split_by=split_by,
    split_groups=split_groups,
    #time_var=time_var,
    block_var=block_var,
    include=include
  )
  
  mmat <- rsa_model_mat(des)
  des$model_mat <- mmat
  class(des) <- c("rsa_design", "list")
  des
}

rsa_model_mat <- function(rsa_des) {

  rvars <- labels(terms(rsa_des$formula))
  denv <- list2env(rsa_des$data)
  vset <- lapply(rvars, function(x) eval(parse(text=x), denv))
  
  vmatlist <- lapply(vset, function(v) {
    if (inherits(v, "dist")) {
      ## an distance matrix of class "dist"
      as.vector(v)
    } else if (isSymmetric(v)) {
      ## a full distance matrix
      v[lower.tri(v)]
    } else {
      as.vector(dist(v))
    }
  })
  
  if (!is.null(rsa_des$include)) {
    vmatlist <- lapply(vmatlist, function(v) v[rsa_des$include])
  }
  
  names(vmatlist) <- sanitize(rvars)
  vmatlist
}


#' rsa_model
#' 
#' @param dataset a \code{mvpa_dataset} instance
#' @param design a \code{rsa_design} instance
rsa_model <- function(dataset,
                      design) {
  
  assert_that(inherits(dataset, "mvpa_dataset"))
  structure(list(dataset=dataset,
                 design=design),
            class=c("rsa_model", "list"))
  
  
}




