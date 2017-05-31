

#' @export
#' 
#' 
nobs.mvpa_design <- function(obj) {
  length(obj$y_train)
}


#' @export
has_test_set.mvpa_design <- function(obj) {
  !is.null(obj$y_test) 
}


#' @export
y_train.mvpa_design <- function(obj) obj$y_train


#' @export
y_test.mvpa_design <- function(obj) if (is.null(obj$y_test)) obj$y_train else obj$y_test

parse_variable <- function(var, design) {
  if (purrr::is_formula(var)) {
    vnames <- all.vars(var[[2]])
    var <- if (length(vnames) > 1) {
      do.call("interaction", c(lapply(vnames, function(vname) as.factor(design[[vname]])), sep=":"))
    } else {
      design[[vnames]]
    }
    var
  } else if (is.character(var) && length(var) == 1) {
    design[[var]]
  } else {
    stop("'var' must be a formula, factor, or character vector")
  }
  
}

#' mvpa_design
#' 
#' @param train_design a \code{data.frame} containing training variables
#' @param y_train a \code{formula}, \code{character} name or \code{factor} designating the training response.
#' @param test_design an optional \code{data.frame} containing test variables
#' @param y_test an optional \code{formula}, \code{character} name or \code{factor} designating the test response.
#' @param block_var an optional \code{formula}, \code{character} name or \code{integer} vector designating the block structure.
#' @param split_by an optional \code{formula} indicating grouping structure for evaluating test performance.
#' 
#' @examples 
#' 
#' @export
mvpa_design <- function(train_design, y_train, test_design=NULL, y_test=NULL, block_var=NULL, split_by=NULL) {
 
  y_train <- if (!purrr::is_formula(y_train) && length(y_train) > 1) {
    y_train
  } else {
    parse_variable(y_train, train_design)
  }
  
  if (is.factor(y_train)) {
   
    if (any(table(y_train) == 0)) {
      flog.warning("y_train: ", table(y_train), capture=TRUE)
      flog.warning("y_train factor has at least one level with zero training instances: dropping unused levels.")
      y_train <- droplevels(y_train)
    }
    
    if (length(table(y_train)) <= 1) {
      flog.error("y_train: ", table(y_train), capture=TRUE)
      stop(paste("error: y_train factor must have at least 2 levels with one or more training instances"))
    }
    
    ytab <- table(levels(y_train))
    
    if (any(ytab == 0)) {
      flog.info("y_train: ", table(y_train), capture=TRUE)
      stop(paste("error: y_train factor must have at least 1 training instance for every factor level"))
    }
  }
  
  if (!is.null(y_test)) {
    assert_that(!is.null(test_design))
    y_test <- if (!purrr::is_formula(y_test) && length(y_test) > 1) {
      y_test
    } else {
      parse_variable(y_test, test_design)
    }
  }
  
  check_split <- function(split_var) {
    minSplits <- min(table(split_var))
    if (minSplits < 3) {
      stop(paste("error: splitting condition results in fewer than 3 observations in at least one set"))
    }
  }
  
  if (!is.null(split_by)) {
    des <- if (!is.null(test_design)) test_design else train_design
    split_var <- parse_variable(split_by, des)
    split_groups <- split(1:nrow(des), split_var)
  } else {
    split_groups=NULL
  }
  
  if (!is.null(block_var)) {
    block_var <- parse_variable(block_var, train_design)
  }
 
  des <- list(
    train_design=tibble::as_tibble(train_design),
    y_train=y_train,
    test_design=tibble::as_tibble(test_design),
    y_test=y_test,
    split_by=split_by,
    split_groups=split_groups,
    block_var=block_var
  )
  
  class(des) <- c("mvpa_design", "list")
  des
}

print.mvpa_design <- function(x) {
  cat("mvpa_design:", "\n")
  cat("  training observations: ", length(x$y_train), "\n")
  if (!is.null(x$y_test)) {
    cat("  test observations: ", length(x$y_test), "\n")
  } else {
    cat("  no test observations. \n")
  }
  cat("  training response: ", capture.output(str(x$y_train)), "\n")
  if (!is.null(x$y_test))
    cat("  test response: ", capture.output(str(x$y_test)), "\n")
  if (!is.null(x$block_var))
    cat("  block var: ", capture.output(str(x$block_var)), "\n")
  if (!is.null(x$split_by))
    cat("  split_by", x$split_by)
  
}



