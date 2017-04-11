
#' @export
has_test_set.mvpa_design <- function(obj) {
  !is.null(obj$y_test) 
}

y_train.mvpa_design <- function(obj) obj$y_train

y_test.mvpa_design <- function(obj) if (is.null(obj$y_test)) obj$y_train else obj$y_test

parse_variable <- function(var, design) {
  if (purrr::is_formula(var)) {
    vnames <- all.vars(var[[2]])
    var <- if (length(vnames) > 1) {
      do.call("interaction", c(lapply(vnames, function(vname) as.factor(design[[vname]])), sep=":"))
    } else {
      as.factor(design[[vnames]])
    }
    var
  } else if (is.character(var)) {
    design[[var]]
  } else {
    stop("'var' must be a formula, factor, or character vector")
  }
  
}

#' mvpa_design
#' 
#' 
#' @export
mvpa_design <- function(train_design, y_train, test_design=NULL, y_test=NULL, block_var=NULL, split_by=NULL) {
  
  y_train <- parse_variable(y_train, train_design)
  
  if (!is.null(y_test)) {
    assert_that(!is.null(test_design))
    y_test <- parse_variable(y_test, test_design)
  }
  
  check_split <- function(split_var) {
    minSplits <- min(table(split_var))
    if (minSplits < 3) {
      stop(paste("error: splitting condition results in fewer than 3 observations in at least one set"))
    }
  }
  
  if (!is.null(split_by)) {
    des <- if (!is.null(test_design)) test_design else train_design
    print(names(des))
    split_var <- parse_variable(split_by, des)
    split_by <- split(1:nrow(des), split_var)
  }
  
  des <- list(
    train_design=train_design,
    y_train=y_train,
    test_design=test_design,
    y_test=y_test,
    split_by=split_by,
    block_var=block_var
  )
  
  class(des) <- c("mvpa_design", "list")
  des
}


