
#' @export
has_test_set.mvpa_design <- function(obj) {
  !is.null(obj$y_test) 
}

y_train.mvpa_design <- function(obj) obj$y_train

y_test.mvpa_design <- function(obj) if (is.null(obj$y_test)) obj$y_train else obj$ytest


#' @export
mvpa_design <- function(train_design, y_train, test_design=NULL, y_test=NULL, block_var=NULL, split_by=NULL) {
  
  if (!is.null(split_by)) {
    assert_that(purrr::is_formula(split_by))
    split_vars <- all.vars(split_by[[2]])
    des <- if (!is.null(test_design)) test_design else train_design
    split_var <- do.call("interaction", lapply(vars, function(vname) as.factor(des[[vname]])))
    minSplits <- min(table(split_var))
    if (minSplits < 3) {
      stop(paste("error: splitting condition results in fewer than 3 observations in at least one set"))
    }
    
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


