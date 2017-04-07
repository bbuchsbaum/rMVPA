
#' @export
has_test_set.mvpa_design <- function(obj) {
  !is.null(obj$y_test) 
}

y_train.mvpa_design <- function(obj) obj$y_train

y_test.mvpa_design <- function(obj) obj$y_test


#' @export
mvpa_design <- function(train_design, y_train, test_design=NULL, y_test=NULL, block_var=NULL, split_formula=NULL) {
  des <- list(
    train_design=train_design,
    y_train=y_train,
    test_design=test_design,
    y_test=y_test,
    split_formula=split_formula,
    block_var=block_var
  )
  
  class(des) <- c("mvpa_design", "list")
  des
}


