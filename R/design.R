
#' @export
has_test_set.mvpa_design <- function(obj) {
  !is.null(obj$testY) 
}



#' @export
mvpa_design <- function(Y, train_design, test_design=NULL, testY=NULL, block_var=NULL, split_formula=NULL) {
  des <- list(
    Y=Y,
    train_design=train_design,
    test_design=test_design,
    testY=testY,
    split_formula=split_formula,
    block_var=block_var
  )
  
  class(des) <- c("mvpa_design", "list")
  des
}

# 
# #' @import formula.tools
# mvpa_external_design <-
#   function(Y,
#            train_design,
#            Ytest,
#            test_design,
#            block_var,
#            split_formula = NULL) {
#     des <- list(
#       Y = Y,
#       train_design = train_design,
#       Ytest = Ytest,
#       test_design = train_design,
#       split_formula = split_formula,
#       block_var = block_var
#     )
#     
#     class(des, "mvpa_external_design", "mvpa_design")
#     des
#   }

