




#' Compute Group Means of a Matrix
#'
#' This function calculates the average vector for each level of a grouping variable in a given matrix.
#'
#' @param X A matrix for which group means should be calculated.
#' @param margin An integer specifying the margin to average over. Use 1 for averaging over rows, and 2 for averaging over columns.
#' @param group A grouping variable, either a factor or an integer vector, that defines the groups to calculate the means for.
#' @return A matrix with the same number of rows or columns (depending on the margin) as the input matrix X, and the number of columns or rows corresponding to the number of unique groups in the grouping variable.
#' @examples
#' # Create a random matrix
#' data <- matrix(rnorm(100 * 100), 100, 100)
#'
#' # Define a grouping variable
#' groups <- factor(rep(1:5, each = 20))
#'
#' # Calculate group means for each row
#' row_means <- group_means(data, margin = 1, group = groups)
#'
#' # Calculate group means for each column
#' col_means <- group_means(data, margin = 2, group = groups)
#' @export
group_means <- function(X, margin, group) {
  if (margin == 1) {
    xsum <- rowsum(X, group)
    sweep(xsum, 1, table(group), "/") 
  } else if (margin == 2) {
    xsum <- rowsum(t(X), group)
    t(sweep(xsum, 1, table(group), "/"))
  } else {
    stop("'margin' must be 1 or 2")
  }
}

#' @noRd
spearman_cor <- function(x, y=NULL, use="everything") {
  cor(x,y,use, method="spearman")
}

#' @noRd
kendall_cor <- function(x, y=NULL, use="everything") {
  cor(x,y,use, method="kendall")
}

#' @noRd
zeroVarianceColumns <- function(M) {
  which(apply(M, 2, sd, na.rm=TRUE) == 0)
}

#' @keywords internal
#' @noRd
zeroVarianceColumns2 <- function(M) {
  apply(M, 2, sd, na.rm=TRUE) == 0
}

#' @keywords internal
#' @noRd
na_cols <- function(M) {
  apply(M, 2, function(x) any(is.na(x)))
}

#' @keywords internal
#' @noRd
nonzeroVarianceColumns <- function(M) {
  which(apply(M, 2, sd, na.rm=TRUE) > 0)
}

#' @keywords internal
#' @noRd
nonzeroVarianceColumns2 <- function(M) {
  ret <- apply(M, 2, sd, na.rm=TRUE) > 0
  ret[is.na(ret)] <- FALSE
  ret
}

#' @noRd
removeZeroVarianceColumns <- function(M) {
  hasVariance <- which(apply(M, 2, sd, na.rm=TRUE) != 0)
  if (length(hasVariance) > 0) {
    M[, hasVariance, drop=FALSE]
  } else {
    M
  }
}

## dfferent version
## https://alistaire.rbind.io/blog/coalescing-joins/
#' Coalesce Join Two Data Frames
#'
#' This function performs a specified type of join on two data frames and then coalesces the joined columns based on their common column names.
#'
#' @param x A data frame to be joined.
#' @param y A second data frame to be joined.
#' @param by A character vector of variables to join by. If NULL (the default), the function will use the common column names in 'x' and 'y'.
#' @param suffix A character vector of length 2, specifying the suffixes to be used for making unique the common column names in 'x' and 'y'. The default is c(".x", ".y").
#' @param join A join function to be used for joining the data frames. The default is dplyr::full_join.
#' @param ... Additional arguments passed on to the join function.
#' @return A data frame resulting from the specified join operation and coalesced columns.
#' @keywords internal
coalesce_join2 <- function(x, y, 
                          by = NULL, suffix = c(".x", ".y"), 
                          join = dplyr::full_join, ...) {
  joined <- join(x, y, by = by, suffix = suffix, ...)
  # names of desired output
  cols <- union(names(x), names(y))
  
  to_coalesce <- names(joined)[!names(joined) %in% cols]
  suffix_used <- suffix[ifelse(endsWith(to_coalesce, suffix[1]), 1, 2)]
  # remove suffixes and deduplicate
  to_coalesce <- unique(substr(
    to_coalesce, 
    1, 
    nchar(to_coalesce) - nchar(suffix_used)
  ))
  
  coalesced <- purrr::map_dfc(to_coalesce, ~dplyr::coalesce(
    joined[[paste0(.x, suffix[1])]], 
    joined[[paste0(.x, suffix[2])]]
  ))
  names(coalesced) <- to_coalesce
  
  dplyr::bind_cols(joined, coalesced)[cols]
}


#' @keywords internal
coalesce_join <- function(x, y, 
                          by = NULL, suffix = c(".x", ".y"), 
                          join = dplyr::full_join, ...) {
  joined <- join(x, y, by = by, suffix = suffix, ...)
  # names of desired output
  cols <- union(names(x), names(y))
  
  to_coalesce <- names(joined)[!names(joined) %in% cols]
  
  if (length(to_coalesce) == 0) {
    ## nothing to coalesce...
    return(joined)
  }
  
  suffix_used <- suffix[ifelse(endsWith(to_coalesce, suffix[1]), 1, 2)]
  # remove suffixes and deduplicate
  to_coalesce <- unique(substr(
    to_coalesce, 
    1, 
    nchar(to_coalesce) - nchar(suffix_used)
  ))
  
  coalesced <- purrr::map_dfc(to_coalesce, ~dplyr::coalesce(
    joined[[paste0(.x, suffix[1])]], 
    joined[[paste0(.x, suffix[2])]]
  ))
  names(coalesced) <- to_coalesce
  
  dplyr::bind_cols(joined, coalesced)[cols]
}