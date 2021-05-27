




#' group_means
#' 
#' Given a matrix, compute the average vector for each level of a grouping variable.
#' 
#' @param X the matrix
#' @param margin the margin to average over (1 for rows, 2 for columns)
#' @param group the grouping variable, a \code{factor} or \code{integer} vector
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

spearman_cor <- function(x, y=NULL, use="everything") {
  cor(x,y,use, method="spearman")
}

kendall_cor <- function(x, y=NULL, use="everything") {
  cor(x,y,use, method="kendall")
}

zeroVarianceColumns <- function(M) {
  which(apply(M, 2, sd, na.rm=TRUE) == 0)
}

#' @keywords internal
zeroVarianceColumns2 <- function(M) {
  apply(M, 2, sd, na.rm=TRUE) == 0
}

#' @keywords internal
na_cols <- function(M) {
  apply(M, 2, function(x) any(is.na(x)))
}

#' @keywords internal
nonzeroVarianceColumns <- function(M) {
  which(apply(M, 2, sd, na.rm=TRUE) > 0)
}

#' @keywords internal
nonzeroVarianceColumns2 <- function(M) {
  ret <- apply(M, 2, sd, na.rm=TRUE) > 0
  ret[is.na(ret)] <- FALSE
  ret
}


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