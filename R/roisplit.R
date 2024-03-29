#' roisplit <- function(data, in_id, out_id) {
#'   if (!inherits(data, "data_sample"))
#'     stop("`data` must inherit from `data_sample`.", call. = FALSE)
#'   
#'   if (!is.integer(in_id) | any(in_id < 1))
#'     stop("`in_id` must be a positive integer vector.", call. = FALSE)
#'   
#'   if(!all(is.na(out_id))) {
#'     if (!is.integer(out_id) | any(out_id < 1))
#'       stop("`out_id` must be a positive integer vector.", call. = FALSE)
#'   }
#'   
#'   if (length(in_id) == 0)
#'     stop("At least one row should be selected for the analysis set.",
#'          call. = FALSE)
#'   
#'   structure(
#'     list(
#'       data = data,
#'       in_id = in_id,
#'       out_id = out_id
#'     ),
#'     class = c("roisplit", "rsplit")
#'   )
#' }
#' 
#' # export
#' # print.roisplit <- function(x, ...) {
#' #   out_char <-
#' #     if (is_missing_out_id(x))
#' #       paste(length(complement(x)))
#' #   else
#' #     paste(length(x$out_id))
#' #   
#' #   cat("<Analysis/Assess/Total>\n")
#' #   cat("<",
#' #       length(x$in_id), "/",
#' #       out_char, "/",
#' #       dim(x$data)[4], ">\n",
#' #       sep = "")
#' # }
#' 
#' 
#' 
#' #' Convert an `roisplit` object to a data frame
#' #'
#' #' The analysis or assessment code can be returned as a data
#' #'   frame (as dictated by the `data` argument) using
#' #'   `as.data.frame.roisplit`. `analysis` and
#' #'   `assessment` are shortcuts.
#' #' @param x An `roisplit` object.
#' #' @param row.names `NULL` or a character vector giving the row names for the data frame. Missing values are not allowed.
#' #' @param optional A logical: should the column names of the data be checked for legality?
#' #' @param data Either "analysis" or "assessment" to specify which data are returned.
#' #' @param ...	Additional arguments to be passed to or from methods. Not currently used.
#' #' @export
#' as.data.frame.roisplit <-
#'   function(x,
#'            row.names = NULL,
#'            optional = FALSE,
#'            data = "analysis",
#'            ...) {
#'     
#'     if (!is.null(row.names))
#'       warning( "`row.names` is kept for consistency with the ",
#'                "underlying class but non-NULL values will be ",
#'                "ignored.", call. = FALSE)
#'     if (optional)
#'       warning( "`optional` is kept for consistency with the ",
#'                "underlying class but TRUE values will be ",
#'                "ignored.", call. = FALSE)
#'     if (!is.null(x$col_id)) {
#'       stop("roisplit does not support permuted samples")
#'     }
#'     
#'     series(x$data, as.integer(x, data = data, ... ))
#'     #x$data[as.integer(x, data = data, ...), , drop = FALSE]
#'   }
#' 
#' 
#' as_roi.roisplit <-
#'   function(x,
#'            row.names = NULL,
#'            optional = FALSE,
#'            data = "analysis",
#'            ...) {
#'     
#'     if (!is.null(row.names))
#'       warning( "`row.names` is kept for consistency with the ",
#'                "underlying class but non-NULL values will be ",
#'                "ignored.", call. = FALSE)
#'     if (optional)
#'       warning( "`optional` is kept for consistency with the ",
#'                "underlying class but TRUE values will be ",
#'                "ignored.", call. = FALSE)
#'     if (!is.null(x$col_id)) {
#'       stop("roisplit does not support permuted samples")
#'     }
#'     
#'     series_roi(x$data, as.integer(x, data = data, ... ))
#'     #x$data[as.integer(x, data = data, ...), , drop = FALSE]
#'   }
#' 
#' 
#' #' @export
#' dim.roisplit <- function(x, ...) {
#'   c(
#'     analysis = length(x$in_id),
#'     assessment = length(complement(x)),
#'     n = dim(x$data)[4],
#'     p = prod(dim(x$data)[1:3])
#'   )
#' }

#' @keywords internal
#' @noRd
roi_volume_matrix <- function(mat, refspace, indices, coords) {
  structure(mat,
            refspace=refspace,
            indices=indices,
            coords=coords,
            class=c("roi_volume_matrix", "matrix"))
  
}

#' @keywords internal
#' @noRd
roi_surface_matrix <- function(mat, refspace, indices, coords) {
  structure(mat,
            refspace=refspace,
            indices=indices,
            coords=coords,
            class=c("roi_surface_matrix", "matrix"))

}
