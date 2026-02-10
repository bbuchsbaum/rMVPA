
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

#' @keywords internal
#' @noRd
center_patterns <- function(X, method = c("none", "stimulus_mean")) {
  method <- match.arg(method)
  X <- as.matrix(X)
  if (method == "none") return(X)
  if (nrow(X) < 1L) return(X)
  mu <- colMeans(X)
  sweep(X, 2, mu, "-")
}

#' @keywords internal
#' @noRd
center_patterns_train_test <- function(X_train, X_test = NULL, method = c("none", "stimulus_mean")) {
  method <- match.arg(method)
  X_train <- as.matrix(X_train)
  X_test <- if (!is.null(X_test)) as.matrix(X_test) else NULL
  if (method == "none") {
    return(list(train = X_train, test = X_test))
  }
  if (!is.null(X_test) && ncol(X_train) != ncol(X_test)) {
    stop("center_patterns_train_test: ncol(X_test) must match ncol(X_train).")
  }
  mu <- if (nrow(X_train) > 0) colMeans(X_train) else rep(0, ncol(X_train))
  train_c <- sweep(X_train, 2, mu, "-")
  test_c <- if (!is.null(X_test)) sweep(X_test, 2, mu, "-") else NULL
  list(train = train_c, test = test_c)
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


#' Set rMVPA Logging Level
#'
#' Convenience helper to control the verbosity of rMVPA's logging output.
#' Internally this wraps \code{futile.logger::flog.threshold()} for the
#' package's logger.
#'
#' @param level Logging level to use. Can be a character string
#'   (\code{"TRACE"}, \code{"DEBUG"}, \code{"INFO"}, \code{"WARN"},
#'   \code{"ERROR"}, \code{"FATAL"}, \code{"OFF"}) or a numeric level
#'   constant from the \pkg{futile.logger} package.
#'
#' @details
#' Typical usage:
#' \itemize{
#'   \item \code{set_log_level("INFO")} - default, hides debug messages.
#'   \item \code{set_log_level("DEBUG")} - show detailed per-batch timing and ROI diagnostics.
#'   \item \code{set_log_level("WARN")} - only warnings and errors.
#' }
#'
#' This affects all rMVPA logging performed via \pkg{futile.logger}.
#'
#' @return Invisibly returns the numeric log level.
#'
#' @examples
#' \dontrun{
#'   rMVPA::set_log_level("DEBUG")
#'   rMVPA::set_log_level("WARN")
#' }
#'
#' @export
set_log_level <- function(level = "INFO") {
  if (!requireNamespace("futile.logger", quietly = TRUE)) {
    stop("Package 'futile.logger' is required to control logging. Please install it.")
  }

  lvl_num <- if (is.character(level)) {
    # Normalize and map string levels to flog constants
    lvl <- toupper(level)
    off_const <- get0("OFF", envir = asNamespace("futile.logger"), inherits = FALSE)
    if (is.null(off_const)) off_const <- Inf
    switch(lvl,
           "TRACE" = futile.logger::TRACE,
           "DEBUG" = futile.logger::DEBUG,
           "INFO"  = futile.logger::INFO,
           "WARN"  = futile.logger::WARN,
           "ERROR" = futile.logger::ERROR,
           "FATAL" = futile.logger::FATAL,
           "OFF"   = off_const,
           stop(sprintf("Unknown log level '%s'. Use one of TRACE, DEBUG, INFO, WARN, ERROR, FATAL, OFF.", level))
    )
  } else {
    level
  }

  futile.logger::flog.threshold(lvl_num)
  invisible(lvl_num)
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

# Declare global variables used in non-standard evaluation (dplyr, etc.)
# to satisfy R CMD check's global variable notes.
utils::globalVariables(c(
  "error",
  "item",
  "observation_index",
  "pmap",
  "res_naive",
  "res_remap",
  "rnum",
  "scores_list",
  "vector_rsa_iterate"
))

#' Report System and Package Information for rMVPA
#'
#' Gathers and displays information about the R session, operating system,
#' rMVPA version, and key dependencies. This information is helpful for
#' debugging, reporting issues, and ensuring reproducibility.
#'
#' @return Invisibly returns a list containing the gathered system and package
#'         information. It is primarily called for its side effect: printing
#'         the formatted information to the console.
#' @export
#' @examples
#' \dontrun{
#' # Display system information in the console
#' mvpa_sysinfo()
#'
#' # Capture the information in a variable
#' sys_info <- mvpa_sysinfo()
#' print(sys_info$r_version)
#' print(sys_info$dependencies$rsample)
#' }
mvpa_sysinfo <- function() {
  info <- list()

  # --- R and System Information ---
  r_version <- R.version
  info$r_version <- r_version$string
  info$platform <- r_version$platform
  sys_info <- Sys.info()
  info$os <- paste(sys_info["sysname"], sys_info["release"])
  info$nodename <- sys_info["nodename"]
  info$user <- sys_info["user"]

  # --- rMVPA Version ---
  tryCatch({
    # Use '::' to avoid issues if utils is masked
    info$rmvpa_version <- as.character(utils::packageVersion("rMVPA"))
  }, error = function(e) {
    info$rmvpa_version <- "Error: rMVPA package not found or version unavailable."
  })

  # --- Key Dependency Versions ---
  # Core dependencies and commonly used ones
  deps <- c("neuroim2", "neurosurf", "rsample", "yardstick", "future", "furrr",
            "dplyr", "tibble", "purrr", "stats", "MASS") # Added stats/MASS

  dep_versions <- list()
  for (dep in deps) {
    tryCatch({
      # Check if namespace exists before getting version
      if (requireNamespace(dep, quietly = TRUE)) {
        dep_versions[[dep]] <- as.character(utils::packageVersion(dep))
      } else {
        dep_versions[[dep]] <- "Not installed"
      }
    }, error = function(e) {
      # Fallback if packageVersion fails for an installed package
      dep_versions[[dep]] <- "Error retrieving version"
    })
  }
  info$dependencies <- dep_versions

  # --- Parallel Backend (Future plan) ---
  tryCatch({
     # Get the current plan structure without changing it
     current_plan <- future::plan("list")

     # Format the plan description
     # Based on future:::as.character.FutureStrategy and internal structure
     plan_to_string <- function(p) {
        cl <- class(p)[1]
        details <- ""
        if (inherits(p, "multicore") || inherits(p, "cluster")) {
            workers <- tryCatch(length(p$workers), error=function(e) NA)
            if (!is.na(workers)) {
              details <- paste0(" (workers=", workers, ")")
            }
        } else if (inherits(p, "sequential")) {
            # No extra details needed
        }
        paste0(cl, details)
     }

     if (is.list(current_plan)) {
         # Handle plan stack (e.g., plan(list(cluster, sequential)))
         plan_desc <- paste(sapply(current_plan, plan_to_string), collapse = " -> ")
     } else {
         # Handle single plan function
         plan_desc <- plan_to_string(current_plan)
     }

     info$parallel_backend <- plan_desc
   }, error = function(e) {
     info$parallel_backend <- paste("Error retrieving plan:", e$message)
   })


  # --- Locale ---
  tryCatch({
    # Get all locale categories if possible
    info$locale <- Sys.getlocale("LC_ALL")
  }, error = function(e) {
    # Fallback for specific category if LC_ALL fails
    tryCatch({
       info$locale <- Sys.getlocale("LC_CTYPE")
    }, error = function(e2) {
        info$locale <- "Error retrieving locale"
    })
  })

  # Assign class for custom printing
  class(info) <- "mvpa_sysinfo"

  # Print the formatted info and return the list invisibly
  print(info)
  invisible(info)
}


#' Print mvpa_sysinfo Object
#'
#' Formats and prints the system information gathered by \code{mvpa_sysinfo}.
#' This method provides a user-friendly display of the system configuration.
#'
#' @param x An object of class `mvpa_sysinfo`.
#' @param ... Ignored.
#' @return Invisibly returns the input object \code{x} (called for side effects).
#' @examples
#' \dontrun{
#'   info <- mvpa_sysinfo()
#'   print(info)
#' }
#' @export
#' @keywords internal
print.mvpa_sysinfo <- function(x, ...) {
  # Helper for formatting lines
  cat_line <- function(label, value) {
    # Use a standard check instead of custom %||%
    display_value <- ifelse(is.null(value) || length(value) == 0 || value == "", "N/A", value)
    cat(sprintf("%-25s: %s\n", label, display_value))
  }

  # Use crayon for optional coloring if available
  has_crayon <- requireNamespace("crayon", quietly = TRUE)
  header_style <- if (has_crayon) crayon::bold$cyan else function(s) s
  label_style <- if (has_crayon) crayon::bold else function(s) s

  cat(header_style("------------------------- rMVPA System Information -------------------------\n"))

  # Apply the helper function
  cat_line("R Version", x$r_version)
  cat_line("Platform", x$platform)
  cat_line("Operating System", x$os)
  cat_line("Node Name", x$nodename)
  cat_line("User", x$user)
  cat_line("Locale", x$locale)
  cat_line("rMVPA Version", x$rmvpa_version)
  cat_line("Parallel Backend (Future)", x$parallel_backend)

  cat(label_style("\nKey Dependencies:\n"))
  if (!is.null(x$dependencies) && length(x$dependencies) > 0) {
    max_name_len <- max(nchar(names(x$dependencies))) %||% 10 # Provide default width
    for (i in seq_along(x$dependencies)) {
      dep_version <- x$dependencies[[i]]
      dep_display <- ifelse(is.null(dep_version) || length(dep_version) == 0 || dep_version == "", "N/A", dep_version)
      cat(sprintf("  - %-*s: %s\n", max_name_len, names(x$dependencies)[i], dep_display))
    }
  } else {
    cat("  (No dependency versions found or rMVPA not fully loaded)\n")
  }

  cat(header_style("--------------------------------------------------------------------------\n"))
  invisible(x)
}


# Helper for print method (borrowed from purrr, avoids dependency)
# No longer strictly needed by print.mvpa_sysinfo after the changes, but keep for now 
# in case it's used elsewhere or intended for future use.
#' @noRd 
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

#' Build a group index for multibasis column layout
#'
#' Maps each column to its physical voxel in basis-major layout:
#' \code{[b1_v1..b1_vV, b2_v1..b2_vV, ..., bk_v1..bk_vV]}.
#'
#' @param V_phys Integer. Number of physical voxels.
#' @param basis_count Integer. Number of basis functions (k).
#' @return Integer vector of length \code{V_phys * basis_count}.
#' @noRd
#' @keywords internal
make_group_idx <- function(V_phys, basis_count) {
  rep(seq_len(V_phys), times = basis_count)
}
