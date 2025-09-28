

#' @export
nobs.mvpa_design <- function(x) {
  length(x$y_train)
}

#' @export
nresponses.mvpa_design <- function(x) {
  if (is.factor(x$y_train)) {
    length(levels(x$y_train))
  } else if (is.vector(x$y_train)) {
    length(x$y_train)
  } else if (is.matrix(x$y_train)) {
    ncol(x$y_train)
  } else {
    stop()
  }
}

#' Does the design include a test set?
#'
#' @inheritParams has_test_set
#' @rdname has_test_set-methods
#' @export
has_test_set.mvpa_design <- function(obj) {
  !is.null(obj$y_test)
}


#' @rdname y_train-methods
#' @export
y_train.mvpa_design <- function(obj) obj$y_train


#' @rdname y_test-methods
#' @export
y_test.mvpa_design <- function(obj) if (is.null(obj$y_test)) obj$y_train else obj$y_test


#' @rdname test_design-methods
#' @export
test_design.mvpa_design <- function(obj) {
  if (is.null(obj$y_test)) obj$train_design else obj$test_design
}



#' @keywords internal
#' @noRd
parse_variable <- function(var, design) {
  ret <- if (purrr::is_formula(var)) {
    vnames <- all.vars(var[[2]])
    ret <- if (length(vnames) > 1) {
      do.call("interaction", c(lapply(vnames, function(vname) as.factor(design[[vname]])), sep=":"))
    } else {
      design[[vnames]]
    }
    
    assertthat::assert_that(!is.null(ret), msg=paste("formula variable", 
                                                     paste(as.character(var), collapse=" "), " not found."))
    ret
    
  } else if (is.character(var) && length(var) == 1) {
    if (is.null(design[[var]])) {
      stop(paste("`design` does not contain variable named: ", var))
    }
    design[[var]]
  } else if (is.factor(var) || is.integer(var)) {
    if (is.data.frame(design)) {
      assertthat::assert_that(nrow(design) == length(var)) 
    }
    var
  } else {
    stop("'var' must be a formula, factor, or character vector")
  }
  
  if (is.factor(ret)) {
    droplevels(ret)
  } else {
    ret
  }
  
}


#' Create an MVPA Design Object
#'
#' Creates a design object for MVPA analysis that encapsulates training and testing designs,
#' response variables, and optional blocking and splitting factors.
#'
#' @param train_design A data frame containing the training design matrix
#' @param test_design Optional data frame containing the test design matrix (default: NULL)
#' @param y_train Formula or vector specifying the training response variable
#' @param y_test Optional formula or vector specifying the test response variable (default: NULL)
#' @param block_var Optional formula or vector specifying the blocking variable for cross-validation
#' @param split_by Optional formula or vector for splitting analyses
#' @param ... Additional arguments (currently unused)
#'
#' @return An \code{mvpa_design} object (S3 class) containing:
#'   \describe{
#'     \item{train_design}{Data frame of training design}
#'     \item{test_design}{Data frame of test design (if provided)}
#'     \item{y_train}{Training response variable}
#'     \item{y_test}{Test response variable (if provided)}
#'     \item{block_var}{Blocking variable for cross-validation (if provided)}
#'     \item{split_by}{Splitting factor (if provided)}
#'   }
#'
#' @details
#' The \code{y_train} and \code{y_test} can be specified either as formulas (e.g., ~ condition) 
#' or as vectors. If formulas are used, they are evaluated within the respective design matrices.
#' 
#' The \code{block_var} and \code{split_by} can also be specified as formulas or vectors. 
#' If formulas, they are evaluated within the training design matrix.
#'
#' @examples
#' # Basic design with only training data
#' train_df <- data.frame(condition = rep(c("A", "B"), each = 50),
#'                        block = rep(1:5, each = 20),
#'                        group = rep(c("Group1", "Group2"), 50))
#' design <- mvpa_design(train_df, y_train = ~ condition)
#'
#' # Design with test data and blocking variable
#' test_df <- data.frame(condition = rep(c("A", "B"), each = 25))
#' design_with_test <- mvpa_design(
#'   train_df, 
#'   test_df, 
#'   y_train = ~ condition, 
#'   y_test = ~ condition,
#'   block_var = ~ block
#' )
#'
#' # Design with split_by factor
#' design_split <- mvpa_design(
#'   train_df, 
#'   y_train = ~ condition,
#'   split_by = ~ group
#' )
#'
#' @seealso 
#' \code{\link{mvpa_dataset}} for creating the corresponding dataset object
#'
#' @importFrom stats as.formula
#' @export
mvpa_design <- function(train_design, test_design=NULL, y_train, y_test=NULL, block_var=NULL, split_by=NULL, ...) {
 
  y_train <- if (!purrr::is_formula(y_train) && length(y_train) > 1) {
    y_train
  } else {
    parse_variable(y_train, train_design)
  }
  
  if (is.factor(y_train) || is.character(y_train)) {
    y_train <- as.factor(y_train)
   
    if (any(table(y_train) == 0)) {
      futile.logger::flog.warn("y_train: ", table(y_train), capture=TRUE)
      futile.logger::flog.warn("y_train factor has at least one level with zero training instances: dropping unused levels.")
      y_train <- droplevels(y_train)
    }
    
    if (length(table(y_train)) <= 1) {
      #futile.logger::flog.error("y_train: ", table(y_train), capture=TRUE)
      stop(paste("error: y_train factor must have at least 2 levels with one or more training instances"))
    }
    
    ytab <- table(levels(y_train))
    
    if (any(ytab == 0)) {
      futile.logger::flog.info("y_train: ", table(y_train), capture=TRUE)
      stop(paste("error: y_train factor must have at least 1 training instance for every factor level"))
    }
  }
  
  if (!is.null(y_test)) {
    
    ## must have a test_design
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
    des <- if (!is.null(test_design) && nrow(test_design) > 0) test_design else train_design
    split_var <- parse_variable(split_by, des)
    split_groups <- split(1:nrow(des), split_var)
  } else {
    split_groups=NULL
  }
  
  if (!is.null(block_var)) {
    block_var <- parse_variable(block_var, train_design)
    assertthat::assert_that(!is.null(block_var))
  }
 
  test_design <- if (!is.null(test_design)) {
    tibble::as_tibble(test_design, .name_repair = .name_repair) %>% mutate(.rownum=1:n()) 
  }
  
  train_design <- tibble::as_tibble(train_design,.name_repair = .name_repair) %>% mutate(.rownum=1:n())
  
  des <- list(
    train_design=tibble::as_tibble(train_design, .name_repair = .name_repair),
    y_train=y_train,
    test_design=if (!is.null(test_design)) tibble::as_tibble(test_design, .name_repair = .name_repair) else NULL,
    y_test=y_test,
    split_by=split_by,
    split_groups=split_groups,
    block_var=block_var
  )
  
  class(des) <- c("mvpa_design", "list")
  des
}

#' @export
#' @method print rsa_design
print.rsa_design <- function(x, ...) {
  # Ensure crayon is available
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty printing. Please install it.")
  }
  
  # Define color scheme
  header_style <- crayon::bold$cyan
  section_style <- crayon::yellow
  info_style <- crayon::white
  formula_style <- crayon::italic$green
  var_style <- crayon::blue
  
  # Print header
  cat("\n", header_style("RSA Design"), "\n\n")
  
  # Formula section
  cat(section_style("- Formula"), "\n")
  cat(info_style("  - "), formula_style(deparse(x$formula)), "\n")
  
  # Variables section
  cat(section_style("- Variables"), "\n")
  cat(info_style("  - "), var_style(paste(names(x$data), collapse=", ")), "\n")
  
  # Block information
  cat(section_style("- Structure"), "\n")
  if (!is.null(x$block_var)) {
    blocks <- table(x$block_var)
    cat(info_style("  - Blocking: "), "Present\n")
    cat(info_style("  - Number of Blocks: "), crayon::green(length(blocks)), "\n")
    cat(info_style("  - Block Sizes: "), 
        crayon::green(paste0(names(blocks), ": ", blocks, collapse=", ")), "\n")
  } else {
    cat(info_style("  - Blocking: "), crayon::red("None"), "\n")
  }
  cat("\n")
}

#' @export
#' @method print mvpa_design
print.mvpa_design <- function(x, ...) {
  # Ensure crayon is available
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty printing. Please install it.")
  }
  
  # Define color scheme
  header_style <- crayon::bold$cyan
  section_style <- crayon::yellow
  info_style <- crayon::white
  number_style <- crayon::green
  level_style <- crayon::blue
  stat_style <- crayon::italic$white
  
  # Print header
  cat("\n", header_style("MVPA Design"), "\n\n")
  
  # Training section
  cat(section_style("- Training Data"), "\n")
  cat(info_style("  - Observations: "), number_style(format(length(x$y_train), big.mark=",")), "\n")
  
  if (is.factor(x$y_train)) {
    cat(info_style("  - Response Type: "), "Factor\n")
    cat(info_style("  - Levels: "), level_style(paste(levels(x$y_train), collapse=", ")), "\n")
    level_counts <- table(x$y_train)
    cat(info_style("  - Class Distribution: "), 
        number_style(paste0(names(level_counts), ": ", level_counts, collapse=", ")), "\n")
  } else {
    cat(info_style("  - Response Type: "), "Numeric\n")
    cat(info_style("  - Range: "), 
        number_style(sprintf("[%.2f, %.2f]", min(x$y_train), max(x$y_train))), "\n")
  }
  
  # Test data section
  cat(section_style("- Test Data"), "\n")
  if (!is.null(x$y_test)) {
    cat(info_style("  - Observations: "), number_style(format(length(x$y_test), big.mark=",")), "\n")
    if (is.factor(x$y_test)) {
      test_counts <- table(x$y_test)
      cat(info_style("  - Class Distribution: "), 
          number_style(paste0(names(test_counts), ": ", test_counts, collapse=", ")), "\n")
    } else {
    cat(info_style("  - Range: "), 
        number_style(sprintf("[%.2f, %.2f]", min(x$y_test), max(x$y_test))), "\n")
  }
  } else {
    cat(info_style("  - "), crayon::red("None"), "\n")
  }
  
  # Structure section
  cat(section_style("- Structure"), "\n")
  
  # Block variable info
  if (!is.null(x$block_var)) {
    blocks <- table(x$block_var)
    cat(info_style("  - Blocking: "), "Present\n")
    cat(info_style("  - Number of Blocks: "), number_style(length(blocks)), "\n")
    cat(info_style("  - Mean Block Size: "), 
        number_style(format(mean(blocks), digits=2)),
        stat_style(" (SD: "),
        number_style(format(sd(blocks), digits=2)),
        stat_style(")"), "\n")
  } else {
    cat(info_style("  - Blocking: "), crayon::red("None"), "\n")
  }
  
  # Split information
  if (!is.null(x$split_by)) {
    split_info <- table(x$split_by)
    cat(info_style("  - Split Groups: "), 
        number_style(paste0(names(split_info), ": ", split_info, collapse=", ")), "\n")
  } else {
    cat(info_style("  - Split Groups: "), crayon::red("None"), "\n")
  }
  cat("\n")
}
