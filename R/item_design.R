#' Construct an ITEM design for trial-wise decoding
#'
#' Creates a design object for ITEM-style decoding where the trial-wise design
#' matrix (`X_t`) and supervised trial targets (`T_target`) are decoupled from
#' TR-level observations used in searchlight/ROI extraction.
#'
#' @param train_design Data frame/tibble describing TR-level observations.
#'   If `NULL`, a minimal design with one row per TR is created.
#' @param X_t Numeric trial-wise design matrix (`n_time x n_trials`).
#' @param T_target Trial-level supervised targets (matrix/data.frame/vector/factor).
#'   Length/rows must equal `n_trials`.
#' @param run_id Trial-level run/session identifiers (length `n_trials`).
#' @param C_transform Optional trial-to-condition transform matrix.
#' @param V Optional temporal covariance/precision passed to `fmrilss::item_compute_u()`.
#' @param v_type One of `"cov"` or `"precision"`.
#' @param trial_id Optional trial identifiers.
#' @param trial_hash Optional trial hash used by alignment guards.
#' @param split_by Optional split variable passed to `mvpa_design()`.
#' @param Z Optional nuisance matrix (`n_time x q`) used by LS-A (`fmrilss::lsa`).
#' @param Nuisance Alias nuisance matrix for LS-A; if both provided, `Z` is used.
#' @param meta Optional metadata list.
#' @param ... Reserved for forward compatibility.
#'
#' @return An object of class `c("item_design", "mvpa_design", "list")`.
#' @export
item_design <- function(train_design = NULL,
                        X_t,
                        T_target,
                        run_id,
                        C_transform = NULL,
                        V = NULL,
                        v_type = c("cov", "precision"),
                        trial_id = NULL,
                        trial_hash = NULL,
                        split_by = NULL,
                        Z = NULL,
                        Nuisance = NULL,
                        meta = list(),
                        ...) {
  .item_require_fmrilss("item_design")

  v_type <- match.arg(v_type)

  X_t <- as.matrix(X_t)
  if (!is.numeric(X_t)) {
    stop("X_t must be numeric.", call. = FALSE)
  }

  n_time <- nrow(X_t)
  n_trials <- ncol(X_t)

  if (is.null(train_design)) {
    train_design <- data.frame(.time_index = seq_len(n_time))
  }
  train_design <- as.data.frame(train_design)

  if (nrow(train_design) != n_time) {
    stop(
      sprintf(
        "train_design must have %d rows to match nrow(X_t); got %d.",
        n_time,
        nrow(train_design)
      ),
      call. = FALSE
    )
  }

  Z <- .item_as_nuisance_optional(Z, n_time, "Z")
  Nuisance <- .item_as_nuisance_optional(Nuisance, n_time, "Nuisance")

  bundle <- fmrilss::item_build_design(
    X_t = X_t,
    T_target = T_target,
    run_id = run_id,
    C_transform = C_transform,
    trial_id = trial_id,
    trial_hash = trial_hash,
    meta = meta,
    diagnostics = list(),
    validate = TRUE
  )

  # Keep cv_labels aligned with TR rows; ITEM folds are handled internally by run_id.
  base <- mvpa_design(
    train_design = train_design,
    cv_labels = seq_len(n_time),
    targets = seq_len(n_time),
    split_by = split_by
  )

  base$item_bundle <- bundle
  base$X_t <- bundle$X_t
  base$T_target <- bundle$T_target
  base$run_id <- bundle$run_id
  base$C_transform <- bundle$C_transform
  base$trial_id <- bundle$trial_id
  base$trial_hash <- bundle$trial_hash
  base$V <- V
  base$v_type <- v_type
  base$Z <- Z
  base$Nuisance <- Nuisance
  base$meta <- bundle$meta

  class(base) <- c("item_design", class(base))
  base
}

#' @export
#' @rdname y_train-methods
y_train.item_design <- function(obj) {
  obj$cv_labels
}

#' @export
#' @method print item_design
print.item_design <- function(x, ...) {
  n_time <- if (!is.null(x$X_t)) nrow(x$X_t) else NA_integer_
  n_trials <- if (!is.null(x$X_t)) ncol(x$X_t) else NA_integer_
  n_runs <- if (!is.null(x$run_id)) length(unique(x$run_id)) else NA_integer_

  cat("item_design\n")
  cat("===========\n")
  cat(sprintf("TR rows: %s\n", n_time))
  cat(sprintf("Trials: %s\n", n_trials))
  cat(sprintf("Runs: %s\n", n_runs))
  cat(sprintf("Targets columns: %s\n", if (!is.null(x$T_target)) ncol(x$T_target) else NA_integer_))
  invisible(x)
}

#' @keywords internal
#' @noRd
.item_require_fmrilss <- function(caller = "ITEM") {
  if (!requireNamespace("fmrilss", quietly = TRUE)) {
    stop(
      sprintf(
        "%s requires package 'fmrilss'. Install it to use ITEM helpers.",
        caller
      ),
      call. = FALSE
    )
  }
}

#' @keywords internal
#' @noRd
.item_as_nuisance_optional <- function(x, n_time, name) {
  if (is.null(x)) return(NULL)

  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop(sprintf("%s must be numeric when provided.", name), call. = FALSE)
  }
  if (nrow(x) != n_time) {
    stop(
      sprintf("%s must have %d rows to match X_t; got %d.", name, n_time, nrow(x)),
      call. = FALSE
    )
  }

  x
}
