#' Construct an hrfdecoder design (explicit mvpa_design extension)
#'
#' Creates a specialized design object for continuous-time decoding that carries
#' event metadata in addition to the fields required by `mvpa_design`.
#'
#' This constructor avoids adding ad-hoc members to a plain `mvpa_design` by
#' returning an object with class `c("hrfdecoder_design", "mvpa_design", "list")`.
#'
#' The hrfdecoder approach operates on continuous TR-level data rather than
#' trial-level beta estimates. `cv_labels` is set to TR indices (1:T) and is used
#' only for fold construction in the cross-validation machinery. Critically, fold
#' assignment is determined by `block_var` (run IDs), NOT by `cv_labels` values.
#' The actual decoding targets come from the `event_model` and `events` fields
#' (passed as `targets`), and `train_model.hrfdecoder_model()` ignores the
#' `cv_labels` parameter entirely.
#'
#' @param event_model An event model object (e.g., from fmridesign::event_model).
#'   This must include a sampling_frame attribute defining the TR structure.
#' @param events A data.frame/tibble of trial/event rows with event-level labels.
#'   Must contain columns: `onset` (event onset times in seconds), `condition`
#'   (factor or character labels).
#' @param block_var Integer/factor vector of length T (TRs) specifying run/block ids.
#' @param split_by Optional formula or vector to define split groups (passed to mvpa_design).
#'
#' @return An object inheriting from mvpa_design with explicit class `hrfdecoder_design`
#'         and additional fields `$event_model` and `$events`.
#'
#' @examples
#' \dontrun{
#' library(fmridesign)
#' library(fmrihrf)
#'
#' # Create event table
#' events_df <- data.frame(
#'   onset = c(5, 15, 35, 45),
#'   condition = factor(c("A", "B", "A", "B")),
#'   run = c(1, 1, 2, 2)
#' )
#'
#' # Define sampling frame (TR structure)
#' sframe <- sampling_frame(blocklens = c(60, 60), TR = 2)
#'
#' # Build event model
#' evmod <- event_model(
#'   onset ~ hrf(condition, basis = "spmg1"),
#'   data = events_df,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#'
#' # Create hrfdecoder design
#' block_var <- rep(1:2, each = 60)  # 60 TRs per run
#' design <- hrfdecoder_design(
#'   event_model = evmod,
#'   events = events_df,
#'   block_var = block_var
#' )
#' }
#' @export
hrfdecoder_design <- function(event_model, events, block_var, split_by = NULL) {
  assertthat::assert_that(!is.null(event_model))
  assertthat::assert_that(is.data.frame(events) || tibble::is_tibble(events))
  assertthat::assert_that(length(block_var) > 0)

  # Validate events table structure
  required_cols <- c("onset", "condition")
  missing_cols <- setdiff(required_cols, names(events))
  if (length(missing_cols) > 0) {
    stop("events table must contain columns: ", paste(missing_cols, collapse = ", "))
  }

  # Validate event_model has sampling_frame and check TR count
  sframe <- attr(event_model, "sampling_frame")
  if (is.null(sframe)) {
    warning("event_model has no sampling_frame attribute; cannot validate TR count")
  } else {
    # Check that block_var length matches expected TR count
    expected_TRs <- sum(attr(sframe, "blocklens"))
    if (length(block_var) != expected_TRs) {
      warning(sprintf(
        "block_var length (%d) does not match sampling_frame total TRs (%d)",
        length(block_var), expected_TRs
      ))
    }

    # Check if any events fall outside TR range
    TR <- attr(sframe, "TR")
    total_time <- length(block_var) * TR
    events_outside <- events$onset > total_time
    if (any(events_outside)) {
      warning(sprintf(
        "%d event(s) have onsets beyond total acquisition time (%.1f seconds)",
        sum(events_outside), total_time
      ))
    }
  }

  # Build a minimal training design for fold construction.
  #
  # cv_labels is set to TR indices (1:T) for the sole purpose of length
  # validation and subsetting in the cross-validation machinery.
  # Fold assignment is determined SOLELY by block_var (run IDs), NOT by cv_labels.
  # In blocked_cross_validation, entire runs are held out based on block_var only.
  #
  # The ACTUAL decoding targets come from $event_model and $events (event-level labels),
  # passed as the `targets` argument so they are stored on the design object.
  # train_model.hrfdecoder_model() ignores cv_labels entirely and reads from targets.
  Tlen <- length(block_var)
  train_df <- data.frame(block = block_var)

  mvdes <- mvpa_design(
    train_design = train_df,
    cv_labels = seq_len(Tlen),                              # TR indices for fold construction
    targets = list(event_model = event_model, events = events),  # Real decoding targets
    block_var = ~ block,                                    # Run IDs - determines fold assignment
    split_by = split_by
  )

  # Also attach as top-level fields for convenient direct access by the decoder
  # adapter and print method (mirrors what is stored in $targets above).
  mvdes$event_model <- event_model  # HRF-convolved event design from fmridesign
  mvdes$events <- events            # Event-level labels and timing

  class(mvdes) <- c("hrfdecoder_design", class(mvdes))
  mvdes
}


#' Print method for hrfdecoder_design
#'
#' @param x An hrfdecoder_design object
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the input object \code{x} (called for side effects).
#' @examples
#' \dontrun{
#'   # Requires hrfdecoder and fmridesign packages
#'   # print(hrfdecoder_design_object)
#' }
#' @export
print.hrfdecoder_design <- function(x, ...) {
  cat("hrfdecoder_design\n")
  cat("=================\n\n")

  # TR information
  Tlen <- nrow(x$train_design)
  cat(sprintf("TRs: %d\n", Tlen))

  # Block/run information
  n_blocks <- length(unique(x$train_design$block))
  cat(sprintf("Runs/Blocks: %d\n", n_blocks))

  # Event information
  if (!is.null(x$events)) {
    n_events <- nrow(x$events)
    conditions <- unique(x$events$condition)
    cat(sprintf("Events: %d (conditions: %s)\n",
                n_events, paste(conditions, collapse = ", ")))
  }

  # Sampling frame info (if available)
  sframe <- attr(x$event_model, "sampling_frame")
  if (!is.null(sframe)) {
    TR <- attr(sframe, "TR")
    blocklens <- attr(sframe, "blocklens")
    cat(sprintf("TR: %.2f seconds\n", TR))
    cat(sprintf("Block lengths: %s TRs\n", paste(blocklens, collapse = ", ")))
  }

  # Split groups
  if (!is.null(x$split_groups)) {
    cat(sprintf("Split groups: %d\n", length(x$split_groups)))
  }

  invisible(x)
}

