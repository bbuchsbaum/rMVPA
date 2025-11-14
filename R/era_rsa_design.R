#' Build item-level ERA-RSA confounds and summaries from mvpa_design
#'
#' Constructs per-item summaries (block/run/time/lag) and item-level confound RDMs
#' that plug directly into \code{era_rsa_model(..., confound_rdms=...)}. Reuses
#' temporal helpers in R/temporal_rdms.R to avoid redundancy.
#'
#' @param design mvpa_design (same object passed to era_rsa_model).
#' @param key_var Column name or formula giving the item ID (e.g., ~ ImageID).
#' @param phase_var Column name or formula indicating phase (e.g., ~ Phase).
#' @param encoding_level Encoding phase label; default = first level of phase_var.
#' @param retrieval_level Retrieval phase label; default = second level of phase_var.
#' @param block_var Optional column giving block/run membership.
#' @param time_var Optional column giving trial index or onset time.
#'
#' @return A list with elements:
#'   \itemize{
#'     \item items: factor of item IDs used (intersection of E/R)
#'     \item item_block: optional factor of per-item block labels
#'     \item item_time_enc, item_time_ret: optional numeric vectors
#'     \item item_lag: optional numeric vector (ret - enc)
#'     \item item_run_enc, item_run_ret: optional per-item run factors
#'     \item confound_rdms: named list of item-level RDMs (block/time/run)
#'   }
#' @export
era_rsa_design <- function(design,
                           key_var,
                           phase_var,
                           encoding_level  = NULL,
                           retrieval_level = NULL,
                           block_var = NULL,
                           time_var  = NULL) {

  stopifnot(inherits(design, "mvpa_design"))
  d <- design$train_design

  key   <- parse_variable(key_var,   d)
  phase <- factor(parse_variable(phase_var, d))
  levs  <- levels(phase)
  if (is.null(encoding_level))  encoding_level  <- levs[1L]
  if (is.null(retrieval_level)) retrieval_level <- levs[2L]

  enc_mask <- phase == encoding_level
  ret_mask <- phase == retrieval_level

  block <- if (!is.null(block_var)) parse_variable(block_var, d) else NULL
  timev <- if (!is.null(time_var))  parse_variable(time_var,  d) else NULL

  # items with at least one enc and one ret trial
  key_enc <- unique(key[enc_mask])
  key_ret <- unique(key[ret_mask])
  common_items <- sort(intersect(key_enc, key_ret))

  Mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

  item_block    <- NULL
  item_time_enc <- NULL
  item_time_ret <- NULL
  item_lag      <- NULL
  item_run_enc  <- NULL
  item_run_ret  <- NULL
  confound_rdms <- list()

  if (!is.null(block)) {
    item_block <- sapply(common_items, function(k) {
      b_enc <- block[enc_mask & key == k]
      if (length(b_enc)) Mode(b_enc) else NA
    })
    names(item_block) <- common_items
    B <- outer(item_block, item_block, FUN = function(a, b) as.numeric(a != b))
    rownames(B) <- colnames(B) <- common_items
    confound_rdms$block <- B
  }

  if (!is.null(timev)) {
    item_time_enc <- sapply(common_items, function(k) {
      t_enc <- timev[enc_mask & key == k]
      if (length(t_enc)) mean(t_enc) else NA_real_
    })
    item_time_ret <- sapply(common_items, function(k) {
      t_ret <- timev[ret_mask & key == k]
      if (length(t_ret)) mean(t_ret) else NA_real_
    })
    names(item_time_enc) <- names(item_time_ret) <- common_items
    item_lag <- item_time_ret - item_time_enc

    # Use temporal_rdm for encoding-time distance RDM (distance-like)
    if (all(is.finite(item_time_enc))) {
      Ten <- temporal_rdm(index = item_time_enc,
                          block = item_block, # may be NULL
                          kernel = "exp",
                          metric = "distance",
                          within_blocks_only = TRUE,
                          normalize = "rank",
                          as_dist = TRUE)
      confound_rdms$time_enc <- Ten
    }
  }

  # Optional: per-item run labels by phase when block_var encodes runs
  if (!is.null(block)) {
    item_run_enc <- sapply(common_items, function(k) {
      b <- block[enc_mask & key == k]
      if (length(b)) Mode(b) else NA
    })
    item_run_ret <- sapply(common_items, function(k) {
      b <- block[ret_mask & key == k]
      if (length(b)) Mode(b) else NA
    })
    names(item_run_enc) <- names(item_run_ret) <- common_items

    # Run confound RDMs for encoding and retrieval separately
    Renc <- outer(item_run_enc, item_run_enc, FUN = function(a, b) as.numeric(a == b))
    Rret <- outer(item_run_ret, item_run_ret, FUN = function(a, b) as.numeric(a == b))
    rownames(Renc) <- colnames(Renc) <- common_items
    rownames(Rret) <- colnames(Rret) <- common_items
    confound_rdms$run_enc <- Renc
    confound_rdms$run_ret <- Rret
  }

  list(
    items          = factor(common_items, levels = common_items),
    item_block     = if (!is.null(item_block)) factor(item_block) else NULL,
    item_time_enc  = item_time_enc,
    item_time_ret  = item_time_ret,
    item_lag       = item_lag,
    item_run_enc   = if (!is.null(item_run_enc)) factor(item_run_enc) else NULL,
    item_run_ret   = if (!is.null(item_run_ret)) factor(item_run_ret) else NULL,
    confound_rdms  = confound_rdms
  )
}

