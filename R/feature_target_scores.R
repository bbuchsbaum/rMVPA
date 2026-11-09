# Scores and permutation nulls for held-out feature-space predictions --------
#
# Shared by every model that predicts, or is summarised over, a structured
# feature space: sets (bands), families of sets, and contrasts over feature
# columns. Scores are computed from one pair of held-out matrices, so family
# and contrast summaries never require a refit.

# Row permutations that respect temporal autocorrelation.
.br_perm_indices <- function(T_rec,
                             strategy = c("circular_shift", "block_shuffle"),
                             block_size = NULL,
                             min_shift = 1L) {
  strategy <- match.arg(strategy)
  idx <- seq_len(T_rec)
  if (T_rec <= 2L) return(idx)

  if (strategy == "circular_shift") {
    # Shifts smaller than min_shift (in either direction) leave the series
    # nearly aligned with itself and are excluded when the series allows it.
    min_shift <- as.integer(max(1L, min(min_shift, T_rec %/% 2L)))
    shift <- min_shift - 1L + sample.int(T_rec - 2L * min_shift + 1L, 1L)
    return(((idx - 1L + shift) %% T_rec) + 1L)
  }

  if (is.null(block_size)) {
    block_size <- max(2L, floor(sqrt(T_rec)))
  }
  block_size <- as.integer(max(1L, min(T_rec, block_size)))

  blocks <- split(idx, ceiling(idx / block_size))
  if (length(blocks) <= 1L) {
    return(sample(idx))
  }

  ord <- sample.int(length(blocks))
  perm <- unlist(blocks[ord], use.names = FALSE)
  if (identical(perm, idx)) {
    ord <- sample.int(length(blocks))
    perm <- unlist(blocks[ord], use.names = FALSE)
  }
  perm
}


# Named column groups: every set, then every family.
.fts_groups <- function(set, families = NULL) {
  set <- as.character(set)
  groups <- split(seq_along(set), factor(set, levels = unique(set)))
  for (label in names(families)) {
    groups[[label]] <- which(set %in% families[[label]])
  }
  groups
}

# Scores of `predicted` against `observed` (rows = observations, columns =
# features). Returns a named vector: `all_r2`/`all_cor` (mean over columns),
# `<group>_r2`/`<group>_cor` for each set and family (mean over its columns,
# so every feature counts equally), and `<contrast>_r2`/`<contrast>_cor` for
# each contrast, scored on the projected series.
.feature_target_scores <- function(observed, predicted, set,
                                   families = NULL, contrasts = NULL) {
  if (is.null(colnames(observed))) {
    colnames(observed) <- paste0("feature_", seq_len(ncol(observed)))
  }
  colnames(predicted) <- colnames(observed)
  columns <- .banded_ridge_score(observed, predicted)
  average <- function(index, field) {
    values <- columns[[field]][index]
    if (all(is.na(values))) NA_real_ else mean(values, na.rm = TRUE)
  }
  groups <- c(list(all = seq_len(ncol(observed))), .fts_groups(set, families))
  out <- unlist(lapply(names(groups), function(label) {
    stats::setNames(
      c(average(groups[[label]], "r2"), average(groups[[label]], "correlation")),
      paste0(label, c("_r2", "_cor"))
    )
  }))
  if (!is.null(contrasts)) {
    projected <- .banded_ridge_score(observed %*% contrasts,
                                     predicted %*% contrasts)
    out <- c(out, unlist(lapply(seq_len(ncol(contrasts)), function(k) {
      stats::setNames(c(projected$r2[[k]], projected$correlation[[k]]),
                      paste0(colnames(contrasts)[[k]], c("_r2", "_cor")))
    })))
  }
  out
}

# Null distribution of the scores: predictions are shifted in time against the
# observations, which breaks their alignment but keeps the autocorrelation of
# both. Returns an n_perm x n_scores matrix.
.feature_target_null <- function(observed, predicted, set, families = NULL,
                                 contrasts = NULL, n_perm = 1000L,
                                 strategy = "circular_shift",
                                 block_size = NULL, min_shift = 1L) {
  n <- nrow(observed)
  t(vapply(seq_len(n_perm), function(i) {
    index <- .br_perm_indices(n, strategy = strategy, block_size = block_size,
                              min_shift = min_shift)
    .feature_target_scores(observed, predicted[index, , drop = FALSE], set,
                           families, contrasts)
  }, .feature_target_scores(observed, predicted, set, families, contrasts)))
}

# One-sided permutation p-values with the observed value counted in the null.
.fts_p_values <- function(observed_scores, null_scores) {
  exceed <- colSums(sweep(null_scores, 2L, observed_scores, ">="), na.rm = TRUE)
  valid <- colSums(!is.na(null_scores))
  stats::setNames((exceed + 1) / (valid + 1), names(observed_scores))
}
