# rmvpa DSL: which nulls the R code can compute, per verb, method and scope.
#
# This whitelist is taken from the package source, not from what a model
# might plausibly support. A cell without a valid null gets `test: none`
# ("descriptive only"); asking for a permutation there is an error that names
# the R piece that would be needed.

.dsl_inference_support <- function(verb, vb, scope_kind) {
  sl <- identical(scope_kind, "searchlight")
  roi <- identical(scope_kind, "regions")
  none <- function(missing) {
    list(kind = "none", missing = missing,
         description = "none: descriptive only (no permutation path in rMVPA for this analysis)")
  }
  engine <- function(extra = "") {
    list(kind = "searchlight_permutation", missing = NULL,
         description = paste0("pooled searchlight permutation null (run_permutation_searchlight), ",
                              "adjusted by searchlight size; not a voxelwise exact test", extra))
  }
  # rMVPA >= the release with run_permutation_regional(): per-region exact
  # permutation tests, sharing permutation_control() with the searchlight.
  regional <- function(extra = "") {
    if (!.dsl_has_regional_permutation()) return(none("run_permutation_regional() (a newer rMVPA)"))
    list(kind = "regional_permutation", missing = NULL,
         description = paste0("per-region exact permutation null (run_permutation_regional): every region is refitted ",
                              "under every relabelling; FDR across regions, and family-wise (minP) p-values", extra))
  }
  model_nperm <- function() {
    list(kind = "model_nperm", missing = NULL,
         description = paste0("model-internal permutations, one p-value per ", if (sl) "centre" else "region",
                              ", uncorrected; smallest attainable p is 1/(permutations + 1)"))
  }

  switch(verb,
    decode = if (sl) engine() else if (roi) regional() else none("a whole-brain permutation path (use regions with a whole-brain mask)"),
    rsa = {
      if (!sl && !roi) return(none("a whole-brain permutation path (use regions with a whole-brain mask)"))
      if (identical(vb$distance, "crossvalidated_euclidean")) {
        return(none("permute_labels() support for crossvalidated_euclidean RSA"))
      }
      k <- length(vb$predictors) + length(vb$nuisance)
      joint <- "; joint null of no association with any predictor (rMVPA refuses per-predictor nulls for regression with more than one predictor), so per-predictor maps have no conditional test"
      wants_joint <- vb$method %in% c("lm", "rfit") && k > 1L
      if (sl) { if (wants_joint) engine(joint) else engine() } else { if (wants_joint) regional(joint) else regional() }
    },
    rsa_vector = , rsa_feature = if (sl || roi) model_nperm() else none("whole-brain permutations for this model"),
    transfer = none(if (identical(vb$method, "naive") && is.null(vb$link)) {
      "a permutation-aware naive_xdec kernel"
    } else "permute_labels() methods that preserve train/test pairings"),
    rsa_era = none("a permute_labels() method that shuffles ERA item pairings"),
    none(sprintf("a permutation method for %s", verb))
  )
}

.dsl_has_regional_permutation <- function() "run_permutation_regional" %in% getNamespaceExports("rMVPA")

.dsl_resolve_inference <- function(support, user, perm_defaults, dlayer) {
  origins <- character()
  if (!is.null(user$test) && identical(user$test, "permutation") && identical(support$kind, "none")) {
    stop(paste0("inference.test: permutation is not available for this analysis. rMVPA has no valid null here; ",
                "it would need ", support$missing, ". Use `inference: { test: none }` or drop the inference section."),
         call. = FALSE)
  }
  tag <- function(v, layer) {
    for (k in names(v)) origins[[paste0("inference.", k)]] <<- layer
    v
  }
  if (identical(support$kind, "none") || identical(user$test, "none")) {
    extra <- setdiff(names(user), "test")
    if (length(extra) && identical(user$test, "none")) {
      stop("inference: with test: none, remove the other inference keys (", paste(extra, collapse = ", "), ").",
           call. = FALSE)
    }
    return(list(value = tag(list(test = "none"), if (is.null(user$test)) dlayer else "analysis"),
                origins = origins))
  }

  if (identical(support$kind, "model_nperm")) {
    allowed <- c("test", "permutations", "correction")
    bad <- setdiff(names(user), allowed)
    if (length(bad)) {
      stop("inference: ", paste(bad, collapse = ", "), " do not apply to model-internal permutations. ",
           "Allowed keys: ", paste(allowed, collapse = ", "), ".", call. = FALSE)
    }
    if (identical(user$correction, "fdr")) {
      stop("inference.correction: fdr is not applied by rMVPA to model-internal permutation p-values. ",
           "Use correction: none and correct the saved p-values downstream.", call. = FALSE)
    }
    base <- tag(list(test = "permutation", permutations = perm_defaults$permutations, correction = "none"), dlayer)
  } else if (identical(support$kind, "regional_permutation")) {
    # Each region has its own exact null: the pooled-null keys do not apply.
    allowed <- c("test", "permutations", "shuffle", "correction")
    bad <- setdiff(names(user), allowed)
    if (length(bad)) {
      stop("inference: ", paste(bad, collapse = ", "), " do not apply to regional permutations (each region has its ",
           "own exact null; nothing is pooled or binned). Allowed keys: ", paste(allowed, collapse = ", "), ".", call. = FALSE)
    }
    base <- tag(perm_defaults[intersect(allowed, names(perm_defaults))], dlayer)
    # Exact per-region nulls resolve p only to 1/(permutations + 1), and FDR or
    # family-wise correction over many regions needs that floor well below
    # .05 / regions; the pooled searchlight null does not. So regions default
    # to 1000 relabellings (cheap on the prepared engine).
    base$permutations <- 1000L
    origins[["inference.permutations"]] <- paste(dlayer, "(regions: exact per-region null)")
  } else {
    base <- tag(perm_defaults, dlayer)
  }
  for (k in names(user)) {
    base[[k]] <- user[[k]]
    origins[[paste0("inference.", k)]] <- "analysis"
  }
  list(value = base, origins = origins)
}

# permutation_control() arguments for the searchlight permutation engine.
.dsl_permutation_control <- function(inf, seed, rsa_joint = FALSE) {
  args <- list(
    n_perm = as.integer(inf$permutations),
    shuffle = switch(inf$shuffle, within_run = "within_block", inf$shuffle),
    null_method = inf$null_method,
    adjust_by = inf$adjust_by,
    subsample = inf$subsample,
    correction = inf$correction,
    seed = seed,
    rsa_null = if (rsa_joint) "joint" else "individual"
  )
  do.call(permutation_control, Filter(Negate(is.null), args))
}
