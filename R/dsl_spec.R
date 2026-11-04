# rmvpa DSL 0.1: vocabulary, closed schema and versioned defaults.
#
# The schema is plain R data. Every key has exactly one shape. Unknown keys
# are errors with a "did you mean" suggestion. Defaults are pinned to the
# `rmvpa:` version written in the project file.

.dsl_versions <- "0.1"

# ---- schema node constructors (internal) ---------------------------------

.n_map <- function(fields, required = character(), exactly_one = NULL, doc = NULL) {
  list(type = "map", fields = fields, required = required,
       exactly_one = exactly_one, doc = doc)
}
.n_dict <- function(value, doc = NULL) list(type = "dict", value = value, doc = doc)
.n_str <- function(doc = NULL) list(type = "str", doc = doc)
.n_path <- function(doc = NULL) list(type = "path", doc = doc)
.n_bool <- function(doc = NULL) list(type = "bool", doc = doc)
.n_int <- function(min = -Inf, max = Inf, doc = NULL) list(type = "int", min = min, max = max, doc = doc)
.n_num <- function(min = -Inf, max = Inf, doc = NULL) list(type = "num", min = min, max = max, doc = doc)
.n_enum <- function(values, doc = NULL) list(type = "enum", values = values, doc = doc)
.n_strs <- function(doc = NULL) list(type = "strs", doc = doc)
.n_enums <- function(values, doc = NULL) list(type = "enums", values = values, doc = doc)
# Filter value: a scalar (equality) or a map with a single `in:` list.
.n_filter <- function(doc = NULL) list(type = "filter", doc = doc)
# Classifier names come from the live registry (register_mvpa_model()).
.n_classifier <- function(doc = NULL) list(type = "classifier", doc = doc)

# ---- shared sections -----------------------------------------------------

.dsl_set_entities <- c("task", "ses", "acq", "desc", "rec", "ce", "dir")

.dsl_set_node <- function() {
  fields <- stats::setNames(
    lapply(.dsl_set_entities, function(e) .n_str(paste0("BIDS entity `", e, "`"))),
    .dsl_set_entities
  )
  fields$keep <- .n_dict(.n_filter(), "Row filters applied before folds: equality or {in: [...]}")
  .n_map(fields, doc = "A named set of beta series, selected by BIDS entities")
}

.dsl_assertions <- c(
  "classes_in_every_fold", "no_item_across_folds",
  "one_run_per_fold", "train_test_disjoint_runs"
)

.dsl_data_node <- function() {
  .n_map(list(
    bids = .n_path("BIDS-derivatives root holding the beta series"),
    files = .n_map(list(
      betas = .n_str("Path template for the 4D beta image, e.g. betas/sub-{subject}_task-{task}_betas.nii.gz"),
      design = .n_str("Path template for the per-volume design table (TSV)"),
      mask = .n_str("Path template for the brain mask")
    ), required = c("betas", "design"), doc = "Explicit file templates instead of BIDS discovery"),
    subjects = .n_strs("Subject labels without the sub- prefix; omitted means all discovered"),
    space = .n_str("BIDS space entity of the betas and mask"),
    suffix = .n_str("BIDS suffix of the beta image"),
    mask = .n_str("Mask desc label (e.g. brain) or a path/template"),
    sets = .n_dict(.dsl_set_node(), "Named sets of beta series")
  ), required = "sets", exactly_one = c("bids", "files"))
}

.dsl_folds_node <- function() {
  .n_map(list(
    by = .n_str("Design column that defines the folds"),
    scheme = .n_enum(c("leave_one_out", "twofold", "bootstrap", "kfold")),
    nfolds = .n_int(2, doc = "Number of folds for scheme: kfold"),
    repeats = .n_int(1, doc = "Repeats for scheme: twofold or bootstrap"),
    assert = .n_enums(.dsl_assertions, "Named fold checks run by `rmvpa check`")
  ))
}

.dsl_scope_node <- function() {
  .n_map(list(
    searchlight = .n_map(list(
      radius = .n_num(0, doc = "Searchlight radius in mm"),
      sampling = .n_enum(c("standard", "randomized", "resampled")),
      iterations = .n_int(1, doc = "Iterations for randomized or resampled sampling")
    ), required = "radius"),
    regions = .n_map(list(
      atlas = .n_path("Integer-labelled atlas image (dseg) in the data space"),
      labels = .n_path("Optional label table (TSV with index and name columns)")
    ), required = "atlas"),
    whole_brain = .n_map(list())
  ), exactly_one = c("searchlight", "regions", "whole_brain"))
}

.dsl_inference_node <- function() {
  .n_map(list(
    test = .n_enum(c("permutation", "none")),
    permutations = .n_int(1),
    shuffle = .n_enum(c("within_run", "circular_shift", "global")),
    null_method = .n_enum(c("adjusted", "global")),
    adjust_by = .n_enum(c("nfeatures", "redundancy", "both")),
    subsample = .n_num(0),
    correction = .n_enum(c("fdr", "none"))
  ))
}

.dsl_run_node <- function() {
  .n_map(list(
    workers = .n_int(1),
    seed = .n_int(),
    backend = .n_enum(c("default", "shard", "auto"))
  ))
}

# ---- verbs ---------------------------------------------------------------

.dsl_rdm_dict <- function(doc) .n_dict(.n_path("Square RDM table (CSV/TSV) with item names in the first column"), doc)

.dsl_verbs <- function() {
  list(
    decode = list(family = "decoding", node = .n_map(list(
      set = .n_str(), target = .n_str("Design column to predict"),
      target_type = .n_enum(c("classification", "regression")),
      method = .n_classifier("Registered classifier"),
      select = .n_map(list(by = .n_enum(c("ftest", "catscore")),
                           top = .n_int(1), top_fraction = .n_num(0, 1)),
                      required = "by", exactly_one = c("top", "top_fraction"))
    ), required = "target")),
    decode_trialwise = list(family = "decoding", node = .n_map(list(
      set = .n_str(), target = .n_str(),
      target_type = .n_enum(c("classification", "regression")),
      events = .n_path("Trial-level event table used to build the trial design")
    ), required = c("target", "events"))),
    transfer = list(family = "decoding", node = .n_map(list(
      train_set = .n_str(), test_set = .n_str(), target = .n_str(),
      method = .n_str("naive, remap_rrr, subspace_alignment, or a registered classifier"),
      link = .n_str("Column pairing rows across sets"),
      rank = .n_int(1, doc = "remap_rrr rank; omitted means chosen automatically"),
      dimensions = .n_int(1, doc = "subspace_alignment dimensionality")
    ), required = c("train_set", "test_set", "target"))),
    rsa = list(family = "rsa", node = .n_map(list(
      set = .n_str(), item = .n_str("Design column indexing the RDM rows"),
      predictors = .dsl_rdm_dict("Model RDMs"),
      nuisance = .dsl_rdm_dict("Nuisance RDMs, always included as predictors"),
      distance = .n_enum(c("pearson", "spearman", "crossvalidated_euclidean")),
      method = .n_enum(c("pearson", "spearman", "lm", "rfit")),
      semipartial = .n_bool(),
      within_run_pairs = .n_bool()
    ), required = c("item", "predictors"))),
    rsa_vector = list(family = "rsa", node = .n_map(list(
      set = .n_str(), item = .n_str(),
      distances = .n_path("One square RDM table over items"),
      method = .n_enum(c("pearson", "spearman"))
    ), required = c("item", "distances"))),
    rsa_contrast = list(family = "rsa", node = .n_map(list(
      set = .n_str(), item = .n_str(),
      contrasts = .n_path("TSV: first column item, one column per contrast"),
      method = .n_enum(c("average", "L2_norm", "crossnobis")),
      output = .n_enum(c("beta_delta", "beta_only", "delta_only", "recon_score",
                         "beta_delta_norm", "beta_delta_reliable", "composite"))
    ), required = c("item", "contrasts"))),
    rsa_feature = list(family = "rsa", node = .n_map(list(
      set = .n_str(), item = .n_str(),
      feature_spaces = .n_dict(.n_path("Feature table (CSV/TSV), item names in the first column")),
      method = .n_enum(c("pls", "pca", "glmnet", "ridge"))
    ), required = c("item", "feature_spaces"))),
    rsa_era = list(family = "rsa", node = .n_map(list(
      train_set = .n_str(), test_set = .n_str(),
      link = .n_str("Column pairing encoding and retrieval items"),
      method = .n_enum(c("similarity", "partition"))
    ), required = c("train_set", "test_set", "link"))),
    rsa_connect = list(family = "rsa", node = .n_map(list(
      set = .n_str(), item = .n_str(),
      seed_region = .n_path("RDM table of the seed region over items")
    ), required = c("item", "seed_region"))),
    rsa_mediate = list(family = "rsa", node = .n_map(list(
      set = .n_str(), item = .n_str(),
      predictors = .dsl_rdm_dict("Exactly one predictor RDM (X)"),
      mediator = .n_path("Mediator RDM table (M)"),
      outcome = .n_path("Outcome RDM table (Y)")
    ), required = c("item", "predictors", "mediator", "outcome"))),
    encode = list(family = "models", node = .n_map(list(
      set = .n_str(), item = .n_str(),
      feature_spaces = .n_dict(.n_path("Feature table (CSV/TSV), item names in the first column")),
      method = .n_enum(c("banded_ridge", "grouped_ridge"))
    ), required = c("item", "feature_spaces"))),
    pattern = list(family = "models", node = .n_map(list(
      set = .n_str(), target = .n_str(),
      max_rank = .n_int(1)
    ), required = "target")),
    manova = list(family = "models", node = .n_map(list(
      set = .n_str(), terms = .n_strs("Design columns entering the MANOVA")
    ), required = "terms")),
    custom = list(family = "escape", node = .n_map(list(
      set = .n_str(), r_file = .n_path("R file defining the builder"),
      `function` = .n_str("Builder: function(dataset, design, config) returning a model_spec")
    ), required = "r_file"))
  )
}

.dsl_verb_names <- function() names(.dsl_verbs())

.dsl_project_node <- function() {
  .n_map(list(
    rmvpa = .n_enum(.dsl_versions),
    data = .dsl_data_node(),
    folds = .dsl_folds_node(),
    run = .dsl_run_node()
  ), required = c("rmvpa", "data"))
}

.dsl_analysis_node <- function() {
  verbs <- .dsl_verbs()
  fields <- c(
    list(extends = .n_path("Parent file, relative to this file"),
         description = .n_str()),
    lapply(verbs, `[[`, "node"),
    list(folds = .dsl_folds_node(), scope = .dsl_scope_node(),
         inference = .dsl_inference_node(), run = .dsl_run_node())
  )
  .n_map(fields, required = c("extends", "scope"), exactly_one = names(verbs))
}

# ---- defaults (pinned per version) ----------------------------------------

.dsl_defaults <- function(version = "0.1") {
  stopifnot(identical(version, "0.1"))
  list(
    data = list(suffix = "betas", mask = "brain"),
    folds = list(by = "run", scheme = "leave_one_out"),
    run = list(workers = 1L, seed = 1L, backend = "default"),
    scope = list(searchlight = list(sampling = "standard", iterations = 4L)),
    verbs = list(
      decode = list(method = "sda_notune"),
      decode_trialwise = list(),
      transfer = list(method = "naive"),
      rsa = list(distance = "spearman", method = "pearson",
                 semipartial = FALSE, within_run_pairs = FALSE),
      rsa_vector = list(method = "pearson"),
      rsa_contrast = list(method = "average", output = "beta_delta"),
      rsa_feature = list(method = "pls"),
      rsa_era = list(method = "similarity"),
      rsa_connect = list(),
      rsa_mediate = list(),
      encode = list(method = "banded_ridge"),
      pattern = list(max_rank = 8L),
      manova = list(),
      custom = list(`function` = "build_spec")
    ),
    permutation = list(test = "permutation", permutations = 100L,
                       shuffle = "within_run", null_method = "adjusted",
                       adjust_by = "nfeatures", subsample = 0.1,
                       correction = "fdr")
  )
}

# ---- validation ----------------------------------------------------------

.dsl_did_you_mean <- function(x, choices) {
  if (length(choices) == 0L) return("")
  d <- utils::adist(x, choices, ignore.case = TRUE)
  best <- which.min(d)
  if (d[best] <= max(2, nchar(x) %/% 3)) sprintf("; did you mean '%s'?", choices[best]) else ""
}

.dsl_is_map <- function(x) is.list(x) && (length(x) == 0L || !is.null(names(x)))

.dsl_validate <- function(x, node, path = "", where = "") {
  errs <- character()
  err <- function(msg) errs <<- c(errs, sprintf("%s%s: %s", where, if (nzchar(path)) path else "(top level)", msg))
  sub <- function(k) if (nzchar(path)) paste0(path, ".", k) else k
  scalar <- function(v) is.atomic(v) && length(v) == 1L && !is.na(v)

  switch(node$type,
    map = {
      if (!.dsl_is_map(x)) { err("must be a map of keys"); return(errs) }
      keys <- names(x)
      unknown <- setdiff(keys, names(node$fields))
      for (k in unknown) err(sprintf("unknown key '%s'%s", k, .dsl_did_you_mean(k, names(node$fields))))
      for (k in setdiff(node$required, keys)) err(sprintf("missing required key '%s'", k))
      if (!is.null(node$exactly_one)) {
        present <- intersect(keys, node$exactly_one)
        if (length(present) != 1L) {
          err(sprintf("needs exactly one of: %s (found %s)",
                      paste(node$exactly_one, collapse = ", "),
                      if (length(present)) paste(present, collapse = ", ") else "none"))
        }
      }
      for (k in intersect(keys, names(node$fields))) {
        if (is.null(x[[k]])) next
        errs <- c(errs, .dsl_validate(x[[k]], node$fields[[k]], sub(k), where))
      }
    },
    dict = {
      if (!.dsl_is_map(x) || length(x) == 0L) { err("must be a non-empty map of names"); return(errs) }
      for (k in names(x)) errs <- c(errs, .dsl_validate(x[[k]], node$value, sub(k), where))
    },
    str = , path = if (!(is.character(x) && scalar(x))) err("must be a single string"),
    bool = if (!(is.logical(x) && scalar(x))) err("must be true or false"),
    int = {
      if (!(is.numeric(x) && scalar(x) && x == round(x))) err("must be a whole number")
      else if (x < node$min || x > node$max) err(sprintf("must be between %s and %s", node$min, node$max))
    },
    num = {
      if (!(is.numeric(x) && scalar(x))) err("must be a number")
      else if (x < node$min || x > node$max) err(sprintf("must be between %s and %s", node$min, node$max))
    },
    enum = {
      if (!(is.character(x) && scalar(x))) err(sprintf("must be one of: %s", paste(node$values, collapse = ", ")))
      else if (!x %in% node$values) err(sprintf("'%s' is not one of: %s%s", x, paste(node$values, collapse = ", "),
                                               .dsl_did_you_mean(x, node$values)))
    },
    strs = if (!(is.character(x) && length(x) >= 1L)) err("must be a list of strings"),
    enums = {
      if (!is.character(x)) { err(sprintf("must be a list drawn from: %s", paste(node$values, collapse = ", "))); return(errs) }
      for (v in setdiff(x, node$values)) err(sprintf("'%s' is not a known value%s", v, .dsl_did_you_mean(v, node$values)))
    },
    filter = {
      ok_scalar <- is.atomic(x) && length(x) == 1L
      ok_in <- .dsl_is_map(x) && identical(names(x), "in") && is.atomic(x[["in"]]) && length(x[["in"]]) >= 1L
      if (!(ok_scalar || ok_in)) err("must be a single value (equality) or {in: [...]}")
    },
    classifier = {
      known <- names(MVPAModels)
      if (!(is.character(x) && scalar(x))) err("must be a classifier name")
      else if (!x %in% known) err(sprintf("unknown classifier '%s'%s (see `rmvpa_classifiers()`)", x,
                                          .dsl_did_you_mean(x, known)))
    },
    stop("internal: unknown schema node type ", node$type)
  )
  errs
}

#' List classifiers usable as `decode.method`
#'
#' @return A character vector of registered classifier names.
#' @export
rmvpa_classifiers <- function() sort(names(MVPAModels))
