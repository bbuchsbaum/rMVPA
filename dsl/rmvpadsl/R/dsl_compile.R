# rmvpa DSL: compiling each verb to an ordinary rMVPA model_spec.
#
# Each compiler receives a context with the resolved verb block (`vb`), the
# loaded train (and test) set, the dataset and the fold settings, and returns
# list(model_spec, info, rsa_joint). Nothing here replaces the model_spec /
# fit_roi / roi_result contract; the specs run through run_searchlight(),
# run_regional() or run_global() like any hand-written analysis.

# Which scopes each verb supports, from the R code: run_global() has methods
# only for mvpa_model and pattern_model, and banded_ridge_model is fitted
# voxelwise by run_banded_ridge().
.dsl_scope_support <- function(verb) {
  switch(verb,
    decode = , pattern = c("searchlight", "regions", "whole_brain"),
    encode = "whole_brain",
    custom = c("searchlight", "regions", "whole_brain"),
    # contrast_rsa_model needs a searchlight centre voxel (train_model.contrast_rsa_model).
    rsa_contrast = "searchlight",
    transfer = , rsa = , rsa_vector = , rsa_feature = , rsa_era = ,
    rsa_connect = , rsa_mediate = , manova = c("searchlight", "regions")
  )
}

# Optional packages a verb needs at fit time. Missing ones are a check error,
# because rMVPA reports them per ROI, which would leave empty outputs.
.dsl_required_packages <- function(verb, vb) {
  pk <- switch(verb,
    decode = .dsl_classifier_packages(vb$method),
    transfer = if (vb$method %in% names(MVPAModels)) .dsl_classifier_packages(vb$method),
    rsa = if (identical(vb$method, "rfit")) "Rfit",
    rsa_feature = switch(vb$method, pls = , pca = "pls", glmnet = "glmnet", NULL),
    manova = "ffmanova",
    NULL
  )
  unique(pk)
}

.dsl_classifier_packages <- function(method) {
  m <- MVPAModels[[method]]
  lib <- if (is.list(m)) m$library else NULL
  if (is.null(lib)) character() else setdiff(as.character(lib), "rMVPA")
}

.dsl_check_scope <- function(verb, scope_kind, where) {
  ok <- .dsl_scope_support(verb)
  if (!scope_kind %in% ok) {
    hint <- if (identical(scope_kind, "whole_brain")) {
      " (rMVPA's run_global() supports only decode and pattern; use scope: { regions: { atlas: <whole-brain mask> } })"
    } else if (identical(verb, "encode")) {
      " (banded ridge is fitted voxelwise over the whole mask)"
    } else ""
    stop(sprintf("%sscope: %s is not available for %s; use one of: %s%s", where, scope_kind, verb,
                 paste(ok, collapse = ", "), hint), call. = FALSE)
  }
}

.dsl_term <- function(col) stats::as.formula(paste0("~`", col, "`"))

.dsl_block <- function(ctx, design = ctx$train$design) {
  by <- ctx$cfg$folds$by
  if (!is.null(by) && by %in% names(design)) design[[by]] else NULL
}

.dsl_cv <- function(ctx) .dsl_crossval(ctx$cfg$folds, ctx$train$design, ctx$where)

.dsl_items <- function(ctx, design = ctx$train$design) {
  .dsl_require_columns(design, ctx$vb$item, ctx$where)
  as.character(design[[ctx$vb$item]])
}

.dsl_nperm <- function(ctx) {
  inf <- ctx$cfg$inference
  if (identical(inf$test, "permutation") && !isTRUE(ctx$draft)) as.integer(inf$permutations) else 0L
}

.dsl_target <- function(ctx, design, type = NULL) {
  .dsl_require_columns(design, ctx$vb$target, ctx$where)
  y <- design[[ctx$vb$target]]
  type <- type %||% ctx$vb$target_type %||% (if (is.numeric(y)) "regression" else "classification")
  y <- if (identical(type, "classification")) factor(as.character(y)) else as.numeric(y)
  list(y = y, type = type)
}

.dsl_syntactic <- function(nms, what, where) {
  bad <- nms[!grepl("^[A-Za-z][A-Za-z0-9_.]*$", nms)]
  if (length(bad)) {
    stop(sprintf("%s %s name '%s' must start with a letter and use only letters, digits, _ or .",
                 where, what, bad[1]), call. = FALSE)
  }
}

.dsl_feature_matrix <- function(ctx, spaces, items) {
  mats <- lapply(names(spaces), function(nm) {
    m <- .dsl_read_features(spaces[[nm]], ctx$root, unique(items), sprintf("%s feature_spaces.%s", ctx$where, nm))
    m <- m[items, , drop = FALSE]
    colnames(m) <- paste0(nm, ".", colnames(m))
    rownames(m) <- NULL
    m
  })
  names(mats) <- names(spaces)
  mats
}

.dsl_compilers <- function() {
  list(
    decode = function(ctx) {
      d <- ctx$train$design
      t <- .dsl_target(ctx, d)
      des <- mvpa_design(d, y_train = t$y, block_var = .dsl_block(ctx))
      fs <- NULL
      if (!is.null(ctx$vb$select)) {
        s <- ctx$vb$select
        fs <- if (!is.null(s$top)) {
          feature_selector(switch(s$by, ftest = "FTest", catscore = "catscore"), "top_k", s$top)
        } else feature_selector(switch(s$by, ftest = "FTest", catscore = "catscore"), "top_p", s$top_fraction)
      }
      ms <- mvpa_model(load_model(ctx$vb$method), ctx$dataset, des, model_type = t$type,
                       crossval = .dsl_cv(ctx), feature_selector = fs)
      list(model_spec = ms, info = list(target_type = t$type, classifier = ctx$vb$method))
    },

    transfer = function(ctx) {
      tr <- ctx$train$design; te <- ctx$test$design
      .dsl_require_columns(te, ctx$vb$target, sprintf("%s test set:", ctx$where))
      t <- .dsl_target(ctx, tr)
      y_te <- if (identical(t$type, "classification")) as.character(te[[ctx$vb$target]]) else as.numeric(te[[ctx$vb$target]])
      if (identical(t$type, "classification")) {
        extra <- setdiff(unique(y_te), levels(t$y))
        if (length(extra)) {
          stop(sprintf("%s test-set %s level(s) never seen in training: %s", ctx$where, ctx$vb$target,
                       paste(extra, collapse = ", ")), call. = FALSE)
        }
        y_te <- factor(y_te, levels = levels(t$y))
      }
      link <- ctx$vb$link
      if (!is.null(link)) {
        .dsl_require_columns(tr, link, sprintf("%s train set:", ctx$where))
        .dsl_require_columns(te, link, sprintf("%s test set:", ctx$where))
      }
      des <- mvpa_design(tr, test_design = te, y_train = t$y, y_test = y_te, block_var = .dsl_block(ctx))
      m <- ctx$vb$method
      ms <- switch(m,
        naive = naive_xdec_model(ctx$dataset, des, link_by = link),
        remap_rrr = remap_rrr_model(ctx$dataset, des, link_by = link, rank = ctx$vb$rank %||% "auto",
                                    leave_one_key_out = TRUE),
        subspace_alignment = subspace_alignment_model(ctx$dataset, des, d = ctx$vb$dimensions %||% 20L),
        mvpa_model(load_model(m), ctx$dataset, des, model_type = t$type, crossval = .dsl_cv(ctx))
      )
      list(model_spec = ms, info = list(target_type = t$type, method = m))
    },

    rsa = function(ctx) {
      d <- ctx$train$design
      items <- .dsl_items(ctx)
      vb <- ctx$vb
      .dsl_syntactic(c(names(vb$predictors), names(vb$nuisance)), "RDM", ctx$where)
      cv_mode <- identical(vb$distance, "crossvalidated_euclidean")
      levels_i <- unique(items)
      read <- function(p, nm) {
        m <- .dsl_read_rdm(p, ctx$root, levels_i, sprintf("%s %s", ctx$where, nm))
        if (cv_mode) m else m[items, items, drop = FALSE]
      }
      preds <- Map(read, vb$predictors, names(vb$predictors))
      nuis <- if (length(vb$nuisance)) Map(read, vb$nuisance, names(vb$nuisance)) else list()
      form <- stats::as.formula(paste("~", paste(names(preds), collapse = " + ")))
      block <- .dsl_block(ctx)
      if (cv_mode) {
        rdes <- rsa_design(form, data = preds, nuisance = nuis, condition_ids = levels_i)
        ms <- rsa_model(ctx$dataset, rdes, distmethod = "crossvalidated_euclidean", regtype = vb$method,
                        semipartial = vb$semipartial, condition_labels = items, crossval = .dsl_cv(ctx))
      } else {
        dat <- c(preds, if (!is.null(block)) list(rmvpa_block = block))
        rdes <- rsa_design(form, data = dat, block_var = if (!is.null(block)) "rmvpa_block",
                           nuisance = nuis, keep_intra_run = vb$within_run_pairs)
        ms <- rsa_model(ctx$dataset, rdes, distmethod = vb$distance, regtype = vb$method,
                        semipartial = vb$semipartial)
      }
      k <- length(preds) + length(nuis)
      list(model_spec = ms, rsa_joint = vb$method %in% c("lm", "rfit") && k > 1L,
           info = list(items = length(levels_i), rdm_level = if (cv_mode) "item" else "trial (expanded from item RDMs)"))
    },

    rsa_vector = function(ctx) {
      items <- .dsl_items(ctx)
      D <- .dsl_read_rdm(ctx$vb$distances, ctx$root, unique(items), paste(ctx$where, "distances"))
      block <- .dsl_block(ctx)
      if (is.null(block)) stop(ctx$where, " rsa_vector needs the folds.by column as its block variable", call. = FALSE)
      rdes <- vector_rsa_design(D = D, labels = items, block_var = block)
      ms <- vector_rsa_model(ctx$dataset, rdes, rsa_simfun = ctx$vb$method, nperm = .dsl_nperm(ctx))
      list(model_spec = ms, info = list(items = length(unique(items)), nperm = .dsl_nperm(ctx)))
    },

    rsa_contrast = function(ctx) {
      d <- ctx$train$design
      items <- .dsl_items(ctx)
      des <- mvpa_design(d, y_train = factor(items), block_var = .dsl_block(ctx))
      levs <- levels(y_train(des))
      f <- .dsl_resolve_path(ctx$vb$contrasts, ctx$root)
      if (!file.exists(f)) stop(ctx$where, " contrasts: file not found: ", f, call. = FALSE)
      t <- .dsl_read_table(f)
      C <- as.matrix(t[, -1, drop = FALSE]); storage.mode(C) <- "double"
      rownames(C) <- as.character(t[[1]])
      miss <- setdiff(levs, rownames(C))
      if (length(miss)) stop(sprintf("%s contrasts: no row for item(s) %s", ctx$where, paste(miss, collapse = ", ")), call. = FALSE)
      C <- C[levs, , drop = FALSE]
      msd <- msreve_design(des, C)
      ms <- contrast_rsa_model(ctx$dataset, msd, estimation_method = ctx$vb$method,
                               output_metric = ctx$vb$output)
      list(model_spec = ms, info = list(items = length(levs), contrasts = colnames(C)))
    },

    rsa_feature = function(ctx) {
      items <- .dsl_items(ctx)
      mats <- .dsl_feature_matrix(ctx, ctx$vb$feature_spaces, items)
      Fm <- do.call(cbind, unname(mats))
      fdes <- feature_rsa_design(F = Fm, labels = paste0("obs", seq_len(nrow(Fm))))
      ms <- feature_rsa_model(ctx$dataset, fdes, method = ctx$vb$method, crossval = .dsl_cv(ctx),
                              nperm = .dsl_nperm(ctx))
      list(model_spec = ms, info = list(features = ncol(Fm), nperm = .dsl_nperm(ctx)))
    },

    rsa_era = function(ctx) {
      tr <- ctx$train$design; te <- ctx$test$design
      link <- ctx$vb$link
      .dsl_require_columns(tr, link, sprintf("%s train set:", ctx$where))
      .dsl_require_columns(te, link, sprintf("%s test set:", ctx$where))
      y_tr <- factor(as.character(tr[[link]]))
      y_te <- factor(as.character(te[[link]]), levels = levels(y_tr))
      if (anyNA(y_te)) {
        stop(sprintf("%s test-set %s values never seen in the train set", ctx$where, link), call. = FALSE)
      }
      des <- mvpa_design(tr, test_design = te, y_train = y_tr, y_test = y_te, block_var = .dsl_block(ctx))
      ms <- switch(ctx$vb$method,
        similarity = era_rsa_model(ctx$dataset, des, key_var = .dsl_term(link)),
        partition = era_partition_model(ctx$dataset, des, key_var = .dsl_term(link))
      )
      list(model_spec = ms, info = list(items = nlevels(y_tr)))
    },

    rsa_connect = function(ctx) {
      d <- ctx$train$design
      items <- .dsl_items(ctx)
      des <- mvpa_design(d, y_train = factor(items), block_var = .dsl_block(ctx))
      seed <- .dsl_read_rdm(ctx$vb$seed_region, ctx$root, unique(items), paste(ctx$where, "seed_region"))
      rn <- repnet_design(des, key_var = .dsl_term(ctx$vb$item), seed_rdm = seed)
      ms <- repnet_model(ctx$dataset, des, rn)
      list(model_spec = ms, info = list(items = length(unique(items))))
    },

    rsa_mediate = function(ctx) {
      d <- ctx$train$design
      items <- .dsl_items(ctx)
      lv <- unique(items)
      des <- mvpa_design(d, y_train = factor(items), block_var = .dsl_block(ctx))
      X <- .dsl_read_rdm(ctx$vb$predictor, ctx$root, lv, paste(ctx$where, "predictor"))
      Y <- .dsl_read_rdm(ctx$vb$outcome, ctx$root, lv, paste(ctx$where, "outcome"))
      md <- repmed_design(items = lv, X_rdm = X, Y_rdm = Y)
      ms <- repmed_model(ctx$dataset, des, md, key_var = .dsl_term(ctx$vb$item))
      list(model_spec = ms, info = list(items = length(lv)))
    },

    encode = function(ctx) {
      items <- .dsl_items(ctx)
      mats <- .dsl_feature_matrix(ctx, ctx$vb$feature_spaces, items)
      .dsl_syntactic(names(mats), "feature space", ctx$where)
      fs <- feature_sets(mats)
      block <- .dsl_block(ctx)
      fd <- feature_sets_design(fs, block_var_train = block)
      ms <- banded_ridge_model(ctx$dataset, fd, outer_crossval = .dsl_cv(ctx))
      list(model_spec = ms, info = list(feature_spaces = names(mats)))
    },

    pattern = function(ctx) {
      d <- ctx$train$design
      t <- .dsl_target(ctx, d)
      des <- if (is.null(ctx$test)) {
        mvpa_design(d, y_train = t$y, block_var = .dsl_block(ctx))
      } else {
        mvpa_design(d, test_design = ctx$test$design, y_train = t$y,
                    y_test = ctx$test$design[[ctx$vb$target]], block_var = .dsl_block(ctx))
      }
      ms <- pattern_model(ctx$dataset, des, crossval = .dsl_cv(ctx), max_rank = as.integer(ctx$vb$max_rank))
      list(model_spec = ms, info = list(target_type = t$type))
    },

    manova = function(ctx) {
      d <- ctx$train$design
      .dsl_require_columns(d, ctx$vb$terms, ctx$where)
      form <- stats::as.formula(paste("~", paste(sprintf("`%s`", ctx$vb$terms), collapse = " + ")))
      md <- manova_design(form, d)
      list(model_spec = manova_model(ctx$dataset, md), info = list(terms = ctx$vb$terms))
    },

    custom = function(ctx) {
      f <- .dsl_resolve_path(ctx$vb$r_file, ctx$root)
      if (!file.exists(f)) stop(ctx$where, " custom.r_file not found: ", f, call. = FALSE)
      env <- new.env(parent = globalenv())
      sys.source(f, envir = env)
      fname <- ctx$vb[["function"]]
      if (!exists(fname, envir = env, inherits = FALSE) || !is.function(env[[fname]])) {
        stop(sprintf("%s custom: %s does not define a function '%s'", ctx$where, f, fname), call. = FALSE)
      }
      config <- list(train_design = ctx$train$design, test_design = ctx$test$design,
                     folds = ctx$cfg$folds, subject = ctx$subject, root = ctx$root)
      ms <- env[[fname]](ctx$dataset, ctx$train$design, config)
      if (!inherits(ms, "model_spec")) {
        stop(sprintf("%s custom: %s() must return a model_spec, got <%s>", ctx$where, fname,
                     paste(class(ms), collapse = "/")), call. = FALSE)
      }
      list(model_spec = ms, info = list(custom = basename(f), builder = fname))
    }
  )
}
