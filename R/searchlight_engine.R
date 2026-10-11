#' @keywords internal
#' @noRd
.match_searchlight_engine <- function(engine = "auto") {
  allowed <- c(
    "auto", "legacy", "swift", "dual_lda_fast", "aggregate_fast",
    "sda_fast", "rsa_fast", "naive_xdec_fast", "era_rsa_fast"
  )
  match.arg(as.character(engine)[1], allowed)
}

#' @keywords internal
#' @noRd
.searchlight_engine_registry <- function() {
  list(
    legacy = list(
      label = "General-purpose iterator (compatibility key: legacy)",
      eligible = function(model_spec, method) TRUE
    ),
    swift = list(
      label = "SWIFT nearest-class-mean estimator (explicit opt-in only)",
      eligible = function(model_spec, method) {
        .swift_searchlight_enabled() && .is_swift_fast_path(model_spec, method)
      }
    ),
    dual_lda_fast = list(
      label = "Dual-LDA incremental fast path",
      eligible = function(model_spec, method) {
        .is_dual_lda_fast_path(model_spec, method)
      }
    ),
    sda_fast = list(
      label = "Exact sda_notune engine (shared per-voxel statistics)",
      eligible = function(model_spec, method) {
        .is_sda_fast_path(model_spec, method)
      }
    ),
    aggregate_fast = list(
      label = "Exact sphere-aggregation engine (corclass, naive_bayes)",
      eligible = function(model_spec, method) {
        .is_aggregate_fast_path(model_spec, method)
      }
    ),
    rsa_fast = list(
      label = "RSA per-sphere engine without iterator overhead",
      eligible = function(model_spec, method) {
        .is_rsa_fast_path(model_spec, method)
      }
    ),
    naive_xdec_fast = list(
      label = "Naive cross-decoding matrix fast path",
      eligible = function(model_spec, method) {
        .is_naive_xdec_fast_path(model_spec, method)
      }
    ),
    era_rsa_fast = list(
      label = "ERA-RSA direct-matrix fast path",
      eligible = function(model_spec, method) {
        .is_era_rsa_fast_path(model_spec, method)
      }
    )
  )
}

#' Summarize Available Searchlight Engines
#'
#' Returns a compact table of registered searchlight engines. When a
#' \code{model_spec} is supplied, the output also includes whether each engine
#' is currently eligible for that analysis.
#'
#' @param model_spec Optional model specification.
#' @param method Searchlight method to audit.
#'
#' @return A data frame.
#' @export
searchlight_engines <- function(model_spec = NULL,
                                method = c("standard", "randomized", "resampled")) {
  method <- match.arg(method)
  registry <- .searchlight_engine_registry()
  nms <- names(registry)

  eligible <- if (is.null(model_spec)) {
    rep(NA, length(nms))
  } else {
    vapply(
      nms,
      function(nm) isTRUE(registry[[nm]]$eligible(model_spec, method)),
      logical(1)
    )
  }

  data.frame(
    engine = nms,
    label = vapply(registry, function(x) x$label, character(1)),
    eligible = eligible,
    stringsAsFactors = FALSE
  )
}

#' Explain Searchlight Engine Selection
#'
#' Reports which searchlight engine would be selected for a given model and
#' method, along with the eligibility status of each registered engine.
#'
#' @param model_spec A model specification.
#' @param method Searchlight method to audit.
#' @param engine Requested engine policy: \code{"auto"}, \code{"legacy"}
#'   (the compatibility key for the general-purpose iterator), \code{"swift"},
#'   \code{"dual_lda_fast"}, \code{"aggregate_fast"}, \code{"sda_fast"},
#'   \code{"rsa_fast"}, \code{"naive_xdec_fast"}, or
#'   \code{"era_rsa_fast"}.
#'
#' @return A data frame with selection metadata.
#' @export
explain_searchlight_engine <- function(model_spec,
                                       method = c("standard", "randomized", "resampled"),
                                       engine = c("auto", "legacy", "swift", "dual_lda_fast", "aggregate_fast", "sda_fast", "rsa_fast", "naive_xdec_fast", "era_rsa_fast")) {
  method <- match.arg(method)
  requested <- .match_searchlight_engine(match.arg(engine))
  registry_tbl <- searchlight_engines(model_spec = model_spec, method = method)
  selected <- .resolve_searchlight_engine(model_spec, method = method, engine = requested)

  registry_tbl$requested <- requested
  registry_tbl$selected <- registry_tbl$engine == selected
  registry_tbl
}

#' @keywords internal
#' @noRd
.run_searchlight_engine <- function(model_spec, radius, method,
                                    engine = "auto",
                                    niter = 4L,
                                    combiner = "average",
                                    drop_probs = FALSE,
                                    fail_fast = FALSE,
                                    backend = c("default", "shard", "auto"),
                                    incremental = TRUE,
                                    gamma = NULL,
                                    verbose = FALSE,
                                    ...) {
  UseMethod(".run_searchlight_engine")
}

#' @keywords internal
#' @noRd
.run_searchlight_engine.default <- function(model_spec, radius, method,
                                            engine = "auto",
                                            niter = 4L,
                                            combiner = "average",
                                            drop_probs = FALSE,
                                            fail_fast = FALSE,
                                            backend = c("default", "shard", "auto"),
                                            incremental = TRUE,
                                            gamma = NULL,
                                            verbose = FALSE,
                                            ...) {
  list(
    handled = FALSE,
    result = NULL,
    engine = "legacy"
  )
}

#' Resolve the searchlight engine for a model specification
#'
#' Single entry point for deciding which engine would run for a given
#' \code{model_spec}/\code{method}/\code{engine} request. Both the diagnostic
#' helper \code{explain_searchlight_engine()} and the mvpa_model runner rely on
#' this so the reported engine matches the engine that actually executes.
#'
#' Dispatch is explicit rather than S3 because these engine helpers are
#' unexported, dotted-name functions with no \code{S3method()} registration;
#' relying on \code{UseMethod()} to reach a \code{.default} fallback is not
#' robust for such names. Model classes with a fast path are disjoint (a
#' \code{naive_xdec_model} does not inherit \code{mvpa_model}), so a small
#' \code{inherits()} ladder is unambiguous and keeps resolution in one place.
#'
#' @keywords internal
#' @noRd
.resolve_searchlight_engine <- function(model_spec, method, engine = "auto", combiner = "average") {
  if (inherits(model_spec, "rsa_model") &&
      identical(model_spec$distmethod, "crossvalidated_euclidean") &&
      !engine %in% c("auto", "legacy")) {
    stop("Crossvalidated Euclidean RSA supports engine = 'auto' or 'legacy' only; 'rsa_fast' and other optimized engines are not implemented for this model.",
         call. = FALSE)
  }
  if (inherits(model_spec, "era_rsa_model")) {
    return(.resolve_searchlight_engine.era_rsa_model(model_spec, method, engine))
  }
  if (inherits(model_spec, "naive_xdec_model")) {
    return(.resolve_searchlight_engine.naive_xdec_model(model_spec, method, engine))
  }
  if (inherits(model_spec, "mvpa_model")) {
    return(.resolve_searchlight_engine.mvpa_model(model_spec, method, engine, combiner = combiner))
  }
  if (inherits(model_spec, "rsa_model")) {
    requested <- .match_searchlight_engine(engine)
    if (requested %in% c("auto", "rsa_fast") && .is_rsa_fast_path(model_spec, method)) {
      return("rsa_fast")
    }
    return("legacy")
  }
  "legacy"
}

#' @keywords internal
#' @noRd
.resolve_searchlight_engine.mvpa_model <- function(model_spec, method, engine = "auto",
                                                   combiner = "average") {
  requested <- .match_searchlight_engine(engine)
  registry <- .searchlight_engine_registry()

  if (identical(requested, "legacy")) {
    return("legacy")
  }

  # Fast engines hard-wire the built-in combiner. A different combiner is a
  # general-path feature, so no fast engine is eligible for it.
  combiner_ok <- .fast_combiner_eligible(combiner, method)
  fast_eligible <- function(nm) {
    combiner_ok && isTRUE(registry[[nm]]$eligible(model_spec, method))
  }

  if (!identical(requested, "auto")) {
    if (fast_eligible(requested)) {
      return(requested)
    }
    return("legacy")
  }

  if (fast_eligible("dual_lda_fast")) {
    return("dual_lda_fast")
  }

  if (fast_eligible("aggregate_fast")) {
    return("aggregate_fast")
  }

  if (fast_eligible("sda_fast")) {
    return("sda_fast")
  }

  # SWIFT is never selected automatically: it computes its own z-scored
  # nearest-class-mean estimator rather than the classifier in model_spec, so
  # auto-selecting it would silently substitute a different model. It runs
  # only on an explicit engine = "swift" request.
  "legacy"
}

#' @keywords internal
#' @noRd
.resolve_searchlight_engine.naive_xdec_model <- function(model_spec, method, engine = "auto") {
  requested <- .match_searchlight_engine(engine)
  if (identical(requested, "legacy")) {
    return("legacy")
  }

  registry <- .searchlight_engine_registry()
  eligible <- isTRUE(registry$naive_xdec_fast$eligible(model_spec, method))

  if (identical(requested, "auto")) {
    return(if (eligible) "naive_xdec_fast" else "legacy")
  }

  if (identical(requested, "naive_xdec_fast") && eligible) {
    return("naive_xdec_fast")
  }

  "legacy"
}

#' @keywords internal
#' @noRd
.execute_searchlight_engine <- function(engine,
                                        model_spec,
                                        radius,
                                        method = "standard",
                                        niter = 4L,
                                        combiner = "average",
                                        drop_probs = FALSE,
                                        fail_fast = FALSE,
                                        backend = c("default", "shard", "auto"),
                                        incremental = TRUE,
                                        gamma = NULL,
                                        verbose = FALSE,
                                        ...) {
  if (identical(engine, "swift")) {
    res <- if (identical(method, "standard")) {
      run_searchlight_swift_fast(
        model_spec = model_spec,
        radius = radius,
        verbose = verbose,
        ...
      )
    } else {
      run_searchlight_swift_sampled_fast(
        model_spec = model_spec,
        radius = radius,
        method = method,
        niter = niter,
        combiner = combiner,
        drop_probs = drop_probs,
        fail_fast = fail_fast,
        backend = backend,
        verbose = verbose,
        ...
      )
    }
    attr(res, "searchlight_engine") <- "swift"
    attr(res, "searchlight_estimator") <- "swift_nearest_mean"
    return(res)
  }

  if (identical(engine, "sda_fast")) {
    res <- run_searchlight_sda_fast(model_spec = model_spec, radius = radius, verbose = verbose)
    attr(res, "searchlight_engine") <- "sda_fast"
    return(res)
  }

  if (identical(engine, "aggregate_fast")) {
    res <- run_searchlight_aggregate_fast(model_spec = model_spec, radius = radius,
                                          verbose = verbose)
    attr(res, "searchlight_engine") <- "aggregate_fast"
    return(res)
  }

  if (identical(engine, "dual_lda_fast")) {
    res <- if (identical(method, "standard")) {
      run_searchlight_dual_lda_fast(
        model_spec = model_spec,
        radius = radius,
        incremental = incremental,
        gamma = gamma,
        verbose = verbose
      )
    } else {
      run_searchlight_dual_lda_sampled_fast(
        model_spec = model_spec,
        radius = radius,
        method = method,
        niter = niter,
        combiner = combiner,
        drop_probs = drop_probs,
        fail_fast = fail_fast,
        backend = backend,
        gamma = gamma,
        verbose = verbose
      )
    }
    attr(res, "searchlight_engine") <- "dual_lda_fast"
    return(res)
  }

  # Defensive: every engine that .resolve_searchlight_engine can return for an
  # mvpa_model must have an executor branch above. A fall-through means a new
  # engine was registered without wiring it here; fail loudly rather than
  # returning NULL (which the caller would treat as a successful empty result).
  stop(
    sprintf("Internal error: no executor for searchlight engine '%s'.", engine),
    call. = FALSE
  )
}

#' @keywords internal
#' @noRd
.run_searchlight_engine.mvpa_model <- function(model_spec, radius, method,
                                               engine = "auto",
                                               niter = 4L,
                                               combiner = "average",
                                               drop_probs = FALSE,
                                               fail_fast = FALSE,
                                               backend = c("default", "shard", "auto"),
                                               incremental = TRUE,
                                               gamma = NULL,
                                               verbose = FALSE,
                                               ...) {
  requested <- .match_searchlight_engine(engine)
  engine <- .resolve_searchlight_engine.mvpa_model(
    model_spec = model_spec,
    method = method,
    engine = requested,
    combiner = combiner
  )

  strict_requested <- !identical(requested, "auto") && !identical(requested, "legacy")

  if (strict_requested && identical(engine, "legacy")) {
    stop(
      sprintf(
        "Requested searchlight engine '%s' is not eligible for method '%s'. Use engine='auto' to allow fallback.",
        requested,
        method
      ),
      call. = FALSE
    )
  }

  if (identical(engine, "swift")) {
    message(
      sprintf(
        paste0(
          "engine = 'swift' computes SWIFT's z-scored nearest-class-mean ",
          "estimator, not the requested '%s' classifier."
        ),
        model_spec$model$label %||% "model"
      )
    )
  }

  if (identical(engine, "legacy")) {
    if (!identical(requested, "legacy")) {
      futile.logger::flog.info(
        "searchlight engine: general-purpose iterator (no eligible fast path)"
      )
    }
    return(list(
      handled = FALSE,
      result = NULL,
      engine = "legacy"
    ))
  }

  fast_res <- try(
    .execute_searchlight_engine(
      engine = engine,
      model_spec = model_spec,
      radius = radius,
      method = method,
      niter = niter,
      combiner = combiner,
      drop_probs = drop_probs,
      fail_fast = fail_fast,
      backend = backend,
      incremental = incremental,
      gamma = gamma,
      verbose = verbose,
      ...
    ),
    silent = TRUE
  )

  if (!inherits(fast_res, "try-error")) {
    futile.logger::flog.info("searchlight engine: %s", engine)
    return(list(
      handled = TRUE,
      result = fast_res,
      engine = engine
    ))
  }

  err_cond <- attr(fast_res, "condition")
  err_msg <- tryCatch(
    conditionMessage(err_cond),
    error = function(...) as.character(fast_res)
  )

  # An engine may refuse data outside its proven regime (e.g. missing values).
  # That is expected behaviour, not a failure: fall back quietly under auto.
  if (inherits(err_cond, "rmvpa_engine_ineligible") && !strict_requested) {
    futile.logger::flog.info(
      "searchlight engine '%s' not applicable (%s); using the general-purpose iterator",
      engine, err_msg
    )
    return(list(handled = FALSE, result = NULL, engine = "legacy"))
  }

  if (strict_requested) {
    stop(
      sprintf(
        "Requested searchlight engine '%s' failed: %s",
        engine,
        err_msg
      ),
      call. = FALSE
    )
  }

  warning(
    sprintf(
      "searchlight engine '%s' failed, falling back to the general-purpose iterator: %s",
      engine,
      err_msg
    ),
    call. = FALSE
  )

  list(
    handled = FALSE,
    result = NULL,
    engine = "legacy"
  )
}
