#' @keywords internal
#' @noRd
.match_searchlight_engine <- function(engine = "auto") {
  allowed <- c("auto", "legacy", "swift", "dual_lda_fast", "naive_xdec_fast")
  match.arg(as.character(engine)[1], allowed)
}

#' @keywords internal
#' @noRd
.searchlight_engine_registry <- function() {
  list(
    legacy = list(
      label = "Legacy iterator",
      eligible = function(model_spec, method) TRUE
    ),
    swift = list(
      label = "SWIFT multiclass fast path",
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
    naive_xdec_fast = list(
      label = "Naive cross-decoding matrix fast path",
      eligible = function(model_spec, method) {
        .is_naive_xdec_fast_path(model_spec, method)
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
#' @param engine Requested engine policy: \code{"auto"}, \code{"legacy"},
#'   \code{"swift"}, \code{"dual_lda_fast"}, or \code{"naive_xdec_fast"}.
#'
#' @return A data frame with selection metadata.
#' @export
explain_searchlight_engine <- function(model_spec,
                                       method = c("standard", "randomized", "resampled"),
                                       engine = c("auto", "legacy", "swift", "dual_lda_fast", "naive_xdec_fast")) {
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
.resolve_searchlight_engine <- function(model_spec, method, engine = "auto") {
  if (inherits(model_spec, "naive_xdec_model")) {
    return(.resolve_searchlight_engine.naive_xdec_model(model_spec, method, engine))
  }
  if (inherits(model_spec, "mvpa_model")) {
    return(.resolve_searchlight_engine.mvpa_model(model_spec, method, engine))
  }
  "legacy"
}

#' @keywords internal
#' @noRd
.resolve_searchlight_engine.mvpa_model <- function(model_spec, method, engine = "auto") {
  requested <- .match_searchlight_engine(engine)
  registry <- .searchlight_engine_registry()

  if (identical(requested, "legacy")) {
    return("legacy")
  }

  if (!identical(requested, "auto")) {
    if (isTRUE(registry[[requested]]$eligible(model_spec, method))) {
      return(requested)
    }
    return("legacy")
  }

  if (isTRUE(registry$dual_lda_fast$eligible(model_spec, method))) {
    return("dual_lda_fast")
  }

  if (isTRUE(registry$swift$eligible(model_spec, method))) {
    return("swift")
  }

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
    engine = requested
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

  if (identical(engine, "legacy")) {
    if (!identical(requested, "legacy")) {
      futile.logger::flog.info("searchlight engine: legacy (no eligible fast path)")
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

  err_msg <- tryCatch(
    conditionMessage(attr(fast_res, "condition")),
    error = function(...) as.character(fast_res)
  )

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
      "searchlight engine '%s' failed, falling back to legacy: %s",
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
