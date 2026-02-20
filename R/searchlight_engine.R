#' @keywords internal
#' @noRd
.searchlight_fail_open_enabled <- function() {
  opt <- getOption("rMVPA.searchlight_fail_open", NULL)
  if (!is.null(opt)) {
    val <- suppressWarnings(as.logical(opt)[1])
    if (!is.na(val)) {
      return(isTRUE(val))
    }
  }
  TRUE
}

#' @keywords internal
#' @noRd
.match_searchlight_engine <- function(engine = "auto") {
  allowed <- c("auto", "legacy", "swift", "dual_lda_fast")
  match.arg(as.character(engine)[1], allowed)
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
  list(handled = FALSE, result = NULL, engine = "legacy")
}

#' @keywords internal
#' @noRd
.resolve_searchlight_engine.mvpa_model <- function(model_spec, method, engine = "auto") {
  requested <- .match_searchlight_engine(engine)

  if (identical(requested, "legacy")) {
    return("legacy")
  }

  if (identical(requested, "swift")) {
    if (.is_swift_fast_path(model_spec, method)) {
      return("swift")
    }
    return("legacy")
  }

  if (identical(requested, "dual_lda_fast")) {
    if (.is_dual_lda_fast_path(model_spec, method)) {
      return("dual_lda_fast")
    }
    return("legacy")
  }

  if (.is_dual_lda_fast_path(model_spec, method)) {
    return("dual_lda_fast")
  }

  if (.swift_searchlight_enabled() && .is_swift_fast_path(model_spec, method)) {
    return("swift")
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

  NULL
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
    return(list(handled = FALSE, result = NULL, engine = "legacy"))
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
    return(list(handled = TRUE, result = fast_res, engine = engine))
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

  if (!.searchlight_fail_open_enabled()) {
    stop(
      sprintf(
        "Searchlight engine '%s' failed and fail-open is disabled: %s",
        engine,
        err_msg
      ),
      call. = FALSE
    )
  }

  futile.logger::flog.warn(
    "searchlight engine '%s' failed; falling back to legacy iterator: %s",
    engine,
    err_msg
  )

  list(handled = FALSE, result = NULL, engine = "legacy")
}
