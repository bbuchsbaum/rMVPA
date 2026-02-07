#' run_searchlight method for remap_rrr_model
#'
#' Defaults to a lean, memory-friendly configuration:
#' - drop_probs = TRUE strips per-ROI probability matrices after metrics.
#' - return_pobserved = FALSE skips voxel×trial probability maps.
#' - combiner = "average" (combine_randomized) for randomized/resampled.
#'
#' @param model_spec A \code{remap_rrr_model} object.
#' @inheritParams run_searchlight_base
#' @export
run_searchlight.remap_rrr_model <- function(model_spec,
                                            radius = 8,
                                            method = c("randomized", "standard", "resampled"),
                                            niter = 4,
                                            combiner = "average",
                                            drop_probs = TRUE,
                                            return_pobserved = FALSE,
                                            fail_fast = FALSE,
                                            ...) {
  method <- match.arg(method)

  if (method == "standard") {
    do_standard(model_spec, radius,
                combiner = combine_standard,
                drop_probs = drop_probs,
                fail_fast = fail_fast,
                ...)
  } else if (method == "randomized") {
    # Force lean defaults
    do_randomized(model_spec, radius,
                  niter        = niter,
                  combiner     = combine_randomized,
                  drop_probs   = drop_probs,
                  return_pobserved = return_pobserved,
                  fail_fast    = fail_fast,
                  ...)
  } else { # resampled
    do_resampled(model_spec, radius,
                 niter        = niter,
                 combiner     = combine_randomized,
                 drop_probs   = drop_probs,
                 return_pobserved = return_pobserved,
                 fail_fast    = fail_fast,
                 ...)
  }
}
