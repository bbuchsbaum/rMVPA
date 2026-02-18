has_hrfdecoder_backend <- function() {
  requireNamespace("hrfdecoder", quietly = TRUE) ||
    requireNamespace("hrfdecode", quietly = TRUE)
}

skip_if_no_hrfdecoder_backend <- function() {
  testthat::skip_if_not(
    has_hrfdecoder_backend(),
    "'hrfdecoder' (or legacy 'hrfdecode') is not installed"
  )
}

skip_if_has_hrfdecoder_backend <- function() {
  testthat::skip_if(
    has_hrfdecoder_backend(),
    "'hrfdecoder'/'hrfdecode' installed; skipping negative path test"
  )
}
