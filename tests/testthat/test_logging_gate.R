test_that(".log_debug skips evaluating its arguments when debug is off", {
  old <- futile.logger::flog.threshold()
  on.exit({ set_log_level(old) }, add = TRUE)

  set_log_level("INFO")
  evaluated <- FALSE
  rMVPA:::.log_debug("%s", { evaluated <- TRUE; "msg" })
  expect_false(evaluated)
})

capture_log <- function(code) {
  msgs <- character()
  out <- utils::capture.output(withCallingHandlers(code, message = function(e) {
    msgs <<- c(msgs, conditionMessage(e))
    invokeRestart("muffleMessage")
  }))
  c(out, msgs)
}

test_that(".log_debug emits when debug is on and set_log_level refreshes the gate", {
  old <- futile.logger::flog.threshold()
  on.exit({ set_log_level(old) }, add = TRUE)

  set_log_level("DEBUG")
  expect_true(any(grepl("gate-check 42", capture_log(rMVPA:::.log_debug("gate-check %d", 42L)))))

  set_log_level("WARN")
  expect_false(any(grepl("gate-check 7", capture_log(rMVPA:::.log_debug("gate-check %d", 7L)))))
})

test_that("setup_mvpa_logger leaves the user's threshold alone", {
  old <- futile.logger::flog.threshold()
  on.exit({ set_log_level(old) }, add = TRUE)

  set_log_level("DEBUG")
  rMVPA:::setup_mvpa_logger()
  expect_identical(as.character(futile.logger::flog.threshold()), "DEBUG")
  expect_true(rMVPA:::.rmvpa_log_state$debug)

  set_log_level("ERROR")
  rMVPA:::setup_mvpa_logger()
  expect_identical(as.character(futile.logger::flog.threshold()), "ERROR")
  expect_false(rMVPA:::.rmvpa_log_state$debug)
})


test_that("the debug gate uses the effective rMVPA logger threshold", {
  old <- futile.logger::flog.threshold()
  old_named <- futile.logger::flog.threshold(name = "rMVPA")
  on.exit({
    futile.logger::flog.threshold(old)
    futile.logger::flog.threshold(old_named, name = "rMVPA")
    rMVPA:::.rmvpa_refresh_log_state()
  }, add = TRUE)
  futile.logger::flog.threshold(futile.logger::DEBUG)
  futile.logger::flog.threshold(futile.logger::ERROR, name = "rMVPA")
  rMVPA:::setup_mvpa_logger()
  evaluated <- FALSE
  rMVPA:::.log_debug("%s", { evaluated <- TRUE; "hidden" })
  expect_false(evaluated)
  expect_false(rMVPA:::.rmvpa_log_state$debug)
})
