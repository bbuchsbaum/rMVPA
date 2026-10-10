test_that(".str_squish matches stringr::str_squish", {
  skip_if_not_installed("stringr")
  x <- c("  a  b   c  ", "a\tb\nc", "", "   ", "single", " lead", "trail ", "a   +   b")
  expect_identical(rMVPA:::.str_squish(x), stringr::str_squish(x))
})

test_that(".str_split_keep matches stringr::str_split on a single string", {
  skip_if_not_installed("stringr")
  pat <- "\\s*\\+\\s*"
  cases <- c("a+b", "a+", "+a", "a", "", "a + b +  c", "a++b", "a+b+", "+")
  for (s in cases) {
    expected <- stringr::str_split(s, pattern = pat, simplify = TRUE)[1, ]
    expect_identical(rMVPA:::.str_split_keep(s, pat), as.character(expected),
                     info = paste0("input: '", s, "'"))
  }
})

test_that(".str_split_keep keeps a trailing empty piece", {
  expect_identical(rMVPA:::.str_split_keep("a+", "\\+"), c("a", ""))
  expect_identical(rMVPA:::.str_split_keep("a+b", "\\+"), c("a", "b"))
})

test_that("parse_spec_to_metadata parses factor definitions with spaces and '+'", {
  labels <- c("a", "b", "c")
  spec <- ~ f(a + b ~ c)
  res <- rMVPA:::parse_spec_to_metadata(labels, spec)
  expect_true(is.list(res) || is.data.frame(res))
  expect_error(rMVPA:::parse_spec_to_metadata(labels, ~ f(a ~ z)),
               "not found in provided labels")
})

test_that("initialize_configuration detects .R config files by extension", {
  tmp <- tempfile(fileext = ".R")
  on.exit(unlink(tmp), add = TRUE)
  writeLines("cfg_marker <- 42", tmp)
  cfg <- rMVPA:::initialize_configuration(list(config = tmp))
  expect_identical(cfg$cfg_marker, 42)
})
