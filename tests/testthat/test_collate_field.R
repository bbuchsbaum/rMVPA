context("DESCRIPTION Collate field")

# The Collate field is maintained by hand: no source file uses an @include tag,
# so roxygen2::update_collate() is a no-op and will not catch drift. A new file
# under R/ that never makes it into Collate breaks R CMD INSTALL for everyone
# installing from source (see issue #64).
test_that("Collate lists every file under R/", {
  desc_path <- testthat::test_path("..", "..", "DESCRIPTION")
  r_dir <- testthat::test_path("..", "..", "R")
  source_files <- if (dir.exists(r_dir)) {
    list.files(r_dir, pattern = "\\.[RrSsQq]$")
  } else {
    character()
  }

  # Skip when running against an installed package rather than a source tree.
  skip_if_not(file.exists(desc_path) && length(source_files) > 0L,
              "not running from the package source tree")

  collate <- read.dcf(desc_path, fields = "Collate")[1, 1]
  skip_if(is.na(collate), "package has no Collate field")

  collated <- regmatches(collate, gregexpr("'[^']+'", collate))[[1]]
  collated <- gsub("'", "", collated, fixed = TRUE)

  expect_setequal(collated, source_files)
})
