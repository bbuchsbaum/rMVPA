context("banded ridge documentation and certification artifacts")

.brdoc_source_root <- function() {
  installed_root <- system.file(package = "rMVPA")
  candidates <- c(
    testthat::test_path("..", ".."),
    if (nzchar(installed_root)) {
      file.path(dirname(installed_root), "00_pkg_src", "rMVPA")
    } else {
      character()
    }
  )
  candidates <- unique(candidates[nzchar(candidates)])
  valid <- candidates[file.exists(file.path(candidates, "DESCRIPTION"))]
  if (length(valid)) return(valid[[1L]])
  candidates[[1L]]
}

.brdoc_root_file <- function(...) {
  file.path(.brdoc_source_root(), ...)
}

test_that("source, Rd, NEWS, and pkgdown agree on the public lifecycle", {
  vignette_path <- .brdoc_root_file("vignettes", "Banded_Ridge_Encoding.Rmd")
  model_rd_path <- .brdoc_root_file("man", "banded_ridge_model.Rd")
  run_rd_path <- .brdoc_root_file("man", "run_banded_ridge.Rd")
  news_path <- .brdoc_root_file("NEWS.md")
  pkgdown_path <- .brdoc_root_file("_pkgdown.yml")
  pkgdown_workflow_path <- .brdoc_root_file(
    ".github", "workflows", "pkgdown.yaml"
  )
  source_artifacts <- c(
    vignette_path, model_rd_path, run_rd_path, news_path, pkgdown_path,
    pkgdown_workflow_path
  )
  skip_if_not(
    all(file.exists(source_artifacts)),
    "source-only NEWS/pkgdown consistency artifacts are excluded from the tarball"
  )

  vignette <- readLines(vignette_path, warn = FALSE)
  model_rd <- readLines(model_rd_path, warn = FALSE)
  run_rd <- readLines(run_rd_path, warn = FALSE)
  news <- readLines(news_path, warn = FALSE)
  pkgdown <- readLines(pkgdown_path, warn = FALSE)
  pkgdown_workflow <- readLines(pkgdown_workflow_path, warn = FALSE)
  combined <- paste(c(vignette, model_rd, run_rd, news), collapse = "\n")

  expect_match(combined, "banded_ridge_model")
  expect_match(combined, "run_banded_ridge")
  expect_match(combined, "delta_sets")
  expect_match(combined, "predictive_leave_one_band_out")
  expect_match(combined, "delta_cv_r2_<band>", fixed = TRUE)
  expect_match(combined, "not (an )?additive", ignore.case = TRUE)
  expect_false(any(grepl("banded_ridge_encoding_model", c(
    vignette, model_rd, run_rd, news
  ), fixed = TRUE)))
  expect_true(any(grepl("Banded_Ridge_Encoding", pkgdown, fixed = TRUE)))
  expect_true(any(grepl("Banded_Ridge_Encoding", news, fixed = TRUE)))
  expect_true(any(grepl(
    'RGL_USE_NULL: "TRUE"', pkgdown_workflow, fixed = TRUE
  )))

  vignette_title <- sub("^title:[[:space:]]*", "", vignette[
    grep("^title:", vignette)[[1L]]
  ])
  vignette_title <- gsub("^['\"]|['\"]$", "", vignette_title)
  index_line <- vignette[grep("VignetteIndexEntry", vignette)[[1L]]]
  index_title <- sub(".*VignetteIndexEntry\\{(.*)\\}.*", "\\1", index_line)
  expect_identical(vignette_title, index_title)
})

test_that("the executable guide exposes the complete user journey without internals", {
  vignette <- readLines(
    .brdoc_root_file("vignettes", "Banded_Ridge_Encoding.Rmd"),
    warn = FALSE
  )
  text <- paste(vignette, collapse = "\n")
  expected_calls <- c(
    "feature_sets(", "feature_sets_design(", "banded_ridge_model(",
    "run_banded_ridge("
  )
  for (call in expected_calls) expect_match(text, call, fixed = TRUE)
  expect_false(grepl("rMVPA:::", text, fixed = TRUE))
  expect_false(grepl("eval = FALSE", text, fixed = TRUE))
  expect_match(text, "inner_score")
  expect_match(text, "outer out-of-fold")
  expect_match(text, "max_retained_mb")
  expect_match(text, "memory_limit_mb")
  expect_match(text, "Changing theta changes the scaled design/kernel")
  expect_match(text, "did \\*\\*not\\*\\* reproduce the earlier 78% claim")
  expect_match(text, "Ordinary ridge")
  expect_match(text, "glmnet")
  expect_match(text, "multiridge")
  expect_match(text, "himalaya")
})

test_that("constructor defaults and generated usage are documented exactly", {
  defaults <- formals(banded_ridge_model)
  expect_identical(eval(defaults$theta_method)[[1L]], "grid")
  expect_null(defaults$delta_sets)
  expect_identical(eval(defaults$weight_retention)[[1L]], "none")
  expect_false(eval(defaults$return_predictions))
  expect_identical(eval(defaults$solver)[[1L]], "auto")

  exports <- getNamespaceExports("rMVPA")
  expect_true(all(c("banded_ridge_model", "run_banded_ridge") %in% exports))
  namespace <- readLines(.brdoc_root_file("NAMESPACE"), warn = FALSE)
  expect_true(any(namespace == "S3method(print,banded_ridge_result)"))
  expect_true(any(namespace == "S3method(run_regional,banded_ridge_model)"))
  expect_true(any(namespace == "S3method(run_searchlight,banded_ridge_model)"))

  model_rd <- paste(readLines(
    .brdoc_root_file("man", "banded_ridge_model.Rd"), warn = FALSE
  ), collapse = "\n")
  expect_match(model_rd, "delta_sets = NULL", fixed = TRUE)
  expect_match(model_rd, "theta_method = c(\"grid\", \"fixed\", \"random\")",
               fixed = TRUE)
  expect_match(model_rd, "weight_retention = c(\"none\", \"dual\", \"primal\")",
               fixed = TRUE)
})

test_that("the hosted dependency chain is declared explicitly", {
  description <- read.dcf(.brdoc_root_file("DESCRIPTION"))
  remotes <- trimws(strsplit(
    description[[1L, "Remotes"]], ",", fixed = TRUE
  )[[1L]])

  expect_true(all(c(
    "bbuchsbaum/fmridesign", "bbuchsbaum/fmrilss"
  ) %in% remotes))
  indirect <- c(
    "github::bbuchsbaum/fmriAR", "github::bbuchsbaum/fmrihrf"
  )
  for (field in c(
    "Config/Needs/check", "Config/Needs/coverage", "Config/Needs/website"
  )) {
    needs <- trimws(strsplit(
      description[[1L, field]], ",", fixed = TRUE
    )[[1L]])
    expect_true(all(indirect %in% needs), info = field)
  }
})

test_that("hosted artifact and Windows linkage gates are explicit", {
  workflow_paths <- c(
    .brdoc_root_file(".github", "workflows", "extended-tests.yaml"),
    .brdoc_root_file(".github", "workflows", "capability-tests.yaml")
  )
  makevars_path <- .brdoc_root_file("src", "Makevars")
  skip_if_not(
    all(file.exists(c(workflow_paths, makevars_path))),
    "source-only CI artifacts are excluded from the tarball"
  )

  workflows <- lapply(workflow_paths, readLines, warn = FALSE)
  expect_true(all(vapply(
    workflows,
    function(lines) any(grepl("local::.", lines, fixed = TRUE)),
    logical(1L)
  )))
  makevars <- paste(readLines(makevars_path, warn = FALSE), collapse = "\n")
  expect_match(makevars, "$(LAPACK_LIBS)", fixed = TRUE)
  expect_match(makevars, "$(BLAS_LIBS)", fixed = TRUE)
  expect_match(makevars, "$(FLIBS)", fixed = TRUE)

  buildignore_path <- .brdoc_root_file(".Rbuildignore")
  buildignore <- readLines(buildignore_path, warn = FALSE)
  expect_false(any(grepl(
    "src/Makevars",
    buildignore,
    fixed = TRUE
  )), info = "src/Makevars must be present in the source tarball")
})

test_that("issue 70 receipt declares seeds, folds, candidates, and uncertainty inputs", {
  results_path <- .brdoc_root_file(
    "inst", "extdata", "banded_ridge_issue70_results.csv"
  )
  script_path <- .brdoc_root_file(
    "inst", "benchmarks", "banded_ridge_issue70.R"
  )
  expect_true(file.exists(results_path))
  expect_true(file.exists(script_path))
  results <- utils::read.csv(results_path, stringsAsFactors = FALSE)

  expect_equal(nrow(results), 30L * 40L)
  expect_equal(length(unique(results$seed)), 30L)
  expect_identical(sort(unique(results$snr)), c(0, 0.4, 1))
  expect_true(all(is.finite(results$shared_oof_correlation)))
  expect_true(all(is.finite(results$response_oof_correlation)))
  expect_true(all(nzchar(results$fold_manifest)))
  expect_true(all(nzchar(results$candidate_manifest)))
  expect_true(all(grepl("test_rows=", results$fold_manifest, fixed = TRUE)))
  expect_true(all(grepl("alpha=", results$candidate_manifest, fixed = TRUE)))

  by_rep <- aggregate(
    absolute_correlation_gain ~ seed + snr, results, mean
  )
  moderate <- by_rep[by_rep$snr == 0.4, "absolute_correlation_gain"]
  expect_equal(mean(moderate), -0.047565622, tolerance = 1e-8)
  expect_equal(
    unname(stats::quantile(moderate, c(0.025, 0.975))),
    c(-0.06868121, -0.02957905), tolerance = 1e-8
  )
  noise_by_rep <- aggregate(
    response_oof_correlation ~ seed + snr, results[results$snr == 0, ], mean
  )
  expect_lt(abs(mean(noise_by_rep$response_oof_correlation)), 0.03)
})

test_that("performance receipt separates timing snapshots from allocation contracts", {
  results_path <- .brdoc_root_file(
    "inst", "extdata", "banded_ridge_performance_results.csv"
  )
  script <- readLines(.brdoc_root_file(
    "inst", "benchmarks", "banded_ridge_performance.R"
  ), warn = FALSE)
  results <- utils::read.csv(results_path, stringsAsFactors = FALSE)

  expect_equal(nrow(results), 6L)
  expect_identical(sort(unique(results$shape)), c("small", "wide"))
  expect_identical(sort(unique(results$solver)),
                   sort(c("direct", "svd_primal", "dual_kernel")))
  expect_true(all(results$elapsed_median_seconds > 0))
  expect_true(all(results$estimated_solver_intermediate_peak_bytes > 0))
  expect_true(all(results$max_chunk_result_bytes > 0))
  expect_lt(max(results$max_abs_prediction_error_vs_direct), 1e-10)
  expect_true(all(nzchar(results$r_version)))
  expect_true(all(nzchar(results$platform)))
  expect_true(all(nzchar(results$hardware_chip)))
  expect_true(any(grepl("not CI pass/fail thresholds", script, fixed = TRUE)))
  expect_true(any(grepl("allocation limits are enforced", script, fixed = TRUE)))
})
