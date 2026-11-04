# The studio's data layer, server guards and static export.

studio_project <- function() {
  dir <- dsl_fixture(subjects = c("01", "02", "03", "04")); dsl_write_project(dir)
  dsl_write_analysis(dir, "sl", c("decode: { set: objects, target: category, method: corclass }",
                                  "scope:", "  searchlight: { radius: 4 }",
                                  "inference: { permutations: 3, subsample: 0.5 }"))
  dsl_write_analysis(dir, "reg", c("decode: { set: objects, target: category, method: corclass }",
                                   "scope:", "  regions:",
                                   "    atlas: atlas/space-MNI_atlas-two_dseg.nii.gz",
                                   "    labels: atlas/space-MNI_atlas-two_dseg.tsv"))
  quiet_dsl(rmvpa_run(dir))
  dir
}

skip_without_app <- function() {
  testthat::skip_if(!nzchar(system.file("app", "index.html", package = "rmvpadsl")), "studio app not built")
}

test_that("permutation runs keep the pooled null as quantiles", {
  dsl_skip_if_missing(); skip_on_cran()
  dir <- studio_project()
  f <- list.files(file.path(dir, "derivatives", "rmvpa", "sl", "sub-01"), "_stat-Accuracy_null\\.tsv$", full.names = TRUE)
  expect_length(f, 1L)
  q <- utils::read.delim(f)
  expect_equal(nrow(q), 1001L)
  expect_false(is.unsorted(q$value))
  expect_true(q$n[1] > 0)
})

test_that("the studio payloads summarise outputs on disk", {
  dsl_skip_if_missing(); skip_on_cran()
  dir <- studio_project()
  src <- rmvpadsl:::.studio_source(dir)
  p <- src$project()
  expect_setequal(vapply(p$analyses, `[[`, "", "name"), c("reg", "sl"))
  expect_true(all(unlist(lapply(p$analyses, function(a) vapply(a$cells, `[[`, "", "status"))) == "done"))
  expect_equal(p$sets[[1]]$runs, 4L)
  expect_true(any(vapply(p$analyses[[2]]$explain, `[[`, "", "path") == "scope.searchlight.radius"))

  tp <- src$template()
  expect_equal(tp$dims, c(6L, 6L, 6L))
  expect_length(jsonlite::base64_dec(tp$template), 216L)

  r <- src$results("sl")
  expect_equal(r$metric, "Accuracy")
  expect_equal(r$chance, 0.5)
  expect_equal(r$null$permutations, 3L)
  expect_equal(r$null$p_resolution, 0.25)
  expect_length(r$subjects, 4L)
  s1 <- r$subjects[[1]]
  expect_length(s1$null_hist, length(r$hist_edges) - 1L)
  expect_true(sum(s1$null_hist) > 0)
  expect_equal(sum(s1$centres_hist), 216L)
  expect_true(all(c("group_mean", "prevalence", "prevalence_all") %in% names(r$maps)))
  expect_length(jsonlite::base64_dec(r$maps$group_mean), 4L * 216L)

  rr <- src$results("reg")
  expect_equal(vapply(rr$regions, `[[`, "", "name"), c("left", "right"))
  expect_length(rr$regions[[1]]$values, 4L)
  expect_equal(rr$chance_source, "1 / 2 classes")
  expect_null(rr$atlas$surface)

  expect_error(src$results("nope"), "No analysis")
  expect_error(src$surface("fsaverage5"), "no surface")
})

test_that("the server handler checks Host, Origin and the session token", {
  dsl_skip_if_missing(); skip_on_cran(); skip_without_app()
  dir <- studio_project()
  h <- rmvpadsl:::.studio_handler(rmvpadsl:::.studio_source(dir), system.file("app", package = "rmvpadsl"),
                                  token = "abc123", port = 7781L)
  req <- function(path, query = "", host = "127.0.0.1:7781", origin = NULL, token = NULL, method = "GET") {
    r <- list(PATH_INFO = path, QUERY_STRING = query, HTTP_HOST = host, REQUEST_METHOD = method)
    if (!is.null(origin)) r$HTTP_ORIGIN <- origin
    if (!is.null(token)) r$HTTP_X_RMVPA_TOKEN <- token
    h(r)
  }
  expect_equal(req("/", "?token=abc123")$status, 200L)
  expect_match(req("/", "?token=abc123")$body, 'name="rmvpa-token" content="abc123"', fixed = TRUE)
  expect_equal(req("/")$status, 403L)
  expect_equal(req("/", "?token=abc124")$status, 403L)
  expect_equal(req("/", "?token=abc123", host = "evil.example:7781")$status, 403L)   # DNS rebinding
  expect_equal(req("/api/project")$status, 403L)                                     # no token header
  expect_equal(req("/api/project", token = "abc123", origin = "https://evil.example")$status, 403L)
  expect_equal(req("/api/project", token = "abc123", method = "POST")$status, 405L)
  ok <- req("/api/project", token = "abc123", origin = "http://127.0.0.1:7781")
  expect_equal(ok$status, 200L)
  expect_equal(jsonlite::fromJSON(ok$body)$name, basename(dir))
  expect_equal(req("/api/results", "?analysis=sl", token = "abc123")$status, 200L)
  expect_equal(req("/api/results", "?analysis=nope", token = "abc123")$status, 500L)
  expect_equal(req("/api/other", token = "abc123")$status, 404L)
  expect_equal(req("/assets/studio.js")$status, 200L)
  expect_equal(req("/../DESCRIPTION")$status, 404L)
  expect_equal(req("/assets/../../DESCRIPTION")$status, 404L)
  expect_identical(ok$headers$`Content-Security-Policy`, "frame-ancestors 'none'")
})

test_that("rmvpa_studio serves the app on 127.0.0.1", {
  dsl_skip_if_missing(); skip_on_cran(); skip_without_app()
  testthat::skip_if_not_installed("httpuv")
  testthat::skip_if_not_installed("curl")
  dir <- studio_project()
  st <- suppressMessages(rmvpa_studio(dir, launch = FALSE, block = FALSE))
  on.exit(rmvpa_studio_stop(st))
  expect_match(st$url, "^http://127\\.0\\.0\\.1:[0-9]+/\\?token=[0-9a-f]{48}$")
  # Read the page and the API from a background R process while this one serves.
  token <- sub("^.*token=", "", st$url)
  out <- tempfile()
  script <- tempfile(fileext = ".R")
  writeLines(c(
    sprintf("u <- '%s'", st$url),
    sprintf("api <- 'http://127.0.0.1:%d/api/project'", st$port),
    "h <- curl::new_handle(); curl::handle_setheaders(h, 'X-Rmvpa-Token' = sub('^.*token=', '', u))",
    "page <- rawToChar(curl::curl_fetch_memory(u)$content)",
    "proj <- curl::curl_fetch_memory(api, handle = h)",
    "bad <- curl::curl_fetch_memory(api)",
    sprintf("saveRDS(list(page = page, status = proj$status_code, body = rawToChar(proj$content), bad = bad$status_code), '%s')", out)
  ), script)
  px <- system2(file.path(R.home("bin"), "Rscript"), script, wait = FALSE, stdout = FALSE, stderr = FALSE)
  t0 <- Sys.time()
  while (!file.exists(out) && Sys.time() - t0 < 30) httpuv::service(100)
  Sys.sleep(0.2)
  res <- readRDS(out)
  expect_match(res$page, "rmvpa-token")
  expect_equal(res$status, 200L)
  expect_equal(jsonlite::fromJSON(res$body)$name, basename(dir))
  expect_equal(res$bad, 403L)
})

test_that("a static export opens without a server", {
  dsl_skip_if_missing(); skip_on_cran(); skip_without_app()
  dir <- studio_project()
  out <- file.path(tempdir(), "studio-export")
  idx <- suppressMessages(rmvpa_export_studio(dir, out, overwrite = TRUE))
  expect_true(file.exists(idx))
  expect_setequal(list.files(file.path(out, "data")), c("project.js", "template.js", "results_sl.js", "results_reg.js"))
  html <- readLines(idx)
  data_tags <- grep('<script src="data/', html)
  app_tag <- grep("<script defer", html)
  expect_length(data_tags, 4L)
  expect_true(all(data_tags < app_tag))
  js <- readLines(file.path(out, "data", "results_sl.js"))
  expect_match(js, '^\\(window.__RMVPA_DATA__ = window.__RMVPA_DATA__ \\|\\| \\{\\}\\)\\["results/sl"\\] = \\{')
  expect_error(rmvpa_export_studio(dir, out), "not empty")
})
