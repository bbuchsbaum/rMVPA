# rmvpa DSL 0.1: schema, inheritance, explain, check, run, cache, status, lock,
# generators and the CLI. Uses the synthetic BIDS project in helper-dsl.R.

atlas_line <- "  regions: { atlas: atlas/space-MNI_atlas-two_dseg.nii.gz, labels: atlas/space-MNI_atlas-two_dseg.tsv }"

dsl_decode_project <- function(...) {
  dir <- dsl_fixture(...)
  dsl_write_project(dir)
  dsl_write_analysis(dir, "decode_category",
                     c("decode:", "  set: objects", "  target: category", "scope:", atlas_line))
  dir
}

# ---- schema and resolution -------------------------------------------------

test_that("unknown keys are rejected with a suggestion", {
  dsl_skip_if_missing()
  dir <- dsl_fixture(subjects = "01"); dsl_write_project(dir)
  dsl_write_analysis(dir, "bad", c("decode:", "  set: objects", "  targte: category", "scope:", atlas_line))
  expect_error(rmvpa_analysis(dir, "bad"), "unknown key 'targte'; did you mean 'target'")
})

test_that("an analysis needs exactly one verb, a scope and extends", {
  dsl_skip_if_missing()
  dir <- dsl_fixture(subjects = "01"); dsl_write_project(dir)
  writeLines(c("decode: { set: objects, target: category }", "scope:", atlas_line),
             file.path(dir, "analyses", "noext.yml"))
  expect_error(rmvpa_analysis(dir, "noext"), "missing required key 'extends'")
  dsl_write_analysis(dir, "two", c("decode: { set: objects, target: category }",
                                   "rsa: { item: stimulus, predictors: { c: rdms/category.csv } }",
                                   "scope:", atlas_line))
  expect_error(rmvpa_analysis(dir, "two"), "exactly one of")
})

test_that("set references are checked and inferred when the project has one set", {
  dsl_skip_if_missing()
  dir <- dsl_fixture(subjects = "01"); dsl_write_project(dir)
  dsl_write_analysis(dir, "noset", c("decode: { target: category }", "scope:", atlas_line))
  expect_error(rmvpa_analysis(dir, "noset"), "set is required because the project defines 3 sets")
  dsl_write_analysis(dir, "typo", c("decode: { set: objcts, target: category }", "scope:", atlas_line))
  expect_error(rmvpa_analysis(dir, "typo"), "did you mean 'objects'")

  dir1 <- dsl_fixture(subjects = "01"); dsl_write_project(dir1, extra_sets = FALSE)
  dsl_write_analysis(dir1, "noset", c("decode: { target: category }", "scope:", atlas_line))
  spec <- rmvpa_analysis(dir1, "noset")
  expect_identical(spec$config$decode$set, "objects")
  expect_match(spec$origins[["decode.set"]], "inferred")
})

test_that("values carry their origin and lists replace with a note", {
  dsl_skip_if_missing()
  dir <- dsl_fixture(subjects = "01"); dsl_write_project(dir)
  dsl_write_analysis(dir, "a", c("decode: { set: objects, target: category }",
                                 "folds: { assert: [one_run_per_fold] }", "run: { backend: ~ }",
                                 "scope:", atlas_line))
  spec <- rmvpa_analysis(dir, "a")
  expect_identical(spec$origins[["decode.method"]], "default@0.1")
  expect_identical(spec$origins[["run.seed"]], "project")
  expect_identical(spec$origins[["folds.assert"]], "analysis")
  expect_identical(spec$config$folds$assert, "one_run_per_fold")
  expect_null(spec$config$run$backend)
  notes <- unlist(spec$notes)
  expect_true(any(grepl("folds.assert: analysis replaced", notes)))
  expect_true(any(grepl("run.backend: removed by analysis", notes)))

  txt <- quiet_dsl(rmvpa_explain(dir, "a"))
  out <- utils::capture.output(rmvpa_explain(dir, "a"))
  expect_true(any(grepl("method: sda_notune +# default@0.1", out)))
  expect_true(any(grepl("assert: \\[one_run_per_fold\\]", out)))
})

test_that("explain output is itself a valid analysis file", {
  dsl_skip_if_missing()
  dir <- dsl_decode_project(subjects = "01")
  out <- utils::capture.output(rmvpa_explain(dir, "decode_category"))
  writeLines(out, file.path(dir, "analyses", "roundtrip.yml"))
  spec2 <- rmvpa_analysis(dir, "roundtrip")
  spec1 <- rmvpa_analysis(dir, "decode_category")
  expect_identical(spec2$config, spec1$config)
})

test_that("extends must point at the project file", {
  dsl_skip_if_missing()
  dir <- dsl_fixture(subjects = "01"); dsl_write_project(dir)
  writeLines(c("extends: ../other.yml", "decode: { set: objects, target: category }", "scope:", atlas_line),
             file.path(dir, "analyses", "x.yml"))
  expect_error(rmvpa_analysis(dir, "x"), "must extend its project file")
})

test_that("unsupported scopes and nulls are refused with the missing R piece", {
  dsl_skip_if_missing()
  dir <- dsl_fixture(subjects = "01"); dsl_write_project(dir)
  dsl_write_analysis(dir, "rsa_wb", c("rsa: { set: objects, item: stimulus, predictors: { c: rdms/category.csv } }",
                                      "scope:", "  whole_brain: {}"))
  expect_error(rmvpa_analysis(dir, "rsa_wb"), "whole_brain is not available for rsa")
  dsl_write_analysis(dir, "dec_perm", c("decode: { set: objects, target: category }", "scope:", atlas_line,
                                        "inference: { test: permutation }"))
  expect_error(rmvpa_analysis(dir, "dec_perm"), "run_permutation_regional")
  dsl_write_analysis(dir, "era_perm", c("rsa_era: { train_set: encoding, test_set: retrieval, link: stimulus }",
                                        "scope:", "  searchlight: { radius: 4 }", "inference: { test: permutation }"))
  expect_error(rmvpa_analysis(dir, "era_perm"), "permute_labels\\(\\) method that shuffles ERA")
})

test_that("searchlight decode defaults to the pooled permutation null", {
  dsl_skip_if_missing()
  dir <- dsl_fixture(subjects = "01"); dsl_write_project(dir)
  dsl_write_analysis(dir, "sl", c("decode: { set: objects, target: category }", "scope:", "  searchlight: { radius: 4 }"))
  spec <- rmvpa_analysis(dir, "sl")
  expect_identical(spec$config$inference$test, "permutation")
  expect_identical(spec$config$inference$permutations, 100L)
  expect_identical(spec$inference_support$kind, "searchlight_permutation")
  dsl_write_analysis(dir, "rsa_lm", c("rsa:", "  set: objects", "  item: stimulus",
                                      "  predictors: { c: rdms/category.csv }", "  nuisance: { n: rdms/noise.csv }",
                                      "  method: lm", "scope:", "  searchlight: { radius: 4 }"))
  expect_match(rmvpa_analysis(dir, "rsa_lm")$inference_support$description, "joint null")
})

# ---- check ---------------------------------------------------------------

test_that("check passes on a valid project and reads no images", {
  dsl_skip_if_missing()
  dir <- dsl_decode_project()
  chk <- rmvpa_check(dir, quiet = TRUE)
  expect_true(attr(chk, "ok"))
  expect_true(any(grepl("classes_in_every_fold", chk$message)))
})

test_that("check warns when a searchlight permutation test has fewer than 1000 permutations", {
  dsl_skip_if_missing()
  dir <- dsl_fixture(subjects = "01"); dsl_write_project(dir)
  dsl_write_analysis(dir, "sl", c("decode: { set: objects, target: category }", "scope:", "  searchlight: { radius: 4 }"))
  chk <- rmvpa_check(dir, quiet = TRUE)
  expect_true(attr(chk, "ok"))
  expect_true(any(chk$level == "warning" & grepl("permutations = 100: .*about 0.0099", chk$message)))
  dsl_write_analysis(dir, "sl", c("decode: { set: objects, target: category }", "scope:", "  searchlight: { radius: 4 }",
                                  "inference: { permutations: 1000 }"))
  expect_false(any(rmvpa_check(dir, quiet = TRUE)$level == "warning"))
})

test_that("check reports missing columns and failing assertions", {
  dsl_skip_if_missing()
  dir <- dsl_fixture(subjects = "01")
  dsl_write_project(dir, assert = "no_item_across_folds")
  dsl_write_analysis(dir, "rsa1", c("rsa: { set: objects, item: stimulus, predictors: { c: rdms/category.csv } }",
                                    "scope:", atlas_line))
  dsl_write_analysis(dir, "badcol", c("rsa: { set: objects, item: stimulis, predictors: { c: rdms/category.csv } }",
                                      "scope:", atlas_line))
  chk <- rmvpa_check(dir, quiet = TRUE)
  expect_false(attr(chk, "ok"))
  expect_true(any(grepl("assert no_item_across_folds: 8 item\\(s\\) appear in more than one run", chk$message)))
  expect_true(any(grepl("no column 'stimulis'; did you mean 'stimulus'", chk$message)))
})

test_that("an assertion that does not apply to the verb is an error", {
  dsl_skip_if_missing()
  dir <- dsl_fixture(subjects = "01"); dsl_write_project(dir, assert = "classes_in_every_fold")
  dsl_write_analysis(dir, "rsa1", c("rsa: { set: objects, item: stimulus, predictors: { c: rdms/category.csv } }",
                                    "scope:", atlas_line))
  chk <- rmvpa_check(dir, quiet = TRUE)
  expect_true(any(chk$level == "error" & grepl("does not apply to rsa", chk$message)))
})

test_that("keep filters are applied before folds and counted", {
  dsl_skip_if_missing()
  dir <- dsl_fixture(subjects = "01")
  dsl_write_project(dir, extra_sets = FALSE)
  txt <- readLines(file.path(dir, "rmvpa.yml"))
  txt <- sub("objects: \\{ task: objects, desc: lss \\}", "objects: { task: objects, keep: { response: hit } }", txt)
  writeLines(txt, file.path(dir, "rmvpa.yml"))
  dsl_write_analysis(dir, "d", c("decode: { target: category }", "scope:", atlas_line))
  cell <- quiet_dsl(rmvpa_compile(dir, "d", "01"))
  expect_identical(cell$train$n_rows, 32L)
  expect_identical(cell$train$n_kept, 24L)
  expect_true(all(cell$train$design$response == "hit"))
  expect_s3_class(cell$model_spec, "mvpa_model")
})

# ---- run, cache, status, lock ---------------------------------------------

test_that("run writes BIDS-named outputs with a sidecar and skips unchanged cells", {
  dsl_skip_if_missing()
  dir <- dsl_decode_project()
  res <- quiet_dsl(rmvpa_run(dir))
  expect_identical(res$status, c("ran", "ran"))
  out <- file.path(dir, "derivatives", "rmvpa", "decode_category", "sub-01")
  stem <- "sub-01_task-objects_space-MNI_desc-decodecategory"
  expect_true(file.exists(file.path(out, paste0(stem, "_stat-Accuracy_statmap.nii.gz"))))
  tsv <- utils::read.delim(file.path(out, paste0(stem, "_regions.tsv")))
  expect_identical(tsv$name, c("left", "right"))
  side <- jsonlite::fromJSON(file.path(out, paste0(stem, "_rmvpa.json")))
  expect_identical(side$inference$test, "none")
  expect_true(side$inference$descriptive_only)
  expect_identical(side$chance$Accuracy, 0.5)
  expect_identical(side$chance$AUC, 0L)
  expect_identical(side$lock, "none")
  expect_identical(side$model_class, "mvpa_model")

  again <- quiet_dsl(rmvpa_run(dir))
  expect_identical(again$status, c("up to date", "up to date"))
  st <- rmvpa_status(dir)
  expect_identical(st$status, c("done", "done"))

  # Changing the analysis makes both cells stale and re-runs them.
  dsl_write_analysis(dir, "decode_category",
                     c("decode:", "  set: objects", "  target: category", "  method: corclass", "scope:", atlas_line))
  expect_identical(rmvpa_status(dir)$status, c("stale", "stale"))
  expect_identical(quiet_dsl(rmvpa_run(dir, subjects = "01"))$status, "ran")
  expect_identical(rmvpa_status(dir)$status, c("done", "stale"))
})

test_that("lock is append-only and stamps match or mismatch", {
  dsl_skip_if_missing()
  dir <- dsl_decode_project(subjects = "01")
  quiet_dsl(rmvpa_lock(dir, note = "prereg"))
  expect_identical(rmvpa_status(dir)$lock, "match")
  quiet_dsl(rmvpa_run(dir))
  side <- jsonlite::fromJSON(list.files(file.path(dir, "derivatives"), "_rmvpa.json$", recursive = TRUE,
                                        full.names = TRUE))
  expect_identical(side$lock, "match")
  dsl_write_analysis(dir, "decode_category",
                     c("decode:", "  set: objects", "  target: category", "  method: corclass", "scope:", atlas_line))
  expect_identical(rmvpa_status(dir)$lock, "mismatch")
  quiet_dsl(rmvpa_lock(dir))
  expect_length(readLines(file.path(dir, "rmvpa.lock")), 2L)
  expect_identical(rmvpa_status(dir)$lock, "match")
})

test_that("run refuses to start when check fails", {
  dsl_skip_if_missing()
  dir <- dsl_fixture(subjects = "01"); dsl_write_project(dir)
  dsl_write_analysis(dir, "d", c("decode: { set: objects, target: nosuch }", "scope:", atlas_line))
  expect_error(quiet_dsl(rmvpa_run(dir)), "nothing was run")
  expect_false(dir.exists(file.path(dir, "derivatives", "rmvpa")))
})

test_that("draft runs one subject without permutations into a separate folder", {
  dsl_skip_if_missing()
  dir <- dsl_fixture(); dsl_write_project(dir)
  dsl_write_analysis(dir, "sl", c("decode: { set: objects, target: category, method: corclass }",
                                  "scope:", "  searchlight: { radius: 4 }"))
  res <- quiet_dsl(rmvpa_run(dir, draft = TRUE))
  expect_identical(res$subject, "01")
  side <- jsonlite::fromJSON(list.files(file.path(dir, "derivatives", "rmvpa-draft"), "_rmvpa.json$",
                                        recursive = TRUE, full.names = TRUE))
  expect_true(side$draft)
  expect_identical(side$inference$test, "none (draft)")
  expect_false(dir.exists(file.path(dir, "derivatives", "rmvpa")))
})

test_that("searchlight permutation runs write p maps", {
  dsl_skip_if_missing()
  skip_on_cran()
  dir <- dsl_fixture(subjects = "01"); dsl_write_project(dir)
  dsl_write_analysis(dir, "sl", c("decode: { set: objects, target: category, method: corclass }",
                                  "scope:", "  searchlight: { radius: 4 }",
                                  "inference: { permutations: 2, subsample: 0.5 }"))
  quiet_dsl(rmvpa_run(dir))
  files <- list.files(file.path(dir, "derivatives", "rmvpa", "sl", "sub-01"))
  expect_true(any(grepl("_stat-Accuracyp_statmap", files)))
  side <- jsonlite::fromJSON(file.path(dir, "derivatives", "rmvpa", "sl", "sub-01",
                                       grep("_rmvpa.json$", files, value = TRUE)))
  expect_equal(side$inference$p_floor, 1 / 3)
  expect_match(side$inference$null, "pooled searchlight permutation")
})

# ---- verbs compile to the expected constructors ---------------------------

test_that("verbs compile to ordinary model_specs", {
  dsl_skip_if_missing()
  dir <- dsl_fixture(subjects = "01"); dsl_write_project(dir)
  stims <- c("cat1", "cat2", "dog1", "dog2", "car1", "car2", "bus1", "bus2")
  dir.create(file.path(dir, "features"))
  set.seed(1)
  utils::write.csv(data.frame(item = stims, f1 = stats::rnorm(8), f2 = stats::rnorm(8)),
                   file.path(dir, "features", "vis.csv"), row.names = FALSE)
  cases <- list(
    decode = list(c("decode: { set: objects, target: category }", "scope:", atlas_line), "mvpa_model"),
    transfer = list(c("transfer: { train_set: encoding, test_set: retrieval, target: category }", "scope:", atlas_line),
                    "naive_xdec_model"),
    rsa = list(c("rsa: { set: objects, item: stimulus, predictors: { c: rdms/category.csv } }", "scope:", atlas_line),
               "rsa_model"),
    rsa_vector = list(c("rsa_vector: { set: objects, item: stimulus, distances: rdms/category.csv }", "scope:", atlas_line),
                      "vector_rsa_model"),
    rsa_era = list(c("rsa_era: { train_set: encoding, test_set: retrieval, link: stimulus }", "scope:", atlas_line),
                   "era_rsa_model"),
    rsa_connect = list(c("rsa_connect: { set: objects, item: stimulus, seed_region: rdms/category.csv }", "scope:", atlas_line),
                       "repnet_model"),
    rsa_mediate = list(c("rsa_mediate: { set: objects, item: stimulus, predictor: rdms/category.csv, outcome: rdms/noise.csv }",
                         "scope:", atlas_line), "repmed_model"),
    encode = list(c("encode: { set: objects, item: stimulus, feature_spaces: { vis: features/vis.csv } }",
                    "scope:", "  whole_brain: {}"), "banded_ridge_model"),
    pattern = list(c("pattern: { set: objects, target: category, max_rank: 1 }", "scope:", "  whole_brain: {}"),
                   "pattern_model")
  )
  for (v in names(cases)) {
    dsl_write_analysis(dir, v, cases[[v]][[1]])
    cell <- quiet_dsl(rmvpa_compile(dir, v, "01"))
    expect_s3_class(cell$model_spec, cases[[v]][[2]])
  }
})

test_that("rsa expands item RDMs to trials and flags the joint null", {
  dsl_skip_if_missing()
  dir <- dsl_fixture(subjects = "01"); dsl_write_project(dir)
  dsl_write_analysis(dir, "r", c("rsa:", "  set: objects", "  item: stimulus",
                                 "  predictors: { category: rdms/category.csv }", "  nuisance: { noise: rdms/noise.csv }",
                                 "  method: lm", "scope:", atlas_line))
  cell <- quiet_dsl(rmvpa_compile(dir, "r", "01"))
  expect_true(cell$rsa_joint)
  expect_identical(cell$model_spec$regtype, "lm")
  expect_match(cell$info$rdm_level, "trial")
})

test_that("custom builders must return a model_spec", {
  dsl_skip_if_missing()
  dir <- dsl_fixture(subjects = "01"); dsl_write_project(dir)
  writeLines("build_spec <- function(dataset, design, config) list(a = 1)", file.path(dir, "bad.R"))
  dsl_write_analysis(dir, "c", c("custom: { set: objects, r_file: bad.R }", "scope:", atlas_line))
  expect_error(quiet_dsl(rmvpa_compile(dir, "c", "01")), "must return a model_spec")
})

# ---- generators and CLI ----------------------------------------------------

test_that("new scaffolds a project that checks cleanly", {
  dsl_skip_if_missing()
  src <- dsl_fixture(subjects = "01")
  proj <- file.path(src, "proj")
  expect_error(rmvpa_new(file.path(src, "derivatives", "lss"), dir = proj), "Choose a target column")
  quiet_dsl(rmvpa_new(file.path(src, "derivatives", "lss"), dir = proj, target = "category"))
  expect_true(file.exists(file.path(proj, "analyses", "decode_category.yml")))
  p <- rmvpa_project(proj)
  expect_setequal(names(p$config$data$sets), c("memory", "objects"))
  spec <- rmvpa_analysis(p, "decode_category")
  expect_identical(spec$config$decode$target, "category")
})

test_that("add writes a template for each verb", {
  dsl_skip_if_missing()
  dir <- dsl_fixture(subjects = "01"); dsl_write_project(dir)
  for (v in .dsl_verb_names()) {
    f <- quiet_dsl(rmvpa_add(dir, v, paste0("t_", v)))
    y <- yaml::read_yaml(f)
    expect_true(v %in% names(y), info = v)
    expect_identical(y$extends, "../rmvpa.yml")
  }
  expect_error(rmvpa_add(dir, "decod", "x"), "did you mean 'decode'")
})

test_that("migrate converts an old single-set CLI config", {
  dsl_skip_if_missing()
  src <- dsl_fixture(subjects = "01")
  f <- file.path(src, "derivatives", "lss", "sub-01", "func")
  old <- file.path(src, "old.yaml")
  writeLines(c(sprintf("train_design: %s/sub-01_task-objects_run-1_space-MNI_desc-lss_betas.tsv", f),
               sprintf("train_data: %s/sub-01_task-objects_run-1_space-MNI_desc-lss_betas.nii.gz", f),
               sprintf("mask: %s/sub-01_space-MNI_desc-brain_mask.nii.gz", f),
               "label_column: category", "block_column: stimulus", "model: corclass", "radius: 4"), old)
  proj <- file.path(src, "migrated")
  quiet_dsl(rmvpa_migrate(old, dir = proj))
  spec <- rmvpa_analysis(proj, "migrated")
  expect_identical(spec$verb, "decode")
  expect_identical(spec$config$decode$method, "corclass")
  expect_identical(spec$config$scope$searchlight$radius, 4L)
})

test_that("the rmvpa CLI dispatches commands and reports errors as status 1", {
  dsl_skip_if_missing()
  dir <- dsl_decode_project(subjects = "01")
  expect_output(st <- rmvpa_cli(c("help")), "usage: rmvpa")
  expect_identical(st, 0L)
  expect_output(rmvpa_cli(c("explain", "decode_category", "--project", dir)), "default@0.1")
  expect_output(rmvpa_cli(c("verbs")), "rsa_era")
  expect_message(st <- rmvpa_cli(c("explain", "nope", "--project", dir)), "No analysis 'nope'")
  expect_identical(st, 1L)
  expect_message(st <- rmvpa_cli(c("rn")), "did you mean 'run'")
  expect_identical(st, 1L)
  quiet_dsl(st <- rmvpa_cli(c("run", "--project", dir, "--yes")))
  expect_identical(st, 0L)
  expect_output(rmvpa_cli(c("status", "--project", dir)), "done")
})
