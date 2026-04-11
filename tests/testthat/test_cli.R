cli_wrapper_path <- function(command) {
  wrapper_name <- unname(c(
    searchlight = "rmvpa-searchlight",
    regional = "rmvpa-regional"
  )[[command]])
  path <- system.file("exec", wrapper_name, package = "rMVPA")
  if (!nzchar(path)) {
    stop("Could not locate packaged CLI wrapper for ", command, call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

test_that("cli_finalize_config applies stable defaults and legacy aliases", {
  cfg <- structure(
    list(
      mode = "searchlight",
      ncores = 3L
    ),
    class = c("rmvpa_config", "list")
  )

  out <- rMVPA:::.cli_finalize_config(cfg, "searchlight")

  expect_equal(out$output, "searchlight_results")
  expect_equal(out$model, "corsim")
  expect_equal(out$model_type, "classification")
  expect_equal(out$data_mode, "image")
  expect_equal(out$label_column, "labels")
  expect_equal(out$class_metrics, TRUE)
  expect_equal(out$pthreads, 3L)
  expect_equal(out$radius, 8)
  expect_equal(out$type, "standard")
  expect_equal(out$niter, 4L)
})

test_that("cli_prepare_output enforces overwrite and skip semantics", {
  out_dir <- tempfile("rmvpa-cli-out-")
  dir.create(out_dir)
  writeLines("existing", file.path(out_dir, "marker.txt"))

  expect_error(
    rMVPA:::.cli_prepare_output(out_dir, overwrite = FALSE, skip_if_exists = FALSE),
    "already exists and is not empty"
  )
  expect_identical(
    rMVPA:::.cli_prepare_output(out_dir, overwrite = FALSE, skip_if_exists = TRUE),
    "skip"
  )
  expect_identical(
    rMVPA:::.cli_prepare_output(out_dir, overwrite = TRUE, skip_if_exists = FALSE),
    "run"
  )
  expect_error(
    rMVPA:::.cli_prepare_output(out_dir, overwrite = TRUE, skip_if_exists = TRUE),
    "either overwrite or skip_if_exists"
  )
})

test_that("install_cli copies packaged wrappers", {
  dest_dir <- tempfile("rmvpa-cli-bin-")
  dir.create(dest_dir)

  paths <- rMVPA::install_cli(dest_dir = dest_dir, overwrite = TRUE)

  expect_setequal(names(paths), c("searchlight", "regional"))
  expect_true(all(file.exists(paths)))
  expect_true(all(file.access(paths, mode = 1) == 0))
})

test_that("cli parser accepts canonical and legacy option spellings", {
  parser <- rMVPA:::.cli_parser("searchlight")
  parsed_new <- optparse::parse_args(
    parser,
    args = c("--train-design", "design.tsv", "--train-data", "betas.nii.gz", "--mask", "mask.nii.gz"),
    positional_arguments = TRUE
  )
  parsed_old <- optparse::parse_args(
    parser,
    args = rMVPA:::.cli_normalize_argv(
      c("--train_design", "design.tsv", "--train_data", "betas.nii.gz", "--mask", "mask.nii.gz")
    ),
    positional_arguments = TRUE
  )
  parsed_new$options <- rMVPA:::.cli_standardize_options(parsed_new$options)
  parsed_old$options <- rMVPA:::.cli_standardize_options(parsed_old$options)

  expect_equal(parsed_new$options$train_design, "design.tsv")
  expect_equal(parsed_old$options$train_design, "design.tsv")
  expect_equal(parsed_new$options$train_data, "betas.nii.gz")
  expect_equal(parsed_old$options$train_data, "betas.nii.gz")
})

test_that("installed-style wrappers support help and version", {
  wrappers <- c(searchlight = "searchlight", regional = "regional")
  rscript <- file.path(R.home("bin"), "Rscript")

  for (wrapper_name in names(wrappers)) {
    wrapper <- cli_wrapper_path(wrappers[[wrapper_name]])

    help <- system2(rscript, c(wrapper, "--help"), stdout = TRUE, stderr = TRUE)
    help_status <- attr(help, "status")
    if (is.null(help_status)) help_status <- 0L
    expect_equal(help_status, 0L, info = wrapper_name)
    expect_true(any(grepl("Usage", help, fixed = TRUE)), info = wrapper_name)

    version <- system2(rscript, c(wrapper, "--version"), stdout = TRUE, stderr = TRUE)
    version_status <- attr(version, "status")
    if (is.null(version_status)) version_status <- 0L
    expect_equal(version_status, 0L, info = wrapper_name)
    expect_true(
      any(grepl(as.character(utils::packageVersion("rMVPA")), version, fixed = TRUE)),
      info = wrapper_name
    )
  }
})

test_that("run_cli_command supports dry-run searchlight workflow", {
  skip_if_not_installed("neuroim2")

  dset <- gen_sample_dataset(c(4, 4, 4), nobs = 18, blocks = 3, nlevels = 2)
  train_path <- tempfile(fileext = ".nii.gz")
  mask_path <- tempfile(fileext = ".nii.gz")
  design_path <- tempfile(fileext = ".tsv")

  neuroim2::write_vec(dset$dataset$train_data, train_path)
  neuroim2::write_vol(dset$dataset$mask, mask_path)
  utils::write.table(
    dset$design$train_design,
    file = design_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  expect_no_error(
    rMVPA:::run_cli_command(
      "searchlight",
      argv = c(
        "--train-design", design_path,
        "--train-data", train_path,
        "--mask", mask_path,
        "--label-column", "Y",
        "--radius", "2",
        "--model", "sda_notune",
        "--workers", "1",
        "--dry-run"
      )
    )
  )
})

test_that("run_cli_command supports dry-run regional workflow", {
  skip_if_not_installed("neuroim2")

  dset <- gen_sample_dataset(c(4, 4, 4), nobs = 18, blocks = 3, nlevels = 2)
  train_path <- tempfile(fileext = ".nii.gz")
  region_mask_path <- tempfile(fileext = ".nii.gz")
  design_path <- tempfile(fileext = ".tsv")

  region_mask <- neuroim2::NeuroVol(
    sample(1:3, size = length(dset$dataset$mask), replace = TRUE),
    neuroim2::space(dset$dataset$mask)
  )

  neuroim2::write_vec(dset$dataset$train_data, train_path)
  neuroim2::write_vol(region_mask, region_mask_path)
  utils::write.table(
    dset$design$train_design,
    file = design_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  expect_no_error(
    rMVPA:::run_cli_command(
      "regional",
      argv = c(
        "--train-design", design_path,
        "--train-data", train_path,
        "--mask", region_mask_path,
        "--label-column", "Y",
        "--model", "sda_notune",
        "--workers", "1",
        "--save-predictors",
        "--pool-predictions", "mean",
        "--coalesce-design-vars",
        "--dry-run"
      )
    )
  )
})

test_that("installed-style regional wrapper supports dry-run workflow", {
  skip_if_not_installed("neuroim2")

  dset <- gen_sample_dataset(c(4, 4, 4), nobs = 18, blocks = 3, nlevels = 2)
  train_path <- tempfile(fileext = ".nii.gz")
  region_mask_path <- tempfile(fileext = ".nii.gz")
  design_path <- tempfile(fileext = ".tsv")
  wrapper <- cli_wrapper_path("regional")
  rscript <- file.path(R.home("bin"), "Rscript")

  region_mask <- neuroim2::NeuroVol(
    sample(1:3, size = length(dset$dataset$mask), replace = TRUE),
    neuroim2::space(dset$dataset$mask)
  )

  neuroim2::write_vec(dset$dataset$train_data, train_path)
  neuroim2::write_vol(region_mask, region_mask_path)
  utils::write.table(
    dset$design$train_design,
    file = design_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  out <- system2(
    rscript,
    c(
      wrapper,
      "--train-design", design_path,
      "--train-data", train_path,
      "--mask", region_mask_path,
      "--label-column", "Y",
      "--model", "sda_notune",
      "--workers", "1",
      "--preflight", "off",
      "--dry-run"
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L

  expect_equal(status, 0L)
  expect_true(any(grepl("Dry run summary", out, fixed = TRUE)))
  expect_true(any(grepl("mode: regional", out, fixed = TRUE)))
})

test_that("run_cli_command executes a tiny regional workflow and writes results", {
  skip_if_not_installed("neuroim2")

  dset <- gen_sample_dataset(c(4, 4, 4), nobs = 18, blocks = 3, nlevels = 2)
  train_path <- tempfile(fileext = ".nii.gz")
  region_mask_path <- tempfile(fileext = ".nii.gz")
  design_path <- tempfile(fileext = ".tsv")
  output_dir <- tempfile("rmvpa-cli-regional-")

  region_mask <- neuroim2::NeuroVol(
    sample(1:3, size = length(dset$dataset$mask), replace = TRUE),
    neuroim2::space(dset$dataset$mask)
  )

  neuroim2::write_vec(dset$dataset$train_data, train_path)
  neuroim2::write_vol(region_mask, region_mask_path)
  utils::write.table(
    dset$design$train_design,
    file = design_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  result <- NULL
  expect_no_error({
    result <- rMVPA:::run_cli_command(
      "regional",
      argv = c(
        "--train-design", design_path,
        "--train-data", train_path,
        "--mask", region_mask_path,
        "--label-column", "Y",
        "--model", "sda_notune",
        "--workers", "1",
        "--preflight", "off",
        "--pool-predictions", "mean",
        "--save-level", "minimal",
        "--output", output_dir
      )
    )
  })

  expect_s3_class(result, "rmvpa_analysis_run")
  expect_true(dir.exists(output_dir))
  expect_true(file.exists(file.path(output_dir, "aux", "analysis_config.rds")))
  expect_equal(length(list.files(output_dir, pattern = "^manifest\\.", full.names = TRUE)), 1L)
  expect_true(file.exists(file.path(output_dir, "analysis", "performance_table.txt")))
  expect_true(file.exists(file.path(output_dir, "analysis", "prediction_table.txt")))
  expect_true(file.exists(file.path(output_dir, "analysis", "pooled_prediction_table.txt")))
  expect_true(file.exists(file.path(output_dir, "analysis", "pooled_performance_table.txt")))
  expect_false(dir.exists(file.path(output_dir, "analysis", "fits")))
})
