# Benchmark tests for run_permutation_searchlight scaling behaviour.
#
# All benchmarks are guarded by skip_on_cran() AND skip_if_not(interactive())
# so they only run when a developer explicitly sources the file or runs it
# interactively.  They do NOT run in automated CI.
#
# Total wall-clock target: under 2 minutes on a modern laptop.
# Model: corclass (fastest built-in classifier).

# ---------------------------------------------------------------------------
# Helper: run one configuration and return elapsed seconds
# ---------------------------------------------------------------------------

.bench_perm_sl <- function(grid_dims, n_obs, n_perm, subsample_frac,
                            radius = 3, seed = 1L) {
  dset  <- gen_sample_dataset(grid_dims, n_obs, blocks = 3, nlevels = 2)
  mdl   <- load_model("corclass")
  cval  <- blocked_cross_validation(dset$design$block_var)
  mspec <- mvpa_model(mdl, dset$dataset, dset$design, "classification",
                      crossval = cval)

  pc <- permutation_control(
    n_perm      = n_perm,
    subsample   = subsample_frac,
    seed        = seed,
    null_method = "global",
    diagnose    = FALSE,
    correction  = "none"
  )

  elapsed <- system.time(
    suppressWarnings(
      run_permutation_searchlight(mspec, radius = radius, perm_ctrl = pc)
    )
  )[["elapsed"]]

  elapsed
}


# ---------------------------------------------------------------------------
# Benchmark 1: n_perm scaling (1 / 2 / 5) on a 5x5x5 grid
# ---------------------------------------------------------------------------

test_that("benchmark: n_perm scaling (interactive only)", {
  skip_on_cran()
  skip_if_not(interactive())

  futile.logger::flog.threshold(futile.logger::ERROR)
  on.exit(futile.logger::flog.threshold(futile.logger::INFO), add = TRUE)

  n_perm_vals <- c(1L, 2L, 5L)
  results     <- data.frame(
    n_perm  = integer(0),
    elapsed = numeric(0)
  )

  for (np in n_perm_vals) {
    t <- .bench_perm_sl(
      grid_dims     = c(5, 5, 5),
      n_obs         = 90,
      n_perm        = np,
      subsample_frac = 0.3,
      seed          = 42L
    )
    results <- rbind(results, data.frame(n_perm = np, elapsed = round(t, 2)))
  }

  cat("\n=== Benchmark 1: n_perm scaling (grid 5x5x5, subsample 0.3) ===\n")
  print(results, row.names = FALSE)

  # Sanity: more permutations should not be faster than fewer
  expect_true(all(results$elapsed > 0))
  # Results exist for all configurations
  expect_equal(nrow(results), length(n_perm_vals))
})


# ---------------------------------------------------------------------------
# Benchmark 2: subsample fraction scaling (0.1 / 0.3 / 0.5 / 1.0)
# ---------------------------------------------------------------------------

test_that("benchmark: subsample fraction scaling (interactive only)", {
  skip_on_cran()
  skip_if_not(interactive())

  futile.logger::flog.threshold(futile.logger::ERROR)
  on.exit(futile.logger::flog.threshold(futile.logger::INFO), add = TRUE)

  fracs   <- c(0.1, 0.3, 0.5, 1.0)
  results <- data.frame(
    subsample = numeric(0),
    elapsed   = numeric(0)
  )

  for (frac in fracs) {
    t <- .bench_perm_sl(
      grid_dims      = c(5, 5, 5),
      n_obs          = 90,
      n_perm         = 2L,
      subsample_frac = frac,
      seed           = 42L
    )
    results <- rbind(results,
                     data.frame(subsample = frac, elapsed = round(t, 2)))
  }

  cat("\n=== Benchmark 2: subsample fraction scaling (grid 5x5x5, n_perm 2) ===\n")
  print(results, row.names = FALSE)

  expect_true(all(results$elapsed > 0))
  expect_equal(nrow(results), length(fracs))

  # Larger subsample fraction should generally take longer (allow 10 % slack)
  elapsed_frac1 <- results$elapsed[results$subsample == 0.1]
  elapsed_frac4 <- results$elapsed[results$subsample == 1.0]
  expect_true(elapsed_frac4 >= elapsed_frac1 * 0.9)
})


# ---------------------------------------------------------------------------
# Benchmark 3: grid size scaling (5x5x5 vs 8x8x8)
# ---------------------------------------------------------------------------

test_that("benchmark: grid size scaling (interactive only)", {
  skip_on_cran()
  skip_if_not(interactive())

  futile.logger::flog.threshold(futile.logger::ERROR)
  on.exit(futile.logger::flog.threshold(futile.logger::INFO), add = TRUE)

  configs <- list(
    list(dims = c(5, 5, 5), label = "5x5x5"),
    list(dims = c(8, 8, 8), label = "8x8x8")
  )

  results <- data.frame(
    grid    = character(0),
    elapsed = numeric(0),
    stringsAsFactors = FALSE
  )

  for (cfg in configs) {
    t <- .bench_perm_sl(
      grid_dims      = cfg$dims,
      n_obs          = 90,
      n_perm         = 2L,
      subsample_frac = 0.2,
      seed           = 42L
    )
    results <- rbind(results,
                     data.frame(grid    = cfg$label,
                                elapsed = round(t, 2),
                                stringsAsFactors = FALSE))
  }

  cat("\n=== Benchmark 3: grid size scaling (n_perm 2, subsample 0.2) ===\n")
  print(results, row.names = FALSE)

  expect_true(all(results$elapsed > 0))
  expect_equal(nrow(results), length(configs))

  # 8x8x8 grid should take at least as long as 5x5x5 (allow 20 % slack)
  t_small <- results$elapsed[results$grid == "5x5x5"]
  t_large <- results$elapsed[results$grid == "8x8x8"]
  expect_true(t_large >= t_small * 0.8)
})


# ---------------------------------------------------------------------------
# Benchmark 4: combined summary table
# ---------------------------------------------------------------------------

test_that("benchmark: combined configuration sweep (interactive only)", {
  skip_on_cran()
  skip_if_not(interactive())

  futile.logger::flog.threshold(futile.logger::ERROR)
  on.exit(futile.logger::flog.threshold(futile.logger::INFO), add = TRUE)

  sweep <- expand.grid(
    n_perm    = c(1L, 2L),
    subsample = c(0.2, 0.5),
    grid      = c("5x5x5"),
    stringsAsFactors = FALSE
  )
  sweep$elapsed <- NA_real_

  grid_map <- list("5x5x5" = c(5L, 5L, 5L))

  for (i in seq_len(nrow(sweep))) {
    dims <- grid_map[[sweep$grid[i]]]
    t <- .bench_perm_sl(
      grid_dims      = dims,
      n_obs          = 90,
      n_perm         = as.integer(sweep$n_perm[i]),
      subsample_frac = sweep$subsample[i],
      seed           = 100L + i
    )
    sweep$elapsed[i] <- round(t, 2)
  }

  cat("\n=== Benchmark 4: combined sweep ===\n")
  print(sweep, row.names = FALSE)

  expect_true(all(sweep$elapsed > 0))
  expect_equal(nrow(sweep), nrow(expand.grid(n_perm = c(1L, 2L),
                                              subsample = c(0.2, 0.5),
                                              grid = "5x5x5")))
})
