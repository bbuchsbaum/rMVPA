test_that("mvpa_design accepts supported block variable specifications", {
  design <- data.frame(
    condition = factor(rep(c("a", "b"), 4)),
    subject = rep(c("s1", "s2"), each = 4),
    run = rep(1:4, each = 2)
  )

  formula_design <- mvpa_design(design, y_train = ~ condition, block_var = ~ run)
  named_design <- mvpa_design(design, y_train = ~ condition, block_var = "run")
  multi_formula <- mvpa_design(
    design,
    y_train = ~ condition,
    block_var = ~ subject + run
  )

  expect_identical(formula_design$block_var, design$run)
  expect_identical(named_design$block_var, design$run)
  expect_identical(
    multi_formula$block_var,
    droplevels(interaction(factor(design$subject), factor(design$run), sep = ":"))
  )
})

test_that("mvpa_design accepts row-aligned atomic block vectors", {
  design <- data.frame(condition = factor(rep(c("a", "b"), 4)))
  inputs <- list(
    factor = factor(rep(c("r1", "r2"), each = 4), levels = c("r1", "r2", "unused")),
    integer = rep(1:2, each = 4),
    double = as.double(rep(1:2, each = 4)),
    character = rep(c("r1", "r2"), each = 4),
    logical = rep(c(FALSE, TRUE), each = 4)
  )

  parsed <- lapply(inputs, function(blocks) {
    mvpa_design(design, y_train = ~ condition, block_var = blocks)$block_var
  })

  expect_identical(parsed$factor, droplevels(inputs$factor))
  expect_identical(parsed$integer, inputs$integer)
  expect_identical(parsed$double, inputs$double)
  expect_identical(parsed$character, inputs$character)
  expect_identical(parsed$logical, inputs$logical)

  split_design <- mvpa_design(
    design,
    y_train = ~ condition,
    split_by = inputs$character
  )
  expect_identical(
    split_design$split_groups,
    split(seq_len(nrow(design)), inputs$character)
  )
})

test_that("mvpa_design rejects ambiguous or misaligned variable inputs clearly", {
  design <- data.frame(
    condition = factor(rep(c("a", "b"), 4)),
    run = rep(1:4, each = 2)
  )

  expect_error(
    mvpa_design(design, y_train = ~ condition, block_var = c(1, 2)),
    "must have length 8.*got 2"
  )
  expect_error(
    mvpa_design(design, y_train = ~ condition, block_var = ~ missing_run),
    "Formula variable not found.*missing_run"
  )
  expect_error(
    mvpa_design(design, y_train = ~ condition, block_var = "missing_run"),
    "does not contain variable named"
  )
  expect_error(
    mvpa_design(design, y_train = ~ condition, block_var = matrix(1:8, ncol = 1)),
    "row-aligned"
  )
})
