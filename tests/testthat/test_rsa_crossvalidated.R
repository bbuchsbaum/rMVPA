cv_rsa_dataset <- function(x) {
  dims <- c(ncol(x), 1L, 1L)
  mvpa_dataset(
    neuroim2::NeuroVec(array(t(x), c(dims, nrow(x))),
                      neuroim2::NeuroSpace(c(dims, nrow(x)))),
    mask = neuroim2::LogicalNeuroVol(array(TRUE, dims), neuroim2::NeuroSpace(dims))
  )
}

cv_rsa_fixture <- function() {
  # Three runs, three conditions, two features. A-B differences are
  # (1,0), (-1,0), (0,0), so the six ordered cross-products sum to -2.
  x <- rbind(c(1,0), c(0,0), c(0,2),
             c(-1,0), c(0,0), c(1,3),
             c(0,0), c(0,0), c(2,1))
  ids <- c("A", "B", "C")
  d <- rsa_design(~ feature, list(feature = dist(c(0,1,3))), condition_ids = ids)
  labels <- rep(ids, 3)
  runs <- rep(1:3, each = 3)
  m <- rsa_model(cv_rsa_dataset(x), d, distmethod = "crossvalidated_euclidean",
                 condition_labels = labels, crossval = blocked_cross_validation(runs))
  list(x=x, ids=ids, labels=labels, runs=runs, model=m)
}

# Slow pairwise oracle, independent of the production Gram-matrix identity.
cv_rsa_oracle <- function(x, labels, runs, ids) {
  ans <- numeric(choose(length(ids), 2))
  pairs <- which(lower.tri(matrix(0, length(ids), length(ids))), arr.ind=TRUE)
  for (p in seq_len(nrow(pairs))) {
    delta <- lapply(unique(runs), function(r) {
      colMeans(x[runs == r & labels == ids[pairs[p,1]], , drop=FALSE]) -
        colMeans(x[runs == r & labels == ids[pairs[p,2]], , drop=FALSE])
    })
    cross <- c()
    for (a in seq_along(delta)) for (b in seq_along(delta)) {
      if (a != b) cross <- c(cross, sum(delta[[a]] * delta[[b]]))
    }
    ans[p] <- mean(cross) / ncol(x)
  }
  ans
}

test_that("ordinary RSA preserves the three-run distance oracle and negative values", {
  f <- cv_rsa_fixture()
  response <- rMVPA:::.rsa_crossvalidated_response(f$model, f$x)
  # Hand calculation: cross-run product sums are -2, 22 and 26;
  # divide each by six ordered run pairs and two features.
  expect_equal(unname(response), c(-1/6, 11/6, 13/6), tolerance=1e-14)
  expect_equal(unname(response), cv_rsa_oracle(f$x, f$labels, f$runs, f$ids), tolerance=1e-14)
  expected <- cor(unname(response), as.vector(dist(c(0,1,3))))
  expect_equal(unname(train_model(f$model, f$x, NULL, NULL)), expected)
  f$model$.fast_kernel <- NULL
  expect_equal(unname(train_model(f$model, f$x, NULL, NULL)), expected)
  for (association in c("pearson", "spearman")) {
    m <- rsa_model(f$model$dataset, f$model$design, distmethod="crossvalidated_euclidean",
                   regtype=association, condition_labels=f$labels,
                   crossval=blocked_cross_validation(f$runs))
    expect_equal(unname(train_model(m, f$x, NULL, NULL)),
                 cor(unname(response), as.vector(dist(c(0,1,3))), method=association))
  }
})

test_that("conditions, observations, run names and feature ordering do not alter the estimand", {
  f <- cv_rsa_fixture()
  expected <- train_model(f$model, f$x, NULL, NULL)
  rows <- c(9,2,5,3,7,1,6,4,8)
  m <- rsa_model(cv_rsa_dataset(f$x[rows,]), f$model$design,
                 distmethod="crossvalidated_euclidean", condition_labels=f$labels[rows],
                 crossval=blocked_cross_validation(c("z","a","q")[f$runs[rows]]))
  expect_equal(train_model(m, f$x[rows,2:1], NULL, NULL), expected)
  ord <- c(3,1,2)
  d <- rsa_design(~ feature, list(feature=dist(c(0,1,3)[ord])), condition_ids=f$ids[ord])
  m <- rsa_model(f$model$dataset, d, distmethod="crossvalidated_euclidean",
                 condition_labels=f$labels, crossval=blocked_cross_validation(f$runs))
  expect_equal(train_model(m, f$x, NULL, NULL), expected)
  # Duplicate observations within a run are means, not extra partitions.
  m <- rsa_model(cv_rsa_dataset(f$x[rep(1:9,each=2),]), f$model$design,
                 distmethod="crossvalidated_euclidean", condition_labels=rep(f$labels,each=2),
                 crossval=blocked_cross_validation(rep(f$runs,each=2)))
  expect_equal(train_model(m, f$x[rep(1:9,each=2),], NULL, NULL), expected)
})

test_that("the same 270 pairs feed responses, model, nuisance and diagnostics", {
  set.seed(125)
  ids <- paste0("item", 1:60)
  recording <- rep(1:6, each=10)
  mask <- outer(recording, recording, `==`)
  diag(mask) <- FALSE
  feat <- matrix(rnorm(60*3),60)
  noise <- matrix(rnorm(60*2),60)
  d <- rsa_design(~ feature, list(feature=dist(feat)), nuisance=list(n=dist(noise)),
                  pair_mask=mask, condition_ids=ids)
  expect_equal(sum(d$include), 270)
  expect_equal(d$model_mat$feature, as.vector(dist(feat))[mask[lower.tri(mask)]])
  expect_equal(d$model_mat$n, as.vector(dist(noise))[mask[lower.tri(mask)]])
  expect_equal(rsa_design_diagnostics(d)$n_pairs, 270)
  x <- matrix(rnorm(180*4),180)
  labels <- rep(ids,3); runs <- rep(1:3,each=60)
  m <- rsa_model(cv_rsa_dataset(x), d, distmethod="crossvalidated_euclidean",
                 condition_labels=labels, crossval=blocked_cross_validation(runs),
                 regtype="lm", statistic="beta", check_collinearity=FALSE)
  y <- cv_rsa_oracle(x, labels, runs, ids)[d$include]
  oracle <- coef(lm(y ~ feature + n, data=as.data.frame(d$model_mat)))["feature"]
  expect_equal(train_model(m, x, NULL, NULL), oracle, tolerance=1e-12)
  # Explicit masking permits missing values only outside the declared pair set.
  bad <- as.matrix(dist(feat)); bad[!mask] <- NA_real_; diag(bad) <- 0
  expect_equal(rsa_design(~ feature, list(feature=bad), pair_mask=mask)$model_mat$feature,
               d$model_mat$feature)
  bad[1,2] <- bad[2,1] <- NA_real_
  expect_error(rsa_design(~ feature, list(feature=bad), pair_mask=mask), "finite")
})

test_that("pair masks are validated and intersect block exclusions", {
  D <- dist(1:4)
  expect_error(rsa_design(~ D, list(D=D), pair_mask=c(TRUE,NA)), "logical")
  expect_error(rsa_design(~ D, list(D=D), pair_mask=rep(TRUE,5)), "choose")
  expect_error(rsa_design(~ D, list(D=D), pair_mask=matrix(TRUE,4,4)), "diagonal")
  mask <- matrix(FALSE,4,4); mask[1,2] <- TRUE
  expect_error(rsa_design(~ D, list(D=D), pair_mask=mask), "symmetric")
  expect_error(rsa_design(~ D, list(D=D), pair_mask=rep(FALSE,6)), "no eligible")
  expect_error(rsa_design(~ D, list(D=D), condition_ids=c("a","a","b","c")), "unique")
  labelled <- D; attr(labelled,"Labels") <- letters[4:1]
  expect_error(rsa_design(~ labelled, list(labelled=labelled), condition_ids=letters[1:4]), "labels")
  labelled <- as.matrix(D); dimnames(labelled) <- list(letters[1:4],letters[4:1])
  expect_error(rsa_design(~ labelled, list(labelled=labelled), condition_ids=letters[1:4]), "labels")
  d <- rsa_design(~ D, list(D=D, block=c(1,1,2,2)), block_var="block", pair_mask=rep(TRUE,6))
  expect_identical(d$include, c(FALSE,TRUE,TRUE,TRUE,TRUE,FALSE))
})

test_that("missing cells, overlapping partitions and invalid controls fail explicitly", {
  f <- cv_rsa_fixture()
  make <- function(labels=f$labels, cv=blocked_cross_validation(f$runs), ...) {
    rsa_model(f$model$dataset, f$model$design, distmethod="crossvalidated_euclidean",
               condition_labels=labels, crossval=cv, ...)
  }
  labels <- f$labels; labels[1] <- "B"
  expect_error(make(labels), "Missing condition/partition")
  labels[1] <- NA
  expect_error(make(labels), "condition_labels")
  expect_error(make(cv=blocked_cross_validation(rep(1,9))), "at least two")
  expect_error(make(cv=NULL), "crossval")
  overlap <- custom_cross_validation(list(list(train=4:9, test=1:6),
                                          list(train=1:3, test=4:9)))
  expect_error(make(cv=overlap), "disjoint")
  incomplete <- custom_cross_validation(list(list(train=4:9, test=1:3),
                                             list(train=c(1:3,7:9), test=4:6)))
  expect_error(make(cv=incomplete), "cover every observation")
  expect_error(make(measure="similarity"), "measure")
  x <- f$x; x[1,1] <- NA
  expect_error(train_model(f$model, x, NULL, NULL), "finite")
  expect_error(train_model(f$model, f$x[-1,], NULL, NULL), "declared rows")
  expect_error(rsa_model(f$model$dataset, f$model$design, crossval=blocked_cross_validation(f$runs)), "require distmethod")
})

test_that("native ordinary searchlights match direct sphere fits and reject rsa_fast", {
  f <- cv_rsa_fixture()
  direct <- train_model(f$model, f$x, NULL, NULL)
  expect_error(run_searchlight(f$model, radius=10, engine="rsa_fast", preflight="off"), "not implemented")
  out <- run_searchlight(f$model, radius=10, method="standard", engine="legacy", fail_fast=TRUE, preflight="off")
  expect_equal(out$metrics, "feature")
  # All voxels are included in every sphere of this tiny image.
  expect_equal(as.numeric(neuroim2::values(out$results$feature)), rep(unname(direct), 2), tolerance=1e-12)
})

test_that("new pair masks also preserve correlation RSA and use explicit exchangeability blocks", {
  set.seed(126)
  n <- 8
  recording <- rep(1:2,each=4)
  mask <- outer(recording,recording,`==`); diag(mask) <- FALSE
  d <- rsa_design(~ feature, list(feature=dist(1:n), recording=recording),
                  pair_mask=mask, block_var="recording", keep_intra_run=TRUE)
  x <- matrix(rnorm(n*5),n)
  m <- rsa_model(cv_rsa_dataset(x), d, distmethod="pearson")
  response <- as.vector(as.dist(1-cor(t(x))))[d$include]
  expect_equal(unname(train_model(m,x,NULL,NULL)), cor(response,d$model_mat$feature))
  shuffled <- permute_labels(d, method="within_block", seed=12)
  expect_true(all(recording[shuffled$item_perm] == recording))
  d$block_var <- NULL
  expect_error(permute_labels(d,method="within_block",seed=12), "preserve.*pair_mask")
})

test_that("condition permutations apply jointly to all partitions and centering cancels", {
  f <- cv_rsa_fixture()
  p <- c(3,1,2)
  f$model$design$item_perm <- p
  wanted <- cv_rsa_oracle(f$x, f$labels, f$runs, f$ids[p])
  expect_equal(unname(rMVPA:::.rsa_crossvalidated_response(f$model,f$x)), wanted)
  response <- train_model(f$model,f$x,NULL,NULL)
  f$model$pattern_center <- "stimulus_mean"
  expect_equal(train_model(f$model,f$x,NULL,NULL),response)
  expect_equal(train_model(f$model,sweep(f$x,2,c(10,-4),`+`),NULL,NULL),response)
  # Input units are retained: multiplying by three multiplies distances by nine.
  expect_equal(unname(rMVPA:::.rsa_crossvalidated_response(f$model,3*f$x)),9*wanted)
})

cv_rsa_retention_case <- function(backend) {
  f <- cv_rsa_fixture()
  # A constant third feature must contribute to P, despite its zero difference.
  x <- cbind(f$x, 7)
  m <- rsa_model(cv_rsa_dataset(x), f$model$design,
                 distmethod="crossvalidated_euclidean", condition_labels=f$labels,
                 crossval=blocked_cross_validation(f$runs),
                 regtype="lm", statistic="beta", check_collinearity=FALSE)
  direct <- unname(train_model(m,x,NULL,NULL))
  wanted <- coef(lm(cv_rsa_oracle(x,f$labels,f$runs,f$ids) ~ as.vector(dist(c(0,1,3)))))[2]
  expect_equal(direct,unname(wanted),tolerance=1e-12)
  out <- run_searchlight(m,radius=10,backend=backend,preflight="off",fail_fast=TRUE)
  expect_equal(as.numeric(neuroim2::values(out$results$feature)), rep(direct,3), tolerance=1e-12)
  # NA at voxel three must fail spheres centered on voxels one and two as well.
  x[1,3] <- NA
  m$dataset <- cv_rsa_dataset(x)
  expect_error(run_searchlight(m,radius=10,backend=backend,preflight="off",fail_fast=TRUE),
               "finite|No valid results")
}

test_that("ordinary sphere extraction preserves constant and missing features", {
  cv_rsa_retention_case("default")
})

test_that("shard sphere extraction preserves constant and missing features", {
  skip_if_not_installed("shard")
  cv_rsa_retention_case("shard")
})

test_that("pair-mask diagnostics count only participating items", {
  mask <- matrix(FALSE,100,100); mask[1:3,1:3] <- TRUE; diag(mask) <- FALSE
  d <- rsa_design(~ feature,list(feature=dist(1:100)),pair_mask=mask)
  diagnostics <- rsa_design_diagnostics(d)
  expect_equal(diagnostics$n_pairs,3)
  expect_equal(diagnostics$n_items,3)
  expect_equal(unname(diagnostics$items_per_predictor),3)
  expect_equal(rMVPA:::.rsa_design_n_items(d),100)
})
