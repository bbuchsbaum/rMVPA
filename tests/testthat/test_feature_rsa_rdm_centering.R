library(rMVPA)

# Deliberately dense reference: a projection matrix defines both operations.
centred_rdm_oracle <- function(x, mode) {
  n <- nrow(x)
  h <- diag(n) - matrix(1 / n, n, n)
  if (mode != "none") x <- h %*% x
  d <- suppressWarnings(1 - stats::cor(t(x)))
  if (mode == "double") d <- h %*% d %*% h
  as.numeric(d[lower.tri(d)])
}

test_that("each RDM centring mode agrees with hand-computed distances", {
  # x has voxel mean (1, 2, 4). Deviations are u, v, -2u, u-v,
  # with u=(1,-1,0), v=(0,1,-1), u.v=-1 and |u|=|v|=sqrt(2).
  x <- rbind(c(2,1,4), c(1,3,3), c(-1,4,4), c(2,0,5))
  a <- sqrt(3) / 2
  g <- (14 + 2*a) / 16
  expected <- list(
    none = 1 - c(6/sqrt(42*24), 15/sqrt(42*150), 69/sqrt(42*114),
                 1, 6/sqrt(24*114), 15/sqrt(150*114)),
    items = c(1.5, 2, 1-a, .5, 1+a, 1+a),
    double = c(1.5-7.5/4+g, g, 1-a-7.5/4+g,
               .5-(6.5+2*a)/4+g, 1+a-(6+2*a)/4+g,
               1+a-(6.5+2*a)/4+g)
  )
  for (mode in names(expected)) {
    for (block in c(1L, 3L, 128L)) {
      withr::local_options(rMVPA.feature_rsa_metric_block_rows = block)
      got <- rMVPA:::.feature_rsa_rdm_vector_blockwise(x, mode)
      expect_equal(got, expected[[mode]], tolerance = 2e-14, info = mode)
    }
  }
})

test_that("centred scores and stored vectors use all items within each fold", {
  set.seed(240901)
  p <- matrix(rnorm(17*13),17,13)
  o <- .6*p + matrix(rnorm(length(p)),17,13)
  # Constant raw rows can have spatial variance after voxel centring.
  p[2,] <- 3
  o[3,] <- 5
  fold <- rep(c(2L,1L),c(8L,9L))
  for (mode in c("items","double","none")) {
    got <- evaluate_model.feature_rsa_model(NULL,p,o,fold_id=fold,
      rdm_centering=mode,compute_rdm_vectors=TRUE)
    pd <- od <- matrix(NA_real_,17,17)
    for (f in unique(fold)) {
      ii <- which(fold==f)
      pd[ii,ii][lower.tri(pd[ii,ii])] <- centred_rdm_oracle(p[ii,],mode)
      od[ii,ii][lower.tri(od[ii,ii])] <- centred_rdm_oracle(o[ii,],mode)
    }
    expect_equal(got$predicted_rdm_vec, pd[lower.tri(pd)], tolerance=2e-12)
    expect_equal(got$observed_rdm_vec, od[lower.tri(od)], tolerance=2e-12)
    expect_equal(got$rdm_correlation,
      suppressWarnings(cor(pd[lower.tri(pd)], od[lower.tri(od)],
                           method="spearman",use="complete.obs")), tolerance=2e-12)
    expect_identical(got$rdm_centering,mode)
  }
  expect_error(evaluate_model.feature_rsa_model(NULL,p,o,rdm_centering="bad"),"arg")
})

test_that("centred geometry is invariant to fold-specific shared patterns", {
  set.seed(240902)
  p <- matrix(rnorm(24*20),24,20)
  o <- p + matrix(rnorm(length(p),sd=.2),24,20)
  fold <- rep(1:3,each=8)
  shift <- matrix(rnorm(3*20,sd=20),3,20)[fold,]
  for (mode in c("items","double")) {
    before <- evaluate_model.feature_rsa_model(NULL,p,o,fold_id=fold,rdm_centering=mode)
    after <- evaluate_model.feature_rsa_model(NULL,p+shift,o-2*shift,
      fold_id=fold,rdm_centering=mode)
    expect_equal(after$rdm_correlation,before$rdm_correlation,tolerance=2e-12)
    expect_gt(before$rdm_correlation,.8)
    expect_gt(abs(after$rdm_correlation_raw-before$rdm_correlation_raw),.01)
  }
})

test_that("centred independent-pattern nulls average near zero", {
  set.seed(240903)
  scores <- replicate(30, {
    p <- matrix(rnorm(60*30),60,30)
    o <- matrix(rnorm(60*30),60,30)
    vapply(c("items","double"),function(mode)
      evaluate_model.feature_rsa_model(NULL,p,o,fold_id=rep(1:3,each=20),
                                      rdm_centering=mode)$rdm_correlation,numeric(1))
  })
  expect_true(all(is.finite(scores)))
  expect_lt(max(abs(rowMeans(scores))), .05)
})

test_that("centred geometry preserves undefined distances and missing pairs", {
  x <- rbind(c(1,2,3,4),c(4,3,2,1),c(2,4,1,3),c(3,1,4,2))
  x <- rbind(x,colMeans(x))
  for (mode in c("items","double")) {
    v <- rMVPA:::.feature_rsa_rdm_vector_blockwise(x,mode)
    expect_true(anyNA(v))
    if (mode=="double") expect_true(all(is.na(v)))
    constant <- matrix(c(1,2,3,4),6,4,byrow=TRUE)
    expect_true(all(is.na(rMVPA:::.feature_rsa_rdm_vector_blockwise(constant,mode))))
  }
  x[1,1] <- NA_real_
  xc <- sweep(x,2,colMeans(x,na.rm=TRUE))
  d <- suppressWarnings(1-cor(t(xc),use="pairwise.complete.obs"))
  expect_equal(suppressWarnings(rMVPA:::.feature_rsa_rdm_vector_blockwise(x,"items")),
               d[lower.tri(d)],tolerance=2e-12)
})

test_that("centred permutation caches match complete within-fold recomputation", {
  set.seed(240904)
  p <- matrix(rnorm(15*10),15,10)
  o <- p/2 + matrix(rnorm(length(p)),15,10)
  fold <- rep(1:3,each=5)
  for (mode in c("items","double","none")) {
    set.seed(240905)
    got <- evaluate_model.feature_rsa_model(NULL,p,o,fold_id=fold,
      rdm_centering=mode,nperm=12,save_distributions=TRUE)
    set.seed(240905)
    expected <- replicate(12, {
      idx <- unlist(lapply(split(seq_len(nrow(p)),fold),sample))
      pp <- unlist(lapply(split(seq_len(nrow(p)),fold),function(ii)
        centred_rdm_oracle(p[idx[ii],],mode)))
      oo <- unlist(lapply(split(seq_len(nrow(o)),fold),function(ii)
        centred_rdm_oracle(o[ii,],mode)))
      cor(pp,oo,method="spearman")
    })
    expect_equal(got$permutation_results$permutation_distributions$rdm_correlation,
                 expected,tolerance=2e-12)
    expect_equal(got$permutation_results$p_values[["rdm_correlation"]],
                 (1+sum(expected>=got$rdm_correlation))/13)
    raw <- evaluate_model.feature_rsa_model(NULL,p,o,fold_id=fold,rdm_centering="none")
    expect_equal(got$rdm_correlation_raw,raw$rdm_correlation)
  }
})

test_that("model, retained payloads and disk batches record RDM centring", {
  set.seed(240906)
  sim <- gen_sample_dataset(c(3,3,3),nobs=24,blocks=3)
  design <- feature_rsa_design(F=matrix(rnorm(24*4),24,4),
    labels=paste0("t",1:24),max_comps=3,block_var=sim$design$block_var)
  model <- feature_rsa_model(sim$dataset,design,method="ridge",lambda=1,
    lambda_selection="fixed",return_rdm_vectors=TRUE,return_predictions=TRUE)
  expect_identical(model$rdm_centering,"items")
  expect_error(feature_rsa_model(sim$dataset,design,rdm_centering="invalid"),"arg")
  mask <- neuroim2::NeuroVol(rep(1,27),neuroim2::space(sim$dataset$mask))
  model$rdm_centering <- "double"
  result <- run_regional(model,mask)
  expect_identical(result$model_spec$rdm_centering,"double")
  expect_identical(feature_rsa_rdm_vectors(result)$rdm_centering,"double")
  pred <- feature_rsa_predictions(result)
  expect_identical(pred$rdm_centering,"double")
  rebuilt <- evaluate_model.feature_rsa_model(model,pred$predicted[[1]],pred$observed[[1]],
    fold_id=pred$fold_id[[1]])
  expect_equal(result$performance_table$rdm_correlation,rebuilt$rdm_correlation)
  expect_equal(result$performance_table$rdm_correlation_raw,rebuilt$rdm_correlation_raw)
  batch <- tibble::tibble(id=1L,result=list(list(predictor=result$fits[[1]])))
  folder <- tempfile("rdm-centering-")
  dir.create(folder)
  on.exit(unlink(folder,recursive=TRUE),add=TRUE)
  rMVPA:::.mvpa_write_feature_rsa_rdm_batch(batch,1L,folder)
  result$rdm_batch_dir <- folder
  expect_equal(feature_rsa_rdm_vectors(result)$rdm_centering,"double")
  # Old retained vectors remain identifiable as raw.
  result$rdm_batch_dir <- NULL
  result$fits[[1]]$rdm_centering <- NULL
  expect_identical(feature_rsa_rdm_vectors(result)$rdm_centering,"none")
})

test_that("double centring controls amplitude-linked item eccentricity", {
  # All items have orthogonal directions. Opposing item amplitudes around a
  # dominant common pattern create anticorrelated raw RDM row means. This is
  # a nuisance stress case, not evidence of predictive validity on real data.
  n <- 40L
  basis <- qr.Q(qr(contr.helmert(n+2L)))
  amplitude <- seq(.2,3,length.out=n)
  directions <- t(basis[,2:(n+1L)])
  common <- matrix(20*basis[,1],n,n+2L,byrow=TRUE)
  p <- sweep(directions,1,amplitude,"*") + common
  o <- sweep(directions,1,1/amplitude,"*") + common
  items <- evaluate_model.feature_rsa_model(NULL,p,o)
  double <- evaluate_model.feature_rsa_model(NULL,p,o,rdm_centering="double")
  expect_lt(items$rdm_correlation_raw,-.7)
  # Voxel centring alone need not eliminate item effects or ensure positivity.
  expect_lt(items$rdm_correlation,0)
  expect_gt(double$rdm_correlation,.9)
})

test_that("constant patterns and two-item folds do not rank roundoff noise", {
  constant <- matrix(rep(c(.1,.2,.7,pi),each=25),25,4)
  set.seed(240907)
  p <- matrix(rnorm(20*10),20,10)
  o <- matrix(rnorm(20*10),20,10)
  for (mode in c("items","double")) {
    expect_true(all(is.na(rMVPA:::.feature_rsa_rdm_vector_blockwise(constant,mode))))
    got <- evaluate_model.feature_rsa_model(NULL,p,o,fold_id=rep(1:10,each=2),
      rdm_centering=mode,compute_rdm_vectors=TRUE,nperm=3)
    expect_true(is.na(got$rdm_correlation))
    expect_true(is.na(got$permutation_results$p_values[["rdm_correlation"]]))
    finite <- got$predicted_rdm_vec[is.finite(got$predicted_rdm_vec)]
    expect_equal(finite,rep(if(mode=="items") 2 else 1,10))
  }
})
