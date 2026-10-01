suppressMessages(pkgload::load_all("/Users/bbuchsbaum/code/rMVPA", quiet=TRUE))
futile.logger::flog.threshold(futile.logger::ERROR)
set.seed(1)
ds <- gen_sample_dataset(c(12,12,12), 40, nlevels=4, blocks=5)
nv <- sum(ds$dataset$mask>0)
D1 <- dist(matrix(rnorm(40*3),40)); D2 <- dist(matrix(rnorm(40*3),40))
rdes <- rsa_design(~ D1 + D2, list(D1 = D1, D2 = D2, block = ds$design$block_var), block_var = "block")
ms <- rsa_model(ds$dataset, rdes, regtype = "lm", distmethod = "pearson", check_collinearity = FALSE)
Rprof("rsa.out", interval=0.005)
t <- system.time(r <- run_searchlight(ms, radius=3))[3]
Rprof(NULL)
cat("RSA lm searchlight 40 obs, r=3:", t, "s; ms/sphere", 1000*t/nv, " engine", attr(r,"searchlight_engine"),"\n")
s <- summaryRprof("rsa.out"); print(head(s$by.self,12))
tt<-s$by.total; print(head(tt[grepl("flog|tibble|fit_roi|merge|train_model|pairwise|cor|run_lm|gc|extract_roi|filter_roi|bind_rows|roi_result",rownames(tt)),c(1,2)],25))
# regional: 10 ROIs
set.seed(2)
ds2 <- gen_sample_dataset(c(20,20,20), 100, nlevels=4, blocks=5)
reg <- neuroim2::NeuroVol(sample(1:10, prod(c(20,20,20)), replace=TRUE), neuroim2::space(ds2$dataset$mask))
cv <- blocked_cross_validation(ds2$design$block_var)
for (m in c("corclass","sda_notune")) {
mm <- mvpa_model(load_model(m), ds2$dataset, ds2$design, "classification", crossval=cv)
t <- system.time(rr <- run_regional(mm, reg))[3]
cat("regional", m, "10 ROIs x ~800 vox, 100 obs:", t, "s\n")
}
