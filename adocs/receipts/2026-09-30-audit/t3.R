suppressMessages(pkgload::load_all("/Users/bbuchsbaum/code/rMVPA", quiet=TRUE))
futile.logger::flog.threshold(futile.logger::ERROR)
f <- function() futile.logger::flog.debug("x %s", paste(1:3, collapse=" "))
environment(f) <- asNamespace("rMVPA")
cat("flog.debug per call (us):", 1e6*system.time(for(i in 1:500) f())[3]/500, "\n")
set.seed(1)
ds <- gen_sample_dataset(c(9,9,9), 100, nlevels=4, blocks=5)
cv <- blocked_cross_validation(ds$design$block_var)
nv <- sum(ds$dataset$mask>0)
for (m in c("corclass","sda_notune")) {
ms <- mvpa_model(load_model(m), ds$dataset, ds$design, "classification", crossval=cv)
t0 <- system.time(run_searchlight(ms, radius=3, engine="legacy", backend="default"))[3]
t0s <- system.time(run_searchlight(ms, radius=3, engine="legacy", backend="shard"))[3]
cat(m, "legacy default-backend ms/sphere", 1000*t0/nv, " shard:",1000*t0s/nv, "\n")
}
# no-op logger in this scratch session only
noop <- function(...) invisible(NULL)
utils::assignInNamespace("flog.debug", noop, "futile.logger")
for (m in c("corclass","sda_notune")) {
ms <- mvpa_model(load_model(m), ds$dataset, ds$design, "classification", crossval=cv)
t1 <- system.time(run_searchlight(ms, radius=3, engine="legacy", backend="default"))[3]
cat(m, "legacy, flog.debug no-op ms/sphere", 1000*t1/nv, "\n")
}
# raw kernel cost: corclass on a 123-voxel sphere, 5 folds
X <- matrix(rnorm(100*123),100); y <- ds$design$y_train; b <- ds$design$block_var
k <- function(){ for (f in 1:5){ tr<-b!=f; M<-rowsum(X[tr,],y[tr])/as.vector(table(y[tr])); Xs<-t(scale(t(X[!tr,]))); Ms<-t(scale(t(M))); p<-tcrossprod(Xs,Ms)}}
cat("raw corclass 5-fold kernel us/sphere:", 1e6*system.time(for(i in 1:500) k())[3]/500, "\n")
k2 <- function(){ for (f in 1:5){ tr<-b!=f; m<-sda::sda(X[tr,],y[tr],verbose=FALSE); predict(m,X[!tr,],verbose=FALSE)}}
cat("raw sda 5-fold us/sphere:", 1e6*system.time(for(i in 1:50) k2())[3]/50, "\n")
