library(futile.logger)
library(neurosurf)
args <- list()
args$config = "~/Dropbox/dual_aud/1009a/config_search.R"

config <- initialize_configuration(args)
config$model <- "sda_notune"
config <- initialize_standard_parameters(config, args, "searchlight")

## Searchlight Specific Params
setArg("niter", config, args, 16)
setArg("radius", config, args, 8)
setArg("type", config, args, "randomized")

config <- initialize_tune_grid(args, config)
config_params <- as.list(config)

config$design <- initialize_design(config)

config$crossval <- initialize_crossval(config, config$design)
surfdat <- load_surface_data(config, "train_data")

lh_dset <- mvpa_surface_dataset(surfdat[[1]], hemisphere="lh")
rh_dset <- mvpa_surface_dataset(surfdat[[2]], hemisphere="rh")

mod1 <- load_mvpa_model(config, lh_dset)
res <- run_searchlight(mod1, method="randomized", radius=16, niter=2)

