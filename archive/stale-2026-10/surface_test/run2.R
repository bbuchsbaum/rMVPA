library(futile.logger)
library(neurosurf)
args <- list()
args$config = "/Users/bbuchsbaum/analysis/mvpa/surface_searchlight2/config_vowel.R"

config <- initialize_configuration(args)
config$model <- "sda_notune"
config$model_type = "classification"
config <- initialize_standard_parameters(config, args, "searchlight")

## Searchlight Specific Params
set_arg("niter", config, args, 2)
set_arg("radius", config, args, 8)
set_arg("type", config, args, "randomized")

config$tune_grid <- initialize_tune_grid(args, config)
config_params <- as.list(config)

design <- initialize_design(config)

crossval <- initialize_crossval(config, design)
surfdat <- initialize_surface_data(config)

mod1 <- load_mvpa_model(config, dataset=surfdat[[1]], design=design, crossval=crossval, feature_selector=NULL)
res <- run_searchlight(mod1, method="randomized", radius=8, niter=1)

