library(futile.logger)
library(neurosurf)
args <- list()
args$config = "~/Dropbox/dual_aud/1009a/config_regional2.R"

config <- initialize_configuration(args)
config$model <- "sda"
config$model_type = "classification"
config <- initialize_standard_parameters(config, args, "regional")


#config$tune_grid <- initialize_tune_grid(args, config)
#config_params <- as.list(config)

design <- initialize_design(config)

crossval <- initialize_crossval(config, design)
surfdat <- initialize_surface_data(config)


mod1 <- load_mvpa_model(config, dataset=surfdat[[1]], design=design, crossval=crossval, feature_selector=NULL)
res <- run_regional(mod1, surfdat[[1]]$mask)

