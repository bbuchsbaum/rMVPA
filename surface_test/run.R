args <- list()
args$config = "~/Dropbox/dual_aud/1009a/config_search.R"

config <- initialize_configuration(args)
config <- initialize_standard_parameters(config, args, "searchlight")

## Searchlight Specific Params
setArg("niter", config, args, 16)
setArg("radius", config, args, 8)
setArg("type", config, args, "randomized")

config <- initialize_tune_grid(args, config)
config_params <- as.list(config)

config <- initialize_design(config)