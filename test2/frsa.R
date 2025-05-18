library(neuroim2)
devtools::load_all()
dat <- readRDS("test2/test_feature_rds_roi_data.rds")
# Create empty list to store results
all_results <- list()

# Define methods to test
methods <- c("pca", "pls", "scca", "glmnet")

# Loop through different max_comps values from 100 to 5 in decrements of 5
for (max_comps in seq(100, 10, by = -10)) {
  print(paste("Running for max_comps =", max_comps))
  # Create design with current max_comps
  rdes <- feature_rsa_design(F=dat$F, labels=dat$labels, max_comps=max_comps)
  
  # Create cross-validation object
  cval <- blocked_cross_validation(dat$block)
  
  # Loop through each method
  for (method in methods) {
    print(method)
    lambda <-0
    # Create model with current method and max_comps
    if (method == "glmnet") {
      lambda <- 1/(max_comps^2)
      model <- feature_rsa_model(dat$dset, rdes, method=method, crossval=cval, alpha=.9, lambda=lambda)
    } else {
      model <- feature_rsa_model(dat$dset, rdes, method=method, crossval=cval)
    }
    
    # Run regional analysis
    out <- run_regional(model, region_mask=dat$atlas)
    
    # Extract performance table and add method and max_comps columns
    perf_table <- out$performance_table
    perf_table$method <- method
    perf_table$max_comps <- max_comps
   
    # Store in results list
    all_results[[paste(method, max_comps, sep="_")]] <- perf_table
  }
}

# Combine all results into a single tibble
combined_results <- do.call(rbind, all_results)

# Convert to tibble if not already
if (!inherits(combined_results, "tbl_df")) {
  combined_results <- tibble::as_tibble(combined_results)
}

# Convert to long format for easier plotting
plot_data_long <- tidyr::pivot_longer(
  combined_results,
  cols = c("mean_correlation", "cor_difference", "mean_rank_percentile", "voxel_correlation", "mse", "r_squared", 
           "cor_temporal_means", "mean_voxelwise_temporal_cor"),
  names_to = "metric",
  values_to = "value"
)

# Create plots for each metric
p = ggplot2::ggplot(plot_data_long, ggplot2::aes(x = max_comps, y = value, color = method, group = method)) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::facet_wrap(~ metric, scales = "free_y") +
  ggplot2::labs(
    title = "Model Performance by Number of Components",
    x = "Maximum Number of Components",
    y = "Performance Value",
    color = "Method"
  ) +
  ggplot2::theme_minimal() 
 
