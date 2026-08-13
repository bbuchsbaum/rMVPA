# Package index

## REMAP Diagnostics

Summaries that quantify how far memory is from perception and how large
the learned correction is.

- [`summarize_remap_roi()`](https://bbuchsbaum.github.io/rMVPA/reference/summarize_remap_roi.md)
  : Summarize REMAP diagnostics at the ROI level
- [`summarize_remap_items()`](https://bbuchsbaum.github.io/rMVPA/reference/summarize_remap_items.md)
  : Summarize REMAP item-level residuals for an ROI

## All Functions

All functions exported by rMVPA

- [`MVPAModels`](https://bbuchsbaum.github.io/rMVPA/reference/MVPAModels.md)
  : Pre-defined MVPA Classification Models

- [`add_interaction_contrasts()`](https://bbuchsbaum.github.io/rMVPA/reference/add_interaction_contrasts.md)
  : Add Interaction Contrasts to an msreve_design

- [`as_roi()`](https://bbuchsbaum.github.io/rMVPA/reference/as_roi.md) :
  Convert object to ROI

- [`average_labels()`](https://bbuchsbaum.github.io/rMVPA/reference/average_labels.md)
  : Average NeuroVec Data by Labels

- [`balance_partitions()`](https://bbuchsbaum.github.io/rMVPA/reference/balance_partitions.md)
  : Balance Cross-Validation Partitions

- [`banded_ridge_da()`](https://bbuchsbaum.github.io/rMVPA/reference/banded_ridge_da.md)
  [`grouped_ridge_da()`](https://bbuchsbaum.github.io/rMVPA/reference/banded_ridge_da.md)
  : Convenience wrapper: build a grouped-ridge domain-adaptation model
  from matrices

- [`banded_ridge_da_model()`](https://bbuchsbaum.github.io/rMVPA/reference/banded_ridge_da_model.md)
  [`grouped_ridge_da_model()`](https://bbuchsbaum.github.io/rMVPA/reference/banded_ridge_da_model.md)
  : Grouped (banded) ridge domain-adaptation model (continuous
  predictors -\> brain)

- [`binary_classification_result()`](https://bbuchsbaum.github.io/rMVPA/reference/binary_classification_result.md)
  : Classification results for binary outcome

- [`blocks()`](https://bbuchsbaum.github.io/rMVPA/reference/blocks.md) :
  Define consecutive column blocks as feature sets

- [`build_analysis()`](https://bbuchsbaum.github.io/rMVPA/reference/build_analysis.md)
  : Build an Analysis Object

- [`by_set()`](https://bbuchsbaum.github.io/rMVPA/reference/by_set.md) :
  Define feature sets by a per-column set label

- [`category_rdm()`](https://bbuchsbaum.github.io/rMVPA/reference/category_rdm.md)
  : Create Hypothesis RDM from Category Structure

- [`classification_result()`](https://bbuchsbaum.github.io/rMVPA/reference/classification_result.md)
  :

  Create a `classification_result` instance

- [`combine_prediction_tables()`](https://bbuchsbaum.github.io/rMVPA/reference/combine_prediction_tables.md)
  : Combine prediction tables

- [`compute_local_redundancy()`](https://bbuchsbaum.github.io/rMVPA/reference/compute_local_redundancy.md)
  : Compute Local Redundancy for Each Searchlight Center

- [`contrast_rsa_model()`](https://bbuchsbaum.github.io/rMVPA/reference/contrast_rsa_model.md)
  : Constructor for contrast_rsa_model

- [`contrasts()`](https://bbuchsbaum.github.io/rMVPA/reference/contrasts.md)
  : Generate Contrast Matrices

- [`create_dist()`](https://bbuchsbaum.github.io/rMVPA/reference/create_dist.md)
  : Create a Distance Function Object

- [`create_model_spec()`](https://bbuchsbaum.github.io/rMVPA/reference/create_model_spec.md)
  : Create a Generic Model Specification

- [`bootstrap_blocked_cross_validation()`](https://bbuchsbaum.github.io/rMVPA/reference/cross_validation.md)
  [`blocked_cross_validation()`](https://bbuchsbaum.github.io/rMVPA/reference/cross_validation.md)
  [`sequential_blocked_cross_validation()`](https://bbuchsbaum.github.io/rMVPA/reference/cross_validation.md)
  [`custom_cross_validation()`](https://bbuchsbaum.github.io/rMVPA/reference/cross_validation.md)
  : bootstrap_blocked_cross_validation

- [`crossv_block()`](https://bbuchsbaum.github.io/rMVPA/reference/crossv_block.md)
  : Block Cross-Validation Data Preparation

- [`crossv_k()`](https://bbuchsbaum.github.io/rMVPA/reference/crossv_k.md)
  : K-fold Cross-Validation Data Preparation

- [`crossval_samples()`](https://bbuchsbaum.github.io/rMVPA/reference/crossval_samples.md)
  : Cross-validation samples

- [`custom_performance()`](https://bbuchsbaum.github.io/rMVPA/reference/custom_performance.md)
  : Apply Custom Performance Metric to Prediction Result

- [`cv_evaluate_roi()`](https://bbuchsbaum.github.io/rMVPA/reference/cv_evaluate_roi.md)
  : Evaluate One ROI with rMVPA Cross-Validation Helpers

- [`data_sample()`](https://bbuchsbaum.github.io/rMVPA/reference/data_sample.md)
  : Extract Sample from Dataset

- [`cordist()`](https://bbuchsbaum.github.io/rMVPA/reference/distance-constructors.md)
  [`mahadist()`](https://bbuchsbaum.github.io/rMVPA/reference/distance-constructors.md)
  [`eucdist()`](https://bbuchsbaum.github.io/rMVPA/reference/distance-constructors.md)
  [`euclidean()`](https://bbuchsbaum.github.io/rMVPA/reference/distance-constructors.md)
  [`robustmahadist()`](https://bbuchsbaum.github.io/rMVPA/reference/distance-constructors.md)
  [`pcadist()`](https://bbuchsbaum.github.io/rMVPA/reference/distance-constructors.md)
  : Distance Function Constructors

- [`era_partition_model()`](https://bbuchsbaum.github.io/rMVPA/reference/era_partition_model.md)
  : ERA Variance-Partition Model

- [`era_rsa_design()`](https://bbuchsbaum.github.io/rMVPA/reference/era_rsa_design.md)
  : Build item-level ERA-RSA confounds and summaries from mvpa_design

- [`era_rsa_model()`](https://bbuchsbaum.github.io/rMVPA/reference/era_rsa_model.md)
  : ERA-RSA: Encoding-Retrieval Similarity and ER Geometry

- [`evaluate_model.feature_rsa_model()`](https://bbuchsbaum.github.io/rMVPA/reference/evaluate_model.feature_rsa_model.md)
  : Evaluate model performance for feature RSA

- [`evaluate_model.vector_rsa_model()`](https://bbuchsbaum.github.io/rMVPA/reference/evaluate_model.vector_rsa_model.md)
  : Evaluate model performance for vector RSA

- [`expected_features()`](https://bbuchsbaum.github.io/rMVPA/reference/expected_features.md)
  : Build expected-domain features from a soft alignment matrix

- [`explain_searchlight_engine()`](https://bbuchsbaum.github.io/rMVPA/reference/explain_searchlight_engine.md)
  : Explain Searchlight Engine Selection

- [`extract_weights()`](https://bbuchsbaum.github.io/rMVPA/reference/extract_weights.md)
  : Extract Raw Model Weights

- [`feature_rsa_connectivity()`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_connectivity.md)
  : Compute ROI-by-ROI Representational Connectivity from Feature RSA
  Predictions

- [`feature_rsa_cross_connectivity()`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_cross_connectivity.md)
  : Compute Cross-Connectivity: Predicted-Observed ROI x ROI Matrix

- [`feature_rsa_da_model()`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_da_model.md)
  : Feature-RSA Domain Adaptation Model

- [`feature_rsa_design()`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_design.md)
  : Create a Feature-Based RSA Design

- [`feature_rsa_model()`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_model.md)
  : Create a Feature-Based RSA Model

- [`feature_rsa_rdm_vectors()`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_rdm_vectors.md)
  : Extract Per-ROI Predicted and Observed RDM Vectors from Feature RSA
  Results

- [`feature_selection`](https://bbuchsbaum.github.io/rMVPA/reference/feature_selection.md)
  : Feature Selection Methods

- [`feature_selector()`](https://bbuchsbaum.github.io/rMVPA/reference/feature_selector.md)
  : Create a feature selection specification

- [`feature_sets()`](https://bbuchsbaum.github.io/rMVPA/reference/feature_sets.md)
  : Feature sets: grouped predictor matrices

- [`feature_sets_design()`](https://bbuchsbaum.github.io/rMVPA/reference/feature_sets_design.md)
  : Feature-sets design (mvpa_design extension for continuous
  regression)

- [`fit_model()`](https://bbuchsbaum.github.io/rMVPA/reference/fit_model-methods.md)
  : Fit Model

- [`fit_roi()`](https://bbuchsbaum.github.io/rMVPA/reference/fit_roi.md)
  : Fit a Model on a Single ROI

- [`fit_roi(`*`<contrast_rsa_model>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/fit_roi.contrast_rsa_model.md)
  : Fit ROI for contrast_rsa_model

- [`format_result()`](https://bbuchsbaum.github.io/rMVPA/reference/format_result.md)
  : Format Result Object

- [`gen_clustered_sample_dataset()`](https://bbuchsbaum.github.io/rMVPA/reference/gen_clustered_sample_dataset.md)
  : Generate a Synthetic Clustered MVPA Dataset

- [`gen_sample_dataset()`](https://bbuchsbaum.github.io/rMVPA/reference/gen_sample_dataset.md)
  : Generate Sample Dataset for MVPA Analysis

- [`generate_folds()`](https://bbuchsbaum.github.io/rMVPA/reference/generate_folds.md)
  : Generate Cross-Validation Folds

- [`get_center_ids()`](https://bbuchsbaum.github.io/rMVPA/reference/get_center_ids.md)
  : Get Center IDs for Searchlight Iteration

- [`get_feature_matrix()`](https://bbuchsbaum.github.io/rMVPA/reference/get_feature_matrix.md)
  : Extract Full Feature Matrix from a Dataset

- [`get_nfolds()`](https://bbuchsbaum.github.io/rMVPA/reference/get_nfolds.md)
  : Get the Number of Folds

- [`get_samples()`](https://bbuchsbaum.github.io/rMVPA/reference/get_samples.md)
  : Get Multiple Data Samples

- [`get_searchlight()`](https://bbuchsbaum.github.io/rMVPA/reference/get_searchlight.md)
  : Generate Searchlight Iterator

- [`global_mvpa_result()`](https://bbuchsbaum.github.io/rMVPA/reference/global_mvpa_result.md)
  : Construct a Global MVPA Result

- [`group_means()`](https://bbuchsbaum.github.io/rMVPA/reference/group_means.md)
  : Compute Group Means of a Matrix

- [`has_crossval()`](https://bbuchsbaum.github.io/rMVPA/reference/has_crossval-methods.md)
  : Cross-Validation Availability

- [`has_test_set()`](https://bbuchsbaum.github.io/rMVPA/reference/has_test_set-methods.md)
  : Test Set Availability

- [`haufe_importance()`](https://bbuchsbaum.github.io/rMVPA/reference/haufe_importance.md)
  : Haufe Feature Importance (Activation Patterns)

- [`install_cli()`](https://bbuchsbaum.github.io/rMVPA/reference/install_cli.md)
  : Install rMVPA Command-Line Wrappers

- [`item_design()`](https://bbuchsbaum.github.io/rMVPA/reference/item_design.md)
  : Construct an ITEM design for trial-wise decoding

- [`item_model()`](https://bbuchsbaum.github.io/rMVPA/reference/item_model.md)
  : ITEM decoding model for ROI/searchlight analysis

- [`kfold_cross_validation()`](https://bbuchsbaum.github.io/rMVPA/reference/kfold_cross_validation.md)
  : kfold_cross_validation

- [`load_model()`](https://bbuchsbaum.github.io/rMVPA/reference/load_model.md)
  : Load a Pre-defined MVPA Model

- [`make_feature_contrasts()`](https://bbuchsbaum.github.io/rMVPA/reference/make_feature_contrasts.md)
  : Generate Contrasts from a Feature Matrix (Optional PCA)

- [`manova_design()`](https://bbuchsbaum.github.io/rMVPA/reference/manova_design.md)
  : Create a MANOVA Design

- [`manova_model()`](https://bbuchsbaum.github.io/rMVPA/reference/manova_model.md)
  : Create a MANOVA Model

- [`merge_classif_results()`](https://bbuchsbaum.github.io/rMVPA/reference/merge_classif_results.md)
  : Merge Multiple Classification/Regression Results

- [`merge_predictions()`](https://bbuchsbaum.github.io/rMVPA/reference/merge_predictions.md)
  : Merge Predictions

- [`merge_results(`*`<manova_model>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/merge_results-methods.md)
  [`merge_results(`*`<mvpa_model>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/merge_results-methods.md)
  [`merge_results(`*`<binary_classification_result>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/merge_results-methods.md)
  [`merge_results(`*`<regression_result>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/merge_results-methods.md)
  [`merge_results(`*`<multiway_classification_result>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/merge_results-methods.md)
  [`merge_results(`*`<regional_mvpa_result>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/merge_results-methods.md)
  [`merge_results(`*`<rsa_model>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/merge_results-methods.md)
  [`merge_results(`*`<vector_rsa_model>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/merge_results-methods.md)
  [`merge_results(`*`<feature_rsa_model>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/merge_results-methods.md)
  : Merge Results for MANOVA Model

- [`merge_results()`](https://bbuchsbaum.github.io/rMVPA/reference/merge_results.md)
  : Merge Multiple Results

- [`mock_context()`](https://bbuchsbaum.github.io/rMVPA/reference/mock_context.md)
  : Build a Mock fit_roi Context

- [`mock_roi_data()`](https://bbuchsbaum.github.io/rMVPA/reference/mock_roi_data.md)
  : Build Mock ROI Data for Plugin Tests

- [`model_importance()`](https://bbuchsbaum.github.io/rMVPA/reference/model_importance.md)
  : Per-Feature Model Importance

- [`model_space_connectivity()`](https://bbuchsbaum.github.io/rMVPA/reference/model_space_connectivity.md)
  : Model-space representational connectivity from a fitted rMVPA result

- [`msreve_design()`](https://bbuchsbaum.github.io/rMVPA/reference/msreve_design.md)
  : Constructor for msreve_design

- [`msreve_temporal_confounds()`](https://bbuchsbaum.github.io/rMVPA/reference/msreve_temporal_confounds.md)
  : MS-ReVE: build temporal nuisance RDMs from a spec

- [`multiway_classification_result()`](https://bbuchsbaum.github.io/rMVPA/reference/multiway_classification_result.md)
  : Create a Multiway Classification Result Object

- [`mvpa_clustered_dataset()`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_clustered_dataset.md)
  : Create an MVPA Dataset for Clustered/Parcellated Data

- [`mvpa_config()`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_config.md)
  : Build a Public Analysis Configuration

- [`mvpa_dataset()`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_dataset.md)
  : Create an MVPA Dataset Object

- [`mvpa_design()`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_design.md)
  : Create an MVPA Design Object

- [`mvpa_iterate()`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_iterate.md)
  : Iterate MVPA Analysis Over Multiple ROIs

- [`mvpa_model()`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_model.md)
  : Create an MVPA Model

- [`mvpa_multibasis_dataset()`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_multibasis_dataset.md)
  : Create a Multibasis MVPA Image Dataset

- [`mvpa_multibasis_image_dataset()`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_multibasis_image_dataset.md)
  :

  Alias for `mvpa_multibasis_dataset`

- [`mvpa_surface_dataset()`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_surface_dataset.md)
  : Create a Surface-Based MVPA Dataset Object

- [`mvpa_sysinfo()`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_sysinfo.md)
  : Report System and Package Information for rMVPA

- [`naive_xdec_model()`](https://bbuchsbaum.github.io/rMVPA/reference/naive_xdec_model.md)
  : Naive Cross-Decoding (correlation to source prototypes)

- [`new-analysis-overview`](https://bbuchsbaum.github.io/rMVPA/reference/new-analysis-overview.md)
  : Extending rMVPA: Plugin Analysis Models

- [`nmf_preprocess_maps()`](https://bbuchsbaum.github.io/rMVPA/reference/nmf_preprocess_maps.md)
  : Preprocess Maps for Spatial NMF

- [`nobs()`](https://bbuchsbaum.github.io/rMVPA/reference/nobs.md) : Get
  Number of Observations

- [`nresponses()`](https://bbuchsbaum.github.io/rMVPA/reference/nresponses.md)
  : Number of Response Categories

- [`orthogonalize_contrasts()`](https://bbuchsbaum.github.io/rMVPA/reference/orthogonalize_contrasts.md)
  : Orthogonalize a Contrast Matrix

- [`output_schema(`*`<contrast_rsa_model>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/output_schema.contrast_rsa_model.md)
  : Output schema for contrast_rsa_model

- [`pair_rsa_design()`](https://bbuchsbaum.github.io/rMVPA/reference/pair_rsa_design.md)
  : Construct a pair-observation RSA design

- [`performance()`](https://bbuchsbaum.github.io/rMVPA/reference/performance-methods.md)
  : Compute Performance Metrics

- [`performance(`*`<regression_result>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/performance.regression_result.md)
  : Calculate Performance Metrics for Regression Result

- [`permutation_control()`](https://bbuchsbaum.github.io/rMVPA/reference/permutation_control.md)
  : Create a Permutation Control Object

- [`permute_labels()`](https://bbuchsbaum.github.io/rMVPA/reference/permute_labels.md)
  : Permute Training Labels in an MVPA Design

- [`predict_model()`](https://bbuchsbaum.github.io/rMVPA/reference/predict_model.md)
  : Predict Model Output

- [`predicted_class()`](https://bbuchsbaum.github.io/rMVPA/reference/predicted_class.md)
  : Calculate the Predicted Class from Probability Matrix

- [`prep_regional()`](https://bbuchsbaum.github.io/rMVPA/reference/prep_regional.md)
  : Prepare regional data for MVPA analysis

- [`print(`*`<feature_rsa_design>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/print.feature_rsa_design.md)
  : Print Method for Feature RSA Design

- [`print(`*`<feature_rsa_model>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/print.feature_rsa_model.md)
  : Print Method for Feature RSA Model

- [`print(`*`<feature_sets>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/print.feature_sets.md)
  : Print method for feature_sets

- [`print(`*`<feature_sets_design>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/print.feature_sets_design.md)
  : Print method for feature_sets_design

- [`print(`*`<vector_rsa_design>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/print.vector_rsa_design.md)
  : Print Method for vector_rsa_design

- [`print(`*`<vector_rsa_model>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/print.vector_rsa_model.md)
  : Print Method for vector_rsa_model

- [`prob_observed()`](https://bbuchsbaum.github.io/rMVPA/reference/prob_observed.md)
  : Probability of Observed Class

- [`process_roi()`](https://bbuchsbaum.github.io/rMVPA/reference/process_roi-methods.md)
  : Process ROI

- [`rdm_decorrelate()`](https://bbuchsbaum.github.io/rMVPA/reference/rdm_decorrelate.md)
  : Decorrelation of correlated model RDMs

- [`rdm_model_space_connectivity()`](https://bbuchsbaum.github.io/rMVPA/reference/rdm_model_space_connectivity.md)
  : ROI-by-ROI similarity through a correlated model-RDM space

- [`region_importance()`](https://bbuchsbaum.github.io/rMVPA/reference/region_importance.md)
  : Region Importance via Random Subset Comparison

- [`region_importance_result()`](https://bbuchsbaum.github.io/rMVPA/reference/region_importance_result.md)
  : Construct a Region Importance Result

- [`regional_mvpa_result()`](https://bbuchsbaum.github.io/rMVPA/reference/regional_mvpa_result.md)
  :

  Create a `regional_mvpa_result` instance

- [`register_mvpa_model()`](https://bbuchsbaum.github.io/rMVPA/reference/register_mvpa_model.md)
  : Register a Custom MVPA Model

- [`regression_result()`](https://bbuchsbaum.github.io/rMVPA/reference/regression_result.md)
  : Create a Regression Result Object

- [`remap_rrr_model()`](https://bbuchsbaum.github.io/rMVPA/reference/remap_rrr_model.md)
  : REMAP-RRR: Residual low-rank, domain-adaptive cross-decoding

- [`repmap_design()`](https://bbuchsbaum.github.io/rMVPA/reference/repmap_design.md)
  : Representational mapping design helper (ReNA-Map)

- [`repmap_model()`](https://bbuchsbaum.github.io/rMVPA/reference/repmap_model.md)
  : Representational mapping model (ReNA-Map)

- [`repmed_design()`](https://bbuchsbaum.github.io/rMVPA/reference/repmed_design.md)
  : Representational mediation design helper (ReNA-RM)

- [`repmed_model()`](https://bbuchsbaum.github.io/rMVPA/reference/repmed_model.md)
  : Representational mediation model (ReNA-RM)

- [`repnet_design()`](https://bbuchsbaum.github.io/rMVPA/reference/repnet_design.md)
  : Representational connectivity design helper (ReNA-RC)

- [`repnet_model()`](https://bbuchsbaum.github.io/rMVPA/reference/repnet_model.md)
  : Representational connectivity model (ReNA-RC)

- [`rmvpa_api_lifecycle()`](https://bbuchsbaum.github.io/rMVPA/reference/rmvpa_api_lifecycle.md)
  : List the rMVPA API Lifecycle Registry

- [`rmvpa_stable_api()`](https://bbuchsbaum.github.io/rMVPA/reference/rmvpa_stable_api.md)
  : Return the Stable rMVPA API Surface

- [`roi_result()`](https://bbuchsbaum.github.io/rMVPA/reference/roi_result.md)
  : Construct a Standardized ROI Result

- [`rsa_design()`](https://bbuchsbaum.github.io/rMVPA/reference/rsa_design.md)
  : Construct a design for an RSA (Representational Similarity Analysis)
  model

- [`rsa_model()`](https://bbuchsbaum.github.io/rMVPA/reference/rsa_model.md)
  : Construct an RSA (Representational Similarity Analysis) model

- [`run_analysis()`](https://bbuchsbaum.github.io/rMVPA/reference/run_analysis.md)
  : Run a Built Analysis

- [`run_custom_regional()`](https://bbuchsbaum.github.io/rMVPA/reference/run_custom_regional.md)
  : Run a Custom Analysis Function Regionally

- [`run_custom_searchlight()`](https://bbuchsbaum.github.io/rMVPA/reference/run_custom_searchlight.md)
  : Run a Custom Analysis Function in a Searchlight

- [`run_global()`](https://bbuchsbaum.github.io/rMVPA/reference/run_global.md)
  : Run Global (Whole-Brain) MVPA Analysis

- [`run_permutation_searchlight()`](https://bbuchsbaum.github.io/rMVPA/reference/run_permutation_searchlight.md)
  : Run Permutation Searchlight Inference

- [`run_regional()`](https://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md)
  [`run_regional_base()`](https://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md)
  : Region of Interest Based MVPA Analysis

- [`run_searchlight()`](https://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md)
  : Run Searchlight Analysis

- [`run_searchlight(`*`<default>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.default.md)
  : Default method for run_searchlight

- [`save_results()`](https://bbuchsbaum.github.io/rMVPA/reference/save_results.md)
  : Save MVPA Results to Disk

- [`save_results(`*`<rmvpa_analysis_run>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/save_results.rmvpa_analysis_run.md)
  : Save a High-Level Analysis Run

- [`schema_scalar()`](https://bbuchsbaum.github.io/rMVPA/reference/schema_metric_spec.md)
  [`schema_vector()`](https://bbuchsbaum.github.io/rMVPA/reference/schema_metric_spec.md)
  : Metric Schema Constructors

- [`searchlight_engines()`](https://bbuchsbaum.github.io/rMVPA/reference/searchlight_engines.md)
  : Summarize Available Searchlight Engines

- [`searchlight_mode()`](https://bbuchsbaum.github.io/rMVPA/reference/searchlight_mode.md)
  : Get Or Set Searchlight Mode

- [`second_order_similarity()`](https://bbuchsbaum.github.io/rMVPA/reference/second_order_similarity.md)
  : Compute Second-Order Similarity Scores

- [`select_features()`](https://bbuchsbaum.github.io/rMVPA/reference/select_features-methods.md)
  : Select Features

- [`set_log_level()`](https://bbuchsbaum.github.io/rMVPA/reference/set_log_level.md)
  : Set rMVPA Logging Level

- [`shard_backend`](https://bbuchsbaum.github.io/rMVPA/reference/shard_backend.md)
  : Experimental Shared-Memory Backend for Parallel MVPA

- [`shard_cleanup()`](https://bbuchsbaum.github.io/rMVPA/reference/shard_cleanup.md)
  : Clean Up Shared-Memory Segments

- [`spatial_nmf()`](https://bbuchsbaum.github.io/rMVPA/reference/spatial_nmf.md)
  : Spatial Non-negative Matrix Factorization

- [`spatial_nmf_component_test()`](https://bbuchsbaum.github.io/rMVPA/reference/spatial_nmf_component_test.md)
  : Component-level Inference for Spatial NMF

- [`spatial_nmf_global_test()`](https://bbuchsbaum.github.io/rMVPA/reference/spatial_nmf_global_test.md)
  : Global Cross-validated Group Test for Spatial NMF

- [`spatial_nmf_maps()`](https://bbuchsbaum.github.io/rMVPA/reference/spatial_nmf_maps.md)
  : Spatial NMF on Map Lists

- [`spatial_nmf_stability()`](https://bbuchsbaum.github.io/rMVPA/reference/spatial_nmf_stability.md)
  : Bootstrap Stability for Spatial NMF Components

- [`spatial_nmf_voxelwise_stats()`](https://bbuchsbaum.github.io/rMVPA/reference/spatial_nmf_voxelwise_stats.md)
  : Voxelwise Statistics from Spatial NMF Stability

- [`strip_dataset()`](https://bbuchsbaum.github.io/rMVPA/reference/strip_dataset-methods.md)
  : Strip Dataset from Model Specification

- [`sub_result(`*`<multiway_classification_result>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/sub_result-methods.md)
  [`sub_result(`*`<binary_classification_result>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/sub_result-methods.md)
  : Subset Multiway Classification Result

- [`sub_result()`](https://bbuchsbaum.github.io/rMVPA/reference/sub_result.md)
  : Extract Row-wise Subset of a Result

- [`subsample_centers()`](https://bbuchsbaum.github.io/rMVPA/reference/subsample_centers.md)
  : Subsample Searchlight Centers

- [`subspace_alignment_model()`](https://bbuchsbaum.github.io/rMVPA/reference/subspace_alignment_model.md)
  : Subspace Alignment cross-decoder

- [`summarize_remap_items()`](https://bbuchsbaum.github.io/rMVPA/reference/summarize_remap_items.md)
  : Summarize REMAP item-level residuals for an ROI

- [`summarize_remap_roi()`](https://bbuchsbaum.github.io/rMVPA/reference/summarize_remap_roi.md)
  : Summarize REMAP diagnostics at the ROI level

- [`summary(`*`<feature_rsa_model>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/summary.feature_rsa_model.md)
  : Summary Method for Feature RSA Model

- [`table_to_rdm()`](https://bbuchsbaum.github.io/rMVPA/reference/table_to_rdm.md)
  : Convert Similarity Table to RDM

- [`temporal()`](https://bbuchsbaum.github.io/rMVPA/reference/temporal.md)
  : Temporal RDM wrapper for formula usage

- [`temporal_confounds()`](https://bbuchsbaum.github.io/rMVPA/reference/temporal_confounds.md)
  : Build multiple temporal confounds from a spec

- [`temporal_from_onsets()`](https://bbuchsbaum.github.io/rMVPA/reference/temporal_from_onsets.md)
  : Temporal RDM from onsets (sugar)

- [`temporal_hrf_overlap()`](https://bbuchsbaum.github.io/rMVPA/reference/temporal_hrf_overlap.md)
  : HRF-overlap temporal confound RDM

- [`temporal_nuisance_for_msreve()`](https://bbuchsbaum.github.io/rMVPA/reference/temporal_nuisance_for_msreve.md)
  : Temporal nuisance RDM at condition level (for MS-ReVE)

- [`temporal_rdm()`](https://bbuchsbaum.github.io/rMVPA/reference/temporal_rdm.md)
  : Temporal/ordinal nuisance RDM (trial-level)

- [`test_design()`](https://bbuchsbaum.github.io/rMVPA/reference/test_design-methods.md)
  : Test Design Extraction

- [`train_indices()`](https://bbuchsbaum.github.io/rMVPA/reference/train_indices.md)
  : Get Training Indices for a Fold

- [`transform_contrasts()`](https://bbuchsbaum.github.io/rMVPA/reference/transform_contrasts.md)
  : Apply Transformations to an Existing Contrast Matrix

- [`tune_grid()`](https://bbuchsbaum.github.io/rMVPA/reference/tune_grid-methods.md)
  : Extract Tuning Grid

- [`twofold_blocked_cross_validation()`](https://bbuchsbaum.github.io/rMVPA/reference/twofold_blocked_cross_validation.md)
  : twofold_blocked_cross_validation

- [`use_shard()`](https://bbuchsbaum.github.io/rMVPA/reference/use_shard.md)
  : Enable Shared-Memory Backend for an MVPA Model Spec

- [`validate_analysis()`](https://bbuchsbaum.github.io/rMVPA/reference/validate_analysis.md)
  : Validate an MVPA Analysis for Common Methodological Issues

- [`validate_cutoff()`](https://bbuchsbaum.github.io/rMVPA/reference/validate_cutoff.md)
  : Validate Feature-Selection Cutoff

- [`validate_model_spec()`](https://bbuchsbaum.github.io/rMVPA/reference/validate_model_spec.md)
  [`print(`*`<model_spec_validation_result>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/validate_model_spec.md)
  : Validate a Model Specification for Plugin Readiness

- [`validate_plugin_model()`](https://bbuchsbaum.github.io/rMVPA/reference/validate_plugin_model.md)
  [`print(`*`<plugin_validation_result>`*`)`](https://bbuchsbaum.github.io/rMVPA/reference/validate_plugin_model.md)
  : Validate a Plugin Model Contract

- [`vector_rsa_design()`](https://bbuchsbaum.github.io/rMVPA/reference/vector_rsa_design.md)
  : Construct a design for a vectorized RSA model

- [`vector_rsa_model()`](https://bbuchsbaum.github.io/rMVPA/reference/vector_rsa_model.md)
  : Create a vectorized RSA model

- [`y_test()`](https://bbuchsbaum.github.io/rMVPA/reference/y_test-methods.md)
  : Test Labels/Response Extraction

- [`y_train()`](https://bbuchsbaum.github.io/rMVPA/reference/y_train-methods.md)
  : Training Labels/Response Extraction
