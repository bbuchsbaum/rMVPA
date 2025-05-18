# rMVPA Roxygen Documentation Plan

## Progress Summary (as of last update)

### Completed Tasks:

1. Updated documentation for `contrast_rsa_model` function:
   - ✅ Updated `output_metric` parameter documentation to show it accepts a character vector
   - ✅ Explained multiple metrics can be requested simultaneously
   - ✅ Documented that duplicates are removed while preserving first-occurrence order

2. Updated documentation for `train_model.contrast_rsa_model`:
   - ✅ Updated return value documentation to explain list structure for multiple metrics
   - ✅ Explained structure and naming conventions for returned elements

3. Fixed dependency handling:
   - ✅ Removed `requireNamespace("crayon")` guard clause in `print.contrast_rsa_model`
   - ✅ Added proper `@importFrom crayon` tags
   - ✅ Added crayon to the DESCRIPTION file's Imports list

### Pending Tasks:

1. Examples:
   - ⬜ Add example for `train_model.contrast_rsa_model()` showing the returned list structure
   - ⬜ Update any other examples that use `output_metric` to reflect it's now a vector

2. S3 methods documentation:
   - ⬜ Ensure all S3 generics have proper documentation (title, description, parameters, return value)
   - ⬜ Link S3 methods to their generics using `@rdname`

3. Dependencies:
   - ⬜ Check for other `requireNamespace()` calls that need fixing
   - ⬜ Review and update DESCRIPTION file for all dependencies

## Package Documentation
- [ ] Add package overview documentation in a separate file (e.g., "R/rMVPA-package.R")
- [ ] Document all imports and dependencies in the package documentation

## Function Documentation

### Primary Focus Areas
1. **contrast_rsa_model** - Document changes to the `output_metric` parameter:
   - [x] Update documentation to show `output_metric` accepts a character vector not just a single string
   - [x] Add examples showing multiple metrics requested simultaneously
   - [x] Update return value documentation to explain how multiple metrics are returned in a list

2. **train_model.contrast_rsa_model**:
   - [x] Update return value documentation to explain list structure for multiple metrics
   - [ ] Add examples showing return structure when multiple metrics are requested

### General Documentation Requirements
- [ ] Ensure all exported functions have:
  - [ ] Title
  - [ ] Description
  - [ ] Parameters fully described
  - [ ] Return value documented
  - [ ] **Runnable** examples (this is critical for CRAN)
  - [ ] Appropriate cross-references to related functions
  
- [ ] Properly document S3 generics and methods:
  - [ ] Document all parameters on the generic function
  - [ ] Use `@rdname` to link methods to their generic
  - [ ] Add method-specific documentation only when necessary

## S3 Methods Documentation

### Key S3 Methods to Prioritize
1. `train_model()`
   - [ ] Document the generic with all parameters fully described
   - [ ] Use `@rdname train_model` for all methods
   - [ ] Add method-specific details for parameters unique to each method

2. `merge_results()`
   - [ ] Document the generic with all parameters
   - [ ] Add method specifics for contrast_rsa_model to handle multiple metrics 

3. `print()` methods
   - [x] Ensure `print.contrast_rsa_model` shows multiple metrics correctly

### Complete List of S3 Generics to Document

1. `get_unique_regions()`
2. `strip_dataset()`
3. `select_features()`
4. `format_result()`
5. `merge_results()`
6. `run_future()`
7. `process_roi()`
8. `train_model()`
9. `y_train()`
10. `y_test()`
11. `test_design()`
12. `fit_model()`
13. `tune_grid()`
14. `has_test_set()`
15. `has_crossval()`
16. `performance()`
17. `compute_performance()`
18. `merge_classif_results()`
19. `get_samples()`
20. `data_sample()`
21. `as_roi()`
22. `get_searchlight()`
23. `wrap_output()`
24. `merge_predictions()`
25. `sub_result()`
26. `nobs()`
27. `prob_observed()`
28. `nresponses()`
29. `predict_model()`
30. `run_searchlight()`
31. `run_regional()`
32. `crossval_samples()`
33. `pairwise_dist()`
34. `filter_roi()`
35. `get_nfolds()`
36. `train_indices()`
37. `balance_partitions()`
38. `print()`

## Documentation Approach for S3 Methods

For each generic above:
1. **On the generic function**:
   - [ ] Document all parameters
   - [ ] Provide a comprehensive description
   - [ ] Document the return value structure
   - [ ] Include a basic example

2. **On each method**:
   - [ ] Use `@rdname generic_name` to link to the generic's documentation
   - [ ] Only document parameters specific to the method
   - [ ] Include method-specific details about return values
   - [ ] Add method-specific examples only if necessary

## Examples Requirements

- [ ] **CRITICAL**: Every exported function MUST have at least one runnable example
- [ ] Examples must NOT use `\dontrun{}` (CRAN will reject)
- [ ] Only use `\donttest{}` for long-running examples
- [ ] Examples must be complete (include necessary setup code)
- [ ] Examples should demonstrate practical usage
- [x] For `contrast_rsa_model`, include an example showing multiple metrics:
  ```r
  model_multi_metric <- contrast_rsa_model(
    dataset = mvpa_dat,
    design = msreve_des,
    output_metric = c("beta_delta", "recon_score", "beta_only")
  )
  ```

## Code Style and Best Practices

- [x] Remove all uses of `requireNamespace()` except for packages in `Suggests` 
- [x] Move package dependencies from `requireNamespace()` checks to proper `Imports` field
- [ ] Ensure consistent parameter naming across related functions
- [ ] Use markdown formatting in roxygen comments where appropriate

## Testing Documentation Code

- [ ] Run `devtools::document()` to generate documentation
- [ ] Check for warnings and fix issues
- [ ] Verify all .Rd files are created correctly
- [ ] Run `R CMD check` to validate documentation against CRAN requirements
