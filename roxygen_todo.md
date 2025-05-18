# rMVPA Roxygen Documentation Plan

## Package-Level Documentation

- [ ] Create package-level documentation using `@docType package`
- [ ] Add package title, description, and authors
- [ ] Include relevant package-level tags: `@import`, `@importFrom`, etc.
- [ ] Include package references and citations

## High Priority Documentation Tasks

### `contrast_rsa_model` Function and Methods

1. **Parameter Documentation Update**
   - [x] Update `@param output_metric` documentation in `contrast_rsa_model()` to clearly state it accepts a character vector of metrics
   - [x] Include the exact list of allowed values: "beta_delta", "beta_only", "delta_only", "recon_score", "beta_delta_norm", "composite"
   - [x] Explain that multiple metrics can be requested simultaneously
   - [x] Document that duplicates are removed while preserving first-occurrence order

2. **Return Value Documentation Update**
   - [x] Update `@return` for `train_model.contrast_rsa_model()` to explain it returns a named list where:
     - [x] Each list element corresponds to a requested metric from `output_metric`
     - [x] For metrics like "beta_delta", "beta_only", "delta_only" - elements are Q-length vectors
     - [x] For metrics like "recon_score", "composite" - elements are single numeric values
     - [x] Attribute "na_reason" may be present if any metric calculation failed

3. **Examples Update**
   - [x] Add example for `contrast_rsa_model()` showing multiple metrics requested:
     ```r
     model_multi_metric <- contrast_rsa_model(
       dataset = mvpa_dat,
       design = msreve_des,
       output_metric = c("beta_delta", "recon_score", "beta_only")
     )
     ```
   - [ ] Add example for `train_model.contrast_rsa_model()` showing the structure of the returned list when multiple metrics are requested
   - [ ] Update any other examples that use `output_metric` to reflect it's now a vector

4. **Print Method Update**
   - [x] Ensure `print.contrast_rsa_model()` correctly shows all requested metrics

## Function Documentation

### Exported Functions
- [ ] Ensure all exported functions (those with `@export` tag) have full documentation
- [ ] Include title, description, parameter details, return value descriptions
- [ ] Add usage examples that run without errors (CRAN requirement)
- [ ] **IMPORTANT**: NEVER use `\dontrun{}` for examples (CRAN will reject these)
- [ ] Only use `\donttest{}` when absolutely necessary (for long-running examples)
- [ ] Ensure all exported functions have at least one runnable example

### Parameters
- [ ] Document all function parameters with `@param`
- [ ] Be consistent in parameter descriptions across similar functions
- [ ] Document inherited parameters in S3 methods appropriately 

### Return Values
- [ ] Document all return values with `@return`
- [ ] For complex return structures, use `\describe{}` or itemize to detail components

## Data Documentation

- [ ] Document all exported datasets using `@format` and `@source`
- [ ] Include dimensions and variable descriptions for data frames
- [ ] Add references to data sources where applicable

## S3/S4 Methods and Classes

- [ ] Document S3 generic functions fully (all parameters, return values)
- [ ] For S3 methods of a generic:
  - [ ] Primary documentation should be on the generic function
  - [ ] Add method-specific documentation only when the method adds parameters or has significantly different behavior
  - [ ] Use `@rdname generic_function` to link method documentation to the generic
  - [ ] Do not repeat parameter documentation that is shared with the generic
- [ ] For class constructors (e.g., `rsa_model()`, `mvpa_dataset()`):
  - [ ] Document all parameters thoroughly
  - [ ] Include complete return value documentation
  - [ ] Add comprehensive, runnable examples
- [ ] Ensure proper S3 method registration with `@export` and appropriate naming (`function.class`)

## Dependencies and requireNamespace()

- [ ] **CRITICAL**: Remove all unnecessary `requireNamespace()` calls:
  - [ ] Keep `requireNamespace()` ONLY for packages listed in `Suggests` field
  - [x] Move essential packages to `Imports` field and remove `requireNamespace()` guard clauses
  - [x] For `rsa_model.R`, ensure dependencies like "crayon" are properly listed in either `Imports` or `Suggests`
  - [ ] For `contrast_rsa_model.R`, ensure packages like "glmnet" (if needed) are properly declared
  
## Vignettes

- [ ] Plan vignettes to demonstrate package functionality
- [ ] Ensure vignettes can be built without errors
- [ ] Consider adding introductory vignette for new users

## NAMESPACE Management

- [ ] Review `@export` tags for consistency
- [ ] Check `@import` vs. `@importFrom` usage (prefer specific imports)
- [ ] Ensure all S3 methods are properly registered

## DESCRIPTION File

- [ ] Ensure all dependencies are properly listed in Imports/Suggests
- [ ] Move unnecessary dependencies to Suggests
- [ ] Add `Collate` field if specific file order is needed
- [ ] Make sure Version, Date, and other fields are up to date

## Code Style and Documentation Readability

- [ ] Use consistent formatting in roxygen comments
- [ ] Keep line lengths reasonable
- [ ] Use markdown formatting where appropriate
- [ ] Check for spelling errors

## Final Checks and Iteration

- [ ] Run `devtools::document()` to generate documentation
- [ ] Review generated .Rd files for issues
- [ ] Run R CMD check to verify documentation completeness
- [ ] Address any documentation warnings or notes

## rMVPA Specific Items

- [ ] Document all model functions (create_, train_, etc.)
- [ ] Document dataset creation and manipulation functions
- [ ] Document analysis utilities and visualization functions
- [ ] Add cross-references between related functions
- [ ] Verify `requireNamespace()` is not used unless the package is in `Suggests` field

## Appendix: S3 Method Families

List of all S3 generic/method families found in the package:

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
38. `print()` (for custom classes)

For each of these generics:
1. Document the generic function fully
2. For methods, use `@rdname` to link to the generic's documentation
3. Only add method-specific documentation if the method has additional parameters or significantly different behavior

## Roxygen Documentation To-Do for CRAN-Positive Status

This list outlines key areas and actionable items to ensure Roxygen documentation meets CRAN standards.

check off completed items when done.

### I. Package-Level Documentation
- [ ] Create/Review main package documentation file (e.g., `packageName-package.R`).
  - [ ] Add `#' @keywords internal` and `"_PACKAGE"` alias if using the common pattern.
  - [ ] Ensure it provides a good overview of the package's purpose and main functionalities.
  - [ ] Link to key functions or vignettes using `\\link{}` or `\\code{\\link{}}`.

### II. Function Documentation (for all EXPORTED functions)
- [ ] **Completeness Check:**
  - [ ] Verify every exported function (generics, standalone functions) has a comprehensive Roxygen block.
  - [ ] For S3/S4 methods, ensure they are correctly associated with their generic for documentation inheritance (see Section IV and Appendix A).
  - [ ] Ensure no undocumented exported functions (`codetools::checkUsageEnv(desc::desc_get_field("Package"))` can help).
- [ ] **Title (`@title`):**
  - [ ] Present for every exported generic function and standalone (non-S3/S4) exported function.
  - [ ] Clear, concise, and in sentence case (capitalize first word, no period at the end unless it's a full sentence).
- [ ] **Description (`@description`):**
  - [ ] Present for every exported generic function and standalone (non-S3/S4) exported function.
  - [ ] Provides a more detailed explanation than the title.
  - [ ] Well-formatted, potentially multiple paragraphs.
- [ ] **Parameters (`@param parameter_name`):**
  - [ ] Document every parameter for generic functions and standalone functions.
  - [ ] For S3/S4 methods: Parameters are generally inherited from the generic. Only add `@param` tags in the method's Roxygen block if:
    - [ ] The method introduces *new* parameters not present in the generic.
    - [ ] The method implements `...` from the generic, and you want to document specific arguments passed through `...` for this particular method.
    - [ ] The behavior of an inherited parameter is *significantly* different or needs clarification for this specific method (use sparingly; prefer documenting general behavior in the generic).
  - [ ] Clearly state the expected type of the parameter (e.g., `numeric vector`, `character string`, `data.frame`).
  - [ ] Describe the purpose and effect of each parameter.
  - [ ] For optional parameters, indicate default values if not obvious from function signature.
- [ ] **Return Value (`@return`):**
  - [ ] Describe the object returned by the generic function or standalone function.
  - [ ] For S3/S4 methods: The return value description is generally inherited. Only add/override `@return` in the method's Roxygen block if the return type or structure is specific to that method and differs from the generic's description.
  - [ ] Specify the class/type of the returned object.
  - [ ] Explain the structure if it's a complex object (e.g., a list with named elements).
- [ ] **Details (`@details`):**
  - [ ] Use for more in-depth explanations, algorithms, or methodological notes if necessary. This is often placed in the generic's documentation.
  - [ ] For S3/S4 methods, `@details` can be used in the method's Roxygen block if there are nuances highly specific to that method's implementation not covered by the generic's details.
- [ ] **Examples (`@examples`):**
  - [ ] **CRITICAL NOTE:** CRAN strictly disallows `@dontrun`. **Never use `@dontrun`**. Use `@donttest` for examples that are computationally intensive, require specific external resources not available on CRAN, or might fail due to stochasticity not controlled by `set.seed`. However, code within `@donttest` blocks **must still be theoretically runnable** (e.g., use loadable data, no placeholder file paths or undefined variables).
  - [ ] Provide runnable examples for *every* exported generic and standalone function.
  - [ ] Examples for generics should ideally demonstrate usage with key S3/S4 methods if applicable.
  - [ ] Examples should be self-contained and demonstrate typical usage.
  - [ ] Use `set.seed()` for examples involving randomness to ensure reproducibility.
  - [ ] Ensure examples run without errors during `R CMD check`.
  - [ ] If example data is needed, consider creating small, self-contained mock data within the example or including minimal datasets in the package (see Data Documentation).
- [ ] **Export Status (`@export` for generics/standalone, `@exportS3Method` or auto-detection for methods):**
  - [ ] Ensure `@export` tag is present for all generic functions and standalone functions intended for users.
  - [ ] For S3 methods, ensure they are correctly registered (e.g., `S3method(generic, class)` in `NAMESPACE`). This is often handled automatically by Roxygen if the generic is exported and methods are named `generic.class`. Explicit `@exportS3Method generic class` can also be used if needed.
  - [ ] Ensure `@noRd` is used for all internal/helper functions not meant for users.
- [ ] **NAMESPACE Management:**
  - [ ] Use `@importFrom package function` for specific function imports.
  - [ ] Alternatively, list essential packages in `DESCRIPTION` Imports field and use `package::function` calls (though `@importFrom` is often preferred for clarity and avoiding masking).
  - [ ] Avoid using `@import package` unless absolutely necessary, as it imports the entire namespace.
  - [ ] Do not use `requireNamespace()` calls to guard code unless the package is in `Suggests` and its absence should not break the function (provide alternative behavior).
- [ ] **Cross-References (`@seealso`):**
  - [ ] Add links to related functions within the package or other relevant documentation using `\\code{\\link{function_name}}` or `\\link[package_name]{function_name}`. Typically placed in the generic's documentation.
- [ ] **Keywords (`@keywords`):**
  - [ ] Add relevant keywords (e.g., `internal`, `datasets`) if applicable, though often less critical if other documentation is thorough.

### III. Data Documentation (if package includes datasets)
- [ ] Ensure every dataset exported by the package has a documentation file (e.g., in `R/data.R`).
- [ ] Use standard Roxygen tags (`@docType data`, `@name data_name`, `@aliases data_name`, `@title`, `@description`, `@format`, `@source`, `@keywords datasets`, `@examples`).
- [ ] `@format`: Describe the structure of the dataset (e.g., "A data frame with X rows and Y variables: ..."). List and describe each variable/column.
- [ ] Ensure `LazyData: true` is in the `DESCRIPTION` file if applicable.

### IV. S3/S4 Methods and Classes Documentation
- [ ] **S3 Methods (See Appendix A for a list of identified S3 families):**
  - [ ] **Generic Function:** Ensure the S3 generic function (e.g., `my_generic <- function(x, ...) UseMethod("my_generic")`) is fully documented with `@title`, `@description`, `@param` for all its arguments (including `...`), `@return`, `@details`, and `@examples` that ideally show usage with important methods.
  - [ ] **S3 Method Files (`R/generic.class.R` or similar):**
    - [ ] Roxygen blocks for S3 methods are often minimal or sometimes not even present if the method introduces no new parameters and its behavior is adequately covered by the generic's documentation.
    - [ ] If an S3 method has parameters *not* in the generic (e.g., if the generic has `...` and the method defines specific arguments passed via `...`), document these new parameters using `@param` in the method's Roxygen block.
    - [ ] If an S3 method has significantly different behavior or return structure for a specific class that warrants detailed explanation beyond the generic's documentation, use `@details` or refine `@return` in the method's Roxygen block.
    - [ ] **Linking to Generic's Documentation:** Ensure S3 methods are linked to the generic's documentation page. This is often automatic if the generic is exported and methods are named `generic.class`. If necessary, use `@rdname generic_function_name` in the method's Roxygen block to explicitly link it to the generic's documentation topic, ensuring all methods appear on the same help page as the generic.
  - [ ] **NAMESPACE:** Ensure Roxygen generates correct `S3method(generic, class)` directives in `NAMESPACE`. This is crucial for method dispatch and is usually handled automatically by `devtools::document()`.
- [ ] **S4 Methods and Classes:**
  - [ ] Document S4 generics similarly to S3 generics.
  - [ ] Document S4 classes using `@slot` for slots.
  - [ ] Document S4 methods, noting that they also inherit documentation but can provide specifics.
  - [ ] Ensure Roxygen generates correct `exportClasses`, `exportMethods` directives in `NAMESPACE`.

### V. Vignettes
- [ ] Ensure at least one vignette demonstrating key package functionality.
- [ ] Vignettes should be listed in the `DESCRIPTION` file under `VignetteBuilder` (e.g., `knitr`) and `Suggests` (e.g., `knitr`, `rmarkdown`).
- [ ] Ensure vignettes build correctly and without errors/warnings during `R CMD check`.
- [ ] Use `\\VignetteIndexEntry{Vignette Title}` in the vignette source.

### VI. NAMESPACE File
- [ ] Primarily, ensure the `NAMESPACE` file is generated by `devtools::document()` (which runs Roxygen) and is not manually edited.
- [ ] Review the generated `NAMESPACE` for correctness (exports, S3 methods, imports).

### VII. DESCRIPTION File
- [ ] **Accuracy and Completeness:**
  - [ ] Title: Sentence case, informative.
  - [ ] Version: Follow `x.y.z` or `x.y.z.9000+` conventions.
  - [ ] Author(s) & Maintainer: Correctly formatted with roles (e.g., `[aut, cre]`).
  - [ ] Description: A more detailed paragraph about the package.
  - [ ] License: CRAN-compatible license (e.g., `GPL-3`, `MIT + file LICENSE`).
- [ ] **Dependencies:**
  - [ ] Imports: List all packages whose functions are essential for your package to work.
  - [ ] Suggests: List packages used in examples, vignettes, or optional functionality.
  - [ ] Depends: Avoid if possible; used for packages whose S4 classes/methods are heavily inherited or if your package relies on a specific version of R.
  - [ ] Specify version requirements if needed (e.g., `dplyr (>= 1.0.0)`).
- [ ] `URL`: Link to package repository (e.g., GitHub).
- [ ] `BugReports`: Link to issue tracker.
- [ ] `Encoding: UTF-8` is recommended.

### VIII. Code and Documentation Style
- [ ] **Readability:** Ensure Roxygen comments are well-formatted and easy to read.
- [ ] **Line Lengths:** Keep Roxygen comment lines to a reasonable length (e.g., <80 characters).
- [ ] **Markdown:** Use markdown for formatting within Roxygen tags where appropriate (e.g., lists in `@description` or `@details`).

### IX. Final Checks and Iteration
- [ ] **Run `devtools::document()`:** Regenerate documentation and `NAMESPACE` frequently.
- [ ] **Run `devtools::check(cran = TRUE)` or `R CMD check --as-cran yourpackage_version.tar.gz`:**
  - [ ] Address ALL errors.
  - [ ] Address ALL warnings.
  - [ ] Address ALL notes (provide explanations for unavoidable notes if necessary).
- [ ] **Spell Check:** Use a spell checker for documentation text.
- [ ] **Proofread:** Read through the generated PDF manual and HTML help pages for clarity, correctness, and formatting issues.
- [ ] **Test on Multiple Platforms:** If possible, test building and checking on different operating systems (Windows, macOS, Linux), especially if system-dependent code or compiled code is involved. Consider using `rhub::check_for_cran()`.

### X. Specific rMVPA Items (Based on Past Interactions & Custom Instructions)
- [ ] Ensure all examples for exported functions are runnable and use `@donttest` appropriately (see note in Section II about `@dontrun` vs `@donttest`).
- [ ] Verify `requireNamespace()` is not used unless the package is in `Suggests` field