library(neuroim2)
library(rMVPA)
##library(bruceR)
# creat a sample dataset
dataset <- gen_sample_dataset(D=c(6,6,6), nobs = 80, blocks = 4, nlevels = 2)
print(dataset)
# indicate the design
block <- dataset$design$block_var
crossval <- blocked_cross_validation(block)
crossval
# ready to run MVPA
svm_model <- load_model('svmLinear')
model <- mvpa_model(model = svm_model, dataset = dataset$dataset, design = dataset$design, crossval = crossval)
model
# run MVPA in a randomized searchilight analysis
result <- run_searchlight(model, radius = 4, method = 'randomized', niter = 2)

###INFO [2021-07-31 10:55:17] 
# model is: svmLinear
# INFO [2021-07-31 10:55:17] running randomized searchlight with 4 radius and 2 iterations
# Error in `$<-`(`*tmp*`, obsLevels, value = levels(y)) : 
#   no method for assigning subsets of this S4 class
# Error in `$<-`(`*tmp*`, obsLevels, value = levels(y)) : 
#   no method for assigning subsets of this S4 class
# Error in `$<-`(`*tmp*`, obsLevels, value = levels(y)) : 
#   no method for assigning subsets of this S4 class
# Error in `$<-`(`*tmp*`, obsLevels, value = levels(y)) : 
#   no method for assigning subsets of this S4 class
# Error in `$<-`(`*tmp*`, obsLevels, value = levels(y)) : 
#   no method for assigning subsets of this S4 class
# Error in `$<-`(`*tmp*`, obsLevels, value = levels(y)) : 
#   no method for assigning subsets of this S4 class
# Error in `$<-`(`*tmp*`, obsLevels, value = levels(y)) : 
#   no method for assigning subsets of this S4 class
# Error in `$<-`(`*tmp*`, obsLevels, value = levels(y)) : 
#   no method for assigning subsets of this S4 class
# Error in train_model.mvpa_model(mspec, tibble::as_tibble(train), ytrain,  : 
#   training data must have more than one valid feature
# Error in train_model.mvpa_model(mspec, tibble::as_tibble(train), ytrain,  : 
#   training data must have more than one valid feature
# Error in train_model.mvpa_model(mspec, tibble::as_tibble(train), ytrain,  : 
#   training data must have more than one valid feature
# Error in train_model.mvpa_model(mspec, tibble::as_tibble(train), ytrain,  : 
#   training data must have more than one valid feature
# INFO [2021-07-31 10:55:17] searchlight iteration: 1
# WARN [2021-07-31 10:55:17] error fitting model 51 : no method for assigning subsets of this S4 class
# WARN [2021-07-31 10:55:17] error fitting model 51 : no method for assigning subsets of this S4 class
# WARN [2021-07-31 10:55:17] error fitting model 51 : no method for assigning subsets of this S4 class
# WARN [2021-07-31 10:55:17] error fitting model 51 : no method for assigning subsets of this S4 class
# INFO [2021-07-31 10:55:17] searchlight iteration: 2
# WARN [2021-07-31 10:55:17] error fitting model 92 : no method for assigning subsets of this S4 class
# WARN [2021-07-31 10:55:17] error fitting model 92 : no method for assigning subsets of this S4 class
# WARN [2021-07-31 10:55:17] error fitting model 92 : no method for assigning subsets of this S4 class
# WARN [2021-07-31 10:55:17] error fitting model 92 : no method for assigning subsets of this S4 class
# WARN [2021-07-31 10:55:17] error fitting model 77 : training data must have more than one valid feature
# WARN [2021-07-31 10:55:17] error fitting model 77 : training data must have more than one valid feature
# WARN [2021-07-31 10:55:17] error fitting model 77 : training data must have more than one valid feature
# WARN [2021-07-31 10:55:17] error fitting model 77 : training data must have more than one valid feature
# INFO [2021-07-31 10:55:17] number of models fit: 3
# INFO [2021-07-31 10:55:17] no method for assigning subsets of this S4 class
# INFO [2021-07-31 10:55:17] no method for assigning subsets of this S4 class
# INFO [2021-07-31 10:55:17] training data must have more than one valid feature
# ERROR [2021-07-31 10:55:17] no valid results for randomized searchlight, exiting.
# Error in combiner(model_spec, good_results) : 
# searchlight: no searchlight samples produced valid results
# In addition: Warning messages:
# 1: UNRELIABLE VALUE: Future (‘<none>’) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced via the L'Ecuyer-CMRG method. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". 
# 2: UNRELIABLE VALUE: Future (‘<none>’) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced via the L'Ecuyer-CMRG method. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". 
# 3: UNRELIABLE VALUE: Future (‘<none>’) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced via the L'Ecuyer-CMRG method. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". 