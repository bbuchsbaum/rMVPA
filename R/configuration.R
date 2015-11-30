#dataset, regionMask, ncores=1, savePredictors=FALSE, autobalance=FALSE, bootstrap=FALSE, 
#featureSelector=NULL, featureParcellation=NULL, classMetrics=FALSE, ensemblePredictor=FALSE) {  


MVPAConfiguration <- function(modelName="corclass", 
                              tuneGrid=NULL, 
                              tuneLength=1, 
                              ncores=1, 
                              savePredictors=FALSE,
                              ensemblePredictor=FALSE,
                              balanceSamples=FALSE, 
                              bootstrapSamples=FALSE,
                              featureSelector=NULL, 
                              featureParcellation=NULL, 
                              saveClassMetrics=FALSE,
                              modelControl=NULL) {
  
  ret <- list(modelName=modelName,
              tuneGrid=tuneGrid,
              tuneLength=tuneLength,
              ncores=ncores,
              savePredictors=savePredictors,
              balanceSamples=balanceSamples,
              bootstrapSamples=bootstrapSamples,
              featureSelector=featureSelector,
              featureParcellation=featureParcellation,
              saveClassMetrics=saveClassMetrics,
              modelControl=modelControl)
  
  class(ret) <- "MVPAConfiguration"
  ret
              
}

MVPASearchlightConfiguration <- function(modelName="sda_notune", 
                              tuneGrid=NULL, 
                              tuneLength=1, 
                              ncores=1, 
                              savePredictors=FALSE,
                              ensemblePredictor=FALSE,
                              balanceSamples=FALSE, 
                              bootstrapSamples=FALSE,
                              featureSelector=NULL, 
                              featureParcellation=NULL, 
                              saveClassMetrics=FALSE,
                              modelControl=NULL,
                              radius=8,
                              iterations=16) {
  
  ret <- list(modelName=modelName,
              tuneGrid=tuneGrid,
              tuneLength=tuneLength,
              ncores=ncores,
              savePredictors=savePredictors,
              balanceSamples=balanceSamples,
              bootstrapSamples=bootstrapSamples,
              featureSelector=featureSelector,
              featureParcellation=featureParcellation,
              saveClassMetrics=saveClassMetrics,
              modelControl=modelControl,
              radius=radius,
              iterations=iterations,
              method="randomized")
  
  class(ret) <- "MVPAConfiguration"
  ret
  
}



