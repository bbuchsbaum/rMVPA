

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


