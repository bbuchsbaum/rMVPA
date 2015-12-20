
BootstrapSampler <- R6Class("BootstrapSampler",
                            public=list(
                              dataset=NA,
              
                              initialize= function(dataset) {
                                dataset=daatset
                              },
                              
                            )
                              
                              
MVPASearchlight <- R6Class("MVPASearchlight",
                        public = list(
                          dataset=NA,
                          radius=8,
                          method="standard",
                          niter=16,
                          
                          initialize = function(dataset, model) {
                            if (!missing(name)) self$name <- name
                            if (!missing(hair)) self$hair <- hair
                            self$greet()
                          },
                          set_hair = function(val) {
                            self$hair <- val
                          },
                          greet = function() {
                            cat(paste0("Hello, my name is ", self$name, ".\n"))
                          }
                        )
)


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


