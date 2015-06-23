setClass("MigModel", contains = "SuperModel")

setMethod(
  f = "applyModel",
  signature = "MigModel",
  definition =  function(object){
    callNextMethod()
    migRates <- meanVal/rowSums(meanVal)
    return(migRates)
  })

setMethod(f="ToStream",
          signature = "MigModel",
          definition = function(object){
            cat("Migration model --------------------------------------- \n")
            callNextMethod()
          }
)

setMethod(f="getParameters",
          signature = "MigModel",
          definition = function(object){
            p <-  sapply(X = object@models,
                         FUN = function(x){paste("Mig.", getParameters(x), sep="")})
            return(as.vector(p))
          }
)