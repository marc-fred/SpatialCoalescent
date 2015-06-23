setClass("MigModel", contains = "SuperModel")

setMethod(
  f="ToStream",
  signature = "MigModel",
  definition = function(object){
    cat("Mig model --------------------------------------- \n")
    callNextMethod(object)
  }
)

setMethod(
  f="getParameters",
  signature = "MigModel",
  definition = function(object){
    p <-  sapply(X = object@models,
                 FUN = function(x){paste("Mig.", getParameters(x), sep="")})
    return(as.vector(p))
  }
)

setMethod(
  f = "applyModel",
  signature = "MigModel",
  definition =  function(object){
    meanVal <- callNextMethod(object)
    migRates <- meanVal/rowSums(meanVal)
    return(migRates)
  }
)
