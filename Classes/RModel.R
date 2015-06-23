setClass("RModel", contains = "SuperModel")

setMethod(
  f="ToStream",
  signature = "RModel",
  definition = function(object){
    cat("R model --------------------------------------- \n")
    callNextMethod(object)
  }
)

setMethod(
  f="getParameters",
  signature = "RModel",
  definition = function(object){
    p <-  sapply(X = object@models,
                 FUN = function(x){paste("R.", getParameters(x), sep="")})
    return(as.vector(p))
  }
)