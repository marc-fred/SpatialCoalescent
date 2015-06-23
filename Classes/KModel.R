setClass("KModel", contains = "SuperModel")

setMethod(
  f="ToStream",
  signature = "KModel",
  definition = function(object){
    cat("K model --------------------------------------- \n")
    callNextMethod(object)
  }
)

setMethod(
  f="getParameters",
  signature = "KModel",
  definition = function(object){
    p <-  sapply(X = object@models,
                 FUN = function(x){paste("K.", getParameters(x), sep="")})
    return(as.vector(p))
  }
)