setClass("Model", contains = "AbstractModel")

setMethod(
  f = "applyModel",
  signature = "Model",
  definition =  function(object){
    val <- getValues(object@varEnv)
    val[] <- applyFunction(object@fun, xval =c(val))
    return(val)
  })