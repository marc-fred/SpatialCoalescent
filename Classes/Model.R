setClass("Model", contains = "AbstractModel")

setMethod(
  f = "applyModel",
  signature = "Model",
  definition =  function(object){
    val <- getValues(object@varEnv)
    val[] <- applyFunction(object@fun, xval =c(val))
    return(val)
  })

model <- function(varName, varEnv, fun){
  if(missing(varName)){
    varName <- deparse(substitute(varEnv))
  }
  new("Model", varName = varName, varEnv = varEnv, fun = fun)
}
