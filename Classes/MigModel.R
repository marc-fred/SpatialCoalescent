setClass("MigModel", contains = "Model")

setMethod(
  f = "applyModel",
  signature = "MigModel",
  definition =  function(object){
    val <- getValues(object@varEnv)
    val[] <- applyFunction(object@fun, xval =c(val))
    migRates <- val/rowSums(val)
    return(migRates)
  })

setMethod("initialize",
          "MigModel",
          function(.Object, varName, varEnv, fun) {
            print(" ----- initialize from MigModel called---- ")
            .Object@varName <- varName
            .Object@varEnv <- varEnv
            .Object@fun <- fun
            # validObject(.Object) ## valide l'objet
            return(.Object)
          })

setMethod(f="getParameters",
          signature = "MigModel",
          definition = function(object){
            p <-  sapply(X = callNextMethod(),
                         FUN = function(x){paste("Mig", x, sep=".")})
            return(as.vector(p))
          }
)