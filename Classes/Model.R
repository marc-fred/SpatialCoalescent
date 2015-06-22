setClass("Model",
         representation(varName = "character",
                        varEnv = "Environment",
                        fun = "Function"),
         prototype(varName = "Undefined",
                   varEnv = new("Environment"),
                   fun = new("Function")),
         validity = function(object) { ## object : nom reserve !
           if (FALSE)
             return(FALSE)
           else
             return(TRUE)
         })

setMethod("show", "Model",
          function(object){
            ToStream(object)
          })

setMethod(f="ToStream",
          signature = "Model",
          definition = function(object){
            cat("###", object@varName, "model\n")
            cat("\t")
            ToStream(object@fun)
            cat("\n")
          }
)

setMethod(f="myWrite",
          signature = "Model",
          definition = function(object, file){
            sink(file = file)
            ToStream(object)
            sink()
          }
)

setMethod(
  f = "getParameters",
  signature = "Model",
  definition =  function(object){
    vapply(X = getParameters(object@fun),
           FUN = function(x) {paste(object@varName, x, sep =".")},
           FUN.VALUE = c("a"),
           USE.NAMES = FALSE
    )
  })

setMethod(
  f = "applyModel",
  signature = "Model",
  definition =  function(object){
    val <- getValues(object@varEnv)
    val[] <- applyFunction(object@fun, xval =c(val))
    return(val)
  })