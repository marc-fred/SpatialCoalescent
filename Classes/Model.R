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

setGeneric(
  name = "ToStream",
  def = function(object, con, description) { return(standardGeneric("ToStream"))})

setGeneric(
  name = "myWrite",
  def = function(object, file) { return(standardGeneric("myWrite"))})

setMethod(f="ToStream",
          signature = "Model",
          definition = function(object){
            cat("MODEL FOR ENVIRONMENTAL VARIABLE", object@varName, ":\n")
            ToStream(object@fun)
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
