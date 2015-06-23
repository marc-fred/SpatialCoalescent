setClass("AbstractModel",
         representation(varName = "character",
                        varEnv = "Lattice",
                        fun = "Function",
                        "VIRTUAL"),
         prototype(varName = "Undefined",
                   varEnv = new("Lattice"),
                   fun = new("Function")),
         validity = function(object) { ## object : nom reserve !
           if (FALSE)
             return(FALSE)
           else
             return(TRUE)
         })

setMethod("show", "AbstractModel",
          function(object){
            ToStream(object)
          })

setMethod(f="ToStream",
          signature = "AbstractModel",
          definition = function(object){
            cat("###", object@varName, "model\n")
            cat("\t")
            ToStream(object@fun)
            cat("\n")
          }
)

setMethod(f="myWrite",
          signature = "AbstractModel",
          definition = function(object, file){
            sink(file = file)
            ToStream(object)
            sink()
          }
)

setMethod(
  f = "getParameters",
  signature = "AbstractModel",
  definition =  function(object){
    vapply(X = getParameters(object@fun),
           FUN = function(x) {paste(object@varName, x, sep =".")},
           FUN.VALUE = c("a"),
           USE.NAMES = FALSE
    )
  })
