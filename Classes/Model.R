setClass("Model",
         representation(varEnv = "Environment", 
                        fun = "Function", 
                        param = "Parameters"),
         prototype(varEnv = new("Environment"), 
                   fun = new("Function"), 
                   param = new("Parameters")),
         validity = function(object) { ## object : nom reserve !
           if (FALSE)
             return(FALSE)
           else
             return(TRUE)
         })

setGeneric("ToStream", function(obj, separator) {
  return(standardGeneric("getFullName"))
})