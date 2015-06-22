
setClass("SurModel",
         representation(models = "list"),
         prototype(models = list(new("Model"), new("Model"))),
         validity = function(object) { ## object : nom reserve !
           if (FALSE)
             return(FALSE)
           else
             return(TRUE)
         }
)

setClass("KModel", contains = "SurModel")

setClass("RModel", contains = "SurModel")


setMethod(f="ToStream",
          signature = "RModel",
          definition = function(object){
            cat("R model --------------------------------------- \n")
            for(i in 1:length(object@models)){
              ToStream(object@models[[i]])
            }
          }
)

setMethod(f="ToStream",
          signature = "KModel",
          definition = function(object){
            cat("K model --------------------------------------- \n")
            for(i in 1:length(object@models)){
              ToStream(object@models[[i]])
            }
          }
)

setMethod("show", "SurModel",
          definition = function(object){
            ToStream(object)
          })

setMethod(f="getParameters",
          signature = "KModel",
          definition = function(object){
            p <-  sapply(X = object@models,
                         FUN = function(x){paste("K.", getParameters(x), sep="")})
            return(as.vector(p))
          }
)

setMethod(f="getParameters",
          signature = "RModel",
          definition = function(object){
            p <-  sapply(X = object@models,
                         FUN = function(x){paste("R.", getParameters(x), sep="")})
            return(as.vector(p))
          }
)