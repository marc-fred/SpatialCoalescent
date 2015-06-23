setClass("SuperModel",
         representation(models = "list"),
         prototype(models = list(new("Model"), new("Model"))),
         validity = function(object) { ## object : nom reserve !
           if (FALSE)
             return(FALSE)
           else
             return(TRUE)
         }
)

setMethod(
  f = "show",
  signature = "SuperModel",
  definition = function(object){
    ToStream(object)
  }
)

setMethod(
  f="ToStream",
  signature = "SuperModel",
  definition = function(object){
    for(i in 1:length(object@models)){
      ToStream(object@models[[i]])
    }}
)

setMethod(
  f = "applyModel", 
  signature = "SuperModel",
  definition = function(object){
    rawList <- lapply(X = object@models,
                      FUN = applyModel)
    meanVal <- Reduce("+", rawList) / length(rawList)
    return(meanVal)
  })












