setClass("Environment",
         representation(values = "matrix"),
         prototype(values = matrix()),
         validity = function(object) { ## object : nom reserve !
           if (FALSE)
             return(FALSE)
           else
             return(TRUE)
         })

setMethod("show", "Environment",
          function(object){
            ToStream(object)
          })

setMethod("myPlot", "Environment",
          function(object){
            m <- list(object@values)
            rl = lapply(m, function(X) raster(X))
            d <- stack(rl)
            spplot(d)
          })

setMethod("ToStream", "Environment",
          function(object){
            show(object@values)
          })

setMethod(f="getValues",
          signature = "Environment",
          definition = function(object){
            return(object@values)
          })