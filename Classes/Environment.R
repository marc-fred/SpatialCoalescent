setClass("Lattice",
         representation(values = "matrix"),
         prototype(values = matrix()),
         validity = function(object) { ## object : nom reserve !
           if (FALSE)
             return(FALSE)
           else
             return(TRUE)
         })

setClass("Environment", contains = "Lattice")

setMethod("show", "Lattice",
          function(object){
            ToStream(object)
          })

setMethod("myPlot", "Lattice",
          function(object){
            m <- list(object@values)
            rl = lapply(m, function(X) raster(X))
            d <- stack(rl)
            spplot(d)
          })

setMethod("ToStream", "Lattice",
          function(object){
            show(object@values)
          })

setMethod(f="getValues",
          signature = "Lattice",
          definition = function(object){
            return(object@values)
          })