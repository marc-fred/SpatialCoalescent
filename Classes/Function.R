setClass("Function",
         representation(name = "character",
                        fun = "function",
                        param = "list"),
         prototype(name = "", fun = function(x){}, param = list()),
                   validity = function(object) { ## object : nom reserve !
                     if ( FALSE)
                       return(FALSE)
                     else
                       return(TRUE)
                   })

setMethod("show", "Function",
          function(object){
            cat("Function", object@name, "with parameters","\t")
            pnames <- lapply(X = 1:length(object@param),
                             FUN = function(x, l){ paste(names(l)[x], "=", l[x], sep ="")},
                             l = object@param
            )
            for(i in 1:length(pnames)){
              cat(pnames[[i]], "\t" )
            }
          })

new("Function", fun = dnorm, param = list(mean = 0, sd = 1))
new("Function", name = "linear", fun = linearTwoParameters, param = list(X0 = 0, slope = 1))
