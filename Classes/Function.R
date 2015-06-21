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
            ToStream(object)
          })

setMethod(f="ToStream",
          signature = "Function",
          definition = function(object){        
            cat("Function", object@name, "with parameters","\t")
            pnames <- lapply(X = 1:length(object@param),
                             FUN = function(x, l){ paste(names(l)[x], "=", l[x], sep ="")},
                             l = object@param)
            
            for(i in 1:length(pnames)) cat(pnames[[i]], "\t" )
          })
