setClass("NicheFunction",
         contains = "Function",
         representation(fun = "NicheFunction"),
         prototype(fun = function(x){}),
         validity = function(object) { ## object : nom reserve !
           if ( FALSE)
             return(FALSE)
           else
             return(TRUE)
           }
         )