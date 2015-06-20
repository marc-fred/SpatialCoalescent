setClass("Function",
         representation(f = "function"),
         prototype(f = function(x){}),
                   validity = function(object) { ## object : nom reserve !
                     if ( FALSE)
                       return(FALSE)
                     else
                       return(TRUE)
                   })