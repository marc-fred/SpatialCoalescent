setClass("Environment",
         representation(values = "matrix"),
         prototype(values = matrix()),
         validity = function(object) { ## object : nom reserve !
           if (FALSE)
             return(FALSE)
           else
             return(TRUE)
         })