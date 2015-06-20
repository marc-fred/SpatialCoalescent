setClass("Parameters",
         representation(name = "character", val = "numeric"),
         prototype(name = NULL, val = NULL),
                   validity = function(object) { ## object : nom reserve !
                     if (FALSE)
                       return(FALSE)
                     else
                       return(TRUE)
                   })
