setwd("/home/arno/Documents/These/SpatialCoalescent/Classes")
source("Environment.R")
source("Function.R")
source("Model.R")

a <- new("Model", varName = "Pluviométrie",  varEnv = new("Environment"),
         fun = new("Function",
                   name = "Linear", 
                   fun = linearTwoParameters, 
                   param = list(X0 = 0, slope = 1)))

b <- new("Model", varName = "Température",  varEnv = new("Environment"),
         fun = new("Function",
                   name = "Linear", 
                   fun = linearTwoParameters, 
                   param = list(X0 = 3, slope = 2)))

show(b)
myWrite(b, file = "toto.txt")

getwd()
new("Function",name= "Gaussian", fun = dnorm, param = list(mean = 0, sd = 1))
show(a)
