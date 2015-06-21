setwd("/home/arno/Documents/These/SpatialCoalescent/Classes")
library(methods)
library(raster)
source("NicheFunction.R")
source("Generics.R")
source("Environment.R")
source("Function.R")
source("Model.R")
source("SurModel.R")

linearTwoParameters <- function(x,X0,slope)
{
  # Computes a linear response within the enveloppe, else returns 0.
  #
  # Args:
  #   x: numeric providing the values of variable to calculate reaction norm
  #   X0: value of the function at x=0
  #   slope: the value of the slope of the function
  #
  # Returns:
  #   The value of the reaction norm
  return(slope*x-slope*X0)
}


pluie <- new("Environment", values= matrix(1:9, 3))
temp <- new("Environment", values= matrix(1:9, 3))
myPlot(pluie)

mk1 <- new("Model", varName = "Pluviométrie",  varEnv = pluie,
         fun = new("Function",
                   name = "Linear", 
                   fun = linearTwoParameters, 
                   param = list(X0 = 0, slope = 1)))
mk1

mk2 <- new("Model", varName = "Température",  varEnv = temp,
         fun = new("Function",
                   name = "Linear", 
                   fun = linearTwoParameters, 
                   param = list(X0 = 3, slope = 2)))
mk2

Kmodel <- new("KModel", models = list(mk1, mk2))
Kmodel
Rmodel <- new("RModel", models = list(mk1, mk2))
Rmodel