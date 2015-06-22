setwd("/home/arno/Documents/These/SpatialCoalescent/Classes")
setwd("/home/arnaudb/Documents/SpatialCoalescent/Classes")
library(methods)
library(raster)
source("Generics.R")
source("Environment.R")
source("Function.R")
source("Model.R")
source("SurModel.R")

pluie <- new("Environment", values= matrix(1:9, 3))
temp <- new("Environment", values= matrix(1:9, 3))
myPlot(pluie)

mk1 <- new("Model", varName = "Pluviométrie",  varEnv = pluie,
         fun = new("Function",
                   name = "Linear", 
                   fun = linearTwoParameters, 
                   param = list(f_0 = 0, slope = 1)))
mk1

mk2 <- new("Model", varName = "Température",  varEnv = temp,
         fun = new("Function",
                   name = "Linear", 
                   fun = linearTwoParameters, 
                   param = list(f_0 = 0, slope = 2)))
mk2

Kmodel <- new("KModel", models = list(mk1, mk2))
Kmodel
getParameters(Kmodel)
applyModel(KModel)

Rmodel <- new("RModel", models = list(mk1, mk2))
Rmodel
getParameters(Rmodel)
applyModel(Rmodel)


