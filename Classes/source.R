# setwd("/home/arno/Documents/These/SpatialCoalescent/Classes")
setwd("/home/arnaudb/Documents/SpatialCoalescent")
library(methods)
library(raster)
source("NicheFunctions.R")
source('DispersionFunctions.R')
source("Classes/Generics.R")
source("Classes/Environment.R")
source("Classes/Function.R")
source("Classes/Model.R")
source("Classes/SurModel.R")
source("Classes/RasterLayer.R")

rasterE1 <- raster(x = matrix(data = sample(1:100, 9), ncol = 3),
                   xmn = 40, xmx = 50, ymn = 0, ymx = 10, crs = "+proj=longlat +datum=WGS84")

rasterE2 <- raster(x = matrix(data = sample(1:100, 9), ncol = 3),
                   xmn = 40, xmx = 50, ymn = 0, ymx = 10, crs = "+proj=longlat +datum=WGS84")

pluie <- new("Environment", values= as.matrix(rasterE1))
temp <- new("Environment", values= as.matrix(rasterE2))
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
K <- applyModel(Kmodel)

Rmodel <- new("RModel", models = list(mk1, mk2))
Rmodel
getParameters(Rmodel)
R <- applyModel(Rmodel)

distances <- new("Lattice", values= computeDistanceMatrix(rasterE1))
migModel <- new("Model", varName = "Distances", varEnv = distances,
              fun = new("Function",
                        name = "Gaussian",
                        fun = gaussianDisp,
                        param = list(mean=0, sd = 10 )))
mig <- applyModel(migModel)

