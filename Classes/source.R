# setwd("/home/arno/Documents/These/SpatialCoalescent/Classes")
setwd("/home/arnaudb/Documents/SpatialCoalescent")
library(methods)
library(raster)
source("Classes/Generics.R")
source("Classes/Environment.R")
source("Classes/Function.R")
source("Classes/AbstractModel.R")
source("Classes/Model.R")
source("Classes/MigModel.R")
source("Classes/SurModel.R")
source("Classes/RasterLayer.R")
source("CoalescentFunctions.R")
source("NicheFunctions.R")
source("DispersionFunctions.R")
source("MutationFunctions.R")
source("PriorFunctions.R")
source("MarkovProcess.R")
source("generalFunctions.R")
source("demographicGrowth.R")

# "Real" Data
rasterE1 <- raster(x = matrix(data = sample(1:100, 9), ncol = 3),
                   xmn = 40, xmx = 50, ymn = 0, ymx = 10, crs = "+proj=longlat +datum=WGS84")

rasterE2 <- raster(x = matrix(data = sample(1:100, 9), ncol = 3),
                   xmn = 40, xmx = 50, ymn = 0, ymx = 10, crs = "+proj=longlat +datum=WGS84")

dataCoord <- xyFromCell(rasterE1, sample(1:ncell(rasterE1), 20, replace = TRUE))

nbLocus <- 10
steps <- sample(1:10, size = nbLocus ) 

# Model Implementation
pluie <- new("Environment", values= as.matrix(rasterE1))
temp <- new("Environment", values= as.matrix(rasterE2))
myPlot(pluie)

mk1 <- new("Model", varName = "Pluviométrie",  varEnv = pluie, fun = new("Function",
                                                                         name = "Linear", 
                                                                         fun = linearTwoParameters, 
                                                                         param = list(f_0 = 0, slope = 1)))

mk2 <- new("Model", varName = "Température",  varEnv = temp, fun = new("Function",
                                                                       name = "Linear", 
                                                                       fun = linearTwoParameters, 
                                                                       param = list(f_0 = 0, slope = 2)))

Kmodel <- new("KModel", models = list(mk1, mk2))
K_m <- applyModel(Kmodel)

Rmodel <- new("RModel", models = list(mk1, mk2))
R_m <- applyModel(Rmodel)

distances <- new("Lattice", values= computeDistanceMatrix(rasterE1))
migFun <- new("Function", name = "Gaussian", fun = gaussianDisp, param = list(mean=0, sd = 250 ))
migModel <- new("MigModel", varName = "Distances", varEnv = distances, fun = migFun)
M_m <- applyModel(migModel)

demoInit <- createInitialDemographicsLandscape(K_m)
myprint(demoInit)

history <- demographicSimulation(numberOfGenerations = 20,
                                 demographicMatrix = demoInit, 
                                 kMatrix = K_m, 
                                 rMatrix = R_m,
                                 migMatrix = M_m)

demoHistory_l <- history$demoHistory
myprint(m = demoHistory_l)
migHistory_l <- history$migHistory
myprint(migHistory_l)

genetValues <- spatialCoalescenceForMultipleLoci(migHistory_l = migHistory_l,
                                                 demoHistory_l = demoHistory_l, 
                                                 localizationData = localizationData, 
                                                 nbLocus = nbLocus, 
                                                 theta_rate = runif(n = nbLocus, min=0.1, max = 0.5),
                                                 steps = steps,
                                                 rasterLandscape = rasterE1)

x <- 1
dataFileName = paste("genetics_", x , ".txt", sep="")
writeDataOutputInFile(Kmodel, Rmodel, MigModel, theta_rate, genetData = genetValues, file = dataFileName)
  
)



