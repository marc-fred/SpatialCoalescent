setwd("/home/arnaudb/Documents/SpatialCoalescent")
source("CoalescentFunctions.R")
source("NicheFunctions.R")
source("DispersionFunctions.R")
source("MutationFunctions.R")
source("PriorFunctions.R")
source("MarkovProcess.R")
source("generalFunctions.R")
source("demographicGrowth.R")
library(raster)
library(parallel)

# "Real" Data
rasterE1 <- raster(x = matrix(data = sample(1:100, 9), ncol = 3),
                   xmn = 40, xmx = 50, ymn = 0, ymx = 10, crs = "+proj=longlat +datum=WGS84")

rasterE2 <- raster(x = matrix(data = sample(1:100, 9), ncol = 3),
                   xmn = 40, xmx = 50, ymn = 0, ymx = 10, crs = "+proj=longlat +datum=WGS84")

dataCoord <- xyFromCell(rasterE1, sample(1:ncell(rasterE1), 20, replace = TRUE))
localizationData <- cellFromXY(rasterE1, dataCoord)

nbLocus <- 10
steps <- sample(1:10, size = nbLocus ) 

# Model Implementation I

E1 <- as.matrix(rasterE1)
myprint(E1)
E2 <- as.matrix(rasterE2)
myprint(E2)

R_m <- demoEnvironmentalVariableConstructor(list(E1, linearTwoParameters, list(X0 = 3, slope = 2)),
                                            list(E2, gaussianNiche, list(mean = 0, sd = 1)))
myprint(R_m)

K_m <- demoEnvironmentalVariableConstructor(list(E1, constant, list(Y = 20)))
myprint(K_m)

dist_m <- distanceMatrixFromRaster(object = rasterE1)/1000
M_m <- migrationMatrixConstructor(list(dist_m , gaussianDisp, list(sd = 1)))

land <- E1 ; land[] <- NA
demoInit <- createInitialDemographicsLandscape(land)
myprint(demoInit)

history <- demographicSimulation(numberOfGenerations = 10,
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
                                                 steps = steps)
