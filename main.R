
source("AskModelsFunctions.R")
source("NicheFunctions.R")
source("DispersionFunctions.R")
source("MutationFunctions.R")
source("PriorFunctions.R")
source("MarkovProcess.R")

library(raster)

###### Environmental data :

environmentVariableRaster <- raster(matrix(data = sample(1:100, 9), ncol = 3),
                              xmn = 40,
                              xmx = 50,
                              ymn = 0,
                              ymx = 10,
                              crs = "+proj=longlat +datum=WGS84")
geoDistMatrix <- distanceMatrixFromRaster(object = environmentVariableRaster)/1000
envMatrix <- as.matrix(environmentVariableRaster)

###### Position data :
dataCoord <- xyFromCell(environmentVariableRaster, sample(1:ncell(environmentVariableRaster), 3))
localizationData <- cellFromXY(environmentVariableRaster, dataCoord)

###### Locus information :
nbLocus <- 2
# assuming we have the step values for each locus
steps <- c(2,3)


###### Model :

# Parameters : 
# dispersion : gaussian(0, sigma)
theta_sigma <- 1000
# niche : constant(Y)
theta_Y_k <- 5
theta_Y_r <- 1
# mutation rate
theta_rate <- 10^-5

##### 
# launch simulations
simulateSpatialCoalescent(theta_sigma, theta_Y, theta_rate, envMatrix, nbLocus, localizationData, steps)
