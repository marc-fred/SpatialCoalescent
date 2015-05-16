
source("AskModelsFunctions.R")
source("NicheFunctions.R")
source("DispersionFunctions.R")
source("MutationFunctions.R")
source("CoalescentFunctions.R")
source("PriorFunctions.R")
source("MarkovProcess.R")

library(raster)

###### Environmental data :

environmentVariableRaster <- raster(matrix(data = sample(1:100, 100), ncol = 10),
                              xmn = 40,
                              xmx = 50,
                              ymn = 0,
                              ymx = 10,
                              crs = "+proj=longlat +datum=WGS84")
geographicDistanceMatrix <- distanceMatrixFromRaster(object = environmentVariableRaster)/1000
environmentVariableMatrix <- as.matrix(environmentVariableRaster)

###### Position data :
dataCoord <- xyFromCell(environmentVariableRaster, sample(1:ncell(environmentVariableRaster), 10))

###### Locus information :
nbLocus <- 2
# assuming we have the step values for each locus
steps <- c(2,3)


###### Model :

# Parameters : 
# dispersion : gaussian(0, sigma)
theta_sigma <- 1
# niche : constant(Y)
theta_Y <- 2
# mutation rate
theta_rate <- 10^-5

##### 
# launch simulations
simulateSpatialCoalescent(theta_sigma, theta_Y, theta_rate, environmentVariableMatrix, nbLocus, dataCoord)
