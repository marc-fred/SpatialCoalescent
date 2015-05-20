source("CoalescentFunctions.R")
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
dataCoord <- xyFromCell(environmentVariableRaster, sample(1:ncell(environmentVariableRaster), 20, replace = TRUE))
localizationData <- cellFromXY(environmentVariableRaster, dataCoord)

###### Locus information :
nbLocus <- 10
# assuming we have the step values for each locus
steps <- sample(1:10, size = nbLocus ) 


###### Model :

# Parameters : 
# dispersion : gaussian(0, sigma)
theta_sigma <- uniform(n=1, min = 1, max = 10000)
# niche : constant(Y)
theta_Y_k <- uniform(n=1, min = 0, max = 100)
theta_Y_r <- uniform(n=1, min = 0, max = 100)
# mutation rate
theta_rate <- uniform(n=1, min=0, max = 10)

##### 
# launch simulations
genetics <- simulateSpatialCoalescent(theta_sigma = theta_sigma, 
                                      theta_Y_r = theta_Y_r,
                                      theta_Y_k = theta_Y_k,
                                      theta_rate = theta_rate,
                                      envMatrix = envMatrix,
                                      nbLocus = nbLocus,
                                      localizationData = localizationData, 
                                      steps = steps,
                                      geoDistMatrix = geoDistMatrix)

# write result
writeDataOutputInFile(theta_sigma,theta_Y_k,theta_Y_r,theta_rate, file="myFile")
