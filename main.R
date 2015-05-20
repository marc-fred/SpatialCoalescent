source("CoalescentFunctions.R")
source("NicheFunctions.R")
source("DispersionFunctions.R")
source("MutationFunctions.R")
source("PriorFunctions.R")
source("MarkovProcess.R")
source("generalFunctions.R")

library(raster)
library(parallel)

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


local({
  
  # open a connection to a temporary file : pipe between master process and child
  f <- fifo(tempfile(), open="w+b", blocking=T)
  
  if (inherits(parallel:::mcfork(), "masterProcess")) {
    # Child
    progress <- 0.0
    
    while (progress < 1 && !isIncomplete(f)) {
      msg <- readBin(f, "double")
      progress <- progress + as.numeric(msg)
      # send a message in C-style
      cat(sprintf("Progress: %.2f%%\n", progress * 100))
    } 
    
    # close the current child process, informing master process
    parallel:::mcexit()
  }
  
  
  numJobs <- 50000
  
  mclapply(X = 1:numJobs, FUN = function(x, geoDistMatrix, envMatrix, localizationData, nbLocus, steps){
    
    # Draw parameters : 
    theta_sigma <- uniform(n=1, min = 1, max = 10000)
    theta_Y_k <- uniform(n=1, min = 0, max = 100)
    theta_Y_r <- uniform(n=1, min = 0, max = 100)
    theta_rate <- uniform(n=1, min=0, max = 10)
    
    # Launch simulation
    genetics <- simulateSpatialCoalescent(theta_sigma = theta_sigma, 
                                          theta_Y_r = theta_Y_r,
                                          theta_Y_k = theta_Y_k,
                                          theta_rate = theta_rate,
                                          envMatrix = envMatrix,
                                          nbLocus = nbLocus,
                                          localizationData = localizationData, 
                                          steps = steps,
                                          geoDistMatrix = geoDistMatrix)
    
    # write results of genetic data 
    fileName = paste("Genetics_", x , ".txt", sep="")
    writeDataOutputInFile(theta_sigma, theta_Y_k, theta_Y_r, theta_rate, data = genetics, file = fileName)
    
    # Send progress update
    writeBin(1/numJobs, f)
    
  }, 
  geoDistMatrix = geoDistMatrix, 
  envMatrix = envMatrix, 
  localizationData = localizationData, 
  nbLocus = nbLocus, 
  steps = steps,
  mc.cores = 20,
  mc.preschedule = FALSE)
  
})
cat("Simulations Done\n")
