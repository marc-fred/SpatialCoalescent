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

x <- 1

set.seed(x)
# Draw parameters : 
theta_sigma <- uniform(n=1, min = 1, max = 10000)
theta_Y_k <- uniform(n=1, min = 0, max = 100)
theta_Y_r <- uniform(n=1, min = 0, max = 100)
theta_rate <- uniform(n=1, min=0, max = 10)

kMatrix <- constructEnvironmentalDemographicMatrix(env = envMatrix, param = theta_Y_k)
rMatrix <- constructEnvironmentalDemographicMatrix(env = envMatrix, param = theta_Y_r)
migMatrix <- constructMigrationMatrix(dist = geoDistMatrix , param = theta_sigma)

createInitialDemographicsLandscape <- function(landscape){
  # Function to create an initial matrix of demographic size : fondator
  #
  # Args : 
  #   landscape : a model of landscape matrix
  #
  # Returns : 
  #   A matrix
  N_m <- landscape
  N_m[] <- 0
  N_m[1,1] <- 1
  return(N_m)
}

reproductionStep <- function(demographicMatrix, kMatrix, rMatrix){
  # Function to create an population matrix of demographic size after reproduction
  #
  # Args : 
  #   demographicMatrix : a initial matrix of population size
  #   kMatrix : a matrix of carrying capacity
  #   rMatrix : a matrix of growth rate
  #
  # Returns : 
  #   A matrix of population after reproduction
  
  N_v <- as.vector(demographicMatrix)
  k_v <- as.vector(kMatrix)
  r_v <- as.vector(rMatrix)
  N_tilde_v <- mapply(FUN=function(N_v, r_v, K_v){ N_v*(1+r_v)/(1+(r_v*N_v)/K_v) }, N_v, r_v, k_v)
  N_tilde_m <- matrix(data = N_tilde_v, ncol = ncol(demographicMatrix))
  return(N_tilde_m)
}

migrationStep <- function(demographicMatrix, migMatrix){
  # Function to create a population matrix of demographic size after dispersion
  #
  # Args : 
  #   demographicMatrix : a initial matrix of population size
  #   migMatrix : a matrix of migration probabilities between cells
  #
  # Returns : 
  #   A matrix of population size after dispersion
  N_tilde_v <- as.vector(demographicMatrix)
  N_v <- N_tilde_v %*% migMatrix
  N_m <- matrix(data = N_v, ncol = ncol(demographicMatrix))
  return(N_m)
}

demographicGeneration <- function(demographicMatrix, kMatrix, rMatrix, migMatrix){
  # Function to create a population matrix of demographic size after reproduction then dispersion
  #
  # Args : 
  #   demographicMatrix : a initial matrix of population size
  #   kMatrix : a matrix of carrying capacity
  #   rMatrix : a matrix of growth rate
  #   migMatrix : a matrix of migration probabilities between cells
  #
  # Returns : 
  #   A list with $demography : matrix of population size and $migration : matrix of pop flux.
  N_tilde_m <- reproductionStep(demographicMatrix, kMatrix, rMatrix)
  
  migHistory_m <- computeMigrationHistory(N_m = N_tilde_m, mig_m = migMatrix)
  N_m <- migrationStep(demographicMatrix = N_tilde_m, migMatrix = migMatrix)
  
  generation_l <- list(demography = N_m, migration = migHistory_m)
  return(generation_l)
}

demographicSimulation <- function(numberOfGenerations, demographicMatrix, kMatrix, rMatrix, migMatrix){
  # Function to create a population matrix of demographic size after reproduction then dispersion
  #
  # Args :
  #   numberOfGenerations : the number of generations along demographic history
  #   demographicMatrix : a initial matrix of population size
  #   kMatrix : a matrix of carrying capacity
  #   rMatrix : a matrix of growth rate
  #   migMatrix : a matrix of migration probabilities between cells
  #
  # Returns : 
  #   A list of matrix of population size, forward in time
  demoHistory_l <- list()
  migHistory_l <- list()

  i_time <- 1
  while(i_time < numberOfGenerations){
    demoHistory_l[[i_time]] <- demographicMatrix
    
    generationList <- demographicGeneration(demographicMatrix, kMatrix, rMatrix, migMatrix)
    
    demographicMatrix <- generationList$demography
    migHistory_l[[i_time]] <- generationList$migration
    
    i_time <- i_time + 1
  }
  
  history_l <- list(demoHistory = demoHistory_l, migHistory = migHistory_l)
  return(history_l)
}

computeMigrationHistory <- function(N_m, mig_m){
  # Function to create a matrix of migration effective number from the migration and the demographic matrix
  #
  # Args :
  #   N_m : a matrix of population size. No need to be integer.
  #   r_m : a matrix of growth rate
  #   mig_m : a matrix of migration probabilities between cells
  #
  # Returns : 
  #   A matrix of effective migrant number : each column i gives the distribution across demes of migrants coming to the deme i.
  N_v <- as.vector(N_m)
  migrants_m <- N_v*mig_m
  return(migrants_m)
}
