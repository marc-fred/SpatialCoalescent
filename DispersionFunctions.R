gaussian <- function(x, sigma){
  # Computes a value for kernel dispersion using a gaussian model, i.e. a simple normal density distribution (sigma, mean=0)
  #
  # Args: 
  #   x: the distance between two points
  #   sigma: the value of the standard deviation
  #
  # Returns:
  #   The value of dispersion kernel for x
  return(dnorm(x, mean = 0, sd = sigma, log = FALSE))
}


distanceMatrixFromRaster <- function(object){
  # Computes a pairwise distance matrix from a raster object
  #
  # Args:
  #   object: a raster object from which computes distances
  #
  # Returns:
  #   A matrix of distances in meters if a coordinate system is precised
  
  # Extract coordinates from raster object
  coords = xyFromCell(object = object, cell = 1:ncell(object), spatial=FALSE)
  
  # Compute distance matrix
  dist <- apply(X = coords,
                MARGIN = 1,
                FUN = function(x){ spDistsN1(as.matrix(coords), x, longlat=TRUE) }
  )
    
  return(dist)
}

constructMigrationMatrix<- function(dist, param){
  # Computes a migration matrix given a distance matrix, gaussian kernel
  #
  # Args : 
  #   dist : a pairwise distance matrix between all cells
  #   param : the sigma parameter of the gaussian dispersion kernel function
  #
  # Returns :
  #   A migration matrix, with all rows summing to 1
  kernel <- apply(dist, c(1,2), gaussian, sigma = param)
  migrationRates <- kernel/rowSums(kernel)
  return(migrationRates)
}

