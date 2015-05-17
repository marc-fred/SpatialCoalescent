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
  kernel <- apply(dist, c(1,2), gaussian, sigma = param)
  migrationRates <- kernel/rowSums(kernel)
  return(migrationRates)
}

dispersionFunctionForValue <- function(dispersionFunction, x, args){
  # Compute a dispersion kernel function over a single value
  #
  # Args:
  #   dispersionFunction: the name of the dispersion kernel function which is called
  #   x: the distance between two points
  #   args : a list of the arguments of the dispersion function
  # 
  # Returns:
  #   The value of dispersion kernel for x
  args <- c(list(x), args)
  return(do.call(dispersionFunction, args))
  
  # Ex : 
  # dispersionFunctionForValue(fatTail1, x=4, args=list(alpha=0.5, beta=0.2))
}

dispersionFunctionForArray <- function(dispersionFunction, array, args){
  # Apply a dispersion kernel function over an array
  #
  # Args:
  #   dispersionFunction: the name of the dispersion kernel function which is called
  #   array: the array of distances used to compute dispersion function.
  #   args : a list of the arguments of the dispersion kernel function
  # 
  # Returns:
  #   An array corresponding to the dispersion kernel values
  apply(X=array, MARGIN=1, FUN=dispersionFunctionForValue, args=args, dispersionFunction=dispersionFunction)
  
  # Ex:
  # dispersionFunctionForArray(dispersionFunction=fatTail2,
  #                            array=array(data= 1:10, dim =10),
  #                            args=list(sigma = 0.2, gamma=0.3 ))
}


computeDispersionKernel <- function(dispersionFunction, distanceMatrix, args){
  # Apply a dispersion kernel function over distances between the cells of a rasterLayer
  #
  # Args:
  #   dispersionFunction: the dispersion kernel function to apply
  #   distanceMatrix: a pairwise distance matrix, computed with distanceMatrixFromRaster function
  #   args: a list of the dispersion function arguments
  #
  # Returns:
  #   A matrix of dispersion kernel values, similar in shape with the distance matrix between cells.
  
  # Apply dispersionFunction
  dispersion <- apply(X=distanceMatrix, 
                      MARGIN=2, 
                      FUN=dispersionFunctionForValue,
                      dispersionFunction=dispersionFunction,
                      args=args)
  return(dispersion)
}

migrationRateMatrix <- function(dispersion){
  # Normalizes a matrix of dispersion kernel between cells to get a migration rate matrix between cells.
  #
  # Args:
  #   dispersion: a matrix representing the values of a specified kernel (function of distances between cells)
  #
  # Returns:
  #   A migration rate matrix (note that rowSums and colSums are not 1: cause of bordure effect, individuals go "out of the world")
  return(dispersion/max(c(colSums(dispersion),rowSums(dispersion))))
}