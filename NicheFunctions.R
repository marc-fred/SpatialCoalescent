constant <- function(x,Y)
{
  # Compute a constant response
  #
  # Args:
  #   x : numeric providing the values of variable to calculate reaction norm
  #   Y : value of the constant response
  #
  # Returns:
  #   The value of the reaction norm
  return(Y)
}

linearTwoParameters <- function(x,X0,slope)
{
  # Computes a linear response within the enveloppe, else returns 0.
  #
  # Args:
  #   x: numeric providing the values of variable to calculate reaction norm
  #   X0: value of the function at x=0
  #   slope: the value of the slope of the function
  #
  # Returns:
  #   The value of the reaction norm
  return(slope*x-slope*X0)
}

gaussianNiche <- function(x, mean, sd){
  y <- dnorm(x, mean, sd)
  return(y)
}

constructEnvironmentalDemographicMatrix <- function(env, param){
  # Compute a demographic variable from environmental data, with constant niche function
  #
  # Args :
  #   env : a matrix of the environmental variable values.
  #   param : the parameter Y of the constant niche function
  #
  # Returns :
  #   A matrix giving the value of the computed demographic variable
  
  matrix <- apply(env, c(1,2), constant, Y = param)
  matrix <- round(matrix)
  
  # Prevent impossible transition :
  if( sum(matrix) <= 0){
    stop(paste("The parameters values for niche models led to null demographic variable over all cells of the landscape :
             coalescence is impossible to simulate"))
  }
  return(matrix)
}

ComputeVariableFunction <- function(...){
  # Apply function with variable arguments over a matrix
  #
  # Args:
  #
  #   ... various lists of form list(matrix, function, list())
  #
  # Returns:
  #   A list of matrix with length equal to the number of lists passed as arguments
  
  ellipsisArgs <- list(...)
  l <- lapply( X = ellipsisArgs,
               FUN = function(model){
                 apply(X = model[[1]],
                       MARGIN = c(1,2),
                       FUN = function(x, model){
                         args <- c(x = x, model[[3]])
                         do.call(what = model[[2]], args)
                       },
                       model = model
                 )
               }
  )
  return(l)
}


arithmeticsMeanOverListOfMatrix <-function(myList){
  # Computes  the arithmetic mean matrix on an element by element basic
  #
  # Args:
  #   myList : a list of matrixof same dimensions
  #
  # Returns:
  #   A matrix
  Reduce("+", myList) / length(myList)  
}

demoEnvironmentalVariableConstructor <-function(...){
  # Compute a response demographic variable according to a niche function
  # environmental variable
  #
  # Args:
  #   ... : various lists of form list(matrix, function, list()) where
  #         the elements of the third elements are named according to the
  #         names of the arguments of the function called
  #
  # Returns:
  #   A matrix
  myList <- ComputeVariableFunction(...)
  meanLandscape <- arithmeticsMeanOverListOfMatrix(myList) 
  return(meanLandscape)
}
