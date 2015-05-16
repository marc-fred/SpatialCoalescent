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