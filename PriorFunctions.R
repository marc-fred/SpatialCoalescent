uniform <- function(x, min, max){
  # A uniform function used to sample parameters in the specified prior distribution
  # 
  # Args:
  #   n: number of observations
  #   min: lower limit of the distribution
  #   max: upper limit of the distribution
  #
  # Returns:
  #   a random sample
  return(runif(n=x, min, max))
}

beta <- function(x, shape1, shape2){
  # random generation for the beta distribution 
  # 
  # Args:
  #   n: number of observations
  #   shape1: non-negative parameter of the Beta distribution
  #   shape2: non-negative parameter of the Beta distribution
  #
  # Returns:
  #   A random sample
  return(rbeta(n=x, shape1, shape2))
}

normal <- function(x, mean, sd){
  # random generation for the normal distribution 
  # 
  # Args:
  #   n: number of observations
  #   mean: vector of means
  #   sd: vector of standard deviations
  #
  # Returns:
  #   A random sample
  return(rnorm(n=x, mean, sd))
}

gamma <- function(x, shape){
  # random generation for the gamma distribution 
  # 
  # Args:
  #   n: number of observations
  #   shape: shape parameter. Must be positive.
  #
  # Returns:
  #   A random sample
  return(rgamma(n=x, shape))
}

