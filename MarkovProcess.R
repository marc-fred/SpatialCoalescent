transitionMatrixBackward <- function(r,K, m){
  # Compute the backward probabilities of migration obtained with an isotropic migration hypothesis
  #
  # Args:
  #   r: a vector a values for growth rate for the various cells
  #   K: a vector of values for carrying capacities of the cells
  #   migration: a migration matrix
  #
  # Returns:
  #   The backward probabilities matrix
  
  kVec <- as.vector(t(K))
  rVec <- as.vector(t(r))
  temp <- apply(t(m) , MARGIN = 1, FUN = function(x, k, r){ k*r*x}, k=kVec, r=rVec)
  transition <- temp/(rowSums(temp))
  return(transition)
}
