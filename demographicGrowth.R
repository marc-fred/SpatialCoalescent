envMatrix
kMatrix
rMatrix
migMatrix

# First event : fondator
demographicMatrix <- envMatrix
demographicMatrix[] <- 0
demographicMatrix[1,1] <- 1

#Reproduction
N <- demographicMatrix
kMatrix
rMatrix
N_tilde <- matrix(mapply(FUN=function(N, r, K){ N*(1+r)/(1+(r*N)/K) }, N, rMatrix, kMatrix), ncol = ncol(N))

#Dispersion
matrix(as.vector(N_tilde) %*% migMatrix, ncol=ncol(N))

for?
apply(X = N, MARGIN = c(1,2), FUN=function(x){
  N[i+1] <- (N[i] *(1+r))/(1+(r*N[i])/K)
   
})
