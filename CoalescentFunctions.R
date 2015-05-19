

simulateSpatialCoalescent <- function(theta_sigma, theta_Y_r, theta_Y_k, theta_rate, EnvMatrix, geoDistMatrix, nbLocus, localizationData,steps){
  # Simulate a genetic dataset given parameters
  #
  # Args :
  #
  #
  #
  # Returns : 
  # a matrix of genetic values for each individual (rows) at each locus (column)
  
  kMatrix <- constructEnvironmentalDemographicMatrix(env = envMatrix, param = theta_Y_k)
  rMatrix <- constructEnvironmentalDemographicMatrix(env = envMatrix, param = theta_Y_r)
  migMatrix <- constructMigrationMatrix(dist = geoDistMatrix , param = theta_sigma)
  transitionBackward <- transitionMatrixBackward(r = rMatrix, K = kMatrix, m = migMatrix)
  
  spatialCoalescenceForMultipleLoci(transitionBackward, localizationData, nbLocus, theta_rate, steps)
  

}
 
spatialCoalescenceForMultipleLoci <- function(transitionBackward, kMatrix, localizationData, nbLocus, theta_rate, steps){
  # Repeat coalescence simulation, once for each locus
  #
  # Args :
  #   backMat: the backward transition matrix
  #   coord: the location of the sampled genes
  #   rep: the number of loci
  #
  # Returns :
  #   A matrix of genetic differenciation
  
  list <- vapply(X = 1:nbLocus, 
                 FUN = spatialCoalescenceForOneLocus(transitionBackward, localizationData, kMatrix, stepValue[x]),
                 FUN.VALUE = matrix(1, nrow = nrow(localizationData), ncol = nbLocus),
                 transitionBackward = transitionBackward, 
                 localizationData = localizationData,
                 stepValue = stepValue)
  return(list)
}


spatialCoalescenceForOneLocus <- function(transitionBackward, localizationData, kMatrix, stepValue){
  # Simulate the genetic values of various individuals at one locus
  #
  # Args: 
  #
  #
  # Returns :
  # A vector of genetic values for each individual.
  
  # coalescent informations : (time of coalescence, Child 1, Child 2, Parent)
  coal <- coalescentCore(tipDemes = localizationData, 
                         transitionBackward = transitionBackward, 
                         N = round(as.vector(t(kMatrix))))
  
  # branches informations : (in columns : Child/Parent/Branch length/Number of mutation/Resultant)
  branch <- computeCoalescentBranchesInformation(coal = coal, stepValue = stepValue)
  
  # compute observed genetic values
  genet <- computePresentGeneticValues(branch)
  
  return(genet)
}

computeCoalescentBranchesInformation <- function(coal, stepValue){
  maxCoalEvent <- nrow(coal)
  # Create a matrice for branches (in columns : Child/Parent/Branch length/Number of mutation/Resultant)
  branchMat <- matrix(NA, nrow = (maxCoalEvent)*2, ncol = 5)
  
  # Fill child -> Parent information (decoupling children nodes)
  branchMat[,c(1,2)] <- rbind(coal[,c(2,4)] , coal[,c(3,4)])
  
  timeC <- vapply(X = branchMat[,1],
                  FUN = function(x, coal){
                    # find position in coalescence table
                    line <- which(coal[, 4] == x)
                    if(length(line) == 1){
                      # it's ok : get time
                      t <- coal[line, 1] 
                    }else if(length(line) == 0) {
                      # node is an initial nod (tip node)
                      t <- 0                          
                    }else{
                      # it's really NOT ok
                      stop("error in filling branch lengths in coalescent : 
                           several times of apparition for one node seem to appear : 
                           please verify code")
                    }
                    return(t)
                    },
                  coal = coal,
                  FUN.VALUE = c(1))
  
  # time of apparition of parent node
  timeP <- vapply(X = branchMat[,2],
                  FUN = function(x, coal){
                    # find position in coalescence table
                    line <- which(coal[, 4] == x)
                    if(length(line) == 1){
                      # it's ok : get time
                      t <- coal[line, 1] 
                    }else if(length(line) == 0) {
                      # node is an initial nod (tip node)
                      t <- 0                          
                    }else{
                      # it's really NOT ok
                      stop("error in filling branch lengths in coalescent : 
                           several times of apparition for one node seem to appear : 
                           please verify code")
                    }
                    return(t)
                  },
                  coal = coal,
                  FUN.VALUE = c(1))
  
  # add branch length
  branchMat[,3] <- timeP - timeC        
  
  # add mutation number
  branchMat[,4] <- vapply(X = branchMat[,3],
                          FUN = function(x){rpois(n = 1, lambda = mutationRate*x)},
                          FUN.VALUE = c(1))
  
  # add resultant 
  branchMat[,5] <- resultantFunction(nbrMutations = branchMat[,4],
                                     stepValue = stepValue,
                                     mutationModel = stepWiseMutationModel,
                                     args = c())
  
  return(branchMat)
}

  
computePresentGeneticValues <- function(branchMat, maxCoalEvent){
  # add genetic values
  values <- rep(NA, times = numNodes + maxCoalEvent)
  values[length(values)] <- initialGenetValue
  for(n in seq(from = length(values)-1, to =1, by = -1 )){
    # find the line of the focal node
    focal <- which(branchMat[,1] == n)
    # find the resultant
    res <- branchMat[focal, 5]
    # find the genetic value of the parent
    values[n] <- values[branchMat[focal, 2]] +res
  }
  return(values)
}     


coalescentCore <- function(tipDemes, transitionBackward, N){
  # Simulate a genealogy backward in the time, accross demes
  # 
  # Args:
  #   tipdDemes: vector of the demes in which each node is found a time 0.
  #   transitionBackward: matrix of transition backward in time
  #   N a vector of population sizes
  #
  # Returns: 
  #   A matrix describing the coalescence events : time/childNode1/childNode2/parentNode
  
  ###### INITIALISATION
  time <- 0
  events <- 0
  headNode <- length(tipDemes)
  maxCoalEvent <- length(tipDemes) - 1
  nodesState <- c(tipDemes, rep(NA, maxCoalEvent))
  
  # coalescent informations : (time of coalescence, Child 1, Child 2, Parent)
  coalescent <- matrix(data = NA, nrow = maxCoalEvent, ncol = 4)
  
  ###### REPEAT UNTIL TOTAL COALESCENCE
  while (is.na(tail(nodesState, n=1))){
    time <- time +1
    
    #### MIGRATION
    nodesState[!is.na(nodesState)] <- vapply(X = nodesState[!is.na(nodesState)],
                                             FUN = function(x, N, transitionBackward)
                                             {sample( length(N), size = 1, prob = c(transitionBackward[x,]) )},
                                             N = N, transitionBackward = transitionBackward,
                                             FUN.VALUE = c(1))
        
    ####### CANDIDATES NODES FOR COALESCENCE
    # for active nodes, i.e which are not coded by NA : 
    activeNodes <- which(!is.na(nodesState))
    activeDemes <- nodesState[activeNodes]
    # gives indices of the demes that are duplicated
    dup <- which(duplicated(activeDemes) | duplicated(activeDemes, fromLast= TRUE))
    # gives the demes in which more than one node exist :
    demes <- unique(activeDemes[dup])
    # gives a list of nodes who can perhaps coalesce 
    candidates <- lapply(X = demes,
                         FUN = function(x, nodesState){which(nodesState == x)},
                         nodesState = nodesState)
    
    ####### COALESCENCE
    if(length(candidates) > 0){
      
      for(x in seq(from = 1, to = length(candidates))){ # x <- 1
        
        focalDeme <- demes[x]
        # /!\ If N=0, the nodes would automatically coalesce (parents n°0 for everyone) -> make sure this does not happen !
        if(N[focalDeme]==0){stop(paste("in coalescentCore you are trying to coalesce in an empty deme : in deme", x ,", N=0"))}
        
        # Attribute parents (among N possible parents) to each node present in the deme
        parents <- sample(N[focalDeme], size = length(candidates[[x]]), replace = TRUE) # parents[1] <- parents[2]
        # Test for equality of parents :
        anonymous <- which(duplicated(parents) | duplicated(parents, fromLast= TRUE))
        
        # If nodes have same parent node
        if(length(anonymous) > 1) {
          
          # sample in the candidates nodes who will coalesce
          children <- sample(x = candidates[[x]], size = length(anonymous), replace = FALSE)
          # number of new coalescent events
          nEvents <- length(children) -1
          
          # Move header node forward, and skip ephemeral ones
          headNode <- headNode + nEvents
          # Precise the deme were araised the new node
          nodesState[headNode] <- focalDeme
          # Shut down children nodes
          nodesState[children] <- NA
          
          lines <- seq(from = events+1, to = events + nEvents)
          parentNodes <- seq(from = headNode - nEvents + 1, to = headNode )
          # Fill time
          coalescent[lines, 1] <- rep(x = time, times = nEvents)
          # Fill Child1
          coalescent[lines, 2] <- c(children[1], parentNodes[-length(parentNodes)])
          # Fill Child2
          coalescent[lines, 3] <- c(children[-1])
          # Fill parents
          coalescent[lines, 4] <- parentNodes
          events <- events + nEvents
          
        } # end of if there are coaelescing nodes
      } # end of for loop over demes
    } # end of if there are co occuring nodes in the same deme
  } # end of while coalescence is not complete
  return(coalescent)
}

timeFinder <- function(x, coal){
  # Find time at which a node appeared for the first time. Function used to compute branch length
  #
  # Args:
  #   x : the ID of the node
  #   coal : the table returned by coalescentCore
  #
  # Returns :
  #   The time at which the node x appeared for the first time
  
  # find position in coalescence table
  line <- which(coal[, 4] == id)
  if(length(line) == 1){
    # it's ok : get time
    t <- coal[line, 1] 
  }else if(length(line) == 0) {
    # node is an initial nod (tip node)
    t <- 0                          
  }else{
    # it's really NOT ok
    stop("error in filling branch lengths in coalescent : 
                                 several times of apparition for one node seem to appear : 
                                 please verify code")
  }
  return(t)
}
