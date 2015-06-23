
spatialCoalescenceForMultipleLoci <- function(migHistory_l, demoHistory_l, localizationData, nbLocus, theta_rate, steps, rasterLandscape){
  # Repeat coalescence simulation, once for each locus
  #
  # Args :
  #
  # Returns :
  #   a matrix of genetic values for each individual (rows) at each locus (column)
  
  # localizationData gives cell number from left to right, then top to bottom
  # raster(demoHistory_l[[length(demoHistory_l)]])[2]
  # c(demoHistory_l[[length(demoHistory_l)]])[2]
  localizationData <- cellFromXY(rasterE1, dataCoord)
  
  # checking for possibility of coalescence : 
  # ACCEPT if sampled individuals in localization data are well simulated by demoHistory_l
  # REJECT if the simulated data did not agree with real data
  r <- raster(demoHistory_l[[length(demoHistory_l)]])
  if( any(sapply(X = unique(localizationData), FUN = function(x) { r[x] == 0}))){
    stop("\n*** Demographic simulation ***\n
         The simulated population did not reach the demes where\n
         real data were sampled. Conclude to aberrant values for\n
         demographicparameters values\n ")
  }

    unique(localizationData)

    
  genetValues <- mapply(FUN = spatialCoalescenceForOneLocus, 
                        steps, theta_rate, 
                        MoreArgs = list(migHistory_l = migHistory_l,
                                        localizationData = localizationData,
                                        demoHistory_l = demoHistory_l),
                        SIMPLIFY = TRUE, USE.NAMES = TRUE)
  
  return(genetValues)
}


spatialCoalescenceForOneLocus <- function(migHistory_l, localizationData, demoHistory_l, stepValue, theta_rate){
  # Simulate the genetic values of various individuals at one locus
  #
  # Args: 
  #
  #
  # Returns :
  #   A vector of genetic values for each individual.
  
  # coalescent informations : (time of coalescence, Child 1, Child 2, Parent)
  coal <- coalescentCore(tipDemes = localizationData, 
                         migHistory_l = migHistory_l, 
                         N_l = demoHistory_l)
  
  # branches informations : (in columns : Child/Parent/Branch length/Number of mutation/Resultant)
  branch <- computeCoalescentBranchesInformation(coal = coal, stepValue = stepValue, mutationRate = theta_rate)
  
  # compute observed genetic values
  genet <- computePresentGeneticValues(branchMat = branch, coal = coal, localizationData = localizationData, initialGeneticValue = 0)
  
  return(genet)
}


coalescentCore <- function(tipDemes, migHistory_l, N_l){
  # Simulate a genealogy backward in the time, accross demes
  # 
  # Args:
  #   tipdDemes: vector of the demes in which each node is found a time 0.
  #   migHistory_l: matrix of transition backward in time
  #   N a vector of population sizes
  #
  # Returns: 
  #   A matrix describing the coalescence events : time/childNode1/childNode2/parentNode
  
  ###### INITIALISATION
  events <- 0
  headNode <- length(tipDemes)
  maxCoalEvent <- length(tipDemes) - 1
  nodesState <- c(tipDemes, rep(NA, maxCoalEvent))
  generationNumber <- length(N_l)
  
  # coalescent informations : (time of coalescence, Child 1, Child 2, Parent)
  coalescent <- matrix(data = NA, nrow = maxCoalEvent, ncol = 4, dimnames = list(c(), c("coalTime", "Child1", "Child2", "Parent")))
  
  #### Loop over all generations backward in time
  for(i_time in generationNumber:1){ # i_time <- 9
    forwardMigHistory_m <- migHistory_l[[i_time]]
    backwardMigHistory_m <- t(forwardMigHistory_m)
    N_m <- N_l[[i_time]]
    N_v <- as.vector(N_m)
    
    nodesState[!is.na(nodesState)] <- vapply(X = nodesState[!is.na(nodesState)],
                                             FUN = sampleParentalDemeInMigrationHistory,
                                             migHistory_m = backwardMigHistory_m,
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
        if(N_v[focalDeme]==0){stop(paste("in coalescentCore you are trying to coalesce in an empty deme : in deme", x ,", N=0"))}
        
        # Attribute parents (among N possible parents) to each node present in the deme
        parents <- sample(N_v[focalDeme], size = length(candidates[[x]]), replace = TRUE) # parents[1] <- parents[2]
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
          coalescent[lines, 1] <- rep(x =  generationNumber - i_time + 1 , times = nEvents)
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

sampleParentalDemeInMigrationHistory <- function(presentDeme, migHistory_m){
  # Function to sample a parental deme knowing the present deme, according to weights
  # given by the migration history
  #
  # Args :
  #   presentDeme : the present deme of the node
  #   migHistory_m : the weight matrix of the migration history, col=present state, row=parent deme
  #
  # Returns :
  #   The parental deme
  #
  parentalDeme <- sample(x = 1:nrow(migHistory_m), size = 1, replace = FALSE, prob = migHistory_m[,presentDeme])
  return(parentalDeme)
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
}


computeCoalescentBranchesInformation <- function(coal, stepValue, mutationRate = theta_rate){
  # Computes the informations along the branches of the coalescent
  #
  # Args :
  #
  #
  # Returns : 
  #   A matrix giving the info along branches : (in columns : Child/Parent/Branch length/Number of mutation/Resultant)
  
  maxCoalEvent <- nrow(coal)
  # Create a matrice for branches (in columns : Child/Parent/Branch length/Number of mutation/Resultant)
  branchMat <- matrix(NA, nrow = (maxCoalEvent)*2, ncol = 5)
  
  # Fill child -> Parent information (decoupling children nodes)
  branchMat[,c(1,2)] <- rbind(coal[,c(2,4)] , coal[,c(3,4)])
  
  timeChild <- mapply(FUN = timeFinder, branchMat[,1], MoreArgs = list(coal = coal))
  
  # time of apparition of parent node
  timeParent <- mapply(FUN = timeFinder, branchMat[,2], MoreArgs = list(coal = coal))
  
  # add branch length
  branchMat[,3] <- timeParent - timeChild        
  
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


computePresentGeneticValues <- function(branchMat, coal, localizationData, initialGeneticValue){
  # Uses branchMat to compute the genetic values of the sample, going down in the coalescent
  #
  # Args :
  #
  #
  # Returns :
  #   a vector of genetic values at the locus for each individuals of the sample.
  
  # add genetic values
  numNodes = length(localizationData)
  maxCoalEvent = numNodes -1
  values <- rep(NA, times = numNodes + maxCoalEvent)
  values[length(values)] <- initialGeneticValue
  for(n in seq(from = length(values)-1, to =1, by = -1 )){
    # find the line of the focal node
    focal <- which(branchMat[,1] == n)
    # find the resultant
    res <- branchMat[focal, 5]
    # find the genetic value of the parent
    values[n] <- values[branchMat[focal, 2]] +res
  }
  return(values[1:numNodes])
}     
