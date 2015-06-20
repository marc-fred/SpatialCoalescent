writeDataOutputInFile <- function(theta_sigma, theta_Y_k, theta_Y_r, theta_rate, data, file){
  # Writes arguments value and genetic data in file
  #
  # Args:
  #
  # Returns:
  # A file written 
  dir.create(path = paste(getwd(),"/Simulations", sep=""), showWarnings = FALSE)
  
  con <- file(paste("Simulations/", file, sep=""), open = "w")
  writeLines(c("theta_sigma : ", as.character(theta_sigma)), con=con)
  writeLines(c("theta_Y_k : ", as.character(theta_Y_k)), con=con)
  writeLines(c("theta_Y_r : ", as.character(theta_Y_r)), con=con)
  writeLines(c("theta_rate : ", as.character(theta_rate)), con=con)
  writeLines(c("", "GENETICS:", ""), con=con)
  write.table(data, file=con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append=TRUE)
  close(con)
}


readGeneticDataFiles <- function(){
  path <- paste(getwd(), "/Simulations", sep = "")
  allFiles <- grep(pattern = "^genetics_\\d*.txt$", x=list.files(path), value = TRUE)
  allPaths <- paste(getwd(), "/Simulations/", allFiles, sep ="")
  allGenetics <- lapply(X = allPaths, FUN = function(x) readGenetics(file = x))
  return(allGenetics)
}

readGenetics <- function(file){
  skipLine <- which(readLines(file)=="GENETICS:")
  genetics <- read.table(file = file, skip = skipLine)
  return(genetics)
}

writeErrorDataOutputFile <- function(cond, x, theta_sigma, theta_Y_k, theta_Y_r,theta_rate){
  # Writes error and parameters values
  #
  # Args:
  #
  # Returns:
  # A file written 
  errorFileName = paste("stderr_", x , ".txt", sep="")
  con <- file(paste("Simulations/", errorFileName, sep=""), open = "w")
  sink(file = con, type = "message", append =TRUE)
  
  message(c("theta_sigma : ", as.character(theta_sigma)), con=con)
  message(c("theta_Y_k : ", as.character(theta_Y_k)), con=con)
  message(c("theta_Y_r : ", as.character(theta_Y_r)), con=con)
  message(c("theta_rate : ", as.character(theta_rate)), con=con)
  message("Here's the original warning message:")
  write(paste("MY ERROR:", cond), file=con, append=TRUE)
}


modifiedMStatisticsExcoffier2005 <- function(mutationMatrix, LocusInColumns = TRUE){ # mutationMatrix <- genetics[[1]]
  if(LocusInColumns == TRUE){
    margin <- 2
  }else if (LocusInColumns == FALSE){
    margin <- 1
  }else{
    stop("LocusInColumn is not a boolean")
  }
  
  k_l <- apply(X = mutationMatrix,
               MARGIN = margin,
               FUN = function(x) length(unique(x))
  )
  
  r_l <- apply(X = mutationMatrix,
               MARGIN = margin,
               FUN = function(x) max(x) - min(x)
  )
  
  stat <- sum(k_l)/sum(1 + r_l)
  return(stat)
}

expectedMeanHeterozygosity <- function(mutationMatrix, LocusInColumns = TRUE){
  
  if(LocusInColumns == TRUE){
    margin <- 2
    n_ind <- nrow(mutationMatrix)
    n_locus <- ncol(mutationMatrix)
    
  }else if (LocusInColumns == FALSE){
    margin <- 1
    n_ind <- ncol(mutationMatrix)
    n_locus <- nrow(mutationMatrix)
    
  }else{
    stop("LocusInColumn is not a boolean")
  }
    
  p_i <- apply(X = mutationMatrix,
                 MARGIN = margin,
                 FUN = function(x) table(x)/n_ind
                 )
  
  p_i_2 <- lapply(X = p_i,
                  FUN = function(x) x^2
                  )
  
  sum_pi_2 <- lapply(X = p_i_2,
                     FUN = sum)
  
  He <- 1 - mean(unlist(sum_pi_2))
  return(He)
}

computeGeographicalClustersForSample <- function(dataCoord, nbCl = 4, max.iter = 5, tolerance = 0.1){
  # construct orig.data
  Latitude <- dataCoord[,"x"]
  Longitude <- dataCoord[,"y"]
  LocationID <- 1:nrow(dataCoord)
  orig.data <- data.frame(Latitude, Longitude, LocationID)
  
  # Constrained clustering
  cl_constrain = dirichletClusters_constrained(orig.data, k=nbCl, max.iter=max.iter, tolerance=tolerance, plot.iter=TRUE)
  # table( cl_constrain$cluster )
  
  # Construct final data
  cluster <- cl_constrain$cluster
  final <- data.frame(orig.data, cluster)
  return(final)
  ### Plot
  # plot(cl_constrain$Longitude, cl_constrain$Latitude, col=cl_constrain$cluster, pch=16, cex=.6,
  #      xlab="Longitude",ylab="Latitude")
  # library(maps)
  # map("state", add=T)
  # points(cl_constrain$centers[,c(2,1)], pch=4, cex=2, col='orange', lwd=4)
  # plot(final$Longitude, final$Latitude, col=final$cluster, pch=16, cex=.6,
  #      xlab="Longitude",ylab="Latitude")
  
  ### Exemple 
  # size <- 80
  # Latitude <- sample(x = seq(42, 45, by = 0.1), size = size, replace = TRUE)
  # Longitude <- sample(x = seq(1, 7, by = 0.1), size = size, replace = TRUE)
  # LocationID <- 1:size
}

dirichletClusters_constrained = function(orig.data, k, max.iter, tolerance, plot.iter=TRUE) {
  # http://statistical-research.com/spatial-clustering-with-equal-sizes/?utm_source=rss&utm_medium=rss&utm_campaign=spatial-clustering-with-equal-sizes
  
  
  # Convert to radian
  as_radians = function(theta=0){
    return(theta * pi / 180)
  }
  
  calc_dist = function(fr, to) {
    lat1 = as_radians(fr$lat)
    lon1 = as_radians(fr$lon)
    lat2 = as_radians(to$lat)
    lon2 = as_radians(to$lon)
    a = 3963.191;
    b = 3949.903;
    numerator = ( a^2 * cos(lat2) )^2 + ( b^2 * sin(lat2) ) ^2
    denominator = ( a * cos(lat2) )^2 + ( b * sin(lat2) )^2
    radiusofearth = sqrt(numerator/denominator) #Accounts for the ellipticity of the earth.
    d = radiusofearth * acos( sin(lat1) * sin(lat2) + cos(lat1)*cos(lat2)*cos(lon2 - lon1) )
    d.return = list(distance_miles=d)
    return(d.return)
  }
  
  fr = to = NULL
  
  r.k.start = sample(seq(1:k))
  n = nrow( orig.data )
  k.size = ceiling(n/k)
  initial.clusters = rep(r.k.start, k.size)
  
  if(n%%length(initial.clusters)!=0){
    exclude.k = length(initial.clusters) - n%%length(initial.clusters)
  } else {
    exclude.k = 0
  }
  orig.data$cluster = initial.clusters[1:(length(initial.clusters)-exclude.k)]
  orig.data$cluster_original = orig.data$cluster
  
  ## Calc centers and merge
  mu = cbind( by(orig.data$Latitude, orig.data$cluster, mean), by(orig.data$Longitude, orig.data$cluster, mean), seq(1:k) )
  tmp1 = matrix( match(orig.data$cluster, mu[,3]) )
  orig.data.centers = cbind(as.matrix(orig.data), mu[tmp1,])[,c(1:2,4:6)]
  
  ## Calc initial distance from centers
  fr$lat = orig.data.centers[,3]; fr$lon = orig.data.centers[,4]
  to$lat = orig.data.centers[,1]; to$lon = orig.data.centers[,2]
  orig.data$distance.from.center = calc_dist(fr, to)$distance_miles
  orig.data$distance.from.center_original = orig.data$distance.from.center
  
  ## Set some initial configuration values
  is.converged = FALSE
  iteration = 0
  error.old = Inf
  error.curr = Inf
  
  while ( !is.converged && iteration < max.iter ) { # Iterate until threshold or maximum iterations
    
    if(plot.iter==TRUE){
      plot(orig.data$Longitude, orig.data$Latitude, col=orig.data$cluster, pch=16, cex=.6,
           xlab="Longitude",ylab="Latitude")
    }
    iteration = iteration + 1
    start.time = as.numeric(Sys.time())
    cat("Iteration ", iteration,sep="")
    for( i in 1:n ) {
      # Iterate over each observation and measure the distance each observation' from its mean center
      # Produces an exchange. It takes the observation closest to it's mean and in return it gives the observation
      # closest to the giver, k, mean
      fr = to = distances = NULL
      for( j in 1:k ){
        # Determine the distance from each k group
        fr$lat = orig.data$Latitude[i]; fr$lon = orig.data$Longitude[i]
        to$lat = mu[j,1]; to$lon = mu[j,2]
        distances[j] = as.numeric( calc_dist(fr, to) )
      }
      
      # Which k cluster is the observation closest.
      which.min.distance = which(distances==min(distances), arr.ind=TRUE)
      previous.cluster = orig.data$cluster[i]
      orig.data$cluster[i] = which.min.distance # Replace cluster with closest cluster
      
      # Trade an observation that is closest to the giving cluster
      if(previous.cluster != which.min.distance){
        new.cluster.group = orig.data[orig.data$cluster==which.min.distance,]
        
        fr$lat = mu[previous.cluster,1]; fr$lon = mu[previous.cluster,2]
        to$lat = new.cluster.group$Latitude; to$lon = new.cluster.group$Longitude
        new.cluster.group$tmp.dist = calc_dist(fr, to)$distance_miles
        
        take.out.new.cluster.group = which(new.cluster.group$tmp.dist==min(new.cluster.group$tmp.dist), arr.ind=TRUE)
        LocationID = new.cluster.group$LocationID[take.out.new.cluster.group]
        orig.data$cluster[orig.data$LocationID == LocationID] = previous.cluster
      }
      
    }
    
    # Calculate new cluster means
    mu = cbind( by(orig.data$Latitude, orig.data$cluster, mean), by(orig.data$Longitude, orig.data$cluster, mean), seq(1:k) )
    tmp1 = matrix( match(orig.data$cluster, mu[,3]) )
    orig.data.centers = cbind(as.matrix(orig.data), mu[tmp1,])[,c(1:2,4:6)]
    mu = cbind( by(orig.data$Latitude, orig.data$cluster, mean), by(orig.data$Longitude, orig.data$cluster, mean), seq(1:k) )
    
    ## Calc initial distance from centers
    fr$lat = orig.data.centers[,3]; fr$lon = orig.data.centers[,4]
    to$lat = orig.data.centers[,1]; to$lon = orig.data.centers[,2]
    orig.data$distance.from.center = calc_dist(fr, to)$distance_miles
    
    # Test for convergence. Is the previous distance within the threshold of the current total distance from center
    error.curr = sum(orig.data$distance.from.center)
    
    error.diff = abs( error.old - error.curr )
    error.old = error.curr
    if( !is.nan( error.diff ) && error.diff < tolerance ) {
      is.converged = TRUE
    }
    
    # Set a time to see how long the process will take is going through all iterations
    stop.time = as.numeric(Sys.time())
    hour.diff = (((stop.time - start.time) * (max.iter - iteration))/60)/60
    cat("\n Error ",error.diff," Hours remain from iterations ",hour.diff,"\n")
    
    # Write out iterations. Can later be used as a starting point if iterations need to pause
    # write.table(orig.data, paste("C:\\optimize_iteration_",iteration,"_instore_data.csv", sep=""), sep=",", row.names=F)
  }
  
  centers = data.frame(mu)
  ret.val = list("centers" = centers, "cluster" = factor(orig.data$cluster), "LocationID" = orig.data$LocationID,
                 "Latitude" = orig.data$Latitude, "Longitude" = orig.data$Longitude,
                 "k" = k, "iterations" = iteration, "error.diff" = error.diff)
  
  return(ret.val)
  

#   
#   raw.og = read.csv("http://statistical-research.com/wp-content/uploads/2013/11/sample_geo.txt", header=T, sep="\t")
#   orig.data = raw.og[,1:3]
#   # Constrained clustering
#   cl_constrain = dirichletClusters_constrained(orig.data, k=4, max.iter=5, tolerance=.0001, plot.iter=TRUE)
#   table( cl_constrain$cluster )
#   plot(cl_constrain$Longitude, cl_constrain$Latitude, col=cl_constrain$cluster, pch=16, cex=.6,
#        xlab="Longitude",ylab="Latitude")
#   library(maps)
#   map("state", add=T)
#   points(cl_constrain$centers[,c(2,1)], pch=4, cex=2, col='orange', lwd=4)
}

meanNumberofAllelesByLocusByLocality <- function(dtfr){
  nbAllelesByLocality <- by(data = dtfr,
                            INDICES = dtfr$cluster,
                            FUN = function(d){
                              # print(d)
                              apply(X = d[, 5:ncol(d)], MARGIN = 2, FUN = function(x) length(unique(x)))
                              # print(means)
                            })
  
  stat <- sapply(nbAllelesByLocality, mean)
  return(stat)
}

myprint <- function(m) {
  # Print a list of numeric matrix in color levels 
  m <- list(m)
  rl = lapply(m, function(X) raster(X))
  d <- stack(rl)
  spplot(d)
}
 
parallelWrapper <- function(func, numJobs, cores = 2, ...) {
  local({
    # open a connection to a temporary file : pipe between master process and child
    f <- fifo(tempfile(), open="w+b", blocking=T)
    if (inherits(parallel:::mcfork(), "masterProcess")) {
      # Child
      progress <- 0.0
      while (progress < 1 && !isIncomplete(f)) {
        msg <- readBin(f, "double")
        progress <- progress + as.numeric(msg)
        # send a message in C-style
        cat(sprintf("Progress: %.2f%%\n", progress * 100))
      } 
      # close the current child process, informing master process
      parallel:::mcexit()
    }
        
    mclapply(X = 1:numJobs, FUN = function(x, ...){ 
      func(...)
    }, ...,
    mc.cores = cores,
    mc.preschedule = FALSE)
    
  })
  cat("Simulations Done\n")
}