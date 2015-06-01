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
  write.table(data, file=con, sep = "\t", quote = FALSE, col.names = FALSE, append=TRUE)
  close(con)
}

readGeneticDataFiles <- function(){
  
  allFiles <- grep(pattern = "^genetics_\\d*.txt$", x=list.files(path), value = TRUE)
  allPaths <- paste(getwd(), "/Simulations/", allFiles, sep ="")
  allGenetics <- lapply(X = allPaths, FUN = readGenetics(x))
  return(genetics)
}

readGenetics <- function(file){
  skipLine <- which(readLines(file)=="GENETICS:")
  genetics <- read.table(file = file, skip = skipLine)
  return(genetics)
}


modifiedMStatisticsExcoffier2005 <- function(mutationMatrix, LocusInColumns = TRUE){
  if(LocusInColumns == TRUE){
    margin <- 2
  }if(LocusInColumns == FALSE){
    margin <- 1
  }else{
    stop("LocusInColumn is not a boolean")
  }
  
  k_l <- apply(X = mutationMatrix,
               MARGIN = margin,
               FUN = length(unique(x))
  )
  
  r_l <- apply(X = mutationMatrix,
               MARGIN = margin,
               FUN = max(x) - min(x)
  )
  
  stat <- sum(k_l)/sum(1 + r_l)
  
}