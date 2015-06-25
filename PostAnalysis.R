source("generalFunctions.R")

readParameters <- function(){
  path <- paste(getwd(), "/Simulations", sep = "")
  allFiles <- grep(pattern = "^genetics_\\d*.txt$", x=list.files(path), value = TRUE)
  allPaths <- paste(getwd(), "/Simulations/", allFiles, sep ="")
  allParameters <- lapply(X = allPaths, 
                          FUN = function(x) {
                            lines <- readLines(x)
                            skipLine <- which(lines =="MODEL") + 1
                            rowsNum <- which(lines == "GENETICS") - 2 - skipLine
                            parameters <- read.table(x, skip = skipLine, nrows = rowsNum)
                          })
  return(allParameters)
}

readGeneticDataFiles <- function(){
  path <- paste(getwd(), "/Simulations", sep = "")
  allFiles <- grep(pattern = "^genetics_\\d*.txt$", x=list.files(path), value = TRUE)
  allPaths <- paste(getwd(), "/Simulations/", allFiles, sep ="")
  allGenetics <- lapply(X = allPaths, 
                        FUN = function(x){
                          skipLine <- which(readLines(x)=="GENETICS")
                          genetics <- read.table(file = x, skip = skipLine +1)
                          return(genetics)
                        })
  return(allGenetics)
}

genet <- readGeneticDataFiles()
param <- readParameters()
dataCoord <- read.table(file = paste(getwd(), "/Simulations/dataCoord.txt", sep = ""), header = TRUE)
dataClusters <- computeGeographicalClustersForSample(dataCoord, nbCl = 4, max.iter = 5, tolerance = 0.1)
clusters <- dataClusters$cluster

sumstat <- sapply(
  X = genet,
  FUN = function(x, clusters){
    M <- modifiedMStatisticsExcoffier2005(x)
    names(M) <- "M"
    
    He <- expectedMeanHeterozygosity(x)
    names(He) <- "He"
    
    # Mean number of Alleles
    nbAlleles <- by(data = x,
                    INDICES = clusters,
                    FUN = function(genet){
                      nbAlleles <- apply(X = genet, MARGIN = 2, FUN = function(x) length(unique(x)))
                      return(nbAlleles)
                    }
    )
    meanLoci <- sapply(X = nbAlleles, FUN = mean, USE.NAMES = TRUE)
    names(meanLoci) <- paste("MeanAllelFr.Locality.", names(meanLoci), sep ="")
    
    # Make the form
    sumstat <- c(M,He, meanLoci)
  },
  clusters = clusters
)