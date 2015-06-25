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

dt <- readGeneticDataFiles()
param <- readParameters()
dataCoord <- read.table(file = paste(getwd(), "/Simulations/dataCoord.txt", sep = ""),
                        header = TRUE)

# Pas besoin des clusters
M <- sapply(dt, modifiedMStatisticsExcoffier2005)
He <- sapply(dt, expectedMeanHeterozygosity)

# Besoin des clusters
dataCl <- computeGeographicalClustersForSample(dataCoord, nbCl = 4, max.iter = 5, tolerance = 0.1)
dtCl <- lapply(dt, FUN = function(x) data.frame(dataCl, x))

mAll <- sapply(X = dtCl, FUN = meanNumberofAllelesByLocusByLocality)

mAll <- meanNumberofAllelesByLocusByLocality(dtCl[[1]])




