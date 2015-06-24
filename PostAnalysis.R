source("generalFunctions.R")

dt <- readGeneticDataFiles()

# Pas besoin des clusters
M <- sapply(dt, modifiedMStatisticsExcoffier2005)
He <- sapply(dt, expectedMeanHeterozygosity)

# Besoin des clusters
dataCl <- computeGeographicalClustersForSample(dataCoord, nbCl = 4, max.iter = 5, tolerance = 0.1)
dtCl <- lapply(dt, FUN = function(x) data.frame(dataCl, x))

mAll <- sapply(X = dtCl, FUN = meanNumberofAllelesByLocusByLocality)

mAll <- meanNumberofAllelesByLocusByLocality(dtCl[[1]])




