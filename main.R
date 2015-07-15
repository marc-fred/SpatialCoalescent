rm(list = ls())
# setwd("/home/arno/Documents/These/SpatialCoalescent/Classes")
#setwd("/home/arnaudb/Documents/SpatialCoalescent")
library(methods)
library(raster)
library(parallel)
source("Classes/Generics.R")
source("Classes/Environment.R")
source("Classes/Function.R")
source("Classes/AbstractModel.R")
source("Classes/Model.R")
source("Classes/SuperModel.R")
source("Classes/KModel.R")
source("Classes/RModel.R")
source("Classes/MigModel.R")
source("Classes/RasterLayer.R")
source("CoalescentFunctions.R")
source("NicheFunctions.R")
source("DispersionFunctions.R")
source("MutationFunctions.R")
source("PriorFunctions.R")
source("MarkovProcess.R")
source("generalFunctions.R")
source("demographicGrowth.R")

# "Real" Data
rasterE1 <- raster(x = matrix(data = sample(1:100, 9), ncol = 3),
                   xmn = 40, xmx = 50, ymn = 0, ymx = 10, crs = "+proj=longlat +datum=WGS84")

rasterE2 <- raster(x = matrix(data = sample(1:100, 9), ncol = 3),
                   xmn = 40, xmx = 50, ymn = 0, ymx = 10, crs = "+proj=longlat +datum=WGS84")

# Environment definition
pluie <- new("Environment", values= as.matrix(rasterE1))
temp <- new("Environment", values= as.matrix(rasterE2))
distances <- new("Lattice", values= computeDistanceMatrix(rasterE1))

# Model Implementation
prior1 <- Function(fun = uniform, param = list(min = 10, max = 50))
prior2 <- Function(fun = uniform, param = list(min = 5, max = 10))
mk1 <- model(varEnv = pluie, fun = Function(fun = linearTwoParameters, param = list(prior1, prior2)))

prior3 <- Function(fun = uniform, param = list(min = 100, max = 500))
prior4 <- Function(fun = uniform, param = list(min = 50, max = 60))
mk2 <- model(varEnv = temp, fun = Function(fun = linearTwoParameters, param = list(prior3, prior4)))

Kmodel <- new("KModel", models = list(mk1, mk2))

prior5 <- Function(fun = uniform, param = list(min = 2, max = 10))
prior6 <- Function(fun = uniform, param = list(min = 10, max = 50))
mr1    <- model(varEnv = pluie, fun = Function(fun = linearTwoParameters, param = list(prior5, prior6)))
Rmodel <- new("RModel", models = list(mr1))

prior7 <- Function(fun = uniform, param = list(min = 0, max = 1000))
prior8 <- Function(fun = uniform, param = list(min = 0, max = 10000))
mig1   <- model(varEnv = distances, fun = Function(fun = gaussianDisp, param = list(prior7, prior8)))
migModel <- new("MigModel", models = list(mig1))

# Simulation
K_m <- applyModel(Kmodel)
R_m <- applyModel(Rmodel)
M_m <- applyModel(migModel)

# Construct the vector of names for transition matrix
lifeLength <- 3
mat <- as.matrix(rasterE1)
demeNames <- paste("D", 1:length(mat), sep = ".")
adultNames <- paste("Adult", 1:lifeLength, sep =".")
lifeNames <- c("Eggs", "Larva", adultNames)
statesNames_m <- sapply(X = demeNames,
                        FUN = function(x, lifeNames){paste(x, lifeNames, sep =".")},
                        lifeNames = lifeNames
)
statesNames <- c(as.vector(statesNames_m), "Cimetery")

# Construct the transition matrix
transition <- matrix(0, ncol = length(statesNames), nrow = length(statesNames), 
                     dimnames = list(statesNames, statesNames))

### Fill the cimetery
nr <- nrow(transition)
# Death of old age
transition[, nr] <- c(rep(c(rep(0, times = length(lifeNames) - 1), 1), times = length(demeNames)), NA)
# No resurrection : stay dead please.
transition[nr,] <- c(rep(0, times = nr -1), 1)

### Fill the age descent with a smaller repetitive matrix
pattern <- matrix(NA, ncol = length(lifeNames), nrow = length(lifeNames))
pattern <- as(pattern, "CsparseMatrix")
# No Fountain of Youth
pattern <- triu(pattern)
# No Premature Ageing
pattern <- tril(pattern, 1)
# Limited life spend
pattern[nrow(pattern), ncol(pattern)] <- 0

# fill for the migration


demoInit <- createInitialDemographicsLandscape(K_m)

history <- demographicSimulation(numberOfGenerations = 20,
                                 demographicMatrix = demoInit, 
                                 kMatrix = K_m, 
                                 rMatrix = R_m,
                                 migMatrix = M_m)

genetValues <- spatialCoalescenceForMultipleLoci(migHistory_l = history$migHistory,
                                                 demoHistory_l = history$demoHistory, 
                                                 localizationData = localizationData, 
                                                 nbLocus = nbLocus, 
                                                 theta_rate = theta_rate,
                                                 steps = steps,
                                                 rasterLandscape = rasterE1)

fileName = paste("genetics_", x , ".txt", sep="")
con <- file(paste("Simulations/", fileName, sep=""), open = "w")
writeLines(text = "MODEL\n", con = con)
writeLines(text = c(getParameters(Kmodel),
                    getParameters(Rmodel),
                    getParameters(migModel)), con = con, sep = "\n")

param_v <- mapply(FUN = function(names, theta_rate){
  paste("Locus.",names,".MutationRate"," ", theta_rate,"\n", sep ="")
},
names = 1:length(theta_rate),
theta_rate = theta_rate)
writeLines(text = param_v, con = con, sep ="")

writeLines("\nGENETICS\n", con=con)
write.table(genetValues, file=con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append=TRUE)
close(con)

)