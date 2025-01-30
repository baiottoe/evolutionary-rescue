# load necessary packages
install.packages(c("raster", "foreach", "doParallel", "cubature"), dependencies = T)
library(raster)
library(foreach)
library(doParallel)

# Set working directory to main folder

####################################################
############# ENV LANDSCAPE GENERATION #############
####################################################

# Read in landscape generation functions
source("script/R/addt_functions/landscape_ac.R")

# Generate 50 random, spatially autocorrelated landscapes
# Periodic along x and y axes
envSeed <- 1
envMaps <- list()
while (length(envMaps) < 50) {
  # generate landscape
  env <- generateLandscape(envSeed, 0.0, 0.0, 1.0, 0.1, minSize=256, periodic_x=T)
  
  # add map to list
  envMaps[[as.character(envSeed)]] <- env
  
  # new seed, next integer
  envSeed = envSeed + 1
}

# get vector of all environment seeds used that meet criteria
envMapSeeds <- names(envMaps)
envMapSeedsInt <- as.integer(envMapSeeds)

for (envMapSeed in envMapSeeds){
  # copy environment
  env <- envMaps[[envMapSeed]]
  
  # Save landscape, COLUMN NAMES STILL BEING SAVED?
  write.csv(env, paste0("outputs/spatial/landscapes/", envMapSeed, "_", "ENV.csv"), row.names = FALSE, col.names = FALSE)
}

####################################################
############# HABITAT LOSS GENERATION ##############
####################################################
HLSampling <- function(envSeed, hlSeed) {
  propLoss <- 0.667
  
  envCopy <- env
  
  # prior conditions
  envBreadthENV <- round(quantile(envCopy, 0.975, na.rm=T) - quantile(envCopy, 0.025, na.rm=T), digits=3)
  
  # select autocorrelation length of HL landscape (patchiness)
  acLen <- round(runif(1, 0.02, 0.18), digits=3)
  
  # generate HL landscape
  hl <- generateLandscape(hlSeed, 0.0, 0.0, 1.0, acLen, minSize=256, periodic_x=T)
  
  # Calculate the 66.7th percentile
  quantileCut <- quantile(hl, propLoss)
  
  # Replace values below 66th percentile with 0 and above with 1 (0 means HL)
  hl[hl >= quantileCut] <- 1
  hl[hl < quantileCut] <- 0
  
  # Set env to NA where hl == 0
  envCopy[hl == 0] <- NA
  
  # Convert env matrix to raster
  envRast <- raster::raster(envCopy)
  
  # Get mean environment in remaining landscape, equivalent to After-prior since mean(prior)=0.0
  envMean <- round(mean(envCopy, na.rm=T), digits=2)
  
  # Calculate environmental breadth (97.5 percentile minus 2.5 percentile)
  envBreadth <- round(quantile(envCopy, 1-(0.025/(1-propLoss)), na.rm=T) - quantile(envCopy, 0.025/(1-propLoss), na.rm=T), digits=3)
  
  breadthLoss <- round((envBreadthENV - envBreadth), digits=3)
  
  result <- list(envSeed=envSeed, hlSeed=hlSeed, envMean=envMean, envBreadthLoss=breadthLoss, acLen=acLen)
  return(result)
}

for (envMapSeed in envMapSeeds) {
  envSeedInt <- as.integer(envMapSeed)
  
  # copy environment
  env <- envMaps[[envMapSeed]]
  
  # get combinations of 
  allSeeds <- list()
  # new seeds for HL; no overlap with previous
  for (hlSeed in (max(envMapSeedsInt)+1):(100000+max(envMapSeedsInt))) {
    myArgs <- list(envSeed = envSeedInt, hlSeed = hlSeed)
      
    # Add the named list to the result list
    allSeeds[[paste0(envSeedInt, "_", hlSeed)]] <- myArgs
  }
  
  # Register a parallel backend
  num_cores <- min(32, detectCores())
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Initialize an empty data frame to store results
  samplingDfFull <- data.frame(envSeed = numeric(0), hlSeed = numeric(0), envMean = numeric(0), envBreadthLoss = numeric(0), acLen = numeric(0))
  
  # Run the function in parallel and append results to the data frame
  samplingDfFull <- foreach(args = allSeeds, .combine = rbind) %dopar% {
    result <- HLSampling(args$envSeed, args$hlSeed)
    data.frame(envSeed = result$envSeed, hlSeed = result$hlSeed, envMean = result$envMean, envBreadthLoss = result$envBreadthLoss, acLen = result$acLen)
  }
  
  # Stop the parallel backend
  stopCluster(cl)
  
  # Optionally, save results to a CSV file
  write.csv(samplingDfFull, paste0("outputs/spatial/landscapes/", envMapSeed, "_", "HL_samples.csv"), row.names = FALSE)
}
# Time for 10000 in parallel, 10 cores: 380 seconds
# time for 100 in sequence, 1 core: 33 seconds