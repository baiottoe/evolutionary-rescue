# load necessary packages
install.packages(c("dplyr", "dbscan"), dependencies = T)
library(dplyr)
library(dbscan)

# Set working directory to main folder

####################################################
############# DENSITY WEIGHTED SAMPLING ############
####################################################

# how many habitat loss scenarios per network property?
n <- 50

# Read in HL sample files
matching_files <- list.files("outputs/spatial/landscapes", pattern = "HL_samples", full.names = TRUE)

all <- data.frame()
for (file in matching_files) {
  # Get all HL scenarios sampled for this environment/landscape
  df <- read.csv(file)
  all <- rbind(all, df)
}
hist(all$envBreadthLoss, na.rm=T)
quantile(all$envBreadthLoss, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
hist(all$envMean, na.rm=T)
quantile(all$envMean, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

HLSeedsMainTrend <- list()

for (file in matching_files) {
  # Get all HL scenarios sampled for this environment/landscape
  HLSamples <- read.csv(file)
  
  HLSeeds <- c()
  envSeed <- HLSamples$envSeed[1]
  
  # Sampling for envBreadth
  s <- HLSamples %>%
    filter(acLen > 0.09 & acLen < 0.11) %>%
    filter(envMean > -0.05 & envMean < 0.05)
  # weighted sampling
  weights = 1/(pointdensity(s$envBreadthLoss, eps=0.1)/nrow(s))
  seeds <- sample(s$hlSeed, size = n, replace = FALSE, prob=weights)
  HLSeeds <- c(HLSeeds, seeds) # add sampled seeds to HLSeeds vector
  HLSamples <- HLSamples %>%
    filter(!(hlSeed %in% seeds)) # remove sampled seeds from HLSamples, so we don't duplicate
  
  # Sampling for envMean
  s <- HLSamples %>%
    filter(acLen > 0.09 & acLen < 0.11) %>%
    filter(envBreadthLoss > 0.75 & envBreadthLoss < 0.85)
  # weighted sampling
  weights = 1/(pointdensity(s$envMean, eps=0.1)/nrow(s))
  seeds <- sample(s$hlSeed, size = n, replace = FALSE, prob=weights)
  HLSeeds <- c(HLSeeds, seeds)
  HLSamples <- HLSamples %>%
    filter(!(hlSeed %in% seeds)) # remove sampled seeds from HLSamples, so we don't duplicate
  
  # Sampling for acLen
  s <- HLSamples %>%
    filter(envMean > -0.05 & envMean < 0.05) %>%
    filter(envBreadthLoss > 0.75 & envBreadthLoss < 0.85)
  # weighted sampling
  weights = 1/(pointdensity(s$acLen, eps=0.01)/nrow(s))
  seeds <- sample(s$hlSeed, size = n, replace = FALSE, prob=weights)
  HLSeeds <- c(HLSeeds, seeds)
  
  # append list of seeds to master list
  HLSeedsMainTrend[[as.character(envSeed)]] <- HLSeeds
  
  print(paste0("done with", file))
}

# Save sampled HL seeds as RDS
saveRDS(HLSeedsMainTrend, file = "outputs/spatial/landscapes/HLSeedsMainTrend.rds")
