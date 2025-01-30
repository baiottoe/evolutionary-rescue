# load necessary packages
install.packages(c("dplyr", "dbscan"), dependencies = T)
library(dplyr)
library(dbscan)

# Set working directory to main folder

####################################################
## INTERACTION SCENARIOS DENSITY WEIGHTED SAMPLING ##
####################################################

# How many habitat loss scenarios per two-way interaction design?
n <- 50

# Get HL sample files
matching_files <- list.files("outputs/spatial/landscapes", pattern = "HL_samples", full.names = TRUE)

all <- data.frame()
for (file in matching_files) {
  # Get all HL scenarios sampled for this environment/landscape
  df <- read.csv(file)
  all <- rbind(all, df)
}

# get quantile values from larger sample pool
hist(all$envBreadthLoss)
quantile(all$envBreadthLoss, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
hist(all$envMean)
quantile(all$envMean, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
hist(all$acLen)
quantile(all$acLen, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

HLSeedsInteractions <- list()

# "Sampling for interactions
for (file in matching_files) {
  # Get all HL scenarios sampled for this environment/landscape
  HLSamples <- read.csv(file)
  
  HLSeeds <- c()
  envSeed <- HLSamples$envSeed[1]
  
  ###### "LOW" END INTERACTIONS ######
  # Sampling for envBreadth - low envMean
  s <- HLSamples %>%
    filter(acLen > 0.09 & acLen < 0.11) %>%
    filter(envMean > -0.17 & envMean < -0.07)
  # weighted sampling
  weights = 1/(pointdensity(s$envBreadthLoss, eps=0.1)/nrow(s))
  seeds <- sample(s$hlSeed, size = n, replace = FALSE, prob=weights)
  HLSeeds <- c(HLSeeds, seeds) # add sampled seeds to HLSeeds vector
  HLSamples <- HLSamples %>%
    filter(!(hlSeed %in% seeds)) # remove sampled seeds from HLSamples, so we don't duplicate
  
  # Sampling for envBreadth - low acLen
  s <- HLSamples %>%
    filter(acLen > 0.05 & acLen < 0.07) %>%
    filter(envMean > -0.05 & envMean < 0.05)
  # weighted sampling
  weights = 1/(pointdensity(s$envBreadthLoss, eps=0.1)/nrow(s))
  seeds <- sample(s$hlSeed, size = n, replace = FALSE, prob=weights)
  HLSeeds <- c(HLSeeds, seeds) # add sampled seeds to HLSeeds vector
  HLSamples <- HLSamples %>%
    filter(!(hlSeed %in% seeds)) # remove sampled seeds from HLSamples, so we don't duplicate
  
  
  # Sampling for envMean - low envBreadthLoss
  s <- HLSamples %>%
    filter(acLen > 0.09 & acLen < 0.11) %>%
    filter(envBreadthLoss > 0.62 & envBreadthLoss < 0.72)
  # weighted sampling
  weights = 1/(pointdensity(s$envMean, eps=0.1)/nrow(s))
  seeds <- sample(s$hlSeed, size = n, replace = FALSE, prob=weights)
  HLSeeds <- c(HLSeeds, seeds)
  HLSamples <- HLSamples %>%
    filter(!(hlSeed %in% seeds)) # remove sampled seeds from HLSamples, so we don't duplicate
  
  # Sampling for envMean - low acLen
  s <- HLSamples %>%
    filter(acLen > 0.05 & acLen < 0.07) %>%
    filter(envBreadthLoss > 0.75 & envBreadthLoss < 0.85)
  # weighted sampling
  weights = 1/(pointdensity(s$envMean, eps=0.1)/nrow(s))
  seeds <- sample(s$hlSeed, size = n, replace = FALSE, prob=weights)
  HLSeeds <- c(HLSeeds, seeds)
  HLSamples <- HLSamples %>%
    filter(!(hlSeed %in% seeds)) # remove sampled seeds from HLSamples, so we don't duplicate
  
  
  # Sampling for acLen - low envBreadthLoss
  s <- HLSamples %>%
    filter(envMean > -0.05 & envMean < 0.05) %>%
    filter(envBreadthLoss > 0.62 & envBreadthLoss < 0.72)
  # weighted sampling
  weights = 1/(pointdensity(s$acLen, eps=0.01)/nrow(s))
  seeds <- sample(s$hlSeed, size = n, replace = FALSE, prob=weights)
  HLSeeds <- c(HLSeeds, seeds)
  HLSamples <- HLSamples %>%
    filter(!(hlSeed %in% seeds)) # remove sampled seeds from HLSamples, so we don't duplicate
  
  # Sampling for acLen - low envMean
  s <- HLSamples %>%
    filter(envMean > -0.17 & envMean < -0.07) %>%
    filter(envBreadthLoss > 0.75 & envBreadthLoss < 0.85)
  # weighted sampling
  weights = 1/(pointdensity(s$acLen, eps=0.01)/nrow(s))
  seeds <- sample(s$hlSeed, size = n, replace = FALSE, prob=weights)
  HLSeeds <- c(HLSeeds, seeds)
  HLSamples <- HLSamples %>%
    filter(!(hlSeed %in% seeds)) # remove sampled seeds from HLSamples, so we don't duplicate
  
  ###### "HIGH" END INTERACTIONS ######
  # Sampling for envBreadth - high envMean
  s <- HLSamples %>%
    filter(acLen > 0.09 & acLen < 0.11) %>%
    filter(envMean > 0.07 & envMean < 0.17)
  # weighted sampling
  weights = 1/(pointdensity(s$envBreadthLoss, eps=0.1)/nrow(s))
  seeds <- sample(s$hlSeed, size = n, replace = FALSE, prob=weights)
  HLSeeds <- c(HLSeeds, seeds) # add sampled seeds to HLSeeds vector
  HLSamples <- HLSamples %>%
    filter(!(hlSeed %in% seeds)) # remove sampled seeds from HLSamples, so we don't duplicate
  
  # Sampling for envBreadth - high acLen
  s <- HLSamples %>%
    filter(acLen > 0.13 & acLen < 0.15) %>%
    filter(envMean > -0.05 & envMean < 0.05)
  # weighted sampling
  weights = 1/(pointdensity(s$envBreadthLoss, eps=0.1)/nrow(s))
  seeds <- sample(s$hlSeed, size = n, replace = FALSE, prob=weights)
  HLSeeds <- c(HLSeeds, seeds) # add sampled seeds to HLSeeds vector
  HLSamples <- HLSamples %>%
    filter(!(hlSeed %in% seeds)) # remove sampled seeds from HLSamples, so we don't duplicate
  
  
  # Sampling for envMean - high envBreadthLoss
  s <- HLSamples %>%
    filter(acLen > 0.09 & acLen < 0.11) %>%
    filter(envBreadthLoss > 0.92 & envBreadthLoss < 1.02)
  # weighted sampling
  weights = 1/(pointdensity(s$envMean, eps=0.1)/nrow(s))
  seeds <- sample(s$hlSeed, size = n, replace = FALSE, prob=weights)
  HLSeeds <- c(HLSeeds, seeds)
  HLSamples <- HLSamples %>%
    filter(!(hlSeed %in% seeds)) # remove sampled seeds from HLSamples, so we don't duplicate
  
  # Sampling for envMean - high acLen
  s <- HLSamples %>%
    filter(acLen > 0.13 & acLen < 0.15) %>%
    filter(envBreadthLoss > 0.75 & envBreadthLoss < 0.85)
  # weighted sampling
  weights = 1/(pointdensity(s$envMean, eps=0.1)/nrow(s))
  seeds <- sample(s$hlSeed, size = n, replace = FALSE, prob=weights)
  HLSeeds <- c(HLSeeds, seeds)
  HLSamples <- HLSamples %>%
    filter(!(hlSeed %in% seeds)) # remove sampled seeds from HLSamples, so we don't duplicate
  
  
  # Sampling for acLen - high envBreadthLoss
  s <- HLSamples %>%
    filter(envMean > -0.05 & envMean < 0.05) %>%
    filter(envBreadthLoss > 0.92 & envBreadthLoss < 1.02)
  # weighted sampling
  weights = 1/(pointdensity(s$acLen, eps=0.01)/nrow(s))
  seeds <- sample(s$hlSeed, size = n, replace = FALSE, prob=weights)
  HLSeeds <- c(HLSeeds, seeds)
  HLSamples <- HLSamples %>%
    filter(!(hlSeed %in% seeds)) # remove sampled seeds from HLSamples, so we don't duplicate
  
  # Sampling for acLen - high envMean
  s <- HLSamples %>%
    filter(envMean > 0.07 & envMean < 0.17) %>%
    filter(envBreadthLoss > 0.75 & envBreadthLoss < 0.85)
  # weighted sampling
  weights = 1/(pointdensity(s$acLen, eps=0.01)/nrow(s))
  seeds <- sample(s$hlSeed, size = n, replace = FALSE, prob=weights)
  HLSeeds <- c(HLSeeds, seeds)
  HLSamples <- HLSamples %>%
    filter(!(hlSeed %in% seeds)) # remove sampled seeds from HLSamples, so we don't duplicate
  
  # append list of seeds to master list
  HLSeedsInteractions[[as.character(envSeed)]] <- HLSeeds
  
  print(paste0("done with", file))
}

# Save HL seeds as RDS
saveRDS(HLSeedsInteractions, file = "outputs/spatial/landscapes/HLSeedsInteractions.rds")
