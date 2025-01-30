# Load necessary packages
install.packages(c("future", "future.apply", "progressr", "dplyr", "stringr"), dependencies = T)
library(future)
library(future.apply)
library(dplyr)
library(progressr)
library(stringr)

# Set working directory to main folder

####################################################
###### GLOBAL SENSITIVITY ANALYSIS SPINUP ######
####################################################

# Source slim functions
source("script/R/addt_functions/slim_functions.R")

# How many processes to run at once?
num_processes <- 120

# Read sample design
seeds <- readRDS("outputs/spatial/landscapes/HLSeedsMainTrend.rds")

# get all combinations of spinup seeds
env_seeds <- as.numeric(names(seeds))
slim_seeds <- as.numeric(1:10)
all_spinup <- expand.grid(env_seeds, slim_seeds, stringsAsFactors = FALSE)
colnames(all_spinup) <- c("ENV_SEED", "SLiM_SEED")

# # add all combinations of sensitivity analysis parameters
global_SA <- read.csv("other/sensitivity_analysis_global.csv")

num_global_sensitivity <- 50

# Create empty dataframe to store global SA values in
global_SA_matrix <- matrix(nrow = num_global_sensitivity, 
                        ncol = nrow(global_SA))
sensitivity_samples <- as.data.frame(global_SA_matrix)
colnames(sensitivity_samples) <- global_SA$parameter

# Fill each column with random draws
for (i in 1:nrow(global_SA)) {
  sensitivity_samples[,i] <- runif(num_global_sensitivity, 
                                   min = global_SA$min[i], 
                                   max = global_SA$max[i])
}

# now standardize values of numbers, K must be integer, and add columns to use in IDing what param combo each global SA is
sensitivity_samples <- round(sensitivity_samples, 3)
sensitivity_samples$K <- round(sensitivity_samples$K)
sensitivity_samples$global_round <- 1:num_global_sensitivity
sensitivity_samples$spinup_suffix <- paste0(sensitivity_samples$global_round, "globalSA")

saveRDS(sensitivity_samples, "other/sensitivity_analysis_global_samples.rds")

# create df of all seeds x parameters for sensitivity analysis
all_spinup_run <- merge(all_spinup, sensitivity_samples, by = NULL)

# If it has been run already, remove from queue
file_list <- list.files(path = "outputs/spatial/spinup_files", pattern = "globalSA", full.names = F)
completed_spinup <- data.frame(
  ENV_SEED = as.numeric(str_extract(file_list, "^\\d+")),             # first number
  SLiM_SEED = as.numeric(str_extract(file_list, "(?<=_)\\d+(?=_)")),   # second number
  spinup_suffix = as.character(str_extract(file_list, "\\d+globalSA"))  # third number
) %>%
  unique()

all_spinup_run <- all_spinup_run %>%
  anti_join(completed_spinup, by=c("ENV_SEED", "SLiM_SEED", "spinup_suffix")) %>%
  # also, no need to keep basic spinups (default parameters, same as main)
  filter(!spinup_suffix == "") %>%
  # rename seeds
  rename(envSeed=ENV_SEED, slimSeed=SLiM_SEED)

# Process the data in batches
process_in_batches(all_spinup_run, batch_size = 500, slim_function = run_slim_spinup, num_processes = num_processes)
