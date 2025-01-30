# Load necessary packages
install.packages(c("future", "future.apply", "progressr", "dplyr", "stringr"), dependencies = T)
library(future)
library(future.apply)
library(dplyr)
library(progressr)
library(stringr)

# Set working directory to main folder

####################################################
###### LOCAL SENSITIVITY ANALYSIS SIMULATIONS ######
####################################################

# Source slim functions
source("script/R/addt_functions/slim_functions.R")

# How many processes to run at once?
num_processes <- 120

# Read sample design
seeds <- readRDS("outputs/spatial/landscapes/HLSeedsMainTrend.rds")

# Get all seed combos
seeds_df <- data.frame(
  envSeed = as.integer(rep(names(seeds), sapply(seeds, length))),
  hlSeed = unlist(seeds))
all_seeds <- seeds_df[rep(row.names(seeds_df), each = 10), ]
all_seeds$slimSeed <- rep(1:10, times = nrow(seeds_df))

# Get properties of post habitat loss landscapes
file_list <- list.files(path = "outputs/spatial/landscapes", pattern = "\\_samples.csv$", full.names = TRUE)
data_list <- lapply(file_list, read.csv)
hl_properties <- do.call(rbind, data_list)

# Combine
all_normal <- merge(all_seeds, hl_properties, by = c("envSeed", "hlSeed"), all.x = TRUE)

# # add all combinations of sensitivity analysis parameters
local_SA <- read.csv("other/sensitivity_analysis.csv")

# create df of all seeds x parameters for sensitivity analysis
all_normal_run <- merge(all_normal, local_SA, by = NULL)

## If it has been run already, remove from queue
file_list <- list.files(path = "outputs/spatial/simulation_output", pattern = "sensitivity.txt", full.names = TRUE)
# Create an empty list to store data frames
data_frames <- list()

# Iterate over each file
for (file_name in file_list) {
  # Read the file
  df <- read.csv(file_name)
  
  # Extract part of the file name
  file_part <- gsub(".*out_summary(.+?)\\.txt.*", "\\1", file_name)
  
  # Create a new column with the extracted part of the file name
  df$output_suffix <- file_part
  
  # Append the data frame to the list
  data_frames[[file_name]] <- df
}

# Combine all data frames into one
all_out <- do.call(rbind, data_frames)
all_out <- all_out %>%
    rename(envSeed="ENV_SEED", hlSeed="HL_SEED", slimSeed="SLIM_SEED")
all_normal_run <- all_normal_run %>%
  anti_join(all_out, by=c("output_suffix", "envSeed", "hlSeed", "slimSeed"))

# Process the data in batches, randomize order when placing into batches
process_in_batches(all_normal_run[sample(nrow(all_normal_run)), ], batch_size = 500, slim_function = run_slim_normal, num_processes = num_processes)

