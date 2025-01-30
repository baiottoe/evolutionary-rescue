# Load necessary packages
install.packages(c("future", "future.apply", "progressr", "dplyr"), dependencies = T)
library(future)
library(future.apply)
library(dplyr)
library(progressr)

# Set working directory to main folder

####################################################
############ INTERACTION SIMULATIONS SPINUP ########
####################################################

# Source slim functions
source("script/R/addt_functions/slim_functions.R")

# How many processes to run at once?
num_processes <- 40

# Read sample design
seeds <- readRDS("outputs/spatial/landscapes/HLSeedsInteractions.rds")

# Get all seed combos
seeds_df <- data.frame(
  envSeed = as.integer(rep(names(seeds), sapply(seeds, length))),
  hlSeed = unlist(seeds))
all_seeds <- seeds_df[rep(row.names(seeds_df), each = 100), ]
all_seeds$slimSeed <- rep(1:100, times = nrow(seeds_df))

# Get properties of post habitat loss landscapes
file_list <- list.files(path = "outputs/spatial/landscapes", pattern = "\\_samples.csv$", full.names = TRUE)
data_list <- lapply(file_list, read.csv)
hl_properties <- do.call(rbind, data_list)

# Combine
all_normal <- merge(all_seeds, hl_properties, by = c("envSeed", "hlSeed"), all.x = TRUE)

# Run SLiM in parallel for each set of seeds, with CC
all_normal_run <- all_normal %>%
 mutate(spinup_suffix = "",
       output_suffix ="interactions",
       delta_e = 3.0,
       sigma_p = 0.1,
       lambda_o = 0.25,
       sigma_f = 0.25,
       sigma_QTL = 0.1,
       K = 1000,
       l_e = 0.1,
       p_QTL = 0.05,
       p_hl = 2/3)


# Check for existing outputs and filter
file_list <- list.files(path = "outputs/spatial/simulation_output", pattern = "\\_summaryinteractions.txt$", full.names = TRUE)
if (length(file_list) > 0) {
  data_list <- lapply(file_list, read.csv)
  all_out <- do.call(rbind, data_list)
  all_out <- all_out %>%
    rename(envSeed = "ENV_SEED", hlSeed = "HL_SEED", slimSeed = "SLIM_SEED")
  all_normal_run <- all_normal_run %>%
    anti_join(all_out, by = c("envSeed", "hlSeed", "slimSeed"))
}

# Process the data in batches, randomize order when placing into batches
process_in_batches(all_normal_run[sample(nrow(all_normal_run)), ], batch_size = 500, slim_function = run_slim_normal)

