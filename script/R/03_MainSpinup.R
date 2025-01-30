# load necessary packages
install.packages(c("future", "future.apply", "progressr", "dplyr"), dependencies = T)
library(future)
library(future.apply)
library(dplyr)
library(progressr)

# Set working directory to main folder

####################################################
############# SIMULATION SPINUP ############
####################################################

# Source slim functions
source("script/R/addt_functions/slim_functions.R")

# How many processes to run at once?
num_processes <- 40

# Read sample design
seeds <- readRDS("outputs/spatial/landscapes/HLSeedsMainTrend.rds")

# get all combinations of spinup seeds
env_seeds <- as.numeric(names(seeds))
slim_seeds <- as.numeric(1:100)
all_spinup <- expand.grid(env_seeds, slim_seeds, stringsAsFactors = FALSE)
colnames(all_spinup) <- c("ENV_SEED", "SLiM_SEED")

# set up all param values to run
all_spinup_run <- all_spinup %>%
  mutate(spinup_suffix = "", sigma_p = 0.1, lambda_o = 0.25, sigma_f=0.25, sigma_QTL=0.1, K=1000, l_e=0.1, p_QTL=0.05)

# Spinup in batches
process_in_batches(all_spinup_run, batch_size = 500, slim_function = run_slim_spinup, num_processes = num_processes)
