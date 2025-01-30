# load necessary packages
install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.5-0.tar.gz", repos=NULL, type="source")
install.packages(c("dplyr", "brms", "modelr", "future", "parallel"), dependencies = T)
library(dplyr) 
library(brms) 
library(modelr)
library(future) 
library(parallel)

# Set working directory to main folder

####################################################
###### STATS MODELS FOR ANALYSIS+VISUALIZATION #####
####################################################

##### Setup #####
# Set up hierarchical parallelization for faster processing
plan(
  list(
    tweak(multisession, workers = 6),
    tweak(multisession, workers = 4)
  )
)

# Define function to run brms model with future - note: temporarily removed this
run_my_models <- function(myModel, myData, myPrior, myFile) {
  brm(myModel,
      data = myData,
      family = bernoulli(),
      cores = 4,
      prior = myPrior,
      control = list(adapt_delta = 0.995, max_treedepth=20),
      iter = 2500,
      future = T,
      file=myFile)
}

# custom function to extract certain output from all
get_specific_output <- function(seeds, output, lower, upper, level) {
  model_seeds <- lapply(seeds, function(x) x[lower:upper])
  model_seeds_df <- data.frame(
    envSeed = as.integer(rep(names(model_seeds), sapply(model_seeds, length))),
    hlSeed = unlist(model_seeds)
  )
  data_interaction <- merge(output, model_seeds_df, by.x = c("ENV_SEED", "HL_SEED"), by.y = c("envSeed", "hlSeed")) %>%
    mutate(level = level)
  return(data_interaction)
}

# Get properties of post habitat loss landscapes
# Normal (HL+CC)
file_list <- list.files(path = "outputs/spatial/simulation_output", pattern = "\\_summaryCC3.txt$", full.names = TRUE)
data_list <- lapply(file_list, read.csv)
all_out <- do.call(rbind, data_list)

# Group outputs by organism type, then get summary statistics for each PA design by looking at results from each seed
grouped_data <- all_out %>%
  dplyr::group_by(ENV_SEED, HL_SEED)

# Convert to binary outcome (persist or not)
result <- grouped_data %>%
  mutate(persist = ifelse(proportion_original_pop_size==0, 0, 1))

# Read sample design
seeds <- readRDS("outputs/spatial/landscapes/HLSeedsMainTrend.rds")

# Get all seed combos
seeds_df <- data.frame(
  envSeed = as.integer(rep(names(seeds), sapply(seeds, length))),
  hlSeed = unlist(seeds))

# Get properties of post habitat loss landscapes
file_list <- list.files(path = "outputs/spatial/landscapes", pattern = "\\_samples.csv$", full.names = TRUE)
data_list <- lapply(file_list, read.csv)
hl_properties <- do.call(rbind, data_list)

# Combine to get properties for all seeds
all_normal <- merge(seeds_df, hl_properties, by = c("envSeed", "hlSeed"), all.x = TRUE)

# Combine results with properties
results_and_properties <- merge(result, all_normal, by.x = c("ENV_SEED", "HL_SEED"), by.y = c("envSeed", "hlSeed"), all.x = TRUE)

# Which seeds for which properties? then extract that data
# ENVIRONMENTAL Breadth LOSS
envBreadthLoss_data <- get_specific_output(seeds, results_and_properties, 1, 50, 1) %>%
  mutate(scaled_envBreadthLoss = scale(envBreadthLoss)) #scale

# ENVIRONMENTAL MEAN
envMean_data <-  get_specific_output(seeds, results_and_properties, 51, 100, 1) %>%
  mutate(scaled_envMean = scale(envMean)) #scale

# FRAGMENTATION
acLen_data <-  get_specific_output(seeds, results_and_properties, 101, 150, 1) %>%
  mutate(scaled_acLen = scale(acLen)) #scale

##### Calculate Scaling Quantities #####
envBreadthLoss_mean <- mean(envBreadthLoss_data$envBreadthLoss)
envBreadthLoss_sd <- sd(envBreadthLoss_data$envBreadthLoss)
envMean_mean <- mean(envMean_data$envMean)
envMean_sd <- sd(envMean_data$envMean)
acLen_mean <- mean(acLen_data$acLen)
acLen_sd <- sd(acLen_data$acLen)


##### Models, Main #####
# fixed effect
bf_envBreadthLoss_fixed <- bf(persist ~ scaled_envBreadthLoss)
bf_envMean_fixed <- bf(persist ~ scaled_envMean)
bf_acLen_fixed <- bf(persist ~ scaled_acLen)

# random effect of environmental landscape
bf_envBreadthLoss_random <- bf(persist ~ scaled_envBreadthLoss + (scaled_envBreadthLoss|ENV_SEED))
bf_envMean_random <- bf(persist ~ scaled_envMean + (scaled_envMean|ENV_SEED))
bf_acLen_random <- bf(persist ~ scaled_acLen + (scaled_acLen|ENV_SEED))

## Priors
# Fixed
envBreadthLoss_fixed_prior <- prior(normal(0,2),class = 'b', coef = scaled_envBreadthLoss) +
  prior(student_t(5, 0, 3), class = 'Intercept')

envMean_fixed_prior <- prior(normal(0,2),class = 'b', coef = scaled_envMean) +
  prior(student_t(5, 0, 3), class = 'Intercept')

acLen_fixed_prior <- prior(normal(0,2),class = 'b', coef = scaled_acLen) +
  prior(student_t(5, 0, 3), class = 'Intercept')

# Random
envBreadthLoss_random_prior <- prior(normal(0,2),class = 'b', coef = scaled_envBreadthLoss) +
  prior(student_t(5, 0, 3), class = 'Intercept') +
  prior(student_t(6, 0, 2.9), class = 'sd')

envMean_random_prior <- prior(normal(0,2),class = 'b', coef = scaled_envMean) +
  prior(student_t(5, 0, 3), class = 'Intercept') +
  prior(student_t(6, 0, 2.9), class = 'sd')

acLen_random_prior <- prior(normal(0,2),class = 'b', coef = scaled_acLen) +
  prior(student_t(5, 0, 3), class = 'Intercept') +
  prior(student_t(6, 0, 2.9), class = 'sd')

## Run!
# fixed
fits1 %<-% run_my_models(bf_envBreadthLoss_fixed, envBreadthLoss_data, envBreadthLoss_fixed_prior, "outputs/stats/envBreadthLoss_fixed_scaled.rds")
fits2 %<-% run_my_models(bf_envMean_fixed, envMean_data, envMean_fixed_prior, "outputs/stats/envMean_fixed_scaled.rds")
fits3 %<-% run_my_models(bf_acLen_fixed, acLen_data, acLen_fixed_prior, "outputs/stats/acLen_fixed_scaled.rds")

# random
fits4 %<-% run_my_models(bf_envBreadthLoss_random, envBreadthLoss_data, envBreadthLoss_random_prior, "outputs/stats/envBreadthLoss_random_scaled.rds")
fits5 %<-% run_my_models(bf_envMean_random, envMean_data, envMean_random_prior, "outputs/stats/envMean_random_scaled.rds")
fits6 %<-% run_my_models(bf_acLen_random, acLen_data, acLen_random_prior, "outputs/stats/acLen_random_scaled.rds")



##### Read in Data, Interactions #####
# Get properties of post habitat loss landscapes
file_list <- list.files(path = "outputs/spatial/simulation_output", pattern = "\\_summaryinteractions.txt$", full.names = TRUE)
data_list <- lapply(file_list, read.csv)
all_out <- do.call(rbind, data_list)


# Group outputs by organism type, then get summary statistics for each PA design by looking at results from each seed
grouped_data_interactions <- all_out %>%
  dplyr::group_by(ENV_SEED, HL_SEED)

result_interactions <- grouped_data_interactions %>%
  mutate(persist = ifelse(proportion_original_pop_size==0, 0, 1))

# Read sample design
seeds_interactions <- readRDS("outputs/spatial/landscapes/HLSeedsInteractions.rds")

# Get all seed combos
seeds_df_interactions <- data.frame(
  envSeed = as.integer(rep(names(seeds_interactions), sapply(seeds_interactions, length))),
  hlSeed = unlist(seeds_interactions))

# Get properties of post habitat loss landscapes
file_list <- list.files(path = "outputs/spatial/landscapes", pattern = "\\_samples.csv$", full.names = TRUE)
data_list <- lapply(file_list, read.csv)
hl_properties <- do.call(rbind, data_list)

# Combine to get properties for all seeds
all_normal_interactions <- merge(seeds_df_interactions, hl_properties, by = c("envSeed", "hlSeed"), all.x = TRUE)

# Combine results with properties
results_and_properties_interactions <- merge(result_interactions, all_normal_interactions, by.x = c("ENV_SEED", "HL_SEED"), by.y = c("envSeed", "hlSeed"), all.x = TRUE)

# Which seeds for which properties? then extract that data
# ENVIRONMENTAL Breadth LOSS
envBreadthLoss_data_interaction_low_envMean <- get_specific_output(seeds_interactions, results_and_properties_interactions, 1, 50, 2) %>%
  mutate(scaled_envBreadthLoss = (envBreadthLoss-envBreadthLoss_mean)/envBreadthLoss_sd, scaled_envMean = (envMean-envMean_mean)/envMean_sd) #scale
envBreadthLoss_data_interaction_high_envMean <- get_specific_output(seeds_interactions, results_and_properties_interactions, 301, 350, 3) %>%
  mutate(scaled_envBreadthLoss = (envBreadthLoss-envBreadthLoss_mean)/envBreadthLoss_sd, scaled_envMean = (envMean-envMean_mean)/envMean_sd) #scale
envBreadthLoss_data_interaction_low_acLen <- get_specific_output(seeds_interactions, results_and_properties_interactions, 51, 100, 2) %>%
  mutate(scaled_envBreadthLoss = (envBreadthLoss-envBreadthLoss_mean)/envBreadthLoss_sd, scaled_acLen = (acLen-acLen_mean)/acLen_sd) #scale
envBreadthLoss_data_interaction_high_acLen <- get_specific_output(seeds_interactions, results_and_properties_interactions, 351, 400, 3) %>%
  mutate(scaled_envBreadthLoss = (envBreadthLoss-envBreadthLoss_mean)/envBreadthLoss_sd, scaled_acLen = (acLen-acLen_mean)/acLen_sd) #scale

envBreadthLoss_data_interaction_envMean <- bind_rows(envBreadthLoss_data, envBreadthLoss_data_interaction_low_envMean, envBreadthLoss_data_interaction_high_envMean)
envBreadthLoss_data_interaction_acLen <- bind_rows(envBreadthLoss_data, envBreadthLoss_data_interaction_low_acLen, envBreadthLoss_data_interaction_high_acLen)


# ENVIRONMENTAL MEAN
envMean_data_interaction_low_envBreadthLoss <- get_specific_output(seeds_interactions, results_and_properties_interactions, 101, 150, 2) %>%
  mutate(scaled_envMean = (envMean-envMean_mean)/envMean_sd, scaled_envBreadthLoss = (envBreadthLoss-envBreadthLoss_mean)/envBreadthLoss_sd) #scale
envMean_data_interaction_high_envBreadthLoss <- get_specific_output(seeds_interactions, results_and_properties_interactions, 401, 450, 3) %>%
  mutate(scaled_envMean = (envMean-envMean_mean)/envMean_sd, scaled_envBreadthLoss = (envBreadthLoss-envBreadthLoss_mean)/envBreadthLoss_sd) #scale
envMean_data_interaction_low_acLen <- get_specific_output(seeds_interactions, results_and_properties_interactions, 151, 200, 2) %>%
  mutate(scaled_envMean = (envMean-envMean_mean)/envMean_sd, scaled_acLen = (acLen-acLen_mean)/acLen_sd) #scale
envMean_data_interaction_high_acLen <- get_specific_output(seeds_interactions, results_and_properties_interactions, 451, 500, 3) %>%
  mutate(scaled_envMean = (envMean-envMean_mean)/envMean_sd, scaled_acLen = (acLen-acLen_mean)/acLen_sd) #scale

envMean_data_interaction_envBreadthLoss <- bind_rows(envMean_data, envMean_data_interaction_low_envBreadthLoss, envMean_data_interaction_high_envBreadthLoss)
envMean_data_interaction_acLen <- bind_rows(envMean_data, envMean_data_interaction_low_acLen, envMean_data_interaction_high_acLen)


# AC LEN
acLen_data_interaction_low_envBreadthLoss <- get_specific_output(seeds_interactions, results_and_properties_interactions, 201, 250, 2) %>%
  mutate(scaled_acLen = (acLen-acLen_mean)/acLen_sd, scaled_envBreadthLoss = (envBreadthLoss-envBreadthLoss_mean)/envBreadthLoss_sd) #scale
acLen_data_interaction_high_envBreadthLoss <- get_specific_output(seeds_interactions, results_and_properties_interactions, 501, 550, 3) %>%
  mutate(scaled_acLen = (acLen-acLen_mean)/acLen_sd, scaled_envBreadthLoss = (envBreadthLoss-envBreadthLoss_mean)/envBreadthLoss_sd) #scale
acLen_data_interaction_low_envMean <- get_specific_output(seeds_interactions, results_and_properties_interactions, 251, 300, 2) %>%
  mutate(scaled_acLen = (acLen-acLen_mean)/acLen_sd, scaled_envMean = (envMean-envMean_mean)/envMean_sd) #scale
acLen_data_interaction_high_envMean <- get_specific_output(seeds_interactions, results_and_properties_interactions, 551, 600, 3) %>%
  mutate(scaled_acLen = (acLen-acLen_mean)/acLen_sd, scaled_envMean = (envMean-envMean_mean)/envMean_sd) #scale

acLen_data_interaction_envBreadthLoss <- bind_rows(acLen_data, acLen_data_interaction_low_envBreadthLoss, acLen_data_interaction_high_envBreadthLoss)
acLen_data_interaction_envMean <- bind_rows(acLen_data, acLen_data_interaction_low_envMean, acLen_data_interaction_high_envMean)

##### Models, Interactions #####
# fixed effect model; random effects model has trouble completing
bf_envBreadthLoss_envMean_fixed <- bf(persist ~ scaled_envBreadthLoss + scaled_envMean + scaled_envBreadthLoss*scaled_envMean)
bf_envBreadthLoss_acLen_fixed <- bf(persist ~ scaled_envBreadthLoss + scaled_acLen + scaled_envBreadthLoss*scaled_acLen)
bf_envMean_envBreadthLoss_fixed <- bf(persist ~ scaled_envMean + scaled_envBreadthLoss + scaled_envMean*scaled_envBreadthLoss)
bf_envMean_acLen_fixed <- bf(persist ~ scaled_envMean + scaled_acLen + scaled_envMean*scaled_acLen)
bf_acLen_envBreadthLoss_fixed <- bf(persist ~ scaled_acLen + scaled_envBreadthLoss + scaled_acLen*scaled_envBreadthLoss)
bf_acLen_envMean_fixed <- bf(persist ~ scaled_acLen + scaled_envMean + scaled_acLen*scaled_envMean)


# Priors
# Fixed
envBreadthLoss_envMean_fixed_prior <- prior(normal(0,2),class = 'b', coef = scaled_envBreadthLoss) +
  prior(normal(0,2),class = 'b', coef = scaled_envMean) +
  prior(student_t(5, 0, 3), class = 'Intercept')
envBreadthLoss_acLen_fixed_prior <- prior(normal(0,2),class = 'b', coef = scaled_envBreadthLoss) +
  prior(normal(0,2),class = 'b', coef = scaled_acLen) +
  prior(student_t(5, 0, 3), class = 'Intercept')
envMean_envBreadthLoss_fixed_prior <- prior(normal(0,2),class = 'b', coef = scaled_envMean) +
  prior(normal(0,2),class = 'b', coef = scaled_envBreadthLoss) +
  prior(student_t(5, 0, 3), class = 'Intercept')
envMean_acLen_fixed_prior <- prior(normal(0,2),class = 'b', coef = scaled_envMean) +
  prior(normal(0,2),class = 'b', coef = scaled_acLen) +
  prior(student_t(5, 0, 3), class = 'Intercept')
acLen_envBreadthLoss_fixed_prior <- prior(normal(0,2),class = 'b', coef = scaled_acLen) +
  prior(normal(0,2),class = 'b', coef = scaled_envBreadthLoss) +
  prior(student_t(5, 0, 3), class = 'Intercept')
acLen_envMean_fixed_prior <- prior(normal(0,2),class = 'b', coef = scaled_acLen) +
  prior(normal(0,2),class = 'b', coef = scaled_envMean) +
  prior(student_t(5, 0, 3), class = 'Intercept')

# Run!
fits7 %<-% run_my_models(bf_envBreadthLoss_envMean_fixed, unique(envBreadthLoss_data_interaction_envMean), envBreadthLoss_envMean_fixed_prior, "outputs/stats/envBreadthLoss_fixed_interaction_envMean_scaled.rds")
fits8 %<-% run_my_models(bf_envBreadthLoss_acLen_fixed, unique(envBreadthLoss_data_interaction_acLen), envBreadthLoss_acLen_fixed_prior, "outputs/stats/envBreadthLoss_fixed_interaction_acLen_scaled.rds")
fits9 %<-% run_my_models(bf_envMean_envBreadthLoss_fixed, unique(envMean_data_interaction_envBreadthLoss), envMean_envBreadthLoss_fixed_prior, "outputs/stats/envMean_fixed_interaction_envBreadthLoss_scaled.rds")
fits10 %<-% run_my_models(bf_envMean_acLen_fixed, unique(envMean_data_interaction_acLen), envMean_acLen_fixed_prior, "outputs/stats/envMean_fixed_interaction_acLen_scaled.rds")
fits11 %<-% run_my_models(bf_acLen_envBreadthLoss_fixed, unique(acLen_data_interaction_envBreadthLoss), acLen_envBreadthLoss_fixed_prior, "outputs/stats/acLen_fixed_interaction_envBreadthLoss_scaled.rds")
fits12 %<-% run_my_models(bf_acLen_envMean_fixed, unique(acLen_data_interaction_envMean), acLen_envMean_fixed_prior, "outputs/stats/aclen_fixed_interaction_envMean_scaled.rds")

##### Read in Data + Models, Local Sensitivity Analysis #####
# Get properties of post habitat loss landscapes
# Read sample design
seeds <- readRDS("outputs/spatial/landscapes/HLSeedsMainTrend.rds")

# Get all seed combos
seeds_df <- data.frame(
  envSeed = as.integer(rep(names(seeds), sapply(seeds, length))),
  hlSeed = unlist(seeds))

# Get properties of post habitat loss landscapes
file_list <- list.files(path = "outputs/spatial/landscapes", pattern = "\\_samples.csv$", full.names = TRUE)
data_list <- lapply(file_list, read.csv)
hl_properties <- do.call(rbind, data_list)

### No need to re-define models or priors, same as in main

SA_prep_function <- function(SA_run){
  file_list <- list.files(path = "outputs/spatial/simulation_output", pattern = paste0(SA_run, "sensitivity.txt$"), full.names = TRUE)
  data_list <- lapply(file_list, read.csv)
  all_out <- do.call(rbind, data_list)

  # Group outputs by organism type, then get summary statistics for each PA design by looking at results from each seed
  grouped_data <- all_out %>%
    dplyr::group_by(ENV_SEED, HL_SEED)

  # Convert to binary outcome (persist or not)
  result <- grouped_data %>%
    mutate(persist = ifelse(proportion_original_pop_size==0, 0, 1))


  # Combine to get properties for all seeds
  all_normal <- merge(seeds_df, hl_properties, by = c("envSeed", "hlSeed"), all.x = TRUE)

  # Combine results with properties
  results_and_properties <- merge(result, all_normal, by.x = c("ENV_SEED", "HL_SEED"), by.y = c("envSeed", "hlSeed"), all.x = TRUE)

  # Which seeds for which properties? then extract that data
  # ENVIRONMENTAL Breadth LOSS
  envBreadthLoss_data <- get_specific_output(seeds, results_and_properties, 1, 50, 1) %>%
    mutate(scaled_envBreadthLoss = (envBreadthLoss-envBreadthLoss_mean)/envBreadthLoss_sd) #scale

  # ENVIRONMENTAL MEAN
  envMean_data <-  get_specific_output(seeds, results_and_properties, 51, 100, 1) %>%
    mutate(scaled_envMean = (envMean-envMean_mean)/envMean_sd) #scale

  # FRAGMENTATION
  acLen_data <-  get_specific_output(seeds, results_and_properties, 101, 150, 1) %>%
    mutate(scaled_acLen = (acLen-acLen_mean)/acLen_sd) #scale

  # Run!
  # random effect of env_seed
  fits1 %<-% run_my_models(bf_envBreadthLoss_random, envBreadthLoss_data, envBreadthLoss_random_prior, paste0("outputs/stats/envBreadthLoss_random_scaled_sensitivity", SA_run,".rds"))
  fits2 %<-% run_my_models(bf_envMean_random, envMean_data, envMean_random_prior, paste0("outputs/stats/envMean_random_scaled", SA_run,".rds"))
  fits3 %<-% run_my_models(bf_acLen_random, acLen_data, acLen_random_prior, paste0("outputs/stats/acLen_random_scaled", SA_run,".rds"))
}

# Run seperate for all local SA combinations
mclapply(1:18, SA_prep_function)

