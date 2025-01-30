install.packages(c("dplyr", "brms", "modelr", "future", "ggplot2", "tidybayes", "tidyr", "forcats", "ggdist", "patchwork", "cowplot", "raster"), dependencies = T)
library(ggplot2)
library(dplyr) 
library(brms) 
library(tidybayes) 
library(tidyr)
library(forcats)
library(modelr)
library(future) 
library(ggdist)
library(patchwork)
library(cowplot)
library(raster)

# Set working directory to main folder

####################################################
################### PLOTTING #######################
####################################################
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

##### Read in main Outputs ######
# Get properties of post habitat loss landscapes
### Normal (HL+CC)
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
# ENVIRONMENTAL BREATH LOSS
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


### HL only (no CC)
file_list_noCC <- list.files(path = "outputs/spatial/simulation_output", pattern = "\\_summaryCC0.txt$", full.names = TRUE)
data_list_noCC <- lapply(file_list_noCC, read.csv)
all_out_noCC <- do.call(rbind, data_list_noCC)

# Group outputs by organism type, then get summary statistics for each PA design by looking at results from each seed
grouped_data_noCC <- all_out_noCC %>%
  dplyr::group_by(ENV_SEED, HL_SEED)

# Convert to binary outcome (persist or not)
result_noCC <- grouped_data_noCC %>%
  mutate(persist = ifelse(proportion_original_pop_size==0, 0, 1))

# Combine results with properties
results_and_properties_noCC <- merge(result_noCC, all_normal, by.x = c("ENV_SEED", "HL_SEED"), by.y = c("envSeed", "hlSeed"), all.x = TRUE)

# Which seeds for which properties? then extract that data
# ENVIRONMENTAL BREATH LOSS
envBreadthLoss_data_noCC <- get_specific_output(seeds, results_and_properties_noCC, 1, 50, 1)

# ENVIRONMENTAL MEAN
envMean_data_noCC <-  get_specific_output(seeds, results_and_properties_noCC, 51, 100, 1)

# FRAGMENTATION
acLen_data_noCC <-  get_specific_output(seeds, results_and_properties_noCC, 101, 150, 1)


##### Set up Figure Plotting, read stats models ######
theme_set(theme_tidybayes() + panel_border())

# Random 
# Read in random effect model outputs
envBreadthLoss_fixed_output <- readRDS("outputs/stats/envBreadthLoss_fixed_scaled.rds")
envMean_fixed_output <- readRDS("outputs/stats/envMean_fixed_scaled.rds")
acLen_fixed_output <- readRDS("outputs/stats/acLen_fixed_scaled.rds")

# Read in random effect model outputs
envBreadthLoss_random_output <- readRDS("outputs/stats/envBreadthLoss_random_scaled.rds")
envMean_random_output <- readRDS("outputs/stats/envMean_random_scaled.rds")
acLen_random_output <- readRDS("outputs/stats/acLen_random_scaled.rds")

# Summary of the models
summary(envMean_random_output)

# What variables?
get_variables(envBreadthLoss_random_output)

##### Figure 2 - Landscape Variation #########
# read in landscape generation function
source("script/R/addt_functions/landscape_ac.R")

# set theme for plotting
theme_set(theme_tidybayes() + panel_border())

# outputs, with climate change
file_list <- list.files(path = "outputs/spatial/simulation_output", pattern = "\\_summaryCC3.txt$", full.names = TRUE)
data_list <- lapply(file_list, read.csv)
all_out <- do.call(rbind, data_list)
persist_prob <- all_out %>%
  mutate(persist = ifelse(proportion_original_pop_size==0, 0, 1)) %>%
  group_by(ENV_SEED) %>%
  summarise(prob_persist = mean(persist))

# outputs, without climate change
file_list <- list.files(path = "outputs/spatial/simulation_output", pattern = "\\_summaryCC0.txt$", full.names = TRUE)
data_list <- lapply(file_list, read.csv)
all_out <- do.call(rbind, data_list)
pop_diversity <- all_out %>%
  group_by(ENV_SEED) %>%
  summarise(phen_diversity = mean(phenotypic_variance), 
            prop_muts = mean(propotion_original_muts), 
            movement = mean(avg_movement),
            pop_size = mean(proportion_original_pop_size))

#merge together and plot
merged <- merge(persist_prob, pop_diversity, by="ENV_SEED")

landscapeVariation_left <- merged %>%
  ggplot(aes(x=phen_diversity, y=prob_persist)) +
  geom_point(size=4) +
  geom_smooth() + 
  labs(x=expression(Var.~Phenotypes~(sigma[P]^2)), y = expression(Probability~of~Persistence~(p[persist])))  +
  scale_y_continuous(limits=c(0,1)) +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))

landscapeVariation_right <- merged %>%
  ggplot(aes(x=movement, y=prob_persist)) +
  geom_point(size=4) +
  geom_smooth(method=) +
  labs(x=expression(Movement~(mu[d])), y = NULL)  +
  scale_y_continuous(limits=c(0,1),labels=NULL, breaks=NULL) +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))

landscapeVariation_bottom <- merged %>%
  ggplot(aes(x=movement, y=phen_diversity)) +
  geom_point(size=4) +
  geom_smooth(method="lm") +
  labs(x=expression(Movement~(mu[d])), y = expression(Var.~Phenotypes~(sigma[P]^2)))  +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))

landscape_combined <- ((landscapeVariation_left + landscapeVariation_right) / landscapeVariation_bottom) + plot_annotation(tag_levels = 'a')

# save plot as pdf
ggsave(landscape_combined, filename = "images/landscapeVariation.pdf", height = 6, width = 10)


## Supplemental: movement, landscape autocorrelation with Geary's C
# Now let's look at the autocorrelation structure of the landscapes
calc_GearyC <- function(ENV_SEED){
  # generate environmental landscape
  env_landscape <- scale(generateLandscape(ENV_SEED, 0.0, 0.0, 1.0, 0.1, minSize=256, periodic_x=T))
  
  # rasterize
  landscape_rast <- raster(env_landscape)
  
  # calculate Geary's C
  gearyC <- Geary(landscape_rast)
  
  return(gearyC)
}


# Convert env matrix to raster
# first, compute Geary's c for each landscape
merged$gearyC <- sapply(merged[, "ENV_SEED"], calc_GearyC)

gearyC_left <- merged %>%
  ggplot(aes(x=gearyC, y=movement)) +
  geom_point(size=4) +
  geom_smooth(method="lm") + 
  labs(x="Geary's C", y = expression(Movement~(mu[d])))  +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))

gearyC_right <- merged %>%
  ggplot(aes(x=gearyC, y=phen_diversity)) +
  geom_point(size=4) +
  geom_smooth(method="lm") +
  labs(x="Geary's C", y = expression(Var.~Phenotypes~(sigma[P]^2)))  +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))

gearyC_combined <- (gearyC_left + gearyC_right) + plot_annotation(tag_levels = 'a')

# save plot as pdf
ggsave(gearyC_combined, filename = "images/gearyC.pdf", height = 4, width = 10)


##### Figure 3 - comparing effects and spaghetti plots #####
# slopes
coef_envBreadthLoss <- coef(envBreadthLoss_random_output, summary=F)
# 1 = iteration, #2 env_seed #3 slopes

# posterior distributiona
# get posterior random slopes
envBreadthLoss_slope <- coef_envBreadthLoss$ENV_SEED[,,2] %>% 
  as.data.frame() %>% 
  pivot_longer(everything(),names_to = "ENV_SEED", values_to = "posterior_samples") %>%
  mutate(ENV_SEED = as.factor(ENV_SEED)) %>%
  group_by(ENV_SEED) %>%
  mutate(mean_slope = mean(posterior_samples)) %>%
  ungroup() 

# reorder factor levels for plotting
envBreadthLoss_slope$ENV_SEED <- fct_reorder(envBreadthLoss_slope$ENV_SEED,envBreadthLoss_slope$mean_slope)

envBreadthLoss_slope <- envBreadthLoss_slope %>%
  ggplot(aes(y = ENV_SEED, x = posterior_samples)) +
  stat_interval(.width = c(.50, .95, .99), linewidth=1.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = mean(envBreadthLoss_slope$posterior_samples)) +
  scale_x_continuous(lim=c(-3.25, 3.25), labels=c('', -2, '', 0, '', 2, ''), breaks=c(-3, -2, -1, 0, 1, 2, 3))  + 
  scale_y_discrete(labels=NULL, breaks=NULL) +
  labs(x=expression(Posterior~Slope~of~Delta*B[e]), y = expression(Landscape)) +
  theme(text = element_text(size = 12)) +
  labs(color = "CI") +
  scale_color_manual(values = c("#C5EEF4", "#9EC6DD", "#558097")) +
  theme(legend.position = "top")
  

coef_envMean <- coef(envMean_random_output, summary=F)

# posterior distributiona
# get posterior random slopes
envMean_slope <- coef_envMean$ENV_SEED[,,2] %>% 
  as.data.frame() %>% 
  pivot_longer(everything(),names_to = "ENV_SEED", values_to = "posterior_samples") %>%
  mutate(ENV_SEED = as.factor(ENV_SEED)) %>%
  group_by(ENV_SEED) %>%
  mutate(mean_slope = mean(posterior_samples)) %>%
  ungroup() 

# reorder factor levels for plotting
envMean_slope$ENV_SEED <- fct_reorder(envMean_slope$ENV_SEED,envMean_slope$mean_slope)

envMean_slope <- envMean_slope %>%
  ggplot(aes(y = ENV_SEED, x = posterior_samples)) +
  stat_interval(.width = c(.50, .95, .99), linewidth=1.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = mean(envMean_slope$posterior_samples)) +
  scale_x_continuous(lim=c(-3.25, 3.25), labels=c('', -2, '', 0, '', 2, ''), breaks=c(-3, -2, -1, 0, 1, 2, 3))  + 
  scale_y_discrete(labels=NULL, breaks=NULL) +
  labs(x=expression(Posterior~Slope~of~Delta*mu[e]), y = NULL) +
  theme(text = element_text(size = 12)) +
  labs(color = "CI") +
  scale_color_manual(values = c("#F77178", "#D90611", "#6A0B12")) +
  theme(legend.position = "top")


coef_acLen <- coef(acLen_random_output, summary=F)
# get estimates
acLen_slope <- coef_acLen$ENV_SEED[,,2] %>% 
  as.data.frame() %>% 
  pivot_longer(everything(),names_to = "ENV_SEED", values_to = "posterior_samples") %>%
  mutate(ENV_SEED = as.factor(ENV_SEED)) %>%
  group_by(ENV_SEED) %>%
  mutate(mean_slope = mean(posterior_samples)) %>%
  ungroup() 

# reorder factor levels for plotting
acLen_slope$ENV_SEED <- fct_reorder(acLen_slope$ENV_SEED,acLen_slope$mean_slope)

acLen_slope <- acLen_slope %>%
  ggplot(aes(y = ENV_SEED, x = posterior_samples)) +
  stat_interval(.width = c(.50, .95, .99), linewidth=1.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = mean(acLen_slope$posterior_samples)) +
  scale_x_continuous(lim=c(-3.25, 3.25), labels=c('', -2, '', 0, '', 2, ''), breaks=c(-3, -2, -1, 0, 1, 2, 3))  + 
  scale_y_discrete(labels=NULL, breaks=NULL) +
  labs(x=expression(Posterior~Slope~of~l[hl]), y = NULL) +
  theme(text = element_text(size = 12)) +
  labs(color = "CI") +
  scale_color_manual(values = c("#B3DEA4", "#6C9E64", "#122916")) +
  theme(legend.position = "top")


# spaghetti plots
# envBreadthLoss
mean_envBreadthLoss <- expand.grid(unique(envBreadthLoss_data$ENV_SEED), seq_range(envBreadthLoss_data$scaled_envBreadthLoss, n = 100)) %>%
  rename(ENV_SEED = Var1, scaled_envBreadthLoss = Var2) %>%
  add_epred_draws(envBreadthLoss_random_output, ndraws = 100) %>%
  dplyr::select(scaled_envBreadthLoss, .epred) %>%
  group_by(scaled_envBreadthLoss) %>%
  summarise(p_persist = mean(.epred))
grouped_envBreadthLoss <- envBreadthLoss_data %>%
  group_by(ENV_SEED) %>%
  data_grid(scaled_envBreadthLoss = seq_range(scaled_envBreadthLoss, n = 101)) %>%
  add_epred_draws(envBreadthLoss_random_output, ndraws = 100) %>%
  summarise(p_persist = mean(.epred)) %>%
  ggplot(aes(x = (scaled_envBreadthLoss*envBreadthLoss_sd + envBreadthLoss_mean), y = p_persist)) + 
  # environmental landscapes (groups)
  geom_line(aes(group = ENV_SEED), color='grey40') +
  geom_line(data = mean_envBreadthLoss, color='black', linewidth=1.5) +
  scale_y_continuous(limits=c(0,1)) +
  labs(x=expression(Environmental~Breadth~Loss~(Delta*B[e])), y = expression(Probability~of~Persistence~(p[persist])))  +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))

mean_envMean <- expand.grid(unique(envMean_data$ENV_SEED), seq_range(envMean_data$scaled_envMean, n = 100)) %>%
  rename(ENV_SEED = Var1, scaled_envMean = Var2) %>%
  add_epred_draws(envMean_random_output, ndraws = 100) %>%
  dplyr::select(scaled_envMean, .epred) %>%
  group_by(scaled_envMean) %>%
  summarise(p_persist = mean(.epred))
grouped_envMean <- envMean_data %>%
  group_by(ENV_SEED) %>%
  data_grid(scaled_envMean = seq_range(scaled_envMean, n = 101)) %>%
  add_epred_draws(envMean_random_output, ndraws = 100) %>%
  summarise(p_persist = mean(.epred)) %>%
  ggplot(aes(x = (scaled_envMean*envMean_sd + envMean_mean), y = p_persist)) + 
  # environmental landscapes (groups)
  geom_line(aes(group = ENV_SEED), color='grey40') +
  geom_line(data = mean_envMean, color='black', linewidth=1.5) +
  scale_y_continuous(limits=c(0,1), labels=NULL, breaks=NULL) +
  xlim(-0.8, 0.8) +
  labs(x=expression(Change~Environmental~Mean~(Delta*mu[e])), y = NULL)  +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))

mean_acLen <- expand.grid(unique(acLen_data$ENV_SEED), seq_range(acLen_data$scaled_acLen, n = 100)) %>%
  rename(ENV_SEED = Var1, scaled_acLen = Var2) %>%
  add_epred_draws(acLen_random_output, ndraws = 100) %>%
  dplyr::select(scaled_acLen, .epred) %>%
  group_by(scaled_acLen) %>%
  summarise(p_persist = mean(.epred)) 
grouped_acLen <- acLen_data %>%
  group_by(ENV_SEED) %>%
  data_grid(scaled_acLen = seq_range(scaled_acLen, n = 101)) %>% #set up grid for extracting draws
  add_epred_draws(acLen_random_output, ndraws = 100) %>% #get draws
  summarise(p_persist = mean(.epred)) %>% #convert binomial draws to probability
  ggplot(aes(x = (scaled_acLen*acLen_sd + acLen_mean), y = p_persist)) +
  # environmental landscapes (groups)
  geom_line(aes(group = ENV_SEED),color='grey40') +  
  geom_line(data = mean_acLen, color='black', linewidth=1.5) +
  scale_y_continuous(limits=c(0,1), labels=NULL, breaks=NULL) +
  labs(x=expression(Habitat~Loss~Autocorrelation~Length~(l[hl])), y = NULL)  +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))


main <- (envBreadthLoss_slope | envMean_slope | acLen_slope) / (grouped_envBreadthLoss | grouped_envMean | grouped_acLen)  + 
  plot_annotation(tag_levels = 'a')# + 
  #plot_layout(heights = c(2, 1))

# save plot as pdf
ggsave(main, filename = "images/main.pdf", height = 8, width = 10)



## Supplemental - what if spaghetti plots were just made using smoothers?
envBreathLoss_main_gam <- envBreadthLoss_data %>%
  dplyr::select(ENV_SEED, HL_SEED, SLIM_SEED, persist, envBreadthLoss) %>%
  group_by(ENV_SEED, HL_SEED) %>%
  summarise(p_persist = mean(persist), envBreadthLoss=mean(envBreadthLoss)) %>% # taking mean of envBreadthLoss OK since it is always the same for every ENV_SEED/HL_SEED combo
  ungroup() %>%
  ggplot(aes(x = envBreadthLoss, y = p_persist)) + 
  geom_smooth(aes(group=ENV_SEED), method="loess", span=1, se=F, color='grey40', linewidth=0.5) +
  geom_smooth(method="loess", span=1, se=F, color='black', linewidth=1.5) +
  scale_y_continuous(limits=c(0,1)) +
  labs(x=expression(Environmental~Breadth~Loss~(Delta*B[e])), y = expression(Probability~of~Persistence~(p[persist])))  +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))

envMean_main_gam <- envMean_data %>%
  dplyr::select(ENV_SEED, HL_SEED, SLIM_SEED, persist, envMean) %>%
  group_by(ENV_SEED, HL_SEED) %>%
  summarise(p_persist = mean(persist), envMean=mean(envMean)) %>% # taking mean of envMean OK since it is always the same for every ENV_SEED/HL_SEED combo
  ungroup() %>%
  ggplot(aes(x = envMean, y = p_persist)) + 
  geom_smooth(aes(group=ENV_SEED), method="loess", span=1, se=F, color='grey40', linewidth=0.5) +
  geom_smooth(method="loess", span=1, se=F, color='black', linewidth=1.5) +
  scale_y_continuous(limits=c(0,1), labels=NULL, breaks=NULL) +
  labs(x=expression(Change~Environmental~Mean~(Delta*mu[e])), y = NULL)  +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))

acLen_main_gam <- acLen_data %>%
  dplyr::select(ENV_SEED, HL_SEED, SLIM_SEED, persist, acLen) %>%
  group_by(ENV_SEED, HL_SEED) %>%
  summarise(p_persist = mean(persist), acLen=mean(acLen)) %>% # taking mean of acLen OK since it is always the same for every ENV_SEED/HL_SEED combo
  ungroup() %>%
  ggplot(aes(x = acLen, y = p_persist)) + 
  geom_smooth(aes(group=ENV_SEED), method="loess", span=1, se=F, color='grey40', linewidth=0.5) +
  geom_smooth(method="loess", span=1, se=F, color='black', linewidth=1.5) +
  scale_y_continuous(limits=c(0,1), labels=NULL, breaks=NULL) +
  labs(x=expression(Habitat~Loss~Autocorrelation~Length~(l[hl])), y = NULL)  +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))
  
main <- (envBreathLoss_main_gam + envMean_main_gam + acLen_main_gam)  + 
  plot_annotation(tag_levels = 'a')# + 
#plot_layout(heights = c(2, 1))

# save plot as pdf
ggsave(main, filename = "images/main_LOESSversion.pdf", height = 4, width = 10)


##### Supplemental: alternate output metrics #####
### First, change in population size for populations that persisted
envBreathLoss_populationsize_gam <- envBreadthLoss_data %>%
  # filter only for populations that persisted
  filter(persist == 1) %>%
  dplyr::select(ENV_SEED, HL_SEED, SLIM_SEED, proportion_original_pop_size, envBreadthLoss) %>%
  group_by(ENV_SEED, HL_SEED) %>%
  summarise(p_persist = mean(proportion_original_pop_size), envBreadthLoss=mean(envBreadthLoss)) %>% # taking mean of envBreadthLoss OK since it is always the same for every ENV_SEED/HL_SEED combo
  ungroup() %>%
  ggplot(aes(x = envBreadthLoss, y = p_persist)) + 
  geom_smooth(aes(group=ENV_SEED), method="loess", span=1, se=F, color='grey40', linewidth=0.5) +
  geom_smooth(method="loess", span=1, se=F, color='black', linewidth=1.5) +
  scale_y_continuous(limits=c(0,0.5)) +
  labs(x=expression(Environmental~Breadth~Loss~(Delta*B[e])), y = expression(Prop.~Pop~Size~(Delta*n[i])))  +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))

envMean_populationsize_gam <- envMean_data %>%
  # filter only for populations that persisted
  filter(persist == 1) %>%
  dplyr::select(ENV_SEED, HL_SEED, SLIM_SEED, proportion_original_pop_size, envMean) %>%
  group_by(ENV_SEED, HL_SEED) %>%
  summarise(p_persist = mean(proportion_original_pop_size), envMean=mean(envMean)) %>% # taking mean of envMean OK since it is always the same for every ENV_SEED/HL_SEED combo
  ungroup() %>%
  ggplot(aes(x = envMean, y = p_persist)) + 
  geom_smooth(aes(group=ENV_SEED), method="loess", span=1, se=F, color='grey40', linewidth=0.5) +
  geom_smooth(method="loess", span=1, se=F, color='black', linewidth=1.5) +
  scale_y_continuous(limits=c(0,0.5), labels=NULL, breaks=NULL) +
  labs(x=expression(Change~Environmental~Mean~(Delta*mu[e])), y = NULL)  +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))

acLen_populationsize_gam <- acLen_data %>%
  # filter only for populations that persisted
  filter(persist == 1) %>%
  dplyr::select(ENV_SEED, HL_SEED, SLIM_SEED, proportion_original_pop_size, acLen) %>%
  group_by(ENV_SEED, HL_SEED) %>%
  summarise(p_persist = mean(proportion_original_pop_size), acLen=mean(acLen)) %>% # taking mean of acLen OK since it is always the same for every ENV_SEED/HL_SEED combo
  ungroup() %>%
  ggplot(aes(x = acLen, y = p_persist)) + 
  geom_smooth(aes(group=ENV_SEED), method="loess", span=1, se=F, color='grey40', linewidth=0.5) +
  geom_smooth(method="loess", span=1, se=F, color='black', linewidth=1.5) +
  scale_y_continuous(limits=c(0,0.5), labels=NULL, breaks=NULL) +
  labs(x=expression(Habitat~Loss~Autocorrelation~Length~(l[hl])), y = NULL)  +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))

main_populationsize <- (envBreathLoss_populationsize_gam + envMean_populationsize_gam + acLen_populationsize_gam)  + 
  plot_annotation(tag_levels = 'a')# + 
#plot_layout(heights = c(2, 1))

# save plot as pdf
ggsave(main_populationsize, filename = "images/proportionPopulation_LOESSversion.pdf", height = 4, width = 10)


### Second, phenotypic variance
envBreathLoss_propmuts_gam <- envBreadthLoss_data %>%
  # filter only for populations that persisted
  filter(persist == 1) %>%
  dplyr::select(ENV_SEED, HL_SEED, SLIM_SEED, propotion_original_muts, envBreadthLoss) %>%
  group_by(ENV_SEED, HL_SEED) %>%
  summarise(p_persist = mean(propotion_original_muts), envBreadthLoss=mean(envBreadthLoss)) %>% # taking mean of envBreadthLoss OK since it is always the same for every ENV_SEED/HL_SEED combo
  ungroup() %>%
  ggplot(aes(x = envBreadthLoss, y = p_persist)) + 
  geom_smooth(aes(group=ENV_SEED), method="loess", span=1, se=F, color='grey40', linewidth=0.5) +
  geom_smooth(method="loess", span=1, se=F, color='black', linewidth=1.5) +
  scale_y_continuous(limits=c(0.1,0.4)) +
  labs(x=expression(Environmental~Breadth~Loss~(Delta*B[e])), y = expression(Prop.~SGV~(Delta*n[a])))  +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))

envMean_propmuts_gam <- envMean_data %>%
  # filter only for populations that persisted
  filter(persist == 1) %>%
  dplyr::select(ENV_SEED, HL_SEED, SLIM_SEED, propotion_original_muts, envMean) %>%
  group_by(ENV_SEED, HL_SEED) %>%
  summarise(p_persist = mean(propotion_original_muts), envMean=mean(envMean)) %>% # taking mean of envMean OK since it is always the same for every ENV_SEED/HL_SEED combo
  ungroup() %>%
  ggplot(aes(x = envMean, y = p_persist)) + 
  geom_smooth(aes(group=ENV_SEED), method="loess", span=1, se=F, color='grey40', linewidth=0.5) +
  geom_smooth(method="loess", span=1, se=F, color='black', linewidth=1.5) +
  scale_y_continuous(limits=c(0.1,0.4), labels=NULL, breaks=NULL) +
  labs(x=expression(Change~Environmental~Mean~(Delta*mu[e])), y = NULL)  +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))

acLen_propmuts_gam <- acLen_data %>%
  # filter only for populations that persisted
  filter(persist == 1) %>%
  dplyr::select(ENV_SEED, HL_SEED, SLIM_SEED, propotion_original_muts, acLen) %>%
  group_by(ENV_SEED, HL_SEED) %>%
  summarise(p_persist = mean(propotion_original_muts), acLen=mean(acLen)) %>% # taking mean of acLen OK since it is always the same for every ENV_SEED/HL_SEED combo
  ungroup() %>%
  ggplot(aes(x = acLen, y = p_persist)) + 
  geom_smooth(aes(group=ENV_SEED), method="loess", span=1, se=F, color='grey40', linewidth=0.5) +
  geom_smooth(method="loess", span=1, se=F, color='black', linewidth=1.5) +
  scale_y_continuous(limits=c(0.1,0.4), labels=NULL, breaks=NULL) +
  labs(x=expression(Habitat~Loss~Autocorrelation~Length~(l[hl])), y = NULL)  +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))

main_proportionmuts <- (envBreathLoss_propmuts_gam + envMean_propmuts_gam + acLen_propmuts_gam)  + 
  plot_annotation(tag_levels = 'a')# + 
#plot_layout(heights = c(2, 1))

# save plot as pdf
ggsave(main_proportionmuts, filename = "images/proportionNonNeutralVariation_LOESSversion.pdf", height = 4, width = 10)


### Lastly, time to extirpation
envBreathLoss_timeextinct_gam <- envBreadthLoss_data %>%
  # filter only for populations that persisted
  filter(persist == 0) %>%
  # Calculate time to extirpation
  mutate(time_to_extirpation = end_tick - 10000) %>%
  dplyr::select(ENV_SEED, HL_SEED, SLIM_SEED, time_to_extirpation, envBreadthLoss) %>%
  group_by(ENV_SEED, HL_SEED) %>%
  summarise(p_persist = mean(time_to_extirpation), envBreadthLoss=mean(envBreadthLoss)) %>% # taking mean of envBreadthLoss OK since it is always the same for every ENV_SEED/HL_SEED combo
  ungroup() %>%
  ggplot(aes(x = envBreadthLoss, y = p_persist)) + 
  geom_smooth(aes(group=ENV_SEED), method="loess", span=1, se=F, color='grey40', linewidth=0.5) +
  geom_smooth(method="loess", span=1, se=F, color='black', linewidth=1.5) +
  scale_y_continuous(limits=c(75,200)) +
  labs(x=expression(Environmental~Breadth~Loss~(Delta*B[e])), y = expression(Time~to~Extinction~(t[extinct])))  +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))

envMean_timeextinct_gam <- envMean_data %>%
  # filter only for populations that persisted
  filter(persist == 0) %>%
  # Calculate time to extirpation
  mutate(time_to_extirpation = end_tick - 10000) %>%
  dplyr::select(ENV_SEED, HL_SEED, SLIM_SEED, time_to_extirpation, envMean) %>%
  group_by(ENV_SEED, HL_SEED) %>%
  summarise(p_persist = mean(time_to_extirpation), envMean=mean(envMean)) %>% # taking mean of envMean OK since it is always the same for every ENV_SEED/HL_SEED combo
  ungroup() %>%
  ggplot(aes(x = envMean, y = p_persist)) + 
  geom_smooth(aes(group=ENV_SEED), method="loess", span=1, se=F, color='grey40', linewidth=0.5) +
  geom_smooth(method="loess", span=1, se=F, color='black', linewidth=1.5) +
  scale_y_continuous(limits=c(75,200), labels=NULL, breaks=NULL) +
  labs(x=expression(Change~Environmental~Mean~(Delta*mu[e])), y = NULL)  +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))

acLen_timeextinct_gam <- acLen_data %>%
  # filter only for populations that persisted
  filter(persist == 0) %>%
  # Calculate time to extirpation
  mutate(time_to_extirpation = end_tick - 10000) %>%
  dplyr::select(ENV_SEED, HL_SEED, SLIM_SEED, time_to_extirpation, acLen) %>%
  group_by(ENV_SEED, HL_SEED) %>%
  summarise(p_persist = mean(time_to_extirpation), acLen=mean(acLen)) %>% # taking mean of acLen OK since it is always the same for every ENV_SEED/HL_SEED combo
  ungroup() %>%
  ggplot(aes(x = acLen, y = p_persist)) + 
  geom_smooth(aes(group=ENV_SEED), method="loess", span=1, se=F, color='grey40', linewidth=0.5) +
  geom_smooth(method="loess", span=1, se=F, color='black', linewidth=1.5) +
  scale_y_continuous(limits=c(75,200), labels=NULL, breaks=NULL) +
  labs(x=expression(Habitat~Loss~Autocorrelation~Length~(l[hl])), y = NULL)  +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))

main_timeextinct <- (envBreathLoss_timeextinct_gam + envMean_timeextinct_gam + acLen_timeextinct_gam)  + 
  plot_annotation(tag_levels = 'a')# + 
#plot_layout(heights = c(2, 1))

# save plot as pdf
ggsave(main_timeextinct, filename = "images/timeExtinct_LOESSversion.pdf", height = 4, width = 10)


##### Figure 4 - mechanisms #####
# Main figure 3
a1 <- ggplot() + 
  geom_smooth(data=envBreadthLoss_data_noCC[envBreadthLoss_data_noCC$end_tick==10200,], aes(x=envBreadthLoss, y=avg_movement, linetype='noCC'), color='black', linewidth=1.5, se=F) +
  geom_smooth(data=envBreadthLoss_data[envBreadthLoss_data$end_tick==10200,], aes(x=envBreadthLoss, y=avg_movement, linetype='CC'), color='black', linewidth=1.5, se=F) +
  labs(x=expression(Environmental~Breadth~Loss~(Delta*B[e])), y=expression(Movement~(mu[d])))  +
  scale_x_continuous(limits=c(0,2.25)) +
  scale_y_continuous(limits=c(0.05,0.075)) +
  scale_linetype_manual(name = "Model Scenario", values = c("CC" = "solid", "noCC" = "twodash"), labels = c("CC" = "With Environmental Change", "noCC" = "Without Environmental Change"))+
  theme_set(theme_tidybayes() + panel_border()) +
  # theme(legend.position = "none") +
  theme(text = element_text(size = 12))
b1 <- ggplot() + 
  geom_smooth(data=envMean_data_noCC[envMean_data_noCC$end_tick==10200,], aes(x=envMean, y=avg_movement, linetype='noCC'), color='black', linewidth=1.5, se=F) +
  geom_smooth(data=envMean_data[envMean_data$end_tick==10200,], aes(x=envMean, y=avg_movement, linetype='CC'), color='black', linewidth=1.5, se=F) +
  labs(x=expression(Change~Environmental~Mean~(Delta*mu[e])), y=NULL)  +
  scale_x_continuous(limits=c(-0.8,0.8)) +
  scale_y_continuous(limits=c(0.05,0.075),labels=NULL, breaks=NULL) +
  scale_linetype_manual(name = "Model Scenario", values = c("CC" = "solid", "noCC" = "twodash"), labels = c("CC" = "With Environmental Change", "noCC" = "Without Environmental Change"))+
  theme_set(theme_tidybayes() + panel_border()) +
  # theme(legend.position = "none") +
  theme(text = element_text(size = 12))
c1 <- ggplot() + 
  geom_smooth(data=acLen_data_noCC[acLen_data_noCC$end_tick==10200,], aes(x=acLen, y=avg_movement, linetype='noCC'), color='black', linewidth=1.5, se=F) +
  geom_smooth(data=acLen_data[acLen_data$end_tick==10200,], aes(x=acLen, y=avg_movement, linetype='CC'), color='black', linewidth=1.5, se=F) +
  labs(x=expression(Habitat~Loss~Autocorrelation~Length~(l[hl])), y=NULL)  +
  scale_x_continuous(limits=c(0.02,0.18)) +
  scale_y_continuous(limits=c(0.05,0.075),labels=NULL, breaks=NULL) +
  scale_linetype_manual(name = "Model Scenario", values = c("CC" = "solid", "noCC" = "twodash"), labels = c("CC" = "With Environmental Change", "noCC" = "Without Environmental Change"))+
  theme_set(theme_tidybayes() + panel_border()) +
  # theme(legend.position = "none") +
  theme(text = element_text(size = 12))

mechanisms_movement <- (a1|b1|c1) + plot_annotation(tag_levels = 'a')  + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

# save plot as pdf
ggsave(mechanisms_movement, filename = "images/Mechanisms_movement_only.pdf", height = 4, width = 10)


# Supplement plot, all
a1 <- ggplot() + 
  geom_smooth(data=envBreadthLoss_data_noCC[envBreadthLoss_data_noCC$end_tick==10200,], aes(x=envBreadthLoss, y=avg_movement, linetype='noCC'), color='black', linewidth=1.5, se=F) +
  geom_smooth(data=envBreadthLoss_data[envBreadthLoss_data$end_tick==10200,], aes(x=envBreadthLoss, y=avg_movement, linetype='CC'), color='black', linewidth=1.5, se=F) +
  labs(x=NULL, y=expression(Movement~(mu[d])))  +
  scale_x_continuous(limits=c(0,2.25), labels=NULL, breaks=NULL) +
  scale_y_continuous(limits=c(0.05,0.075)) +
  scale_linetype_manual(name = "Model Scenario", values = c("CC" = "solid", "noCC" = "twodash"), labels = c("CC" = "With Environmental Change", "noCC" = "Without Environmental Change"))+
  theme_set(theme_tidybayes() + panel_border()) +
  # theme(legend.position = "none") +
  theme(text = element_text(size = 12))
a2 <- ggplot() + 
  geom_smooth(data=envBreadthLoss_data_noCC[envBreadthLoss_data_noCC$end_tick==10200,], aes(x=envBreadthLoss, y=propotion_original_muts, linetype='noCC'), color='black', linewidth=1.5, se=F) +
  geom_smooth(data=envBreadthLoss_data[envBreadthLoss_data$end_tick==10200,], aes(x=envBreadthLoss, y=propotion_original_muts, linetype='CC'), color='black', linewidth=1.5, se=F) +
  labs(x=NULL, y=expression(Prop.~SGV~(Delta*n[a])))  +
  scale_x_continuous(limits=c(0,2.25), labels=NULL, breaks=NULL) +
  scale_y_continuous(limits=c(0,0.7)) +
  scale_linetype_manual(name = "Model Scenario", values = c("CC" = "solid", "noCC" = "twodash"), labels = c("CC" = "With Environmental Change", "noCC" = "Without Environmental Change"))+
  theme_set(theme_tidybayes() + panel_border()) +
  # theme(legend.position = "none") +
  theme(text = element_text(size = 12))
a3 <- ggplot() + 
  geom_smooth(data=envBreadthLoss_data_noCC[envBreadthLoss_data_noCC$end_tick==10200,], aes(x=envBreadthLoss, y=phenotypic_variance, linetype='noCC'), color='black', linewidth=1.5, se=F) +
  geom_smooth(data=envBreadthLoss_data[envBreadthLoss_data$end_tick==10200,], aes(x=envBreadthLoss, y=phenotypic_variance, linetype='CC'), color='black', linewidth=1.5, se=F) +
  labs(x=NULL, y=expression(Var.~Phenotypes~(sigma[P]^2)))  +
  scale_x_continuous(limits=c(0,2.25), labels=NULL, breaks=NULL) +
  scale_y_continuous(limits=c(0,1)) +
  scale_linetype_manual(name = "Model Scenario", values = c("CC" = "solid", "noCC" = "twodash"), labels = c("CC" = "With Environmental Change", "noCC" = "Without Environmental Change"))+
  theme_set(theme_tidybayes() + panel_border()) +
  # theme(legend.position = "none") +
  theme(text = element_text(size = 12))
a4 <- ggplot() + 
  geom_smooth(data=envBreadthLoss_data_noCC[envBreadthLoss_data_noCC$end_tick==10200,], aes(x=envBreadthLoss, y=proportion_original_pop_size, linetype='noCC'), color='black', linewidth=1.5, se=F) +
  geom_smooth(data=envBreadthLoss_data[envBreadthLoss_data$end_tick==10200,], aes(x=envBreadthLoss, y=proportion_original_pop_size, linetype='CC'), color='black', linewidth=1.5, se=F) +
  labs(x=expression(Environmental~Breadth~Loss~(Delta*B[e])), y=expression(Prop.~Pop~Size~(Delta*n[i])))  +
  scale_x_continuous(limits=c(0,2.25)) +
  scale_y_continuous(limits=c(0,0.75)) +
  scale_linetype_manual(name = "Model Scenario", values = c("CC" = "solid", "noCC" = "twodash"), labels = c("CC" = "With Environmental Change", "noCC" = "Without Environmental Change"))+
  theme_set(theme_tidybayes() + panel_border()) +
  # theme(legend.position = "none") +
  theme(text = element_text(size = 12))
b1 <- ggplot() + 
  geom_smooth(data=envMean_data_noCC[envMean_data_noCC$end_tick==10200,], aes(x=envMean, y=avg_movement, linetype='noCC'), color='black', linewidth=1.5, se=F) +
  geom_smooth(data=envMean_data[envMean_data$end_tick==10200,], aes(x=envMean, y=avg_movement, linetype='CC'), color='black', linewidth=1.5, se=F) +
  labs(x=NULL, y=NULL)  +
  scale_x_continuous(limits=c(-0.8,0.8), labels=NULL, breaks=NULL) +
  scale_y_continuous(limits=c(0.05,0.075),labels=NULL, breaks=NULL) +
  scale_linetype_manual(name = "Model Scenario", values = c("CC" = "solid", "noCC" = "twodash"), labels = c("CC" = "With Environmental Change", "noCC" = "Without Environmental Change"))+
  theme_set(theme_tidybayes() + panel_border()) +
  # theme(legend.position = "none") +
  theme(text = element_text(size = 12))
b2 <- ggplot() + 
  geom_smooth(data=envMean_data_noCC[envMean_data_noCC$end_tick==10200,], aes(x=envMean, y=propotion_original_muts, linetype='noCC'), color='black', linewidth=1.5, se=F) +
  geom_smooth(data=envMean_data[envMean_data$end_tick==10200,], aes(x=envMean, y=propotion_original_muts, linetype='CC'), color='black', linewidth=1.5, se=F) +
  labs(x=NULL, y=NULL)  +
  scale_x_continuous(limits=c(-0.8,0.8), labels=NULL, breaks=NULL) +
  scale_y_continuous(limits=c(0,0.7),labels=NULL, breaks=NULL) +
  scale_linetype_manual(name = "Model Scenario", values = c("CC" = "solid", "noCC" = "twodash"), labels = c("CC" = "With Environmental Change", "noCC" = "Without Environmental Change"))+
  theme_set(theme_tidybayes() + panel_border()) +
  # theme(legend.position = "none") +
  theme(text = element_text(size = 12))
b3 <- ggplot() + 
  geom_smooth(data=envMean_data_noCC[envMean_data_noCC$end_tick==10200,], aes(x=envMean, y=phenotypic_variance, linetype='noCC'), color='black', linewidth=1.5, se=F) +
  geom_smooth(data=envMean_data[envMean_data$end_tick==10200,], aes(x=envMean, y=phenotypic_variance, linetype='CC'), color='black', linewidth=1.5, se=F) +
  labs(x=NULL, y=NULL)  +
  scale_x_continuous(limits=c(-0.8,0.8), labels=NULL, breaks=NULL) +
  scale_y_continuous(limits=c(0,1),labels=NULL, breaks=NULL) +
  scale_linetype_manual(name = "Model Scenario", values = c("CC" = "solid", "noCC" = "twodash"), labels = c("CC" = "With Environmental Change", "noCC" = "Without Environmental Change"))+
  theme_set(theme_tidybayes() + panel_border()) +
  # theme(legend.position = "none") +
  theme(text = element_text(size = 12))
b4 <- ggplot() + 
  geom_smooth(data=envMean_data_noCC[envMean_data_noCC$end_tick==10200,], aes(x=envMean, y=proportion_original_pop_size, linetype='noCC'), color='black', linewidth=1.5, se=F) +
  geom_smooth(data=envMean_data[envMean_data$end_tick==10200,], aes(x=envMean, y=proportion_original_pop_size, linetype='CC'), color='black', linewidth=1.5, se=F) +
  labs(x=expression(Change~Environmental~Mean~(Delta*mu[e])), y=NULL)  +
  scale_x_continuous(limits=c(-0.8,0.8)) +
  scale_y_continuous(limits=c(0,0.75),labels=NULL, breaks=NULL) +
  scale_linetype_manual(name = "Model Scenario", values = c("CC" = "solid", "noCC" = "twodash"), labels = c("CC" = "With Environmental Change", "noCC" = "Without Environmental Change"))+
  theme_set(theme_tidybayes() + panel_border()) +
  # theme(legend.position = "none") +
  theme(text = element_text(size = 12))
c1 <- ggplot() + 
  geom_smooth(data=acLen_data_noCC[acLen_data_noCC$end_tick==10200,], aes(x=acLen, y=avg_movement, linetype='noCC'), color='black', linewidth=1.5, se=F) +
  geom_smooth(data=acLen_data[acLen_data$end_tick==10200,], aes(x=acLen, y=avg_movement, linetype='CC'), color='black', linewidth=1.5, se=F) +
  labs(x=NULL, y=NULL)  +
  scale_x_continuous(limits=c(0.02,0.18),labels=NULL, breaks=NULL) +
  scale_y_continuous(limits=c(0.05,0.075),labels=NULL, breaks=NULL) +
  scale_linetype_manual(name = "Model Scenario", values = c("CC" = "solid", "noCC" = "twodash"), labels = c("CC" = "With Environmental Change", "noCC" = "Without Environmental Change"))+
  theme_set(theme_tidybayes() + panel_border()) +
  # theme(legend.position = "none") +
  theme(text = element_text(size = 12))
c2 <- ggplot() + 
  geom_smooth(data=acLen_data_noCC[acLen_data_noCC$end_tick==10200,], aes(x=acLen, y=propotion_original_muts, linetype='noCC'), color='black', linewidth=1.5, se=F) +
  geom_smooth(data=acLen_data[acLen_data$end_tick==10200,], aes(x=acLen, y=propotion_original_muts, linetype='CC'), color='black', linewidth=1.5, se=F) +
  labs(x=NULL, y=NULL)  +
  scale_x_continuous(limits=c(0.02,0.18),labels=NULL, breaks=NULL) +
  scale_y_continuous(limits=c(0,0.7),labels=NULL, breaks=NULL) +
  scale_linetype_manual(name = "Model Scenario", values = c("CC" = "solid", "noCC" = "twodash"), labels = c("CC" = "With Environmental Change", "noCC" = "Without Environmental Change"))+
  theme_set(theme_tidybayes() + panel_border()) +
  # theme(legend.position = "none") +
  theme(text = element_text(size = 12))
c3 <- ggplot() + 
  geom_smooth(data=acLen_data_noCC[acLen_data_noCC$end_tick==10200,], aes(x=acLen, y=phenotypic_variance, linetype='noCC'), color='black', linewidth=1.5, se=F) +
  geom_smooth(data=acLen_data[acLen_data$end_tick==10200,], aes(x=acLen, y=phenotypic_variance, linetype='CC'), color='black', linewidth=1.5, se=F) +
  labs(x=NULL, y=NULL)  +
  scale_x_continuous(limits=c(0.02,0.18),labels=NULL, breaks=NULL) +
  scale_y_continuous(limits=c(0,1),labels=NULL, breaks=NULL) +
  scale_linetype_manual(name = "Model Scenario", values = c("CC" = "solid", "noCC" = "twodash"), labels = c("CC" = "With Environmental Change", "noCC" = "Without Environmental Change"))+
  theme_set(theme_tidybayes() + panel_border()) +
  # theme(legend.position = "none") +
  theme(text = element_text(size = 12))
c4 <- ggplot() + 
  geom_smooth(data=acLen_data_noCC[acLen_data_noCC$end_tick==10200,], aes(x=acLen, y=proportion_original_pop_size, linetype='noCC'), color='black', linewidth=1.5, se=F) +
  geom_smooth(data=acLen_data[acLen_data$end_tick==10200,], aes(x=acLen, y=proportion_original_pop_size, linetype='CC'), color='black', linewidth=1.5, se=F) +
  labs(x=expression(Habitat~Loss~Autocorrelation~Length~(l[hl])), y=NULL)  +
  scale_x_continuous(limits=c(0.02,0.18)) +
  scale_y_continuous(limits=c(0,0.75),labels=NULL, breaks=NULL) +
  scale_linetype_manual(name = "Model Scenario", values = c("CC" = "solid", "noCC" = "twodash"), labels = c("CC" = "With Environmental Change", "noCC" = "Without Environmental Change"))+
  theme_set(theme_tidybayes() + panel_border()) +
  # theme(legend.position = "none") +
  theme(text = element_text(size = 12))

all_mechanisms <- ((a1|b1|c1) / (a2|b2|c2) / (a3|b3|c3) / (a4|b4|c4)) +
  plot_annotation(tag_levels = 'a')  + 
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom')

# save plot as pdf
ggsave(all_mechanisms, filename = "images/Mechanisms_All.pdf", height = 10, width = 10)


##### Figure 5 - interactions #####
# Read in data again
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

envBreadthLoss_data_interaction_envMean <- bind_rows(envBreadthLoss_data_interaction_low_envMean, envBreadthLoss_data_interaction_high_envMean)
envBreadthLoss_data_interaction_acLen <- bind_rows(envBreadthLoss_data_interaction_low_acLen, envBreadthLoss_data_interaction_high_acLen)


# ENVIRONMENTAL MEAN
envMean_data_interaction_low_envBreadthLoss <- get_specific_output(seeds_interactions, results_and_properties_interactions, 101, 150, 2) %>%
  mutate(scaled_envMean = (envMean-envMean_mean)/envMean_sd, scaled_envBreadthLoss = (envBreadthLoss-envBreadthLoss_mean)/envBreadthLoss_sd) #scale
envMean_data_interaction_high_envBreadthLoss <- get_specific_output(seeds_interactions, results_and_properties_interactions, 401, 450, 3) %>%
  mutate(scaled_envMean = (envMean-envMean_mean)/envMean_sd, scaled_envBreadthLoss = (envBreadthLoss-envBreadthLoss_mean)/envBreadthLoss_sd) #scale
envMean_data_interaction_low_acLen <- get_specific_output(seeds_interactions, results_and_properties_interactions, 151, 200, 2) %>%
  mutate(scaled_envMean = (envMean-envMean_mean)/envMean_sd, scaled_acLen = (acLen-acLen_mean)/acLen_sd) #scale
envMean_data_interaction_high_acLen <- get_specific_output(seeds_interactions, results_and_properties_interactions, 451, 500, 3) %>%
  mutate(scaled_envMean = (envMean-envMean_mean)/envMean_sd, scaled_acLen = (acLen-acLen_mean)/acLen_sd) #scale

envMean_data_interaction_envBreadthLoss <- bind_rows(envMean_data_interaction_low_envBreadthLoss, envMean_data_interaction_high_envBreadthLoss)
envMean_data_interaction_acLen <- bind_rows(envMean_data_interaction_low_acLen, envMean_data_interaction_high_acLen)


# AC LEN
acLen_data_interaction_low_envBreadthLoss <- get_specific_output(seeds_interactions, results_and_properties_interactions, 201, 250, 2) %>%
  mutate(scaled_acLen = (acLen-acLen_mean)/acLen_sd, scaled_envBreadthLoss = (envBreadthLoss-envBreadthLoss_mean)/envBreadthLoss_sd) #scale
acLen_data_interaction_high_envBreadthLoss <- get_specific_output(seeds_interactions, results_and_properties_interactions, 501, 550, 3) %>%
  mutate(scaled_acLen = (acLen-acLen_mean)/acLen_sd, scaled_envBreadthLoss = (envBreadthLoss-envBreadthLoss_mean)/envBreadthLoss_sd) #scale
acLen_data_interaction_low_envMean <- get_specific_output(seeds_interactions, results_and_properties_interactions, 251, 300, 2) %>%
  mutate(scaled_acLen = (acLen-acLen_mean)/acLen_sd, scaled_envMean = (envMean-envMean_mean)/envMean_sd) #scale
acLen_data_interaction_high_envMean <- get_specific_output(seeds_interactions, results_and_properties_interactions, 551, 600, 3) %>%
  mutate(scaled_acLen = (acLen-acLen_mean)/acLen_sd, scaled_envMean = (envMean-envMean_mean)/envMean_sd) #scale

acLen_data_interaction_envBreadthLoss <- bind_rows(acLen_data_interaction_low_envBreadthLoss, acLen_data_interaction_high_envBreadthLoss)
acLen_data_interaction_envMean <- bind_rows(acLen_data_interaction_low_envMean, acLen_data_interaction_high_envMean)

# envBreadthLoss
envBreadthLoss_interaction_envMean_fits <- readRDS("outputs/stats/envBreadthLoss_fixed_interaction_envMean_scaled.rds")
envBreadthLoss_interaction_acLen_fits <- readRDS("outputs/stats/envBreadthLoss_fixed_interaction_acLen_scaled.rds")

envBreadthLoss_interaction_envMean <- expand.grid(scaled_envBreadthLoss = seq_range(envBreadthLoss_data_interaction_envMean$scaled_envBreadthLoss, n = 100), scaled_envMean = c(-0.12,0.12)) %>% 
  add_epred_draws(envBreadthLoss_interaction_envMean_fits, ndraws = 100) %>%
  dplyr::select(scaled_envBreadthLoss, scaled_envMean, .epred) %>%
  group_by(scaled_envBreadthLoss, scaled_envMean) %>%
  summarise(p_persist = mean(.epred)) %>%
  mutate(quantile=ifelse(scaled_envMean==-0.12, '25%', '75%')) %>%
  mutate(interaction_variable = 'envMean') %>%
  dplyr::select(-scaled_envMean)
envBreadthLoss_interaction_acLen <- expand.grid(scaled_envBreadthLoss = seq_range(envBreadthLoss_data_interaction_acLen$scaled_envBreadthLoss, n = 100), scaled_acLen = c(0.06,0.14)) %>% 
  add_epred_draws(envBreadthLoss_interaction_acLen_fits, ndraws = 100) %>%
  dplyr::select(scaled_envBreadthLoss, scaled_acLen, .epred) %>%
  group_by(scaled_envBreadthLoss, scaled_acLen) %>%
  summarise(p_persist = mean(.epred)) %>%
  mutate(quantile=ifelse(scaled_acLen==0.06, '25%', '75%')) %>%
  mutate(interaction_variable = 'acLen')%>%
  dplyr::select(-scaled_acLen)
envBreadthLoss_interaction_mean <- expand.grid(scaled_envBreadthLoss = seq_range(envBreadthLoss_data$scaled_envBreadthLoss, n = 100)) %>% 
  add_epred_draws(envBreadthLoss_fixed_output, ndraws = 100) %>%
  dplyr::select(scaled_envBreadthLoss, .epred) %>%
  group_by(scaled_envBreadthLoss) %>%
  summarise(p_persist = mean(.epred)) %>%
  mutate(quantile='50%', interaction_variable = c("none"))

# envBreadthLoss_interaction_mean <- mean_envBreadthLoss %>%
#   mutate(quantile='50%', interaction_variable = c("none"))

# combine all for plotting
envBreadthLoss_interactions <- rbind(envBreadthLoss_interaction_mean, envBreadthLoss_interaction_envMean, envBreadthLoss_interaction_acLen)

a <- ggplot() +
  geom_line(data=envBreadthLoss_interactions, aes(x=(scaled_envBreadthLoss*envBreadthLoss_sd + envBreadthLoss_mean), y=p_persist, linetype=quantile, color=interaction_variable), linewidth=1.5, show.legend=F) +
  scale_y_continuous(limits=c(0,1)) +
  labs(x=expression(Environmental~Breadth~Loss~(Delta*B[e])), y=expression(Probability~of~Persistence~(p[persist]))) +
  scale_colour_manual(name = "Interaction Variable", values = c("envBreadthLoss" = "#558097", "envMean" = "#9F1D22", "acLen" = "#6C9E64", "none" = "black"), labels = c("envBreadthLoss" = expression(Delta*B[e]), "envMean" = expression(Delta*N[e]), "acLen" = expression(l[hl])))+
  scale_linetype_manual(name = "Approximate Percentile", values = c("50%" = "solid", "25%" = "dotted", "75%" = "twodash")) +
  theme(text = element_text(size = 12))



# envMean
envMean_interaction_envBreadthLoss_fits <- readRDS("outputs/stats/envMean_fixed_interaction_envBreadthLoss_scaled.rds")
envMean_interaction_acLen_fits <- readRDS("outputs/stats/envMean_fixed_interaction_acLen_scaled.rds")

envMean_interaction_envBreadthLoss <- expand.grid(scaled_envMean = seq_range(envMean_data$scaled_envMean, n = 100), scaled_envBreadthLoss = c(0.67,0.97)) %>% 
  add_epred_draws(envMean_interaction_envBreadthLoss_fits, ndraws = 100) %>%
  dplyr::select(scaled_envMean, scaled_envBreadthLoss, .epred) %>%
  group_by(scaled_envMean, scaled_envBreadthLoss) %>%
  summarise(p_persist = mean(.epred)) %>%
  mutate(quantile=ifelse(scaled_envBreadthLoss==0.67, '25%', '75%')) %>%
  mutate(interaction_variable = 'envBreadthLoss') %>%
  dplyr::select(-scaled_envBreadthLoss)
envMean_interaction_acLen <- expand.grid(scaled_envMean = seq_range(envMean_data_interaction_acLen$scaled_envMean, n = 100), scaled_acLen = c(0.06,0.14)) %>% 
  add_epred_draws(envMean_interaction_acLen_fits, ndraws = 100) %>%
  dplyr::select(scaled_envMean, scaled_acLen, .epred) %>%
  group_by(scaled_envMean, scaled_acLen) %>%
  summarise(p_persist = mean(.epred)) %>%
  mutate(quantile=ifelse(scaled_acLen==0.06, '25%', '75%')) %>%
  mutate(interaction_variable = 'acLen')%>%
  dplyr::select(-scaled_acLen)
envMean_interaction_mean <- expand.grid(scaled_envMean = seq_range(envMean_data_interaction_acLen$scaled_envMean, n = 100)) %>% 
  add_epred_draws(envMean_fixed_output, ndraws = 100) %>%
  dplyr::select(scaled_envMean, .epred) %>%
  group_by(scaled_envMean) %>%
  summarise(p_persist = mean(.epred)) %>%
  mutate(quantile='50%', interaction_variable = c("none"))

# envMean_interaction_mean <- mean_envMean %>%
#   mutate(quantile='50%', interaction_variable = c("none"))

# combine all for plotting
envMean_interactions <- rbind(envMean_interaction_mean, envMean_interaction_envBreadthLoss, envMean_interaction_acLen)

b <- ggplot() +
  geom_line(data=envMean_interactions, aes(x=(scaled_envMean*envMean_sd + envMean_mean), y=p_persist, linetype=quantile, color=interaction_variable), linewidth=1.5, show.legend=F) +
  scale_y_continuous(limits=c(0,1)) +
  labs(x=expression(Change~Environmental~Mean~(Delta*mu[e])), y=NULL) +
  scale_colour_manual(name = "Interaction Variable", values = c("envBreadthLoss" = "#558097", "envMean" = "#9F1D22", "acLen" = "#6C9E64", "none" = "black"), labels = c("envBreadthLoss" = expression(Delta*B[e]), "envMean" = expression(Delta*N[e]), "acLen" = expression(l[hl])))+
  scale_linetype_manual(name = "Approximate Percentile", values = c("50%" = "solid", "25%" = "dotted", "75%" = "twodash")) +
  theme(text = element_text(size = 12))



# acLen
acLen_interaction_envBreadthLoss_fits <- readRDS("outputs/stats/acLen_fixed_interaction_envBreadthLoss_scaled.rds")
acLen_interaction_envMean_fits <- readRDS("outputs/stats/aclen_fixed_interaction_envMean_scaled.rds")

acLen_interaction_envBreadthLoss <- expand.grid(scaled_acLen = seq_range(acLen_data_interaction_envBreadthLoss$scaled_acLen, n = 100), scaled_envBreadthLoss = c(0.67,0.97)) %>% 
  add_epred_draws(acLen_interaction_envBreadthLoss_fits, ndraws = 100) %>%
  dplyr::select(scaled_acLen, scaled_envBreadthLoss, .epred) %>%
  group_by(scaled_acLen, scaled_envBreadthLoss) %>%
  summarise(p_persist = mean(.epred)) %>%
  mutate(quantile=ifelse(scaled_envBreadthLoss==0.67, '25%', '75%')) %>%
  mutate(interaction_variable = 'envBreadthLoss') %>%
  dplyr::select(-scaled_envBreadthLoss)
acLen_interaction_envMean <- expand.grid(scaled_acLen = seq_range(acLen_data_interaction_envMean$scaled_acLen, n = 100), scaled_envMean = c(-0.12,0.12)) %>% 
  add_epred_draws(acLen_interaction_envMean_fits, ndraws = 100) %>%
  dplyr::select(scaled_acLen, scaled_envMean, .epred) %>%
  group_by(scaled_acLen, scaled_envMean) %>%
  summarise(p_persist = mean(.epred)) %>%
  mutate(quantile=ifelse(scaled_envMean==-0.12, '25%', '75%')) %>%
  mutate(interaction_variable = 'envMean')%>%
  dplyr::select(-scaled_envMean)
acLen_interaction_mean <- expand.grid(scaled_acLen = seq_range(acLen_data$scaled_acLen, n = 100)) %>% 
  add_epred_draws(acLen_fixed_output, ndraws = 100) %>%
  dplyr::select(scaled_acLen, .epred) %>%
  group_by(scaled_acLen) %>%
  summarise(p_persist = mean(.epred)) %>%
  mutate(quantile='50%', interaction_variable = c("none"))

# acLen_interaction_mean <- mean_acLen %>%
#   mutate(quantile='50%', interaction_variable = c("none"))

# combine all for plotting
acLen_interactions <- rbind(acLen_interaction_envBreadthLoss, acLen_interaction_envMean)

c <- ggplot() +
  geom_line(data = acLen_interaction_mean, aes(x=(scaled_acLen*acLen_sd + acLen_mean), y=p_persist, linetype='50%'), colour="black", linewidth=1.5) + # main
  geom_line(data=acLen_interactions, aes(x=(scaled_acLen*acLen_sd + acLen_mean), y=p_persist, linetype=quantile, color=interaction_variable), linewidth=1.5) +
  geom_line(aes(x=0, y=0, linetype='25%', colour = 'acLen'), linewidth=1.5) + # dummy plot for envBreadthLoss so we can include it in the legend below
  scale_y_continuous(limits=c(0,1)) +
  labs(x=expression(Habitat~Loss~Autocorrelation~Length~(l[hl])), y=NULL) +
  scale_colour_manual(name = "Interaction Variable", values = c("envBreadthLoss" = "#558097", "envMean" = "#9F1D22", "acLen" = "#6C9E64", "none" = "black"), labels = c("envBreadthLoss" = expression(Delta*B[e]), "envMean" = expression(Delta*N[e]), "acLen" = expression(l[hl])))+
  scale_linetype_manual(name = "Approximate Percentile", values = c("50%" = "solid", "25%" = "dotted", "75%" = "twodash")) +
  theme(text = element_text(size = 12))

interactions <- (a + b + c) +
  plot_annotation(tag_levels = 'a')  + 
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom')


# save plot as pdf
ggsave(interactions, filename = "images/interactions.pdf", height = 4.5, width = 10)


##### Supplemental: Local Sensitivity analysis #####
### Read in outputs again ###
# Get properties of post habitat loss landscapes
### Normal (HL+CC)
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
# ENVIRONMENTAL BREATH LOSS
envBreadthLoss_data <- get_specific_output(seeds, results_and_properties, 1, 50, 1) %>%
  mutate(scaled_envBreadthLoss = scale(envBreadthLoss)) #scale

# ENVIRONMENTAL MEAN
envMean_data <-  get_specific_output(seeds, results_and_properties, 51, 100, 1) %>%
  mutate(scaled_envMean = scale(envMean)) #scale

# FRAGMENTATION
acLen_data <-  get_specific_output(seeds, results_and_properties, 101, 150, 1) %>%
  mutate(scaled_acLen = scale(acLen)) #scale


# Create an empty list to store plots
plot_list <- list()

sensitivity_plots <- function(SA_num1, SA_num2, plot_list, param_name, param_val1, param_val2, is_bottom_row = FALSE) {
  ### Labels
  # Create legend labels
  legend_labels <- c(param_val1, param_val2)
  names(legend_labels) <- c(param_val1, param_val2)
  
  # For bottom row, add these x-axis labels
  x_lab_left <- if(is_bottom_row) expression(Environmental~Breadth~Loss~(Delta*B[e])) else NULL
  x_lab_middle <- if(is_bottom_row) expression(Change~Environmental~Mean~(Delta*mu[e])) else NULL
  x_lab_right <- if(is_bottom_row) expression(Habitat~Loss~Autocorrelation~Length~(l[hl])) else NULL
  
  # Scale settings for x-axis
  x_scale_settings <- if(is_bottom_row) {
    list(
      breaks = scales::pretty_breaks(n = 4),
      labels = scales::number_format(accuracy = 0.1)
    )
  } else {
    list(breaks = NULL, labels = NULL)
  }
  
  ### Read in model outputs
  # Read in random effect model outputs for first parameter value
  envBreadthLoss_random_output_1 <- readRDS(paste0("outputs/stats/envBreadthLoss_random_scaled_sensitivity", as.character(SA_num1),".rds"))
  envMean_random_output_1 <- readRDS(paste0("outputs/stats/envMean_random_scaled", as.character(SA_num1),".rds"))
  acLen_random_output_1 <- readRDS(paste0("outputs/stats/acLen_random_scaled", as.character(SA_num1),".rds"))
  
  # Read in random effect model outputs for second parameter value
  envBreadthLoss_random_output_2 <- readRDS(paste0("outputs/stats/envBreadthLoss_random_scaled_sensitivity", as.character(SA_num2),".rds"))
  envMean_random_output_2 <- readRDS(paste0("outputs/stats/envMean_random_scaled", as.character(SA_num2),".rds"))
  acLen_random_output_2 <- readRDS(paste0("outputs/stats/acLen_random_scaled", as.character(SA_num2),".rds"))
  
  ### Prep and subplots
  # Environmental breadth loss
  mean_envBreadthLoss_1 <- expand.grid(unique(envBreadthLoss_data$ENV_SEED), 
                                       seq_range(envBreadthLoss_data$scaled_envBreadthLoss, n = 100)) %>%
    rename(ENV_SEED = Var1, scaled_envBreadthLoss = Var2) %>%
    add_epred_draws(envBreadthLoss_random_output_1, ndraws = 100) %>%
    dplyr::select(scaled_envBreadthLoss, .epred) %>%
    group_by(scaled_envBreadthLoss) %>%
    summarise(p_persist = mean(.epred)) %>%
    mutate(param_value = param_val1)
  
  mean_envBreadthLoss_2 <- expand.grid(unique(envBreadthLoss_data$ENV_SEED), 
                                       seq_range(envBreadthLoss_data$scaled_envBreadthLoss, n = 100)) %>%
    rename(ENV_SEED = Var1, scaled_envBreadthLoss = Var2) %>%
    add_epred_draws(envBreadthLoss_random_output_2, ndraws = 100) %>%
    dplyr::select(scaled_envBreadthLoss, .epred) %>%
    group_by(scaled_envBreadthLoss) %>%
    summarise(p_persist = mean(.epred)) %>%
    mutate(param_value = param_val2)
  
  mean_envBreadthLoss_combined <- rbind(mean_envBreadthLoss_1, mean_envBreadthLoss_2) %>%
    ggplot(aes(x = (scaled_envBreadthLoss*envBreadthLoss_sd + envBreadthLoss_mean), 
               y = p_persist,
               linetype = param_value)) +
    geom_line(linewidth=1.2) +
    scale_linetype_manual(values = c("solid", "dashed"),
                          labels = legend_labels,
                          name = param_name) +
    scale_x_continuous(breaks = x_scale_settings$breaks,
                       labels = x_scale_settings$labels) +
    scale_y_continuous(limits=c(0,1)) +
    labs(x=x_lab_left, y=expression(p["persist"])) +
    theme_set(theme_tidybayes() + panel_border()) +
    theme(text = element_text(size = 12),
          legend.position = "none",
          axis.title.x = element_text(size = 10, margin = margin(t = 10)),
          axis.text.x = if(is_bottom_row) element_text(size = 10) else element_blank())
  
  plot_list[[paste0(SA_num1, "l")]] <- mean_envBreadthLoss_combined
  
  # Environmental mean
  mean_envMean_1 <- expand.grid(unique(envMean_data$ENV_SEED), 
                                seq_range(envMean_data$scaled_envMean, n = 100)) %>%
    rename(ENV_SEED = Var1, scaled_envMean = Var2) %>%
    add_epred_draws(envMean_random_output_1, ndraws = 100) %>%
    dplyr::select(scaled_envMean, .epred) %>%
    group_by(scaled_envMean) %>%
    summarise(p_persist = mean(.epred)) %>%
    mutate(param_value = param_val1)
  
  mean_envMean_2 <- expand.grid(unique(envMean_data$ENV_SEED), 
                                seq_range(envMean_data$scaled_envMean, n = 100)) %>%
    rename(ENV_SEED = Var1, scaled_envMean = Var2) %>%
    add_epred_draws(envMean_random_output_2, ndraws = 100) %>%
    dplyr::select(scaled_envMean, .epred) %>%
    group_by(scaled_envMean) %>%
    summarise(p_persist = mean(.epred)) %>%
    mutate(param_value = param_val2)
  
  mean_envMean_combined <- rbind(mean_envMean_1, mean_envMean_2) %>%
    ggplot(aes(x = (scaled_envMean*envMean_sd + envMean_mean), 
               y = p_persist,
               linetype = param_value)) +
    geom_line(linewidth=1.2) +
    scale_linetype_manual(values = c("solid", "dashed"),
                          labels = legend_labels,
                          name = param_name) +
    scale_x_continuous(breaks = x_scale_settings$breaks,
                       labels = x_scale_settings$labels) +
    scale_y_continuous(limits=c(0,1), breaks=NULL, labels=NULL) +
    labs(x=x_lab_middle, y=NULL) +
    theme_set(theme_tidybayes() + panel_border()) +
    theme(text = element_text(size = 12),
          legend.position = "none",
          axis.title.x = element_text(size = 10, margin = margin(t = 10)),
          axis.text.x = if(is_bottom_row) element_text(size = 10) else element_blank())
  
  plot_list[[paste0(SA_num1, "m")]] <- mean_envMean_combined
  
  # Autocorrelation length
  mean_acLen_1 <- expand.grid(unique(acLen_data$ENV_SEED), 
                              seq_range(acLen_data$scaled_acLen, n = 100)) %>%
    rename(ENV_SEED = Var1, scaled_acLen = Var2) %>%
    add_epred_draws(acLen_random_output_1, ndraws = 100) %>%
    dplyr::select(scaled_acLen, .epred) %>%
    group_by(scaled_acLen) %>%
    summarise(p_persist = mean(.epred)) %>%
    mutate(param_value = param_val1)
  
  mean_acLen_2 <- expand.grid(unique(acLen_data$ENV_SEED), 
                              seq_range(acLen_data$scaled_acLen, n = 100)) %>%
    rename(ENV_SEED = Var1, scaled_acLen = Var2) %>%
    add_epred_draws(acLen_random_output_2, ndraws = 100) %>%
    dplyr::select(scaled_acLen, .epred) %>%
    group_by(scaled_acLen) %>%
    summarise(p_persist = mean(.epred)) %>%
    mutate(param_value = param_val2)
  
  mean_acLen_combined <- rbind(mean_acLen_1, mean_acLen_2) %>%
    ggplot(aes(x = (scaled_acLen*acLen_sd + acLen_mean), 
               y = p_persist,
               linetype = param_value)) +
    geom_line(linewidth=1.2) +
    scale_linetype_manual(values = c("solid", "dashed"),
                          labels = legend_labels,
                          name = param_name) +
    scale_x_continuous(breaks = x_scale_settings$breaks,
                       labels = x_scale_settings$labels) +
    scale_y_continuous(limits=c(0,1), breaks=NULL, labels=NULL) +
    labs(x=x_lab_right, y=NULL) +
    theme_set(theme_tidybayes() + panel_border()) +
    theme(text = element_text(size = 12),
          legend.position = "right",
          legend.box.margin = margin(0, 0, 0, -10),
          legend.margin = margin(0, 0, 0, 0),
          legend.key.size = unit(1.5, "lines"),
          legend.key.width = unit(2, "lines"),
          axis.title.x = element_text(size = 10, margin = margin(t = 10)),
          axis.text.x = if(is_bottom_row) element_text(size = 10) else element_blank())
  
  plot_list[[paste0(SA_num1, "r")]] <- mean_acLen_combined
  
  return(plot_list)
}

# Create parameter pairs with their values
param_pairs <- list(
  list(1, 2, expression(sigma["p"]), "0.05", "0.025"),
  list(3, 4, expression(lambda["0"]), "0.125", "0.5"),
  list(5, 6, expression(sigma["f"]), "0.125", "0.5"),
  list(7, 8, expression(sigma["QTL"]), "0.05", "0.2"),
  list(9, 10, "K", "500", "250"),
  list(11, 12, expression(l["e"]), "0.05", "0.2"),
  list(13, 14, expression(p["QTL"]), "0.025", "0.1"),
  list(15, 16, expression(delta["e"]), "2.4", "3.75"),
  list(17, 18, expression(p["hl"]), "8/15", "5/6")
)

# Create plots with bottom row flag
plot_list <- list()
for(i in seq_along(param_pairs)) {
  pair <- param_pairs[[i]]
  is_bottom <- (i == length(param_pairs))
  plot_list <- sensitivity_plots(pair[[1]], pair[[2]], plot_list, pair[[3]], 
                                 pair[[4]], pair[[5]], is_bottom)
}

# Combine plots with reduced spacing
SA_plot <- (plot_list$'1l' + plot_list$'1m' + plot_list$'1r') /
  (plot_list$'3l' + plot_list$'3m' + plot_list$'3r') /
  (plot_list$'5l' + plot_list$'5m' + plot_list$'5r') /
  (plot_list$'7l' + plot_list$'7m' + plot_list$'7r') /
  (plot_list$'9l' + plot_list$'9m' + plot_list$'9r') /
  (plot_list$'11l' + plot_list$'11m' + plot_list$'11r') /
  (plot_list$'13l' + plot_list$'13m' + plot_list$'13r') /
  (plot_list$'15l' + plot_list$'15m' + plot_list$'15r') /
  (plot_list$'17l' + plot_list$'17m' + plot_list$'17r') +
  plot_layout(heights = unit(rep(1, 9), "null")) +  # Equal heights for all rows
  plot_annotation(tag_levels = "a") &
  theme(plot.margin = margin(2, 2, 2, 2))  # Reduce overall plot margins

# Save plot with adjusted height
ggsave(SA_plot, filename = "images/SA_Combined.pdf", height = 10, width = 10)  # Reduced height

##### Supplemental: Global Sensitivity analysis #########
# Read data again
### Normal (HL+CC)
file_list <- list.files(path = "outputs/spatial/simulation_output", pattern = "\\globalSA.txt$", full.names = TRUE)
data_list <- lapply(file_list, read.csv)
all_out <- do.call(rbind, data_list)

# Get the list of files
file_list <- list.files(path = "outputs/spatial/simulation_output", pattern = "\\globalSA.txt$", full.names = TRUE)

# Extract batch numbers (refers to which group of global SA this is from) from filenames using regular expressions
batch_numbers <- as.numeric(gsub(".*?(\\d+)globalSA\\.txt$", "\\1", file_list))

# Read the files and add batch numbers
data_list <- mapply(function(file, batch) {
  df <- read.csv(file)
  df$batch <- batch
  return(df)
}, file_list, batch_numbers, SIMPLIFY = FALSE)

# Combine all data frames
all_out <- do.call(rbind, data_list)

# Group outputs by organism type, then get summary statistics for each PA design by looking at results from each seed
grouped_data <- all_out %>%
  dplyr::group_by(batch, ENV_SEED, HL_SEED)

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
# ENVIRONMENTAL BREATH LOSS
envBreadthLoss_data <- get_specific_output(seeds, results_and_properties, 1, 50, 1) %>%
  mutate(scaled_envBreadthLoss = scale(envBreadthLoss)) #scale

# ENVIRONMENTAL MEAN
envMean_data <-  get_specific_output(seeds, results_and_properties, 51, 100, 1) %>%
  mutate(scaled_envMean = scale(envMean)) #scale

# FRAGMENTATION
acLen_data <-  get_specific_output(seeds, results_and_properties, 101, 150, 1) %>%
  mutate(scaled_acLen = scale(acLen)) #scale

# Supplemental - what if spaghetti plots were just made using smoothers?
envBreathLoss_main_gam <- envBreadthLoss_data %>%
  dplyr::select(batch, ENV_SEED, HL_SEED, SLIM_SEED, persist, envBreadthLoss) %>%
  group_by(batch, ENV_SEED, HL_SEED) %>%
  summarise(p_persist = mean(persist), envBreadthLoss=mean(envBreadthLoss)) %>% # taking mean of envBreadthLoss OK since it is always the same for every ENV_SEED/HL_SEED combo
  ungroup() %>%
  ggplot(aes(x = envBreadthLoss, y = p_persist)) + 
  geom_smooth(aes(group=batch), method="loess", span=1, se=F, color='grey40', linewidth=0.5) +
  geom_smooth(method="loess", span=1, se=F, color='black', linewidth=1.5) +
  scale_y_continuous(limits=c(0,1)) +
  labs(x=expression(Environmental~Breadth~Loss~(Delta*B[e])), y = expression(Probability~of~Persistence~(p[persist])))  +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))

envMean_main_gam <- envMean_data %>%
  dplyr::select(batch, ENV_SEED, HL_SEED, SLIM_SEED, persist, envMean) %>%
  group_by(batch, ENV_SEED, HL_SEED) %>%
  summarise(p_persist = mean(persist), envMean=mean(envMean)) %>% # taking mean of envMean OK since it is always the same for every ENV_SEED/HL_SEED combo
  ungroup() %>%
  ggplot(aes(x = envMean, y = p_persist)) + 
  geom_smooth(aes(group=batch), method="loess", span=1, se=F, color='grey40', linewidth=0.5) +
  geom_smooth(method="loess", span=1, se=F, color='black', linewidth=1.5) +
  scale_y_continuous(limits=c(0,1), labels=NULL, breaks=NULL) +
  labs(x=expression(Change~Environmental~Mean~(Delta*mu[e])), y = NULL)  +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))

acLen_main_gam <- acLen_data %>%
  dplyr::select(batch, ENV_SEED, HL_SEED, SLIM_SEED, persist, acLen) %>%
  group_by(batch, ENV_SEED, HL_SEED) %>%
  summarise(p_persist = mean(persist), acLen=mean(acLen)) %>% # taking mean of acLen OK since it is always the same for every ENV_SEED/HL_SEED combo
  ungroup() %>%
  ggplot(aes(x = acLen, y = p_persist)) + 
  geom_smooth(aes(group=batch), method="loess", span=1, se=F, color='grey40', linewidth=0.5) +
  geom_smooth(method="loess", span=1, se=F, color='black', linewidth=1.5) +
  scale_y_continuous(limits=c(0,1), labels=NULL, breaks=NULL) +
  labs(x=expression(Habitat~Loss~Autocorrelation~Length~(l[hl])), y = NULL)  +
  theme_set(theme_tidybayes() + panel_border())+
  theme(text = element_text(size = 12))

main <- (envBreathLoss_main_gam + envMean_main_gam + acLen_main_gam)  + 
  plot_annotation(tag_levels = 'a')# + 
#plot_layout(heights = c(2, 1))

# save plot as pdf
ggsave(main, filename = "images/globalSA_LOESSversion.pdf", height = 4, width = 10)



##### Extracting Slope+Intercept Estimates ######
# SA_num <- 18
# be <- readRDS(paste0("outputs/stats/envBreadthLoss_random_scaled_sensitivity", as.character(SA_num),".rds"))
# me <- readRDS(paste0("outputs/stats/envMean_random_scaled", as.character(SA_num),".rds"))
# ac <- readRDS(paste0("outputs/stats/acLen_random_scaled", as.character(SA_num),".rds"))
