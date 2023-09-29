# This script runs the full pipeline for data from fixed images of cells with only KHSV terminal repeats (Figure 2)

#####
#  Setup 
#####

# load libraries 
library(tidyverse)
library(here)
# For plotting
library(patchwork)
library(RColorBrewer)
library(ggExtra)
library(bayesplot)
# for markov chain convergence statistics
library(rstan)
# for fitting a poisson distribution
library(MASS)
# for running in parallel
library(furrr)
library(purrr)
library(foreach)
library(parallel)

select <- dplyr::select
theme_set(theme_bw())

# # # create a cluster
# n.cores <- parallel::detectCores() - 3
# my.cluster <- parallel::makeCluster(
#   n.cores,
#   type = "PSOCK",
#   rscript_args = c("--no-init-file", "--no-site-file", "--no-environ")
# )
# doParallel::registerDoParallel(cl = my.cluster)
# foreach::getDoParRegistered()


# Check wd
here()

# Read in helper functions 
source(here("scripts", "likelihood_functions.R"))

# Folder for results
results_folder <- here("results", "Figure_3")
dir.create(results_folder)

#####
# Read in data
#####

daughter_cell_file <- "Fig3 dividing cells.xlsx"
mother_cell_file <- "Fig3 non dividing cells.xlsx"

# Daughter cell data
episome_data <- readxl::read_excel(here("data", "derived", daughter_cell_file))
# Mother cell data
mother_cell_data <- readxl::read_excel(here("data", "derived", mother_cell_file))

# reformat the daughter cell data
names(episome_data) <- gsub(" ", "_", tolower(names(episome_data)))
episome_data <- rbind(
  select(episome_data, mother_cell_id, contains("daughter1")) %>% rename_with(~gsub("_daughter1",  "", .x)) %>% mutate(daughter_cell = 1),
  select(episome_data, mother_cell_id, contains("daughter2")) %>% rename_with(~gsub("_daughter2",  "", .x)) %>% mutate(daughter_cell = 2)
) %>% 
  select(mother_cell_id, daughter_cell, cluster, everything()) %>% 
  arrange(mother_cell_id, daughter_cell, cluster) %>% 
  filter(!is.na(total_cluster_intensity), total_cluster_intensity !=0) %>% 
  mutate(cluster_id  = paste0("n", 1:nrow(.)),
         cell_id = paste(mother_cell_id, daughter_cell, sep = "_"))

# reformat the mother cell data
names(mother_cell_data) <- gsub(" ", "_", tolower(names(mother_cell_data)))
mother_cell_data <- mother_cell_data %>% 
  # remove clusters with no intensity data
  filter(!is.na(total_cluster_intensity)) %>%  
  mutate(cell_id = paste(image, cell, sep = "_"),
         cluster_id = paste0("n", 1:nrow(.))) 


#####
# Gibbs sampling to infer the number of episomes per cluster in each daughter cell
#####

#  set initial conditions:
# We will assume that the value of mu is bounded between the ranges observed in the data. 
# Thus, we will select initial conditions within this range. We will only try one value of tau that is small enough 
# (meaning sigma is large enough) that the chain will converge based on chains run in simulated data. 
# We will determine the initial conditions for nk based on the observed data and the initial condition for mu.
n_iterations <- 100000
burn_in <- 5000
tau0 <- 1e-5 
ICs <- episome_data %>% mutate(ratio = total_cluster_intensity/min_episome_in_cluster) %>% filter(!is.na(ratio)) %>% pull(ratio) %>% summary
mu0 <- ICs[-3]

cat("tau0:", tau0, "\nmu0:", mu0)

# run 5 chains
chain1 <- run_gibbs(tau0, mu0[1], episome_data$total_cluster_intensity, n_iterations)
chain2 <- run_gibbs(tau0, mu0[2], episome_data$total_cluster_intensity, n_iterations)
chain3 <- run_gibbs(tau0, mu0[3], episome_data$total_cluster_intensity, n_iterations)
chain4 <- run_gibbs(tau0, mu0[4], episome_data$total_cluster_intensity, n_iterations)
chain5 <- run_gibbs(tau0, mu0[5], episome_data$total_cluster_intensity, n_iterations)

save(chain1, chain2, chain3, chain4, chain5, file = here(results_folder, "daughter_cell_cluster_Gibbs_samples.RData"))

# # TEST running Gibbs in parallel
# gibbs_samples <- foreach(i = 1:length(mu0), .packages = "tidyverse") %dopar% {
#   run_gibbs(tau0, mu0[i], episome_data$total_cluster_intensity, n_iterations)
# }

# combine all chains together and remove burn-in period
all_chains <- rbind(chain1, chain2, chain3, chain4, chain5) %>%
  mutate(chain = as.factor(rep(c("chain1", "chain2", "chain3", "chain4", "chain5"), each = n_iterations)))  %>%
  filter(iteration > burn_in)

#####
# Analyze results of Gibbs
#####


p1 <- all_chains %>% 
  ggplot(aes(iteration, mu, color = chain)) + 
  geom_line() + 
  geom_point(shape = NA) + 
  theme(legend.position = "bottom") + 
  ylim(c(0, max(all_chains$mu))) +
  labs(title = "Trace of mu")

(p1 <- ggMarginal(p1, margins = "y", groupColour = T))

p2 <- all_chains %>% 
  ggplot(aes(iteration, 1/tau, color = chain)) + 
  geom_line() + 
  geom_point(shape = NA) + 
  theme(legend.position = "bottom") +
  ylim(c(0, max(1/all_chains$tau))) + 
  labs(title = "Trace of sigma^2")

(p2 <- ggMarginal(p2, margins = "y", groupColour = T))


# Based on the median mu and sigma above, the inferred distribution of intensities for a single cluster is:
# plot the distribution of cluster intensities with inferred parameters:
x <- seq(0, 4000, length.out = 500)
inferred_mu <- median(all_chains$mu)
inferred_sigma2 <- median(1/all_chains$tau)
cat('mean:', inferred_mu, '\nvariance:', inferred_sigma2, '\nstandard deviation:', sqrt(inferred_sigma2))
(p3 <- ggplot(tibble(x, y = dnorm(x, inferred_mu, sqrt(inferred_sigma2))), aes(x, y)) + 
    geom_line() + 
    labs(x = "I_k", y = "probability density", title = "Inferred Distribution of intensity for a single cluster"))


## infer nk for each cluster according to the mode of the posterior distribution
all_chain_ns <- all_chains %>% select(-mu, -tau)  %>% pivot_longer(!c(iteration, chain), names_to = "cluster")
modes <- all_chain_ns %>% add_count(chain, cluster, value) %>% group_by(chain, cluster) %>% mutate(mode = value[which.max(n)]) %>% ungroup()

distinct_modes <- modes %>% distinct(chain, cluster, mode)

p4 <- distinct_modes %>% count(mode, chain) %>% mutate(freq = n/sum(n)) %>% 
  ggplot(aes(mode, freq)) + geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.1) + 
  facet_wrap(~chain) +
  labs(x = "mode of nk", title = "histogram of inferred nk given by mode of the posterior distribution")

# How do these inferences relate to the observed intensity data?
p5  <- episome_data %>% 
  mutate(cluster = cluster_id) %>% 
  merge(distinct_modes) %>% 
  ggplot(aes(total_cluster_intensity, color = as.factor(mode))) + 
  geom_density() + 
  labs(title = "Distribution of cluster intensities for clusters with each nk", color  = "mode")

p6 <- episome_data %>% 
  mutate(cluster = cluster_id) %>% 
  merge(distinct_modes) %>% 
  ggplot(aes(total_cluster_intensity, as.factor(min_episome_in_cluster))) + 
  geom_jitter(aes(color = as.factor(mode), shape = chain)) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  labs(title = "Distribution of cluster intensities stratified by minimum number of episomes",
       subtitle = "colored by inferred nk",
       color = "mode", y = "Min # episome in cluster", x = "Total  Cluster Intensity")


#####
# Calculate number of episomes per cell at each iteration
#####

# We are actually interested in the number of episomes per cell, not per cluster. 
# To find this, we will sum samples of the posterior distribution for each cluster in a cell.

# function to get number of episomes in cell:
get_n_cell <- function(cellID, all_chains, data){
  # pull the clusters in the cell
  coi <- data %>% filter(cell_id == cellID) %>% pull(cluster_id)  
  # Sum the sampled number of episomes in these clusters for each iteration
  n_cell <- apply(all_chains[,coi,drop = F], 1, sum)
  all_chains %>% select(iteration, chain) %>% mutate(cell_id = cellID, number_of_episomes = n_cell)
}

cell_samples_long <- as.list(unique(episome_data$cell_id)) %>% lapply(get_n_cell, all_chains, episome_data) %>% bind_rows() 
cell_samples <- cell_samples_long %>% pivot_wider(names_from = cell_id, values_from = number_of_episomes)

save(cell_samples_long, cell_samples, file = here(results_folder, "daughter_cell_Gibbs_samples.RData"))

# look at posterior probabilities
temp_df <- cell_samples_long %>% filter(chain == "chain1") %>% 
  group_by(cell_id) %>% add_count(number_of_episomes) %>% 
  mutate(mode = number_of_episomes[which.max(n)]) %>% 
  ungroup 

p7 <- temp_df %>% 
  ggplot(aes(number_of_episomes, group = cell_id)) + 
  geom_density() + facet_wrap(~mode, scales = "free_y", labeller = "label_both") + 
  labs(x = "number of episomes per cell",  title = "Posterior probabilities for number of episomes per cell")

# look at distribution of inferred episomes per cell according to  mode
cell_modes <- cell_samples_long %>% count(chain, cell_id, number_of_episomes) %>%
  group_by(cell_id, chain) %>% filter(n == max(n)) %>% rename(mode = number_of_episomes) %>% ungroup

p8  <- cell_modes %>% count(chain, mode) %>% mutate(freq = n/sum(n)) %>% 
  ggplot(aes(mode, freq)) + geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.1) + 
  facet_wrap(~chain) +
  labs(x = "mode of number of episomes per cell", title = "histogram of inferred number per cell given by mode of the posterior distribution")


#####
# Convergence Statistics for Daughter Cell Gibbs Sampling
#####

convergence <- convergence_results(all_chains)

# Rhat values near 1 indicate equilibrium
p9 <- mcmc_rhat(convergence$Rhat) +
  ggtitle("Rhat")

p10 <- mcmc_neff(convergence$ESS_bulk/(n_iterations- burn_in)) + 
  ggtitle("Bulk Effective Sample Size Ratio")
p11 <- mcmc_neff(convergence$ESS_tail/(n_iterations- burn_in)) + 
  ggtitle("Tail Effective Sample Size Ratio")


plot_list <- list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11)

ggsave(here(results_folder, "test.pdf"), 
       gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1),
       width = 7, height = 7)

# convergence of mu and tau:
convergence %>% slice(1:2)

#####
# Gibbs sampling to determine the distribution of initial number of episomes
#####

# set initial conditions (use same number of iterations and burn-in as daughter cells)
tau0 <- 1e-5 
ICs <- mother_cell_data %>% mutate(ratio = total_cluster_intensity/min_episome_in_cluster) %>% filter(!is.na(ratio)) %>% pull(ratio) %>% summary
mu0 <- ICs[-3]

# Run Gibs
intensities <- mother_cell_data %>% pull(total_cluster_intensity)
chain1_m <- run_gibbs(tau0, mu0[1], intensities, n_iterations)
chain2_m <- run_gibbs(tau0, mu0[2], intensities, n_iterations)
chain3_m <- run_gibbs(tau0, mu0[3], intensities, n_iterations)
chain4_m <- run_gibbs(tau0, mu0[4], intensities, n_iterations)
chain5_m <- run_gibbs(tau0, mu0[5], intensities, n_iterations)

save(chain1_m, chain2_m, chain3_m, chain4_m, chain5_m, file = here(results_folder, "mother_cell_cluster_Gibbs_samples.RData"))

all_chains_initial <- rbind(chain1_m %>% mutate(chain = "chain 1"),
                            chain2_m %>% mutate(chain = "chain 2"),
                            chain3_m %>% mutate(chain = "chain 3"),
                            chain4_m %>% mutate(chain = "chain 4"),
                            chain5_m %>% mutate(chain = "chain 5")) %>% 
  filter(iteration > burn_in)

# report the inferred mean and variance for the mother cells
all_chains_initial %>% 
  group_by(chain) %>% summarise(m = median(mu), sigma2 = median(1/tau), sd = sd(mu))

#####
# Analyze the result of Gibbs on the mother cells
#####

p1 <- all_chains_initial %>% 
  ggplot(aes(iteration, mu, color = chain)) + 
  geom_line() + 
  geom_point(shape = NA) + 
  theme(legend.position = "bottom") + 
  ylim(c(0, max(all_chains_initial$mu))) +
  labs(title = "Trace of mu")

p1 <- ggMarginal(p1, margins = "y", groupColour = T)

p2 <- all_chains_initial %>% 
  ggplot(aes(iteration, 1/tau, color = chain)) + 
  geom_line() + 
  geom_point(shape = NA) + 
  theme(legend.position = "bottom") +
  ylim(c(0, max(1/all_chains_initial$tau))) + 
  labs(title = "Trace of sigma^2")

p2 <- ggMarginal(p2, margins = "y", groupColour = T)

#####
# Calculate number of episomes per cell at each iteration
#####

## Sum sampled number per cluster at each iteration to get number per cell (IS THERE A FASTER WAY?)
mother_cell_samples_long <- as.list(unique(mother_cell_data$cell_id)) %>% lapply(get_n_cell, all_chains = all_chains_initial, data = mother_cell_data) %>% bind_rows() 
mother_cell_samples <- cell_samples_long %>% pivot_wider(names_from = cell_id, values_from = number_of_episomes)

save(mother_cell_samples, mother_cell_samples_long, file = here(results_folder, "mother_cell_Gibbs_samples.RData"))

# Calculate the mode number of episomes per cell 
p3 <- mother_cell_samples_long %>% 
  count(chain, cell_id, number_of_episomes) %>% 
  group_by(chain, cell_id) %>% filter(n == max(n)) %>% ungroup %>% 
  count(mode=number_of_episomes, chain) %>% mutate(freq = n/sum(n)) %>% 
  ggplot(aes(mode, freq)) + geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.1) + 
  facet_wrap(~chain) +
  labs(x = "mode of nk", title = "histogram of inferred nk given by mode of the posterior distribution")


#####
# Convergence Statistics for Mother Cell Gibbs Sampling
#####

convergence_initial <- convergence_results(all_chains_initial)

# Rhat values near 1 indicate equilibrium
p4 <- mcmc_rhat(convergence_initial$Rhat) + 
  ggtitle("Rhat")

p5 <- mcmc_neff(convergence_initial$ESS_bulk/(n_iterations- burn_in)) + 
  ggtitle("Bulk Effective Sample Size Ratio")
p6 <- mcmc_neff(convergence_initial$ESS_tail/(n_iterations- burn_in)) + 
  ggtitle("Tail Effective Sample Size Ratio")

# convergence of mu and tau:
convergence_initial %>% slice(1:2)


#####
# Fitting distribution of initial number of episomes
#####

# we will construct the distribution of initial number of episomes using the full set of samples from 
# the posterior distribution
# using the mode of the distribution for each cell provides the same estimate for the median with a slighlty 
# larger mean, and the distribution is not a smooth.

# # get modes of each cell
# temp <- mother_cell_samples_long %>% count(cell_id,  number_of_episomes) %>%
#   group_by(cell_id) %>% filter(n == max(n))
# 
# mother_cell_samples_long %>% ggplot(aes(number_of_episomes, y = after_stat(density))) +
#   geom_histogram(aes(color = "whole posterior"), fill = NA) +
#   geom_histogram(data = temp,
# aes(color = "modes"), fill = NA)

# fitting a binomial distribution to the whole posterior results in a very large shape parameter, meaning
# a poisson will fit will
# fitdistr(mother_cell_samples_long$number_of_episomes, "negative binomial")$estimate

(lambda <- fitdistr(mother_cell_samples_long$number_of_episomes, "Poisson")$estimate)

test_lambda <- mother_cell_samples_long %>% filter(chain == 'chain 1') %>% 
  group_by(iteration) %>% 
  summarise(lambda = fitdistr(number_of_episomes, "Poisson")$estimate) %>% 
  ungroup()

p7  <- mother_cell_samples_long %>% 
  count(number_of_episomes) %>% mutate(freq = n/sum(n)) %>% 
  ggplot(aes(number_of_episomes, freq)) + 
  geom_bar(aes(color = "whole posterior"), fill = NA, stat = "identity") +
  geom_point(data = tibble(number_of_episomes = 1:max(mother_cell_samples_long$number_of_episomes), 
                           y= dpois(number_of_episomes, lambda = lambda)),
             aes(y  = y)) + 
  geom_point(data = tibble(number_of_episomes = 1:max(mother_cell_samples_long$number_of_episomes), 
                           y= dpois(number_of_episomes, lambda = 2.28)),
             aes(y  = y, color = "test")) + 
  labs(title = "Sampled distribtion of initial number of episomes",
       subtitle = str_interp("poisson(lambda = ${round(lambda, 2)}), median = ${median(mother_cell_samples_long$number_of_episomes)}"))

# # The distribution fit to the modes of each posterior distribution is very similar to the full distribution,
# # with a slightly smaller mean 
# fitdistr(temp$number_of_episomes, "Poisson")$estimate
# 
# temp %>% ggplot(aes(number_of_episomes)) + 
#   geom_histogram(aes(color = "modes", y = after_stat(density)), fill = NA) +
#   geom_point(data = tibble(number_of_episomes = 1:max(temp$number_of_episomes), 
#                            y= dpois(number_of_episomes, lambda = 4.317073 )),
#              aes(y  = y))


plot_list2 <- list(p1, p2, p3, p4, p5, p6, p7)
ggsave(here(results_folder, "mother_cell_gibbs_results.pdf"), 
       gridExtra::marrangeGrob(plot_list2, nrow=1, ncol=1),
       width = 7, height = 7)


#####
# Make sure inferred distribution of initial episomes makes sense with the inferred episomes per daughter cell
#####

sampled_daughter_cells <- cell_samples_long %>% filter(chain == "chain1") %>%
  left_join(distinct(episome_data, cell_id, mother_cell_id, daughter_cell), by = c("cell_id")) %>% 
  mutate(daughter_cell = paste0("X", daughter_cell)) %>% 
  # make one row per mother cell per sample (iteration:chain)
  select(-cell_id) %>% 
  pivot_wider(names_from = daughter_cell, values_from = number_of_episomes, values_fill = 0) 

initial_episome_constraints <- sampled_daughter_cells %>% mutate(max = X1 + X2, min = ceiling(max/2)) %>% 
  pivot_longer(c(min, max)) %>% 
  select(iteration, chain, cell_id = mother_cell_id, number_of_episomes = value, name) %>% 
  rbind(mother_cell_samples_long %>% mutate(name = "sampled")) %>% 
  count(name, number_of_episomes) %>% 
  group_by(name) %>% mutate(freq = n/sum(n)) %>% 
  ggplot(aes(number_of_episomes, freq, color = name))  + 
  geom_bar(stat = "identity", fill = NA) + 
  labs(x = "Initial number of episomes", y= "Frequency", color = "constraint") + 
  facet_wrap(~name, ncol = 1) + 
  theme(legend.position = "none") + 
  scale_x_continuous(breaks =0:15)

# The inferred distribution of mother cells does not fall bewteen the constraints of the daughter cells
# fit distributions to them

sampled_daughter_cells %>% mutate(max = X1 + X2) %>% pull(max) %>% fitdistr( "Poisson")
sampled_daughter_cells %>% mutate(min = ceiling((X1 + X2)/2)) %>% pull(min) %>% fitdistr( "Poisson")


compare_daughter_cell_posteriors <- sampled_daughter_cells %>% 
  count(X1, X2) %>% mutate(p = n/sum(n)) %>% 
  ggplot(aes(X1, X2, fill = p)) + 
  geom_tile() + 
  scale_x_continuous(breaks = 0:max(sampled_daughter_cells$X1)) + 
  scale_y_continuous(breaks = 0:max(sampled_daughter_cells$X2)) + 
  scale_fill_viridis_c() + 
  labs(title = "Fraction of sampled X1 and X2 combinations across all cells",
       fill = "fraction") 

## Histograms of episomes per cell:
plot_histograms <- function(sampled_daughter_cells, total_only = T){
  temp_df <- sampled_daughter_cells %>% 
    filter(chain == "chain1") %>% 
    pivot_longer(starts_with("X")) %>% 
    count(mother_cell_id, name, value) %>% 
    group_by(mother_cell_id, name) %>% filter(n == max(n))  %>% 
    select(-n) %>% 
    pivot_wider() %>% 
    ungroup %>% 
    mutate(temp_X1 = ifelse(X1 > X2, X1, X2),
           temp_X2 = ifelse(X2 < X1, X2, X1),
           total = X1 + X2) %>% 
    select(-X1, -X2, X1 = temp_X1, X2 = temp_X2, total) 
  
  pair_histogram <- temp_df %>% mutate(pair = paste0("(", X1, ",", X2, ")")) %>% 
    count(pair) %>% 
    arrange(desc(n)) %>% 
    slice(1:10) %>% 
    ggplot(aes(pair, n)) + 
    geom_bar(stat = "identity") + 
    labs(x = "Number of of episomes in each daughter cell (X1, X2)",
         y = "Number of daughter cell pairs")
  
  print(temp_df %>% count(X1 == X2))
  
  if(total_only){
    distribution_df <- rbind(
      temp_df %>% 
        select(!starts_with("X")) %>% 
        pivot_longer(!mother_cell_id, names_to = "type", values_to = "number_of_episomes") %>% 
        rename(cell_id = mother_cell_id),
      
      mother_cell_samples_long %>% filter(chain == "chain 1") %>% 
        count(cell_id, number_of_episomes) %>% 
        group_by(cell_id) %>% filter(n == max(n)) %>% 
        ungroup() %>% select(-n)  %>% 
        mutate(type = "mother cell")
    ) %>% 
      mutate(type = factor(type, levels= c("mother cell", "total"),
                           labels = paste("Frequency in\n", c("non-dividing cells", "daughter cell pairs")))) %>% 
      count(type, number_of_episomes) %>% 
      group_by(type) %>% 
      mutate(freq = n/sum(n))
    
  }else{
    distribution_df <- rbind(
      temp_df %>% 
        pivot_longer(!mother_cell_id, names_to = "type", values_to = "number_of_episomes") %>% 
        rename(cell_id = mother_cell_id),
      
      mother_cell_samples_long %>% filter(chain == "chain 1") %>% 
        count(cell_id, number_of_episomes) %>% 
        group_by(cell_id) %>% filter(n == max(n)) %>% 
        ungroup() %>% select(-n)  %>% 
        mutate(type = "mother cell")
    ) %>% 
      mutate(type = factor(type, levels= c("mother cell", "total", "X1", "X2"),
                           labels = paste("Frequency in\n", c("non-dividing cells", "daughter cell pairs",
                                                              "daughter cell with\nmore episomes", 
                                                              "daughter cell with\nfewer episomes")))) %>% 
      count(type, number_of_episomes) %>% 
      group_by(type) %>% 
      mutate(freq = n/sum(n))
    
  }
  
  all_dist <- distribution_df %>%
    ggplot(aes(number_of_episomes, freq)) +
    geom_bar(stat = "identity") +
    facet_wrap(~type, switch = "y", ncol = 1) +
    labs(x = "number of episomes") +
    theme(legend.position = "none", strip.placement = "outside",
          strip.background = element_blank(),
          axis.title.y = element_blank(), panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = c(0:100))
  
  # all_dist <- distribution_df %>% filter(grepl("mother", type)) %>% 
  #   ggplot(aes(number_of_episomes, freq)) + 
  #   geom_bar(stat = "identity") + 
  #   # facet_wrap(~type, switch = "y", ncol = 1) + 
  #   labs(x = "number of episomes", y = "frequency") +
  #   theme(legend.position = "none", strip.placement = "outside",
  #         strip.background = element_blank(),
  #         panel.grid.minor = element_blank()) + 
  #   scale_x_continuous(breaks = c(0:13), limits = c(0,max(distribution_df$number_of_episomes)+0.5)) + 
  #   
  #   distribution_df %>% filter(grepl("both", type)) %>% 
  #   ggplot(aes(number_of_episomes, freq)) + 
  #   geom_bar(stat = "identity") + 
  #   # facet_wrap(~type, switch = "y", ncol = 1) + 
  #   labs(x = "number of episomes", y = "frequency") +
  #   theme(legend.position = "none", strip.placement = "outside",
  #         strip.background = element_blank(),
  #         panel.grid.minor = element_blank()) + 
  #   scale_x_continuous(breaks = c(0:100), limits = c(0,max(distribution_df$number_of_episomes)+0.5), expand = expansion()) + 
  #   
  #   plot_layout(ncol = 1)
  # 
  
  return(list(pair_histogram, all_dist))
}

figure_plots <- plot_histograms(sampled_daughter_cells )

ggsave(here(results_folder, "pairs_histogram.png"), figure_plots[[1]], width = 5, height = 3)
ggsave(here(results_folder, "all_histograms.png"), figure_plots[[2]], width = 5, height = 5)


#####
# Maximum likelihood estimation with grid search
#####

# randomly choose 100 samples from the markov chain to use as inferred value of number of episomes per cell:
n_samples <- 100
samples <- sample_n(cell_samples, n_samples) %>% mutate(sample_id = 1:n_samples)

# reformat samples so that we can apply the grid search function: 
# need a column for mother cell id, number of cells in daughter cell 1, and number of cells in daughter cecll 2

sampled_experimental_data <- samples %>% 
  pivot_longer(!c(iteration, chain, sample_id), names_to = "cell_id", values_to = "n_per_cell") %>% 
  # pull in mother and daughter cell ids
  left_join(distinct(episome_data, cell_id, mother_cell_id, daughter_cell), by = c("cell_id")) %>% 
  mutate(daughter_cell = paste0("X", daughter_cell)) %>% 
  # make one row per mother cell per sample (iteration:chain)
  select(-cell_id) %>% 
  pivot_wider(names_from = daughter_cell, values_from = n_per_cell, values_fill = 0)

# iterate through the sampled data and apply a grid search to each one
MLE_with_uncertainty <- foreach(i = 1:n_samples, .packages = "tidyverse") %dopar% {
  temp <- sampled_experimental_data %>% 
    filter(sample_id == i) %>% 
    select(id =  mother_cell_id,  X1, X2) %>% 
    #  run a  grid search with a grid ranging by 0.01 and a PMF inferred from the 
    #  Experimental data
    run_grid_search(viz = F, increment = 0.01, 
                    PMF = tibble(X0 = 1:100, prob = dpois(X0, lambda = lambda)), 
                    known_X0 = F, CI = F)
  
  # scale the log  loglikelihood according to the maximum  value and then convert to a probability
  temp <- temp$grid_search %>% mutate(likelihood = exp(log_likelihood - max(log_likelihood)),
                                      probability = likelihood/sum(likelihood))
  temp
}

MLE_with_uncertainty_df <- MLE_with_uncertainty %>% bind_rows() %>% 
  select(Pr, Ps, probability) %>% 
  group_by(Pr, Ps) %>% summarise(probability = sum(probability)) %>% 
  ungroup %>% 
  mutate(probability = probability/sum(probability),
         log_likelihood = log(probability))

MLE_with_uncertainty_CI <- calculate_CI(MLE_with_uncertainty_df)

MLE_uncertainty_grid <- list(grid_search = MLE_with_uncertainty_df, 
                             estimates = MLE_with_uncertainty_CI$estimates,
                             top_95 = MLE_with_uncertainty_CI$probs)

save(MLE_with_uncertainty, MLE_uncertainty_grid, file = here(results_folder, "MLE_with_uncertainty.RData"))

write_csv(MLE_with_uncertainty_df, file = here(results_folder, "grid_search_with_uncertainty.RData"))

ggsave(here(results_folder, "grid_search_with_uncertainty.pdf"), 
       plot_grid_search(MLE_uncertainty_grid, prob = T, simulation = F),
       width = 7, height = 7)



#####
# Try again with different X0 based on my suspician that we are assuuming it is too large

# iterate through the sampled data and apply a grid search to each one
MLE_with_uncertainty_newX0 <- foreach(i = 1:n_samples, .packages = "tidyverse") %dopar% {
  temp <- sampled_experimental_data %>% 
    filter(sample_id == i) %>% 
    select(id =  mother_cell_id,  X1, X2) %>% 
    #  run a  grid search with a grid ranging by 0.01 and a PMF inferred from the 
    #  Experimental data
    run_grid_search(viz = F, increment = 0.01, 
                    PMF = tibble(X0 = 1:100, prob = dpois(X0, lambda = 2.39)), 
                    known_X0 = F, CI = F)
  
  # scale the log  loglikelihood according to the maximum  value and then convert to a probability
  temp <- temp$grid_search %>% mutate(likelihood = exp(log_likelihood - max(log_likelihood)),
                                      probability = likelihood/sum(likelihood))
  temp
}

MLE_with_uncertainty_df_newX0 <- MLE_with_uncertainty_newX0 %>% bind_rows() %>% 
  select(Pr, Ps, probability) %>% 
  group_by(Pr, Ps) %>% summarise(probability = sum(probability)) %>% 
  ungroup %>% 
  mutate(probability = probability/sum(probability),
         log_likelihood = log(probability))

MLE_with_uncertainty_CI_newX0 <- calculate_CI(MLE_with_uncertainty_df_newX0)

MLE_uncertainty_grid_newX0 <- list(grid_search = MLE_with_uncertainty_df_newX0, 
                                   estimates = MLE_with_uncertainty_CI_newX0$estimates,
                                   top_95 = MLE_with_uncertainty_CI_newX0$probs)

plot_grid_search(MLE_uncertainty_grid_newX0, prob = T, simulation = F)

MLE_with_uncertainty_CI_newX0$estimates


mother_cell_data %>% 
  ggplot(aes(total_cluster_intensity, y = after_stat(density))) +
  geom_histogram(aes(fill = "mother intensity"), alpha = 0.5) + 
  geom_histogram(aes(x = 2*total_cluster_intensity, fill = "perfect replication"), alpha = 0.5) + 
  geom_histogram(data = episome_data %>% group_by(mother_cell_id) %>% 
                   mutate(int = sum(total_cluster_intensity)), 
                 aes(int, fill = "total daughter intensity"), alpha = 0.5)


median(all_chains$mu)
median(all_chains_initial$mu)



### FIgure out issue with data

mother_cell_data %>% 
  group_by(cell_id) %>% summarise(total_cluster_intensity = sum(total_cluster_intensity)) %>% 
  ggplot(aes(total_cluster_intensity, y = after_stat(count / sum(count)))) +
  geom_histogram(aes(fill = "mother intensity", color = "mother intensity"), alpha = 0.5) + 
  geom_histogram(aes(x = 2*total_cluster_intensity, fill = "perfect replication\nfrom mother cells",
                     color = "perfect replication\nfrom mother cells"), alpha = 0.5) + 
  geom_histogram(data = episome_data %>% group_by(mother_cell_id) %>% 
                   mutate(int = sum(total_cluster_intensity)), 
                 aes(int, fill = "total daughter intensity", color = "total daughter intensity"), alpha = 0.5) + 
  labs(x = "cluster intensity", fill = "cell set", y= "frequency") + 
  scale_y_continuous(labels = scales::percent) +
  # theme(legend.position = c(1,1), legend.justification = c(1.1,1.1))  +
  guides(color = "none") +
  labs(title = "Distribution of whole cell intensities",
       subtitle = "Live 8TR") +
  
  mother_cell_data %>%
  group_by(cell_id) %>% summarise(total_cluster_intensity = sum(total_cluster_intensity)) %>% 
  ggplot() +
  geom_boxplot(aes(x = total_cluster_intensity, fill = "mother intensity", y = "mother intensity"), alpha = 0.5) + 
  geom_boxplot(aes(x = 2*total_cluster_intensity, fill = "perfect replication from mother cells", y = "perfect replication\nfrom mother cells"), alpha = 0.5) + 
  geom_boxplot(data = episome_data %>% group_by(mother_cell_id) %>% 
                 mutate(int = sum(total_cluster_intensity)), 
               aes(int, fill = "total daughter intensity", y = "total daughter intensity"), alpha = 0.5) + 
  labs(x = "cluster intensity", fill = "cell set") + 
  theme(legend.position = "none", axis.title.y = element_blank()) + 
  
  plot_layout(ncol = 1, heights = c(3,1)) & 
  labs(x = "Intensity in cell(s)")

mother_cell_data %>%
  ggplot(aes(x = total_cluster_intensity)) +
  geom_boxplot(aes(y = "Mother Cells"), outlier.shape = NA) + 
  geom_jitter(aes(y = "Mother Cells")) + 
  geom_boxplot(data = episome_data, aes(y = "Daughter Cells"), outlier.shape = NA) + 
  geom_jitter(data = episome_data, aes( y = "Daughter Cells")) + 
  # coord_flip() + 
  labs(x = "Cluster Intensity", y = "Cell Set", 
       title = "Distribution of cluster intensities",
       subtitle = "Live 8TR") 

daughter_modes <- sampled_daughter_cells %>% 
  filter(chain == "chain1") %>% 
  pivot_longer(starts_with("X")) %>% 
  count(mother_cell_id, name, value) %>% 
  group_by(mother_cell_id, name) %>% filter(n == max(n))  %>% 
  select(-n) %>% 
  pivot_wider() %>% 
  ungroup %>% 
  mutate(temp_X1 = ifelse(X1 > X2, X1, X2),
         temp_X2 = ifelse(X2 < X1, X2, X1),
         total = X1 + X2) %>% 
  select(-X1, -X2, X1 = temp_X1, X2 = temp_X2, total) 


mother_modes <- mother_cell_samples_long %>% filter(chain == "chain 1") %>% 
  count(cell_id, number_of_episomes) %>% 
  group_by(cell_id) %>% filter(n == max(n)) %>% 
  ungroup() %>% select(-n)


mother_modes %>% 
  ggplot(aes(number_of_episomes, y = after_stat(count / sum(count)))) +
  geom_bar(aes(fill = "mother cells", color = "mother cells"), alpha = 0.5, width = 0.9) + 
  geom_bar(aes(x = 2*number_of_episomes, fill = "perfect replication\nfrom mother cells", color = "perfect replication\nfrom mother cells"), alpha = 0.5, width = 0.9) + 
  geom_bar(data = daughter_modes, 
           aes(total, fill = "total daughter cells", color = "total daughter cells"), alpha = 0.5, width = 0.9) + 
  labs(x = "number of episomes", fill = "cell set", y= "frequency") + 
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = c(0:20)) +
  theme(legend.position = c(1,1), legend.justification = c(1.1,1.1)) + 
  guides(color = "none") +
  
  
  mother_modes %>% 
  ggplot(aes(number_of_episomes, y = after_stat(count / sum(count)))) +
  geom_boxplot(aes(fill = "mother cells", y= "mother cells"), alpha = 0.5) + 
  geom_boxplot(aes(x = 2*number_of_episomes, fill = "perfect replication\nfrom mother cells", y ="perfect replication\nfrom mother cells"), alpha = 0.5) + 
  geom_boxplot(data = daughter_modes, 
               aes(total, fill = "total daughter cells", y = "total daughter cells"), alpha = 0.5) + 
  labs(x = "number of episomes", fill = "cell set", y= "frequency") + 
  theme(legend.position = "none") + 
  
  plot_layout(ncol = 1, heights = c(3,1))




