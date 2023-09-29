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
results_folder <- here("results", "Figure_2_subset")
dir.create(results_folder)

#####
# Read in data
#####

daughter_cell_file <- "Figure 2 dots new version.xlsx"
mother_cell_file <- "Figure 2 non dividing cells 2.xlsx"

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
# Subset the daughter cell data so it's the same size as live cell images
#####
set.seed(22)
mother_cells <- sample(unique(episome_data$mother_cell_id), 10)
episome_data <- episome_data %>% filter(mother_cell_id %in% mother_cells) %>% sample_n(27)

# the data are more evenly distributed than the live cell images
episome_data %>% ggplot(aes(x = "", y=total_cluster_intensity)) + geom_boxplot(outlier.shape = NA) + geom_jitter()

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

names(all_chains)[-c(1:3, 31)] <- episome_data$cluster_id

#####
# Analyze results of Gibbs
#####

inferred_mu <- median(all_chains$mu)
inferred_sigma2 <- median(1/all_chains$tau)
inferred_sigma <- median(sqrt(1/all_chains$tau))
cat('mean:', inferred_mu, '\nvariance:', inferred_sigma2, '\nstandard deviation:', inferred_sigma)


p1 <- all_chains %>% 
  ggplot(aes(iteration, mu, color = chain)) + 
  geom_line() +
  geom_point(shape = NA) + 
  geom_hline(yintercept = inferred_mu, lty = "dashed") + 
  geom_text(data = data.frame(NA), aes(5000, inferred_mu, label = str_interp("${round(inferred_mu,2)}")),
            inherit.aes = F, vjust = -0.2, hjust = -0.2) + 
  theme(legend.position = "bottom") + 
  ylim(c(0, max(all_chains$mu))) +
  labs(title = "Trace of the mean (\u00b5)", 
       y = bquote("mean of the distribution of\nfluorescence intensity for one episome (\u00b5)"))

(p1 <- ggMarginal(p1, margins = "y", groupColour = T))

p2 <- all_chains %>% 
  ggplot(aes(iteration, 1/tau, color = chain)) + 
  geom_line() +
  geom_point(shape = NA) +
  theme(legend.position = "bottom") +
  ylim(c(0, max(1/all_chains$tau))) + 
  labs(title = "Trace of the variance (\u03C3\u00b2)", 
       y = bquote("variance of the distribution of\nfluorescence intensity for one episome (\u03C3\u00b2)"))

(p2 <- ggMarginal(p2, margins = "y", groupColour = T))

p2.2 <- all_chains %>% 
  ggplot(aes(iteration, sqrt(1/tau), color = chain)) + 
  geom_line() +
  geom_point(shape = NA) +
  geom_hline(yintercept = inferred_sigma, lty = "dashed") + 
  geom_text(data = data.frame(NA), aes(5000, inferred_sigma, label = str_interp("${round(inferred_sigma,2)}")),
            inherit.aes = F, vjust = -0.2, hjust = -0.2) + 
  theme(legend.position = "bottom") +
  ylim(c(0, max(sqrt(1/all_chains$tau)))) + 
  labs(title = "Trace of the standard deviation (\u03C3)", 
       y = bquote("standard deviation of the distribution of\nfluorescence intensity for one episome (\u03C3)"))

(p2.2 <- ggMarginal(p2.2, margins = "y", groupColour = T))


# Based on the median mu and sigma above, the inferred distribution of intensities for a single cluster is:
# plot the distribution of cluster intensities with inferred parameters:
x <- seq(0, max(episome_data$total_cluster_intensity), length.out = 500)
(p3 <- ggplot(tibble(x, y = dnorm(x, inferred_mu, inferred_sigma)), aes(x, y)) + 
    geom_line() + 
    labs(x = "Fluorescence Intensity", 
         y = "probability density", 
         title = "Inferred distribution of fluorescence intensity\nof a single cluster",
         subtitle = str_interp("mean = ${round(inferred_mu,2)}, standard deviation = ${round(inferred_sigma, 2)}")))


## infer nk for each cluster according to the mode of the posterior distribution
all_chain_ns <- all_chains %>% select(-mu, -tau)  %>% pivot_longer(!c(iteration, chain), names_to = "cluster")
modes <- all_chain_ns %>% add_count(chain, cluster, value) %>% group_by(chain, cluster) %>% mutate(mode = value[which.max(n)]) %>% ungroup()

modes %>% filter(chain == "chain1") %>% 
  ggplot(aes(value, group = cluster)) + geom_density() + 
  facet_wrap(~mode, labeller = "label_both", scales = 'free_y') + 
  labs(x = "number of episomes in cluster",
       y = "density",
       title = "Posterior distribution of the number of episomes per cluster")


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
  filter(chain == "chain1") %>% 
  ungroup %>% 
  ggplot(aes(total_cluster_intensity, color = as.factor(mode))) + 
  geom_density() + 
  labs(title = "Distribution of cluster intensities for clusters with each nk", color  = "mode")

p6 <- episome_data %>% 
  mutate(cluster = cluster_id) %>% 
  merge(distinct_modes) %>% 
  filter(chain == "chain1") %>% 
  ggplot(aes(total_cluster_intensity, as.factor(min_episome_in_cluster))) + 
  geom_jitter(aes(color = as.factor(mode))) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  labs(title = "Distribution of cluster intensities stratified by minimum number of episomes",
       subtitle = "colored by inferred nk",
       color = "mode", y = "Min # episome in cluster", x = "Total  Cluster Intensity")

p6.2 <- episome_data %>% 
  mutate(cluster = cluster_id) %>% 
  merge(distinct_modes) %>% 
  filter(chain == "chain1") %>% 
  ggplot(aes(total_cluster_intensity, as.factor(mode))) + 
  geom_jitter(aes(color = as.factor(min_episome_in_cluster))) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  labs(title = "Distribution of cluster intensities",
       y = "Mode of the posterior distribution", 
       color = "Min # episome in cluster",
       x = "Total  Cluster Intensity") + 
  theme(legend.position = c(1,0), legend.justification = c(1,0), legend.background = element_blank())


merged_df <- all_chains %>% filter(chain == "chain1") %>% 
  pivot_longer(starts_with("n"), names_to = "cluster_id", values_to = "n") %>% 
  left_join(episome_data %>% select(cluster_id, total_cluster_intensity, min_episome_in_cluster)) 



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
# Make sure inferred distribution of initial episomes makes sense with the inferred episomes per daughter cell
#####

lambda <- 4.4176 

sampled_daughter_cells <- cell_samples_long %>% filter(chain == "chain1") %>% 
  left_join(distinct(episome_data, cell_id, mother_cell_id, daughter_cell), by = c("cell_id")) %>% 
  mutate(daughter_cell = paste0("X", daughter_cell)) %>% 
  # make one row per mother cell per sample (iteration:chain)
  select(-cell_id) %>% 
  pivot_wider(names_from = daughter_cell, values_from = number_of_episomes, values_fill = 0) 

initial_episome_constraints <- sampled_daughter_cells %>% mutate(max = X1 + X2, min = ceiling(max/2)) %>% 
  pivot_longer(c(min, max)) %>% 
  count(name, value) %>% 
  group_by(name) %>% mutate(freq = n/sum(n)) %>% 
  ggplot(aes(value, freq, color = name))  + 
  geom_bar(stat = "identity", fill = NA) + 
  geom_point(data = tibble(value = 1:15, freq = dpois(value, lambda), name = "assumed distribution")) +
  labs(x = "Initial number of episomes", y= "Frequency", color = "contraint")

compare_daughter_cell_posteriors <- sampled_daughter_cells %>% 
  count(X1, X2) %>% mutate(p = n/sum(n)) %>% 
  ggplot(aes(X1, X2, fill = p)) + 
  geom_tile() + 
  scale_x_continuous(breaks = 0:max(sampled_daughter_cells$X1)) + 
  scale_y_continuous(breaks = 0:max(sampled_daughter_cells$X2)) + 
  scale_fill_viridis_c() + 
  labs(title = "Fraction of sampled X1 and X2 combinations across all cells",
       fill = "fraction",
       x = "number of epsiomes in daughter cell 1 (X1)",
       y = "number of epsiomes in daughter cell 2 (X2)") 



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
# Compare Gibbs results to full dataset:
#####

full_data <- new.env()
load(here("results", "Figure_2", "daughter_cell_cluster_Gibbs_samples.RData"), envir = full_data)


all_chains_full <- rbind(full_data$chain1, full_data$chain2, full_data$chain3, full_data$chain4, full_data$chain5) %>%
  mutate(chain = as.factor(rep(c("chain1", "chain2", "chain3", "chain4", "chain5"), each = n_iterations)))  %>%
  filter(iteration > burn_in)


compare_runs <- rbind(
  all_chains_full %>% select(all_of(names(all_chains))) %>%
    filter(chain == "chain1") %>% 
    pivot_longer(starts_with("n"), names_to = "cluster", values_to = "nk") %>% 
    mutate(run = "full_data"), 
  
  all_chains %>% 
    filter(chain == "chain1") %>% 
    pivot_longer(starts_with("n"), names_to = "cluster", values_to = "nk") %>% 
    mutate(run = "subset")
)


run_modes <- compare_runs %>% 
  count(cluster, run, nk) %>% 
  group_by(cluster, run) %>% filter(n == max(n)) %>% 
  rename(mode = nk) %>% 
  select(-n) %>% 
  ungroup

p <- compare_runs %>% 
  # filter(cluster %in% c("n57", "n82")) %>%
  count(cluster, run, nk) %>% 
  mutate(frac = n/sum(n)) %>% 
  ggplot(aes(x = nk, y = run, fill = frac)) + 
  geom_tile() + 
  geom_text(data = run_modes , aes(Inf, run, label = paste("mode:", mode)),
            hjust = 1.1, inherit.aes = F) +
  facet_wrap(~cluster) 

ggsave(here(results_folder, "compare_full_to_subset_Gibbs.png"), p, width = 15, height = 10)

run_modes %>% 
  pivot_wider(names_from = run, values_from = mode) %>% 
  count(full_data, subset) %>% 
  ggplot(aes(full_data, subset)) + 
  geom_tile(aes(fill=  n))

# The mode number of episomes per cluster agrees between the two methods for 55% of clusters:
run_modes %>% 
  pivot_wider(names_from = run, values_from = mode) %>% 
  count(full_data == subset)


#####
# What if we run MLE with samples from Gibbs with full dataset?
#####


cell_samples_long_test <- as.list(unique(episome_data$cell_id)) %>% 
  lapply(get_n_cell, all_chains_full, episome_data) %>% bind_rows() 
cell_samples_test <- cell_samples_long_test %>% pivot_wider(names_from = cell_id, values_from = number_of_episomes)


samples_test <- sample_n(cell_samples_test, n_samples) %>% mutate(sample_id = 1:n_samples)

# reformat samples so that we can apply the grid search function: 
# need a column for mother cell id, number of cells in daughter cell 1, and number of cells in daughter cecll 2

sampled_experimental_data_test <- samples_test %>% 
  pivot_longer(!c(iteration, chain, sample_id), names_to = "cell_id", values_to = "n_per_cell") %>% 
  # pull in mother and daughter cell ids
  left_join(distinct(episome_data, cell_id, mother_cell_id, daughter_cell), by = c("cell_id")) %>% 
  mutate(daughter_cell = paste0("X", daughter_cell)) %>% 
  # make one row per mother cell per sample (iteration:chain)
  select(-cell_id) %>% 
  pivot_wider(names_from = daughter_cell, values_from = n_per_cell, values_fill = 0)

# iterate through the sampled data and apply a grid search to each one
MLE_with_uncertainty_test <- foreach(i = 1:n_samples, .packages = "tidyverse") %dopar% {
  temp <- sampled_experimental_data_test %>% 
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

MLE_with_uncertainty_df_test <- MLE_with_uncertainty_test %>% bind_rows() %>% 
  select(Pr, Ps, probability) %>% 
  group_by(Pr, Ps) %>% summarise(probability = sum(probability)) %>% 
  ungroup %>% 
  mutate(probability = probability/sum(probability),
         log_likelihood = log(probability))

MLE_with_uncertainty_CI_test <- calculate_CI(MLE_with_uncertainty_df_test)

MLE_uncertainty_grid_test <- list(grid_search = MLE_with_uncertainty_df_test, 
                             estimates = MLE_with_uncertainty_CI_test$estimates,
                             top_95 = MLE_with_uncertainty_CI_test$probs)

plot_grid_search(MLE_uncertainty_grid_test, prob = T, simulation = F)


sampled_experimental_data_test %>% 
  mutate(max = X1 + X2, min = ceiling(max/2)) %>% 
  pivot_longer(c(min, max)) %>% 
  count(name, value) %>% 
  group_by(name) %>% 
  mutate(freq= n/sum(n)) %>% 
  ungroup %>% 
  ggplot() + 
  geom_histogram(aes(freq, color = name)) + 
  geom_point(data = tibbsle(X0 = 1:13, prob = dpois(X0, lambda = lambda)),
             aes(X0, prob, color  = "assumed_probability"))


### Try with exact subset (not losing clustersr)
load(here("results", "Figure_2", "daughter_cell_Gibbs_samples.RData"), envir = full_data)
cell_samples_test2 <- full_data$cell_samples %>% select(iteration, chain, all_of(episome_data$cell_id))
samples_test2 <- sample_n(cell_samples_test2, n_samples) %>% mutate(sample_id = 1:n_samples)

# reformat samples so that we can apply the grid search function: 
# need a column for mother cell id, number of cells in daughter cell 1, and number of cells in daughter cecll 2

sampled_experimental_data_test2 <- samples_test2 %>% 
  pivot_longer(!c(iteration, chain, sample_id), names_to = "cell_id", values_to = "n_per_cell") %>% 
  # pull in mother and daughter cell ids
  left_join(distinct(episome_data, cell_id, mother_cell_id, daughter_cell), by = c("cell_id")) %>% 
  mutate(daughter_cell = paste0("X", daughter_cell)) %>% 
  # make one row per mother cell per sample (iteration:chain)
  select(-cell_id) %>% 
  pivot_wider(names_from = daughter_cell, values_from = n_per_cell, values_fill = 0)

# iterate through the sampled data and apply a grid search to each one
MLE_with_uncertainty_test2 <- foreach(i = 1:n_samples, .packages = "tidyverse") %dopar% {
  temp <- sampled_experimental_data_test2 %>% 
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

MLE_with_uncertainty_df_test2 <- MLE_with_uncertainty_test2 %>% bind_rows() %>% 
  select(Pr, Ps, probability) %>% 
  group_by(Pr, Ps) %>% summarise(probability = sum(probability)) %>% 
  ungroup %>% 
  mutate(probability = probability/sum(probability),
         log_likelihood = log(probability))

MLE_with_uncertainty_CI_test2 <- calculate_CI(MLE_with_uncertainty_df_test2)

MLE_uncertainty_grid_test2 <- list(grid_search = MLE_with_uncertainty_df_test2, 
                                   estimates = MLE_with_uncertainty_CI_test2$estimates,
                                   top_95 = MLE_with_uncertainty_CI_test2$probs)

plot_grid_search(MLE_uncertainty_grid_test2, prob = T, simulation = F)


sampled_experimental_data_test2 %>% 
  mutate(max = X1 + X2, min = ceiling(max/2)) %>% 
  pivot_longer(c(min, max)) %>% 
  count(name, value) %>% 
  group_by(name) %>% 
  mutate(freq= n/sum(n)) %>% 
  ungroup %>% 
  ggplot() + 
  geom_histogram(aes(freq, color = name)) + 
  geom_point(data = tibbsle(X0 = 1:13, prob = dpois(X0, lambda = lambda)),
             aes(X0, prob, color  = "assumed_probability"))
