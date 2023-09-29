### Propogate Error in MLE Estimates of real data

#####
# create a cluster
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()

library(tidyverse)
library(patchwork)
library(purrr)
library(foreach)

#####
# read in experimental data
episome_data <- readxl::read_excel("Figure 2 dots new version_MGedit.xlsx", skip = 1,  .name_repair = "none")

colnames <- names(episome_data)
names(episome_data) <- c(colnames[1:2], paste(colnames[3:5], "daughter1", sep = "_"), paste(colnames[6:8], "daughter2", sep = "_"), colnames[9:10])
episome_data <- episome_data[,-9]

episome_data <- rbind(episome_data %>% select( `Mother cell id`, contains("daughter1")) %>% mutate(daughter_cell = 1) %>% 
                        rename_with(~gsub("_daughter1", "", .x)),
                      episome_data %>% select( `Mother cell id`, contains("daughter2")) %>% mutate(daughter_cell =2) %>% 
                        rename_with(~gsub("_daughter2", "", .x))
) %>% 
  select(`Mother cell id`, daughter_cell, Cluster, everything()) %>% 
  arrange(`Mother cell id`, daughter_cell, Cluster) %>% 
  filter(!is.na(`Total cluster intensity`), `Total cluster intensity` != 0) %>% 
  mutate(cluster_id = paste0("n", 1:nrow(.)),
         cell_id = paste0("c", `Mother cell id`, "_", daughter_cell))

# load number per cell sampled at each iteration
load("gibbs_cell_samples.RData")

#####
# randomly choose 100 samples to use as inferred value of number of episomes per cell:

set.seed(100)
samples <- sample_n(cell_samples, 100) %>% mutate(sample_id = 1:100)

sampled_n <- episome_data %>% 
  distinct(`Mother cell id`, daughter_cell, cell_id) %>% 
  merge(pivot_longer(samples, !c(iteration, chain, sample_id), names_to = "cell_id")) %>% 
  select(-cell_id) %>% mutate(daughter_cell = paste0("X", daughter_cell)) %>% 
  pivot_wider(names_from = daughter_cell, values_from = value, values_fill = 0) %>% 
  rename(id = `Mother cell id`)

#####
# for each sampled combination, run a grid search and then sum all probabilities at each Pr and Ps combination

temp <- sampled_n %>% filter(sample_id <= 2) %>% 
  mutate(X0 = 4) %>% 
  nest(data = c(id, X1, X2, X0)) %>% 
  mutate(l = map(data, run_grid_search, viz = F)) %>% 
  pull(l)


temp <- sampled_n %>% filter(sample_id == 1) %>% 
  # mutate(X0 = 4) %>% 
  select(id,  X1, X2) %>% 
  run_grid_search(viz = F, increment = 0.01, PMF = tibble(X0 = 1:100, prob = dnbinom(X0, 100, mu = 4.275)), known_X0 = F, CI = F)
temp <- temp$grid_search %>% mutate(likelihood = exp(log_likelihood - max(log_likelihood)),
                                    probability = likelihood/sum(likelihood))

collect <- as.matrix(temp %>% select(Pr, Ps, likelihood))

MLE_with_uncertainty <- foreach(i = 1:100, .packages = "tidyverse") %dopar% {
  temp <- sampled_n %>% filter(sample_id == i) %>% 
    # mutate(X0 = 4) %>% 
    select(id,  X1, X2) %>% 
    run_grid_search(viz = F, increment = 0.01, PMF = tibble(X0 = 1:100, prob = dnbinom(X0, 100, mu = 4.275)), known_X0 = F, CI = F)
  temp <- temp$grid_search %>% mutate(likelihood = exp(log_likelihood - max(log_likelihood)),
                                      probability = likelihood/sum(likelihood))
  # collect[,3] <- collect[,3] + temp$likelihood
  temp
}

save(MLE_with_uncertainty, file = "MLE_with_uncertainty.RData")

MLE_with_uncertainty_df <- MLE_with_uncertainty %>% bind_rows() %>% 
  select(Pr, Ps, probability) %>% 
  group_by(Pr, Ps) %>% summarise(probability = sum(probability)) %>% 
  ungroup %>% 
  mutate(probability = probability/sum(probability),
         log_likelihood = log(probability))

MLE_with_uncertainty_CI <- calculate_CI(MLE_with_uncertainty_df)

MLE_with_uncertainty_df %>% ggplot(aes(Pr, Ps, fill = probability)) + geom_tile()

MLE_uncertainty_grid <- list(grid_search = MLE_with_uncertainty_df, 
                             estimates = MLE_with_uncertainty_CI$estimates,
                             top_95 = MLE_with_uncertainty_CI$probs)
plot_grid_search(MLE_uncertainty_grid, prob = T, simulation = F)
 