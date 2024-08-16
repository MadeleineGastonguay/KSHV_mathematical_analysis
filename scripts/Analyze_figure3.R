# This script runs the full pipeline for data from fixed images of cells with only KHSV terminal repeats (Figure 2)

#####
#  Setup 
#####

# load libraries 
library(tidyverse)
library(here)

theme_set(theme_bw())

# Check wd
here()

# Read in helper functions 
source(here("scripts", "likelihood_functions.R"))
source(here("scripts", "run_pipeline.R"))


#####
# Script inputs

# Folder for results
results_folder <- here("results", "Figure_3_updated_pdf_n_prior")

# Daughter cell intensity data
daughter_cell_file <- "Fig3 dividing cells.xlsx"

# Mother cell intensity data
mother_cell_file <- "Fig3 Non dividing cells.xlsx"
#####

Figure3_data <- load_data(mother_cell_file, daughter_cell_file)

daughter_cell_data <- Figure3_data$daughter_cell_data
mother_cell_data <- Figure3_data$mother_cell_data

Figure3_results <- run_pipeline(daughter_cell_data, mother_cell_data, results_folder,
                                n_prior = list("geom", 0.5), parallel = T, just_Pr = T)

Figure3_results$MLE_grid$grid_search %>% group_by(Pr) %>% summarise(probability = sum(probability)) %>% 
  ggplot(aes(Pr, probability)) + geom_line(aes(color = "joint estimation")) +
  geom_line(data = Figure3_results$MLE_Pr_grid$grid_search, aes(color = "single estimation"))

Figure3_results$MLE_grid$estimates
Figure3_results$MLE_Pr_grid$estimates

make_plots(Figure3_results, daughter_cell_data, mother_cell_data, results_folder)

daughter_cell_samples <- Figure3_results$daughter_cell_samples
figures <- figures(daughter_cell_data, mother_cell_data, daughter_cell_samples, results_folder)

## How many cases did we infer fewer episomes than observed visually?
Figure3_results$all_chains %>% 
  pivot_longer(!c(chain, iteration, mu, tau), names_to = "cluster_id", values_to = "n_epi") %>% 
  count(cluster_id, n_epi) %>% 
  group_by(cluster_id) %>% 
  filter(n == max(n)) %>% 
  ungroup() %>% 
  left_join(rbind(select(Figure3_data$daughter_cell_data, cluster_id, min_episome_in_cluster),
                  select(Figure3_data$mother_cell_data, cluster_id, min_episome_in_cluster))
  ) %>% 
  count(n_epi < min_episome_in_cluster)

## Adjust confidence intervals to be marginal:

Figure3_results$MLE_grid$grid_search %>%
  group_by(Pr) %>%
  summarise(probability = sum(probability)) %>% 
  arrange(desc(probability)) %>% 
  mutate(cum_sum = cumsum(probability)) %>% 
  filter(cum_sum <= 0.95) %>% 
  pull(Pr) %>% range

Figure3_results$MLE_grid$grid_search %>%
  group_by(Ps) %>%
  summarise(probability = sum(probability)) %>% 
  arrange(desc(probability)) %>% 
  mutate(cum_sum = cumsum(probability)) %>%
  filter(cum_sum <= 0.95) %>% 
  pull(Ps) %>% range
