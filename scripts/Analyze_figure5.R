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
results_folder <- here("results", "Figure_5_updated_pdf_n_prior")

# Daughter cell intensity data
daughter_cell_file <- "Fig5 dividing cells.xlsx"

# Mother cell intensity data
mother_cell_file <- "Fig5 Non dividing cells.xlsx"
#####

Figure5_data <- load_data(mother_cell_file, daughter_cell_file)

daughter_cell_data <- Figure5_data$daughter_cell_data
mother_cell_data <- Figure5_data$mother_cell_data

Figure5_results <- run_pipeline(daughter_cell_data, mother_cell_data, results_folder,
                                n_prior = list("geom", 0.5), parallel = T, just_Pr = F)

# Figure5_results$MLE_grid$grid_search %>% group_by(Pr) %>% summarise(probability = sum(probability)) %>% 
#   ggplot(aes(Pr, probability)) + geom_line(aes(color = "joint estimation")) +
#   geom_line(data = Figure5_results$MLE_Pr_grid$grid_search, aes(color = "single estimation"))

Figure5_results$MLE_grid$estimates


## Adjust confidence intervals to be marginal:

Figure5_results$MLE_grid$grid_search %>%
  group_by(Pr) %>%
  summarise(probability = sum(probability)) %>% 
  arrange(desc(probability)) %>% 
  mutate(cum_sum = cumsum(probability)) %>% 
  filter(cum_sum <= 0.95) %>% 
  pull(Pr) %>% range

Figure5_results$MLE_grid$grid_search %>%
  group_by(Ps) %>%
  summarise(probability = sum(probability)) %>% 
  arrange(desc(probability)) %>% 
  mutate(cum_sum = cumsum(probability)) %>%
  filter(cum_sum <= 0.95) %>% 
  pull(Ps) %>% range

make_plots(Figure5_results, daughter_cell_data, mother_cell_data, results_folder)

daughter_cell_samples <- Figure5_results$daughter_cell_samples
figures <- figures(daughter_cell_data, mother_cell_data, daughter_cell_samples, results_folder)

## How many cases did we infer fewer episomes than observed visually?
Figure5_results$all_chains %>% 
  pivot_longer(!c(chain, iteration, mu, tau), names_to = "cluster_id", values_to = "n_epi") %>% 
  count(cluster_id, n_epi) %>% 
  group_by(cluster_id) %>% 
  filter(n == max(n)) %>% 
  ungroup() %>% 
  left_join(rbind(select(Figure5_data$daughter_cell_data, cluster_id, min_episome_in_cluster),
                  select(Figure5_data$mother_cell_data, cluster_id, min_episome_in_cluster))
            ) %>% 
  count(n_epi < min_episome_in_cluster)



## Show inference results for example images:
# 12, 9, #10

example_inference <- Figure5_results$daughter_cell_samples %>% 
  filter(chain == "chain1") %>% 
  select(starts_with("Image 12_") | starts_with("Image 9_") | starts_with("Image #10_") ) %>%  
  pivot_longer(everything(), names_to = "cell_id", values_to = "n_epi") %>% 
  rbind(data.frame(cell_id = "Image 9_2", n_epi = 0)) %>% 
  count(cell_id, n_epi) %>% 
  group_by(cell_id) %>% 
  mutate(prob = n/sum(n),
         cell_id = factor(cell_id, levels = c("Image 12_1", "Image 12_2", "Image 9_2", "Image 9_1", "Image #10_1", "Image #10_2"))) %>%
  ggplot(aes(n_epi, prob)) + 
  geom_col() + 
  facet_wrap(~cell_id, nrow = 1) + 
  labs(x = "Estimated number of episomes in cell", y = "Probability") + 
  theme_classic() +
  scale_x_continuous(breaks = seq(0,10, by = 2), limits = c(0,6)) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1), breaks = c(0, 1))  +
  theme(strip.background = element_blank(), strip.text = element_blank()) 

ggsave(here(results_folder, "example_inference.png"), width = 5, height = 1.2)


Figure5_results$daughter_cell_samples %>% 
  filter(chain == "chain1") %>% 
  select(starts_with("Image 12_") ) %>%  
  pivot_longer(everything(), names_to = "cell_id", values_to = "n_epi") %>% 
  # rbind(data.frame(cell_id = "Image 9_2", n_epi = 0)) %>% 
  count(cell_id, n_epi) %>% 
  group_by(cell_id) %>% 
  mutate(prob = n/sum(n),
         cell_id = factor(cell_id, levels = c("Image 12_1", "Image 12_2", "Image 9_2", "Image 9_1", "Image #10_1", "Image #10_2"))) %>%
  ggplot(aes(n_epi, prob)) + 
  geom_col() + 
  facet_wrap(~cell_id, nrow = 1) + 
  labs(x = "Estimated number of episomes in cell", y = "Probability") + 
  theme_classic() +
  scale_x_continuous(breaks = seq(0,10, by = 2), limits = c(0,4)) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1), breaks = c(0, 1))  +
  theme(strip.background = element_blank(), strip.text = element_blank()) 

Figure5_results$daughter_cell_samples %>% 
  filter(chain == "chain1") %>% 
  select(starts_with("Image 9_") ) %>%  
  pivot_longer(everything(), names_to = "cell_id", values_to = "n_epi") %>% 
  rbind(data.frame(cell_id = "Image 9_2", n_epi = 0)) %>%
  count(cell_id, n_epi) %>% 
  group_by(cell_id) %>% 
  mutate(prob = n/sum(n),
         cell_id = factor(cell_id, levels = c("Image 12_1", "Image 12_2", "Image 9_2", "Image 9_1", "Image #10_1", "Image #10_2"))) %>%
  ggplot(aes(n_epi, prob)) + 
  geom_col() + 
  facet_wrap(~cell_id, nrow = 1) + 
  labs(x = "Estimated number of episomes in cell", y = "Probability") + 
  theme_classic() +
  scale_x_continuous(breaks = seq(0,10, by = 2), limits = c(0,4)) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1), breaks = c(0, 1))  +
  theme(strip.background = element_blank(), strip.text = element_blank()) 



## Mean values of intensity per cell and number of episomes per cell
Figure5_results$daughter_cell_samples %>% 
  filter(chain == "chain1") %>% 
  select(-c(chain, iteration)) %>% 
  as.matrix %>% 
  mean

Figure5_results$mother_cell_samples %>% 
  filter(chain == "chain1") %>% 
  select(-c(chain, iteration)) %>% 
  as.matrix %>% 
  mean

Figure5_data$daughter_cell_data %>% 
  group_by(cell_id) %>% 
  summarise(intensity = sum(total_cluster_intensity)) %>% 
  pull(intensity) %>% mean

Figure5_data$mother_cell_data %>% 
  group_by(cell_id) %>% 
  summarise(intensity = sum(total_cluster_intensity)) %>% 
  pull(intensity) %>% mean

## Get standard deviation of inferred mean and sigma:

Figure5_results$all_chains %>% 
  filter(chain == "chain1") %>% 
  select(mu, tau) %>% 
  mutate(sd = sqrt(1/tau)) %>% 
  apply(2, sd)
