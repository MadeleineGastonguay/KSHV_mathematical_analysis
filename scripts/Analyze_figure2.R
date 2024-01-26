# This script runs the full pipeline for data from fixed images of cells with only KHSV terminal repeats (Figure 2)

#####
#  Setup 
#####

# load libraries 
library(tidyverse)
library(here)
library(patchwork)

theme_set(theme_bw())

# Check wd
here()

# Read in helper functions 
source(here("scripts", "likelihood_functions.R"))
source(here("scripts", "run_pipeline.R"))


#####
# Script inputs

# Folder for results
results_folder <- here("results", "Figure_2_updated_pdf_n_prior")

# Daughter cell intensity data
daughter_cell_file <- "Figure 2 dots new version.xlsx"

# Mother cell intensity data
mother_cell_file <- "Figure 2 non dividing cells 2.xlsx"
#####

Figure2_data <- load_data(mother_cell_file, daughter_cell_file)

daughter_cell_data <- Figure2_data$daughter_cell_data
mother_cell_data <- Figure2_data$mother_cell_data

cat(length(unique(daughter_cell_data$mother_cell_id)), "daughter cell pairs with", nrow(daughter_cell_data), "clusters\n")
cat(length(unique(mother_cell_data$cell_id)), "non-dividing cells with", nrow(mother_cell_data), "clusters")

Figure2_results <- run_pipeline(daughter_cell_data, mother_cell_data, results_folder, 
                                n_prior = list("geom", 0.5), parallel = T, just_Pr = F)

Figure2_results$MLE_grid$estimates

make_plots(Figure2_results, daughter_cell_data, mother_cell_data, results_folder)

daughter_cell_samples <- Figure2_results$daughter_cell_samples
figures <- figures(daughter_cell_data, mother_cell_data, daughter_cell_samples, results_folder)

## How many cases did we infer fewer episomes than observed visually?
Figure2_results$all_chains %>% 
  pivot_longer(!c(chain, iteration, mu, tau), names_to = "cluster_id", values_to = "n_epi") %>% 
  count(cluster_id, n_epi) %>% 
  group_by(cluster_id) %>% 
  filter(n == max(n)) %>% 
  ungroup() %>% 
  left_join(rbind(select(Figure2_data$daughter_cell_data, cluster_id, min_episome_in_cluster),
                  select(Figure2_data$mother_cell_data, cluster_id, min_episome_in_cluster))
  ) %>% 
  count(n_epi < min_episome_in_cluster)

## Show inference results for example images:
# 30, 31, 3, 4

example_inference <- Figure2_results$daughter_cell_samples %>% 
  filter(chain == "chain1") %>% 
  select(starts_with("30") | starts_with("31") | starts_with("3_") | starts_with("4_")) %>% 
  pivot_longer(everything(), names_to = "cell_id", values_to = "n_epi") %>% 
  count(cell_id, n_epi) %>% 
  group_by(cell_id) %>% 
  mutate(prob = n/sum(n),
         cell_id = factor(cell_id, levels = c("30_2", "30_1", "31_1", "31_2", "3_2", "3_1", "4_1"))) %>%
  ggplot(aes(n_epi, prob)) + 
  geom_col() + 
  facet_wrap(~cell_id, nrow = 1) + 
  labs(x = "Estimated number of episomes in cell", y = "Probability") + 
  theme_classic() +
  scale_x_continuous(breaks = seq(0,8, by = 2), limits = c(0,8)) +
  # ylim(c(0,1)) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1), breaks = c(0, 1))  +
  theme(strip.background = element_blank(), strip.text = element_blank()) 

ggsave(here(results_folder, "example_inference.png"), width = 5, height = 1.2)


## Adjust confidence intervals to be marginal:

Figure2_results$MLE_grid$grid_search %>%
  group_by(Pr) %>%
  summarise(probability = sum(probability)) %>% 
  arrange(desc(probability)) %>% 
  mutate(cum_sum = cumsum(probability)) %>% 
  filter(cum_sum <= 0.95) %>% 
  pull(Pr) %>% range

Figure2_results$MLE_grid$grid_search %>%
  group_by(Ps) %>%
  summarise(probability = sum(probability)) %>% 
  arrange(desc(probability)) %>% 
  mutate(cum_sum = cumsum(probability)) %>%
  filter(cum_sum <= 0.95) %>% 
  pull(Ps) %>% range

## Look at the intensity distribution for dots inferred to have a single episome:
single_epi_dist <- function(results, data, set = "all", title = ""){
  temp_df <- results$all_chains %>% 
    filter(chain == "chain1") %>% 
    pivot_longer(grep("[0-9]", names(results$all_chains)), names_to = "cluster", values_to = "n_epi") %>% 
    count(cluster, n_epi) %>% 
    group_by(cluster) %>% 
    mutate(p = n/sum(n)) 
  
  single_episome_clusters <- temp_df %>% filter(n == max(n), n_epi == 1) %>% pull(cluster)
  
  print(temp_df %>% 
          filter(cluster %in% single_episome_clusters) %>% 
          ggplot(aes(as.factor(n_epi), p)) + geom_boxplot(outlier.shape = NA) + geom_point() + 
          geom_line(aes(group = cluster), alpha = 0.1))
  
  if(set == "daughter"){
    results$all_chains <- results$all_chains %>% rename(mu = mu_d, tau = tau_d)
    data_subset <- data$daughter_cell_data %>% select(cluster_id, total_cluster_intensity) %>% 
      filter(cluster_id %in% single_episome_clusters) 
  }else if(set == "mother"){
    results$all_chains <- results$all_chains %>% rename(mu = mu_m, tau = tau_m)
    data_subset <- data$mother_cell_data %>% select(cluster_id, total_cluster_intensity) %>% 
      filter(cluster_id %in% single_episome_clusters) 
  }else{
    data_subset <- rbind(
      data$mother_cell_data %>% select(cluster_id, total_cluster_intensity),
      data$daughter_cell_data %>% select(cluster_id, total_cluster_intensity)
    ) %>% 
      filter(cluster_id %in% single_episome_clusters) 
  }
  
  
  qqplot <- data_subset %>% ggplot(aes(sample= total_cluster_intensity)) + stat_qq() + stat_qq_line() + 
    labs(title = "qqplot", x = "theoretical quantiles", y ="sample quantiles")
  
  mean <- results$all_chains %>% pull(mu) %>% median
  sd <- results$all_chains %>% mutate(sd = sqrt(1/tau)) %>% pull(sd) %>% median  
  
  plot <- data_subset %>% 
    ggplot(aes(total_cluster_intensity)) +
    geom_density(aes(lty = "data")) +
    stat_function(fun = dnorm , args = list(mean = mean, sd = sd), aes(lty = "fit normal distribution")) + 
    labs(y = "density", title = "density plot") + 
    theme(legend.title = element_blank(), legend.position = c(1,0), legend.justification = c(1.1,-0.1))
  
  title <- ifelse(set == "all", title, paste(title, set, "cells"))
  print(qqplot + plot + plot_annotation(title = title))
  
}

single_epi_dist(Figure2_results, Figure2_data, "Fixed 8TR")
single_epi_dist(Figure3_results, Figure3_data, "Live 8TR")
single_epi_dist(Figure5_results, Figure5_data, "Fixed KSHV")
single_epi_dist(Figure6_results, Figure6_data, "daughter", "Live KSHV")
single_epi_dist(Figure6_results, Figure6_data, "mother", "Live KSHV")

