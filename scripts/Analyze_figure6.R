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
results_folder <- here("results", "Figure_6_different_mean")

# Daughter cell intensity data
daughter_cell_file <- "Fig6 dividing cells.xlsx"

# Mother cell intensity data
mother_cell_file <- "Fig6 Non dividing cells.xlsx"
#####

Figure6_data <- load_data(mother_cell_file, daughter_cell_file)

daughter_cell_data <- Figure6_data$daughter_cell_data
mother_cell_data <- Figure6_data$mother_cell_data

Figure6_results <- run_pipeline(daughter_cell_data, mother_cell_data, results_folder, same_mu = F)

Figure6_results$MLE_grid$estimates

make_plots(Figure6_results, daughter_cell_data, mother_cell_data, results_folder)

daughter_cell_samples <- Figure6_results$daughter_cell_samples
figures <- figures(daughter_cell_data, mother_cell_data, daughter_cell_samples, results_folder)

all_chains <- Figure6_results$all_chains

convergence <- Figure6_results$convergence_metrics

convergence %>% 
  filter(!name %in% c("mu", "tau")) %>% 
  mutate(set = ifelse(grepl("d", name), "daughter", "mother")) %>% 
  pivot_longer(contains("ESS"), names_to = "metric") %>% 
  ggplot(aes(value, y = metric, color = set)) + 
  geom_boxplot() 


convergence %>% 
  filter(!name %in% c("mu", "tau")) %>% 
  mutate(set = ifelse(grepl("d", name), "daughter", "mother")) %>% 
  pivot_longer(contains("ESS"), names_to = "metric") %>% 
  filter(value < 50000) %>% 
  ggplot(aes(value, y = metric, color = set)) + 
  geom_boxplot() 


daughter_cell_data %>% left_join(modes) %>% 
  ggplot(aes( total_cluster_intensity, y = "")) + 
  geom_boxplot() + geom_jitter(aes(color = as.factor(mode), shape = mother_cell_id == "Cell 16"), size = 3) + 
  labs(x = "Total Intensity of LANA dot", y = "Live KSHV Daughter Cells", color = "Estimated Number of Episomes in LANA dot") + 
  theme(legend.position = "bottom") + scale_color_manual(values=safe_colorblind_palette) + 
  scale_shape_manual(values = c(19, 8)) + 
  guides(shape = "none")
  