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
results_folder <- here("results", "Figure_2_updated")

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

Figure2_results <- run_pipeline(daughter_cell_data, mother_cell_data, results_folder)

Figure2_results$MLE_grid$estimates

make_plots(Figure2_results, daughter_cell_data, mother_cell_data, results_folder)

daughter_cell_samples <- Figure2_results$daughter_cell_samples
figures <- figures(daughter_cell_data, mother_cell_data, daughter_cell_samples, results_folder)


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
  scale_x_continuous(breaks = seq(0,10, by = 2)) +
  # ylim(c(0,1)) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1), breaks = c(0, 1))  +
  theme(strip.background = element_blank(), strip.text = element_blank()) 

ggsave(here(results_folder, "example_inference.png"), width = 5, height = 1.2)
