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
results_folder <- here("results", "Figure_5_updated_pdf")

# Daughter cell intensity data
daughter_cell_file <- "Fig5 dividing cells.xlsx"

# Mother cell intensity data
mother_cell_file <- "Fig5 Non dividing cells.xlsx"
#####

Figure5_data <- load_data(mother_cell_file, daughter_cell_file)

daughter_cell_data <- Figure5_data$daughter_cell_data
mother_cell_data <- Figure5_data$mother_cell_data

Figure5_results <- run_pipeline(daughter_cell_data, mother_cell_data, results_folder)

Figure5_results$MLE_grid$estimates

make_plots(Figure5_results, daughter_cell_data, mother_cell_data, results_folder)

daughter_cell_samples <- Figure5_results$daughter_cell_samples
figures <- figures(daughter_cell_data, mother_cell_data, daughter_cell_samples, results_folder)


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
  scale_x_continuous(breaks = seq(0,10, by = 2)) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1), breaks = c(0, 1))  +
  theme(strip.background = element_blank(), strip.text = element_blank()) 

ggsave(here(results_folder, "example_inference.png"), width = 5, height = 1.2)

