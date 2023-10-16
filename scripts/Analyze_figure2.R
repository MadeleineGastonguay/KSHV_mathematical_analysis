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

data <- load_data(mother_cell_file, daughter_cell_file)

daughter_cell_data <- data$daughter_cell_data
mother_cell_data <- data$mother_cell_data

cat(length(unique(daughter_cell_data$mother_cell_id)), "daughter cell pairs with", nrow(daughter_cell_data), "clusters\n")
cat(length(unique(mother_cell_data$cell_id)), "non-dividing cells with", nrow(mother_cell_data), "clusters")

Figure2_results <- run_pipeline(daughter_cell_data, mother_cell_data, results_folder)

Figure2_results$MLE_grid$estimates

make_plots(Figure2_results, daughter_cell_data, mother_cell_data, results_folder)

daughter_cell_samples <- Figure2_results$daughter_cell_samples
figures <- figures(daughter_cell_data, mother_cell_data, daughter_cell_samples, results_folder)

