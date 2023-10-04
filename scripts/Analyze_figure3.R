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
results_folder <- here("results", "Figure_3_updated")

# Daughter cell intensity data
daughter_cell_file <- "Fig3 dividing cells.xlsx"

# Mother cell intensity data
mother_cell_file <- "Fig3 Non dividing cells.xlsx"
#####

data <- load_data(mother_cell_file, daughter_cell_file)

daughter_cell_data <- data$daughter_cell_data
mother_cell_data <- data$mother_cell_data

Figure3_results <- run_pipeline(daughter_cell_data, mother_cell_data, results_folder, overwrite = T)

Figure3_results$MLE_grid$estimates

make_plots(Figure3_results, daughter_cell_data, mother_cell_data, results_folder)
