## Understand why n=10 doesn't work

library(tidyverse)
library(here)

# Check wd
here()

# Read in helper functions 
source(here("scripts", "likelihood_functions.R"))

#####
# Simulated data with 10 cells, using MLE for Pr and Ps
#####

n_cells <- 10
Pr <- 0.80
Ps <- 0.89

X0s <- sample(1:100, n_cells, prob = dpois(1:100, 4.42)/sum(dpois(1:100, 4.42)), replace = T)
while(median(X0s) != 4){
  X0s <- sample(1:100, n_cells, prob = dpois(1:100, 4.42)/sum(dpois(1:100, 4.42)), replace = T)
}

simulated_data <- simulate_multiple_cells(X0s,  Pr, Ps, n_cells)

head(simulated_data)

hist(c(simulated_data$X1, simulated_data$X2))

# Apply algorithm with known X0
MLE_estimates_known_X0 <- run_grid_search(simulated_data)

plot_grid_search(MLE_estimates_known_X0, prob = T)


# Apply algorithm with unknown X0
MLE_estimates_unknown_X0 <- run_grid_search(simulated_data, known_X0 = F, 
                                            PMF = tibble(X0 = 1:50, prob = dpois(X0, lambda = 4.42)),
                                            viz = F
                                            )

plot_grid_search(MLE_estimates_unknown_X0, prob = T)



#####
# Test out larger X0 than assumed:
#####
X0s <- sample(1:100, n_cells, prob = dpois(1:100, 6)/sum(dpois(1:100, 6)), replace = T)
simulated_data_big <- simulate_multiple_cells(X0s,  Pr, Ps, n_cells)

MLE_estimates_known_X0 <- run_grid_search(simulated_data_big)
plot_grid_search(MLE_estimates_known_X0, prob = T)


MLE_estimates_unknown_X0 <- run_grid_search(simulated_data_big, known_X0 = F, 
                                            PMF = tibble(X0 = 1:20, prob = dpois(X0, lambda = 4.12)),
                                            viz = F
)
plot_grid_search(MLE_estimates_unknown_X0, prob = T)


#####
# Test out smaller X0 than assumed:
#####
X0s <- sample(1:100, n_cells, prob = dpois(1:100, 2)/sum(dpois(1:100, 2)), replace = T)
simulated_data_small <- simulate_multiple_cells(X0s,  Pr, Ps, n_cells)

MLE_estimates_known_X0 <- run_grid_search(simulated_data_small)
plot_grid_search(MLE_estimates_known_X0, prob = T)


MLE_estimates_unknown_X0 <- run_grid_search(simulated_data_small, known_X0 = F, 
                                            PMF = tibble(X0 = 1:20, prob = dpois(X0, lambda = 4.12)),
                                            viz = F
)
plot_grid_search(MLE_estimates_unknown_X0, prob = T)


# with larger sample size?
X0s <- sample(1:100, 40, prob = dpois(1:100, 2)/sum(dpois(1:100, 2)), replace = T)
simulated_data_small2 <- simulate_multiple_cells(X0s,  Pr, Ps, 40)

MLE_estimates_known_X0 <- run_grid_search(simulated_data_small2)
plot_grid_search(MLE_estimates_known_X0, prob = T)


MLE_estimates_unknown_X0 <- run_grid_search(simulated_data_small2, known_X0 = F, 
                                            PMF = tibble(X0 = 1:20, prob = dpois(X0, lambda = 4.12)),
                                            viz = F
)
plot_grid_search(MLE_estimates_unknown_X0, prob = T)
