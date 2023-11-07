# This script benchmarks the likelihood in simulated data

#####
#  Setup 
#####

# load libraries 
library(tidyverse)
library(here)
# library(furrr)

# future::plan(multisession, workers = 4)

theme_set(theme_bw())

# Check wd
here()

# Read in helper functions 
source(here("scripts", "likelihood_functions.R"))
source(here("scripts", "run_pipeline.R"))

#####
# Simulate data
#####

set.seed(400)

# test a variety of Pr and Ps values and sample sizes:
values <- c(0.1, 0.3, 0.5, 0.8, 0.9)
n <- c(10, 40, 100)

sim_parameters <- expand_grid(Pr = values, Ps = values, n = n)

sim_data <- sim_parameters %>% 
  pmap(function(Pr, Ps, n){
    X0s <- sample(1:100, n, replace = T, prob = dpois(1:100, 3))
    simulate_multiple_cells(Pr = Pr, Ps = Ps, n_cells = n, X0s = X0s )
  } ) %>% 
  list_rbind(names_to = "simulation")


glimpse(sim_data)

## Apply MLE to each dataset, assuming X0 is known:
MLE_with_known_X0 <- sim_data %>% 
  group_by(simulation) %>% 
  nest() %>% 
  pull(data) %>% 
  map(run_grid_search, viz = F, known_X0 = T)

MLE_with_known_X0_estimates <- MLE_with_known_X0 %>% lapply("[[", "estimates") %>% 
  list_rbind(names_to = "simulation") %>% remove_rownames() %>% 
  mutate(rep = 1)
  # merge(sim_parameters %>% mutate(simulation = 1:nrow(.)), by = "simulation") 

MLE_with_known_X0_grid <- MLE_with_known_X0 %>% lapply("[[", "grid_search") %>% 
  list_rbind(names_to = "simulation") %>% remove_rownames() %>% mutate(rep = 1)

#####
# Plots
#####

# MLE_with_known_X0_estimates %>% 
#   ggplot(aes(Pr, MLE_Pr)) + 
#   geom_point() + 
#   geom_errorbar(aes(ymin = min_Pr, ymax = max_Pr), width = 0, alpha = 0.7) +
#   geom_abline() +
#   facet_grid(Ps~n)
# 
# MLE_with_known_X0_estimates %>% 
#   ggplot(aes(Ps, MLE_Ps)) + 
#   geom_point() + 
#   geom_errorbar(aes(ymin = min_Ps, ymax = max_Ps), width = 0, alpha = 0.7) +
#   geom_abline() +
#   facet_grid(Pr~n)
# 
# MLE_with_known_X0_estimates %>% 
#   mutate(Pr_error = MLE_Pr - Pr) %>% 
#   ggplot(aes(Pr, Ps, fill = Pr_error)) + 
#   geom_tile() + 
#   facet_wrap(~n) + 
#   scale_fill_distiller(type = "div")

#####
# Apply with unknown X0
#####

MLE_with_unknown_X0 <- sim_data %>%
  group_by(simulation) %>%
  nest() %>%
  pull(data) %>%
  map(run_grid_search, viz = F, known_X0 = F, lambda = 3)

MLE_with_unknown_X0_estimates <- MLE_with_unknown_X0 %>% lapply("[[", "estimates") %>% 
  list_rbind(names_to = "simulation") %>% remove_rownames() %>% mutate(rep = 1)

MLE_with_unknown_X0_grid <- MLE_with_unknown_X0 %>% lapply("[[", "grid_search") %>% 
  list_rbind(names_to = "simulation") %>% remove_rownames() %>% mutate(rep = 1)

data <- list(sim_data)

## Repeat 9 more times:
for(i in 2:10){
  sim_data.2 <- sim_parameters %>% 
    pmap(function(Pr, Ps, n){
      X0s <- sample(1:100, n, replace = T, prob = dpois(1:100, 3))
      simulate_multiple_cells(Pr = Pr, Ps = Ps, n_cells = n, X0s = X0s )
    } ) %>% 
    list_rbind(names_to = "simulation")
  
  data[[i]] <- sim_data.2
  
  MLE_with_unknown_X0.2 <- sim_data.2 %>% 
    group_by(simulation) %>%
    nest() %>%
    pull(data) %>%
    map(run_grid_search, viz = F, known_X0 = F, lambda = 3)
  
  MLE_with_unknown_X0_estimates <- rbind(MLE_with_unknown_X0_estimates, 
                                         MLE_with_unknown_X0.2 %>% lapply("[[", "estimates") %>% 
                                           list_rbind(names_to = "simulation") %>% remove_rownames() %>% mutate(rep = i)   )
  
  MLE_with_unknown_X0_grid <- rbind(MLE_with_unknown_X0_grid, 
                                         MLE_with_unknown_X0.2 %>% lapply("[[", "grid_search") %>% 
                                           list_rbind(names_to = "simulation") %>% remove_rownames() %>% mutate(rep = i)   )
  
  MLE_with_known_X0.2 <- sim_data.2 %>% 
    group_by(simulation) %>%
    nest() %>%
    pull(data) %>%
    map(run_grid_search, viz = F, known_X0 = T)
  
  MLE_with_known_X0_estimates <- rbind(MLE_with_known_X0_estimates, 
                                       MLE_with_known_X0.2 %>% lapply("[[", "estimates") %>% 
                                         list_rbind(names_to = "simulation") %>% remove_rownames() %>% mutate(rep = i)   )
  
  MLE_with_known_X0_grid <- rbind(MLE_with_known_X0_grid, 
                                    MLE_with_known_X0.2 %>% lapply("[[", "grid_search") %>% 
                                      list_rbind(names_to = "simulation") %>% remove_rownames() %>% mutate(rep = i)   )

}

MLE_with_unknown_X0_estimates <- MLE_with_unknown_X0_estimates %>%  
  left_join(sim_parameters %>% mutate(simulation = 1:nrow(.)),  by = "simulation") 

MLE_with_known_X0_estimates <- MLE_with_known_X0_estimates %>%  
  left_join(sim_parameters %>% mutate(simulation = 1:nrow(.)),  by = "simulation") 


save(MLE_with_known_X0_estimates,  data,  MLE_with_known_X0_grid,
     file = here("results", "benchmarking", "MLE_results_knownX0.RData"))

save(MLE_with_unknown_X0_estimates, data, MLE_with_unknown_X0_grid,
     file = here("results", "benchmarking", "MLE_results_unknownX0.RData"))

#####
# Plots
#####

# MLE_with_unknown_X0_estimates %>% 
#   ggplot(aes(Pr, MLE_Pr)) + 
#   geom_point() + 
#   geom_errorbar(aes(ymin = min_Pr, ymax = max_Pr), width = 0, alpha = 0.7) +
#   geom_abline() +
#   facet_grid(Ps~n)
# 
# MLE_with_unknown_X0_estimates %>% 
#   ggplot(aes(Ps, MLE_Ps)) + 
#   geom_point() + 
#   geom_errorbar(aes(ymin = min_Ps, ymax = max_Ps), width = 0, alpha = 0.7) +
#   geom_abline() +
#   facet_grid(Pr~n)
# 
# rbind(
#   MLE_with_unknown_X0_estimates %>% 
#     mutate(analysis = "unknown X0"),
#   MLE_with_known_X0_estimates %>% 
#     mutate(analysis = "known X0")
# ) %>% 
#   select(MLE_Pr, analysis, simulation, n, Ps) %>% 
#   pivot_wider(names_from = analysis, values_from = MLE_Pr) %>% 
#   ggplot(aes(`known X0`, `unknown X0`, color = Ps)) + 
#   geom_point() + 
#   geom_abline() + 
#   facet_wrap(~n)
# 
# rbind(
#   MLE_with_unknown_X0_estimates %>% 
#     mutate(analysis = "unknown X0"),
#   MLE_with_known_X0_estimates %>% 
#     mutate(analysis = "known X0")
# ) %>% 
#   select(MLE_Ps, analysis, simulation, n, Ps, Pr) %>% 
#   pivot_wider(names_from = analysis, values_from = MLE_Ps) %>% 
#   ggplot(aes(`known X0`, `unknown X0`, color = Ps)) + 
#   geom_point() + 
#   geom_abline() + 
#   facet_wrap(~n)
# 
# 
# rbind(
#   MLE_with_unknown_X0_estimates %>% 
#     mutate(analysis = "unknown X0"),
#   MLE_with_known_X0_estimates %>% 
#     mutate(analysis = "known X0")
# ) %>% 
#   mutate(Pr_width = max_Pr - min_Pr) %>% 
#   select(Pr_width, analysis, simulation, n, Ps, Pr) %>% 
#   pivot_wider(names_from = analysis, values_from = Pr_width) %>% 
#   ggplot(aes(`known X0`, `unknown X0`, color = Pr)) + 
#   geom_point() + 
#   geom_abline() + 
#   facet_wrap(~n)
# 
# rbind(
#   MLE_with_unknown_X0_estimates %>% 
#     mutate(analysis = "unknown X0"),
#   MLE_with_known_X0_estimates %>% 
#     mutate(analysis = "known X0")
# ) %>% 
#   mutate(Ps_width = max_Ps - min_Ps) %>% 
#   select(Ps_width, analysis, simulation, n, Ps, Pr) %>% 
#   pivot_wider(names_from = analysis, values_from = Ps_width) %>% 
#   ggplot(aes(`known X0`, `unknown X0`, color = Pr)) + 
#   geom_point() + 
#   geom_abline() + 
#   facet_wrap(~n)
# 
# 
# plot_grid_search(MLE_with_known_X0[[74]], prob = T) |
#   plot_grid_search(MLE_with_unknown_X0[[74]], prob = T) 
#####
# Supplemental Figures
#####

Pr_estimates <- MLE_with_unknown_X0_estimates %>% 
  mutate(Ps = paste0(Ps*100, "%")) %>% 
  ggplot(aes(MLE_Pr - Pr, as.factor(Pr), group = rep, color = min_Pr <= Pr & max_Pr >= Pr)) + 
  geom_vline(xintercept = 0, lty = "dashed", color = "lightgray")  +
  geom_point(position = position_dodge2(width = 0.7)) + 
  geom_linerange(aes(xmin = min_Pr - Pr, xmax=  max_Pr - Pr ),
                 position = position_dodge2(width = 0.7)) +
  facet_grid(n~Ps, labeller = "label_both") +
  scale_color_manual(values = c("black", "darkgray")) +
  guides(color ="none") + 
  labs(y = "Pr", x = "Difference between MLE and parameter value (MLE_Pr - Pr)") +
  scale_x_continuous(limits = c(-1,1)) + 
  theme(panel.grid.major.y =  element_blank(), panel.grid.minor.x = element_blank())

Ps_estimates <- MLE_with_unknown_X0_estimates %>% 
  mutate(Pr = paste0(Pr*100, "%")) %>% 
  ggplot(aes(MLE_Ps - Ps, as.factor(Ps), group = rep, color = min_Ps <= Ps & max_Ps >= Ps)) + 
  geom_vline(xintercept = 0, lty = "dashed", color = "lightgray")  +
  geom_point(position = position_dodge2(width = 0.7)) + 
  geom_linerange(aes(xmin = min_Ps - Ps, xmax=  max_Ps - Ps ),
                 position = position_dodge2(width = 0.7)) +
  facet_grid(n~Pr, labeller = "label_both") +
  scale_color_manual(values = c("black", "darkgray")) +
  guides(color ="none") + 
  labs(y = "Ps", x = "Difference between MLE and parameter value (MLE_Ps - Ps)") + 
  scale_x_continuous(limits = c(-1,1)) + 
  theme(panel.grid.major.y =  element_blank(), panel.grid.minor.x = element_blank())

ggsave(here("results", "benchmarking", "Pr_estimates.png"), Pr_estimates, width = 7, height = 6)
ggsave(here("results", "benchmarking", "Ps_estimates.png"), Ps_estimates, width = 7, height = 6)

Pr_estimates_knownX0 <- MLE_with_known_X0_estimates %>% 
  mutate(Ps = paste0(Ps*100, "%")) %>% 
  ggplot(aes(MLE_Pr - Pr, as.factor(Pr), group = rep, color = min_Pr <= Pr & max_Pr >= Pr)) + 
  geom_vline(xintercept = 0, lty = "dashed", color = "lightgray")  +
  geom_point(position = position_dodge2(width = 0.7)) + 
  geom_linerange(aes(xmin = min_Pr - Pr, xmax=  max_Pr - Pr ),
                 position = position_dodge2(width = 0.7)) +
  facet_grid(n~Ps, labeller = "label_both") +
  scale_color_manual(values = c("black", "darkgray")) +
  guides(color ="none") + 
  labs(y = "Pr", x = "Difference between MLE and parameter value (MLE_Pr - Pr)") +
  scale_x_continuous(limits = c(-1,1)) + 
  theme(panel.grid.major.y =  element_blank(), panel.grid.minor.x = element_blank())

Ps_estimates_knownX0 <- MLE_with_known_X0_estimates %>% 
  mutate(Pr = paste0(Pr*100, "%")) %>% 
  ggplot(aes(MLE_Ps - Ps, as.factor(Ps), group = rep, color = min_Ps <= Ps & max_Ps >= Ps)) + 
  geom_vline(xintercept = 0, lty = "dashed", color = "lightgray")  +
  geom_point(position = position_dodge2(width = 0.7)) + 
  geom_linerange(aes(xmin = min_Ps - Ps, xmax=  max_Ps - Ps ),
                 position = position_dodge2(width = 0.7)) +
  facet_grid(n~Pr, labeller = "label_both") +
  scale_color_manual(values = c("black", "darkgray")) +
  guides(color ="none") + 
  labs(y = "Ps", x = "Difference between MLE and parameter value (MLE_Ps - Ps)") + 
  scale_x_continuous(limits = c(-1,1)) + 
  theme(panel.grid.major.y =  element_blank(), panel.grid.minor.x = element_blank())

ggsave(here("results", "benchmarking", "Pr_estimates_knownX0.png"), Pr_estimates_knownX0, width = 7, height = 6)
ggsave(here("results", "benchmarking", "Ps_estimates_knownX0.png"), Ps_estimates_knownX0, width = 7, height = 6)

merge(
  MLE_with_known_X0_estimates %>% select(simulation, MLE_Pr, MLE_Ps, rep, Pr, Ps, n),
  MLE_with_unknown_X0_estimates %>% select(simulation, MLE_Pr, MLE_Ps, rep),
  by = c("simulation", "rep")
) %>% 
  ggplot() + 
  geom_point(data = . %>% mutate(param = "Pr"), aes(MLE_Pr.x, MLE_Pr.y, color = Ps))  + 
  geom_point(data = . %>% mutate(param = "Ps"), aes(MLE_Ps.x, MLE_Ps.y, color = Pr)) +
  geom_abline() + 
  facet_grid(n~param)

