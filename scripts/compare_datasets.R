# This script compares the results of our statistical framework applied to datasets from four different experimental designs

#####
#  Setup 
#####

# load libraries 
library(tidyverse)
library(here)
# For plotting
library(patchwork)
library(RColorBrewer)
library(ggExtra)

select <- dplyr::select
theme_set(theme_bw())

# Check wd
here()

#####
# Read in results
#####

# fixed 8TR
load(here("results", "Figure_2", "mother_cell_Gibbs_samples.RData"))
load(here("results", "Figure_2", "daughter_cell_Gibbs_samples.RData"))

mother_cells_figure_2 <- mother_cell_samples_long
daughter_cells_figure_2 <- cell_samples_long

# live 8TR
load(here("results", "Figure_3", "mother_cell_Gibbs_samples.RData"))
load(here("results", "Figure_3", "daughter_cell_Gibbs_samples.RData"))

mother_cells_figure_3 <- mother_cell_samples_long
daughter_cells_figure_3 <- cell_samples_long

# fixed KSHV
load(here("results", "Figure_5", "mother_cell_Gibbs_samples.RData"))
load(here("results", "Figure_5", "daughter_cell_Gibbs_samples.RData"))

mother_cells_figure_5 <- mother_cell_samples_long
daughter_cells_figure_5 <- cell_samples_long

# live KSHV
load(here("results", "Figure_6", "mother_cell_Gibbs_samples.RData"))
load(here("results", "Figure_6", "daughter_cell_Gibbs_samples.RData"))

mother_cells_figure_6 <- mother_cell_samples_long
daughter_cells_figure_6 <- cell_samples_long


#####
# Format data
#####

# find mode of mother distribution for each cell:

get_mode <- function(data){
  data %>% filter(chain %in% c("chain 1", "chain1")) %>% 
    count(cell_id, number_of_episomes)  %>% 
    group_by(cell_id) %>% 
    filter(n == max(n)) %>%
    ungroup %>% 
    select(cell_id, mode = number_of_episomes)
}

mother_cell_distributions <- rbind(get_mode(mother_cells_figure_2) %>% mutate(data = "Fixed 8TR"),
                                   get_mode(mother_cells_figure_3) %>% mutate(data = "Live 8TR"),
                                   get_mode(mother_cells_figure_5) %>% mutate(data = "Fixed KSHV"),
                                   get_mode(mother_cells_figure_6) %>% mutate(data = "Live KSHV")
)

daughter_cell_distributions <- rbind(get_mode(daughter_cells_figure_2) %>% mutate(data = "Fixed 8TR"),
                                     get_mode(daughter_cells_figure_3) %>% mutate(data = "Live 8TR"),
                                     get_mode(daughter_cells_figure_5) %>% mutate(data = "Fixed KSHV"),
                                     get_mode(daughter_cells_figure_6) %>% mutate(data = "Live KSHV")
)

daughter_cell_breakdown <- daughter_cell_distributions %>% 
  separate(cell_id, into = c("mother_cell", "daughter_cell"), sep = "_") %>% 
  group_by(mother_cell, data) %>% 
  summarise(larger_daughter = max(mode), smaller_daughter = ifelse(n() == 2, min(mode), 0 ), total = sum(mode)) %>% 
  ungroup

# 
# mother_cell_distributions %>% 
#   ggplot(aes(mode)) + 
#   geom_histogram(aes(color = "mother cells"), alpha = 0.5) + 
#   # geom_histogram(data = daughter_cell_breakdown, 
#   #          aes(larger_daughter, color = "larger daughter cell"), alpha = 0.5) +
#   # geom_histogram(data = daughter_cell_breakdown, 
#   # aes(smaller_daughter, color = "smaller daughter cell"), alpha = 0.5) +
#   geom_histogram(data = daughter_cell_breakdown, 
#                  aes(total, color = "total episomes"), alpha = 0.5) +
#   facet_wrap(~data) + 
#   scale_x_continuous(breaks = 1:100)
# 
# mother_cell_distributions %>% 
#   filter(data == "Fixed 8TR") %>% 
#   ggplot(aes(mode)) + 
#   geom_bar(aes(color = "mother cells"), alpha = 0.5) +
#   geom_bar(data = daughter_cell_breakdown,
#            aes(larger_daughter, color = "daughter cell\nwith more episomes"), alpha = 0.5) +
#   geom_bar(data = daughter_cell_breakdown,
#            aes(smaller_daughter, color = "daughter cell\nwith fewer episomes"), alpha = 0.5) +
#   geom_bar(data = daughter_cell_breakdown,
#            aes(total, color = "both daughter cells"), alpha = 0.5) +
#   # facet_wrap(~data) + 
#   scale_x_continuous(breaks = 0:100) + 
#   labs(x = "number of episomes per cell", color = "") + 
#   theme(legend.position = c(1,1), legend.justification = c(1,1),
#         legend.background = element_blank())
# 
# 
# 
# daughter_cell_breakdown %>% 
#   filter(data == "Fixed 8TR") %>% 
#   pivot_longer(3:5, names_to = "cell_category", values_to = "n_epi") %>% 
#   ggplot() + 
#   geom_boxplot(data = mother_cell_distributions %>% 
#                  filter(data == "Fixed 8TR") ,
#                aes(mode, y = "mother cells"), fill = NA) +
#   geom_boxplot(data = mother_cell_distributions %>% 
#                  filter(data == "Fixed 8TR") ,
#                aes(mode*2, y = "ideal total"), fill = NA) +
#   geom_boxplot(aes(n_epi,  y = cell_category), fill = NA) + 
#   labs(x = "number of episomes per cell", color = "",
#        title = "Distribution of episomes per cell in fixed 8TR images") + 
#   theme(legend.position = c(1,1), legend.justification = c(1,1),
#         legend.background = element_blank())


rbind(
  daughter_cell_breakdown %>% 
    pivot_longer(3:5, names_to = "cell_category", values_to = "n_epi") %>% select(-mother_cell),
  mother_cell_distributions %>% mutate(cell_category = "mother cells") %>% rename(n_epi = mode ) %>% select(-cell_id)
) %>% 
  count(data, cell_category, n_epi) %>%
  group_by(data, cell_category) %>% mutate(freq = n/sum(n)) %>%
  mutate(cell_category = factor(cell_category, levels = c("mother cells", "total", "larger_daughter", "smaller_daughter"),
                                labels = paste0("Frequency in\n", c("Mother cells", "Both daughter cells", "Daughter cell with\nmore episomes", 
                                                                    "Daughter cell with\nfewer epsiomes")))) %>% 
  # filter(data == "Fixed 8TR") %>% 
  ggplot() + 
  geom_bar(aes(n_epi, freq, color = cell_category), fill = NA, stat = "identity", position = "identity") +
  # facet_wrap(~data) +
  facet_grid(cell_category~data, switch = "y") +
  scale_x_continuous(breaks = 0:100) + 
  labs(x = "number of episomes per cell", color = "", y = "",
       title = "Distribution of episomes per cell in fixed 8TR images") + 
  theme(#legend.position = c(1,1), legend.justification = c(1,1),
    legend.background = element_blank(), strip.placement = "outside",
    axis.title.y = element_blank(), strip.background =  element_blank(), legend.position = "none") + 
  scale_color_brewer(palette = "Set2")

rbind(
  daughter_cell_breakdown %>% 
    pivot_longer(3:5, names_to = "cell_category", values_to = "n_epi") %>% select(-mother_cell),
  mother_cell_distributions %>% mutate(cell_category = "mother cells") %>% rename(n_epi = mode ) %>% select(-cell_id)
) %>% 
  count(data, cell_category, n_epi) %>%
  group_by(data, cell_category) %>% mutate(freq = n/sum(n)) %>%
  mutate(cell_category = factor(cell_category, levels = c("mother cells", "total", "larger_daughter", "smaller_daughter"),
                                labels = paste0("Frequency in\n", c("Mother cells", "Both daughter cells", "Daughter cell with\nmore episomes", 
                                                                    "Daughter cell with\nfewer epsiomes")))) %>% 
  filter(data == "Fixed 8TR") %>%
  ggplot() + 
  geom_bar(aes(n_epi, freq, color = cell_category), fill = NA, stat = "identity", position = "dodge") +
  # facet_wrap(~data) +
  # facet_grid(cell_category~data, switch = "y") +
  scale_x_continuous(breaks = 0:100) + 
  labs(x = "number of episomes per cell", color = "", y = "",
       title = "Distribution of episomes per cell in fixed 8TR images") + 
  theme(#legend.position = c(1,1), legend.justification = c(1,1),
    legend.background = element_blank(), strip.placement = "outside",
    axis.title.y = element_blank(), strip.background =  element_blank(), legend.position = "none") + 
  scale_color_brewer(palette = "Set2")

rbind(
  daughter_cell_breakdown %>% 
    pivot_longer(3:5, names_to = "cell_category", values_to = "n_epi") %>% select(-mother_cell),
  mother_cell_distributions %>% mutate(cell_category = "mother cells") %>% rename(n_epi = mode ) %>% select(-cell_id)
) %>% 
  count(data, cell_category, n_epi) %>%
  group_by(data, cell_category) %>% mutate(freq = n/sum(n)) %>%
  mutate(cell_category = factor(cell_category, levels = c("mother cells", "total", "larger_daughter", "smaller_daughter"),
                                labels =c("Mother cells", "Both daughter cells", "Daughter cell with\nmore episomes", 
                                                                    "Daughter cell with\nfewer epsiomes"))) %>% 
  filter(data == "Fixed 8TR") %>%
  ggplot() + 
  geom_bar(aes(n_epi, freq, fill = cell_category), stat = "identity", position = "dodge") +
  # facet_wrap(~data) +
  # facet_grid(cell_category~data, switch = "y") +
  scale_x_continuous(breaks = 0:100) + 
  labs(x = "number of episomes per cell", fill = "", y = "",
       title = "Distribution of episomes per cell in fixed 8TR images") + 
  theme(legend.position = c(1,1), legend.justification = c(1,1),
    legend.background = element_blank(), strip.placement = "outside",
     strip.background =  element_blank()) + 
  scale_fill_brewer(palette = "Set2")

daughter_cell_breakdown %>% 
  filter(data == "Fixed 8TR") %>% 
  count(data, larger_daughter, smaller_daughter) %>% arrange(data, desc(n)) %>% 
  group_by(data) %>% mutate(freq = n/sum(n)) %>% slice(1:10) %>% ungroup %>%
  mutate(pair = paste0("(", larger_daughter, ",", smaller_daughter, ")")) %>% 
  ggplot(aes(pair, freq)) + 
  geom_bar(stat = "identity") + 
  # facet_wrap(~data, scales = "free_x")
  labs(x = "Number of episomes in each daughter cell (cell 1, cell 2)", y = "frequency",
       title = "Frequency of daughter cell episome-splits") 
