
MCMC_vs_data <-  function(pipeline_output, data, results_folder, dataset, same_mu = T){
  
  all_chains <- pipeline_output$all_chains
  daughter_cell_samples <- pipeline_output$daughter_cell_samples
  mother_cell_samples <- pipeline_output$mother_cell_samples
  lambda <- pipeline_output$X0_lambda
  convergence <- pipeline_output$convergence_metrics
  MLE_grid <- pipeline_output$MLE_grid
  
  mother_cell_data <- data$mother_cell_data
  daughter_cell_data <- data$daughter_cell_data
  
  intensity_data <- rbind(daughter_cell_data %>% select(cell_id, cluster_id, total_cluster_intensity, min_episome_in_cluster) %>% mutate(set = "daughter"), 
                          mother_cell_data  %>% select(cell_id, cluster_id, total_cluster_intensity, min_episome_in_cluster) %>% mutate(set = "mother"))
  
  ## Color-blind friendly color palette
  safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  
  # How do these inferences relate to the observed intensity data?
  cluster_modes <- all_chains %>% 
    filter(chain == "chain1") %>% 
    select(grep("[0-9]", names(all_chains))) %>% 
    pivot_longer(everything(), names_to = "cluster_id", values_to = "n_episomes") %>% 
    group_by(cluster_id) %>% 
    mutate(sd = sd(n_episomes)) %>% 
    count(cluster_id, n_episomes, sd) %>% 
    group_by(cluster_id) %>%
    mutate(p = n/sum(n)) %>% 
    filter(n == max(n)) %>% 
    ungroup %>% 
    select(-n) %>% 
    rename(mode = n_episomes) %>% 
    mutate(set = ifelse(grepl("d", cluster_id), "daughter", "mother"))
  
  if(!same_mu){
    dat <- intensity_data %>% 
      merge(cluster_modes) %>%
      mutate(total_cluster_intensity = ifelse(set == "daughter", total_cluster_intensity/DescTools::Mode(round(all_chains$mu_d)),
                                              total_cluster_intensity/DescTools::Mode(round(all_chains$mu_m))),
             set = str_to_title(paste(set, "cell")),
             min_episome_in_cluster = factor(min_episome_in_cluster, levels = 1:4))
  }else{
    dat <- intensity_data %>% 
      merge(cluster_modes) %>%
      mutate(total_cluster_intensity = total_cluster_intensity/DescTools::Mode(round(all_chains$mu)),
             set = str_to_title(paste(set, "cell")),
             min_episome_in_cluster = factor(min_episome_in_cluster, levels = 1:4))
  }
  
  MCMC_vs_intensity_data <- dat %>% 
    ggplot(aes(total_cluster_intensity, as.factor(mode))) + 
    geom_jitter(aes(shape = set, color = sd)) +
    geom_boxplot(outlier.shape = NA, fill = NA) +
    labs(y = "Mode of the posterior distribution", 
         color = "Standard deviation\nof the\nposterior distribution",
         shape = "Cell set",
         x = "Total intensity of LANA dot normalized by mean intensity of a single episome") + 
    theme(legend.position = c(0,1), legend.justification = c(-0.1,1.1)) + 
    scale_x_continuous(breaks = 0:20) + 
    scale_shape_manual(values = c(2,19)) + 
    scale_color_viridis_c(limits = c(0, 0.9)) + 
    guides(shape = "none")
  
  intensity_data_plot <- intensity_data %>% 
    mutate(set = str_to_title(paste(set, "cell"))) %>% 
    ggplot(aes(total_cluster_intensity, shape = set, y = set)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter() + 
    scale_shape_manual(values = c(2,19)) + 
    labs(x = "Total intensity of LANA dot",
         title = dataset) +
    theme(legend.position = "none", axis.title.y = element_blank())
  
  # mode_plot <- cluster_modes %>%
  #   mutate(set = str_to_title(paste(set, "cell"))) %>%
  #   ggplot(aes(mode, shape = set, x = set, color= sd)) +
  #   geom_boxplot(outlier.shape = NA) +
  #   geom_jitter() +
  #   scale_shape_manual(values = c(2,19)) +
  #   labs(y = "Mode of the posterior distribution",
  #        color = "Standard deviation of\nthe posteiror distribution") +
  #   theme(legend.position = "none", axis.title.y = element_blank()) +
  #   scale_y_continuous(breaks = 0:10)
  
  plot <- intensity_data_plot + #plot_spacer() +
    MCMC_vs_intensity_data  +
    plot_layout(ncol = 1, heights = c(1,2)) & 
    theme(legend.position = "right")
  
  # mode_plot +
  # plot_layout(ncol = 2, heights = c(1,3), widths = c(3,1), guides = "collect")

  file <- switch(tolower(dataset),
                 "fixed 8tr" = "Figure2",
                 "live 8tr" = "Figure3",
                 "fixed kshv" = "Figure5",
                 "live kshv" = "Figure6")
  file <- paste0(file, "_MCMC_inference.png")

  ggsave(here(results_folder, file), plot, width = 8, height = 5)
  
}

file_path <- here("results", "supplemental_figures")

MCMC_vs_data(Figure2_results, Figure2_data, file_path, "Fixed 8TR")  
MCMC_vs_data(Figure3_results, Figure3_data, file_path, "Live 8TR")  
MCMC_vs_data(Figure5_results, Figure5_data, file_path, "Fixed KSHV")  
MCMC_vs_data(Figure6_results, Figure6_data, file_path, "Live KSHV", same_mu = F)  

certain_example <- Figure2_results$all_chains %>% filter(chain == "chain1") %>% 
  ggplot(aes(d107, after_stat(count/sum(count)))) + geom_bar(aes(fill = sd(d107)), show.legend = F) + 
  labs(x = "Number of episomes in LANA dot", y = "Posterior probability") +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank()) +
  ylim(c(0,1)) + 
  scale_fill_viridis_c(limits = c(0, 0.9)) 

uncertain_example <- Figure2_results$all_chains %>% filter(chain == "chain1") %>% 
  ggplot(aes(m78, after_stat(count/sum(count)))) + geom_bar(aes(fill = sd(m78)), show.legend = F) + 
  labs(x = "Number of episomes in LANA dot", y = "Posterior probability") +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank()) +
  ylim(c(0,1)) +
  scale_fill_viridis_c(limits = c(0, 0.9)) 

ggsave(here("results", "supplemental_figures", "certain_example.png"), certain_example, width = 3, height = 2)
ggsave(here("results", "supplemental_figures", "uncertain_example.png"), uncertain_example, width = 3, height = 2)

# Plot MLE grids for all conditions
plot_grid <- function(results, dataset){
  plot <- plot_grid_search(results$MLE_grid, simulation = F, prob = T, error_bars = F) & 
    plot_annotation(title = element_blank()) & 
    theme(plot.caption = element_blank(), plot.background = element_blank())
  
  file <- switch(tolower(dataset),
                 "fixed 8tr" = "Figure2",
                 "live 8tr" = "Figure3",
                 "fixed kshv" = "Figure5",
                 "live kshv" = "Figure6")
  file <- paste0(file, "_MLE_grid.png")
  
  
  ggsave(here("results", "supplemental_figures", file), plot, width = 4.5, height = 3.5, bg = "white") 
  
}

plot_grid(Figure2_results, "Fixed 8TR")
plot_grid(Figure3_results, "Live 8TR")
plot_grid(Figure5_results, "Fixed KSHV")
plot_grid(Figure6_results, "Live KSHV")
 
Figure2_results$all_chains %>% filter(chain == "chain1") %>% select(grep("[0-9]", names(.), perl = T)) %>% apply(2, sd) %>% which.min

