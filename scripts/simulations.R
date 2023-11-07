## This script runs and plots simulations of a cell population with constant size

library(tidyverse)
library(here)
library(patchwork)
library(scico)
library(cowplot)
theme_set(theme_bw())
#####
# Simulation functions
#####

makeChildren<- function(pRep, pSeg, numEpisomes){
  if(numEpisomes == 0){
    return(c(0, 0))
  }
  else{
    outcomes = c(1, 0, 0, 1, 1, 1, 2, 0, 0, 2)
    refMat = matrix(outcomes, nrow= 5, byrow=T) #possible episome fates
    probVec = c(1/2*(1-pRep), 1/2*(1-pRep), pRep*(pSeg+1/2*(1-pSeg)), pRep*1/4*(1-pSeg), pRep*1/4*(1-pSeg))
    result = rmultinom(n=1, prob=probVec, size=numEpisomes) #simulate each independently and add the results in the return statement
    return(t(result)%*%refMat)
  }
}

simStepFlex <- function(pRep, pSeg, cells, birthVec, deathVec, selectAgainstZero = T, max_epi){
  #all FUN arguments are functions
  #cells is a compressed vector of the number of cells with 0, 1, 2, ... episomes
  if(sum(cells) == 0){
    return(cells)
  }
  else{
    type <- sample(length(cells), size=1, replace=TRUE, prob = cells*(birthVec+deathVec))
    b = birthVec[type]
    d = deathVec[type]
    timeAdvance=rexp(n=1, rate=sum(cells*(birthVec+deathVec)))
    test = runif(n=1)
    if(test < d/(b+d)){
      #print("DEATH")
      cells[type] = cells[type] - 1
      if(selectAgainstZero == TRUE){
        cells[1] <- 0
      }
      return(list(timeAdvance, cells))
    }
    else{
      #print("BIRTH")
      children <- makeChildren(pRep, pSeg, type-1)
      if(children[1] > max_epi){
        #print("OVERFLOW")
        children[1] = max_epi
      }
      if(children[2] > max_epi){
        #print("OVERFLOW")
        children[2] = max_epi
      }
      #print(children)
      #print(type)
      cells[type] = cells[type]-1
      cells[children[1]+1] = cells[children[1]+1] + 1
      cells[children[2]+1] = cells[children[2]+1] + 1
      if (selectAgainstZero==TRUE){
        cells[1] <- 0
      }
      return(list(timeAdvance, cells))
    }
  }
}

extinction <- function(pRep, pSeg, nTrials, n_epi){
  
  results = 1:nTrials*0
  times = c()
  totals = c()
  trials = c()
  zero = c()
  one = c()
  two = c()
  three = c()
  four = c()
  five = c()
  six = c()
  seven = c()
  eight = c()
  nine = c()
  
  i = 1
  j = 1
  for(z in 1:nTrials){
    print(z)
    time = 0
    deathVec = rep(1, 10)
    cells = rep(0, 10)
    cells[n_epi + 1] <- 1000
    total = sum((0:9)*cells)
    while(total > 0){
      if(time > 40 & pRep == 1 & pSeg == 1) break
      birthVec = rep(2*(1-sum(cells)/1000) + 1, 10)
      result = simStepFlex(pRep, pSeg, cells, birthVec, deathVec, FALSE, max_epi = 9)
      cells = result[[2]]
      time=time + result[[1]]
      #print(cells)
      total = sum((0:9)*cells)
      if(j %% 100 == 0 | j == 1){ # report out every 100 iterations
        times[i] <- time
        totals[i] <- total
        trials[i] <- z
        zero[i] <- cells[1]
        one[i] <- cells[2]
        two[i] <- cells[3]
        three[i] <- cells[4]
        four[i] <- cells[5]
        five[i] <- cells[6]
        six[i] <- cells[7]
        seven[i] <- cells[8]
        eight[i] <- cells[9]
        nine[i] <- cells[10]
        i <- i + 1
      }
      j <- j + 1
    }
    results[z]= time
  }
  total_df <- data.frame(trial = trials, time = times, total = totals,
                         zero, one, two, three, four, five, six, seven, eight, nine)
  return(list(ExtinctionTime = tibble(ExtinctionTime = results), Totals = total_df))
}

#####
# Plotting functions
#####

pivot_extinction <- function(extinction_out){
  if(!is.data.frame(extinction_out)) extinction_out <- extinction_out$Totals
  extinction_out %>% 
    pivot_longer(c("zero","one", "two", "three", "four", "five", "six", "seven", "eight", "nine"), 
                 names_to = "episomes_per_cell", values_to = "number_of_cells") %>% 
    mutate(episomes_per_cell = factor(fct_recode(episomes_per_cell, 
                                                 !!!c(`0` = "zero", `1`= "one", `2` = "two", `3`= "three", `4` = "four", 
                                                      `5` = "five", `6`= "six", `7` = "seven", `8` = "eight", `9` = "nine")
    ), levels = 0:9)
    )
}

number_over_time <- function(extinction_long){
  extinction_long %>% ggplot(aes(time, number_of_cells, color = episomes_per_cell)) + 
    geom_point(alpha = 0.25) +  
    guides(color = guide_legend(override.aes = list(alpha = 1))) + 
    facet_wrap(~episomes_per_cell)
}

fraction_over_time <- function(extinction_long, multiple = F){
  if(!multiple){
    df <- extinction_long %>% 
      group_by(time, trial) %>% 
      mutate(total_cells = sum(number_of_cells)) %>% 
      ungroup %>% mutate(frac = number_of_cells/total_cells) 
  }else{
    df <- extinction_long %>% 
      group_by(time, trial, pRep, pSeg) %>% 
      mutate(total_cells = sum(number_of_cells)) %>% 
      ungroup %>% mutate(frac = number_of_cells/total_cells) 
  }
  df %>% 
    ggplot(aes(time, frac, color = episomes_per_cell)) + 
    geom_point(alpha = 0.25) +  
    guides(color = guide_legend(override.aes = list(alpha = 1))) #+ 
}

#####
# Run 100 trials of simulations with varying parameters
#####
ntrials <- 100
n_epi <- 3 

# Indicator if simuilations should be rerun (TRUE) or loaded from previous run (FALSE)
rerun <- FALSE

if(!rerun){
  load(here("results", "simulations", "constant_pop_3epi.RData"))
  
  extinction_3epi_df_long <- pivot_extinction(extinction_3epi_df) %>% 
    group_by(time, trial, Pr = pRep, Ps = pSeg) %>% 
    mutate(total_cells = sum(number_of_cells)) %>% 
    ungroup %>% mutate(frac = number_of_cells/total_cells) %>% 
    arrange(desc(Pr), Ps) %>% 
    mutate(Pr = paste0(Pr*100, "%"),
           Ps = paste0(Ps*100, "%"),
           Pr = fct_inorder(factor(Pr)),
           Ps = fct_inorder(factor(Ps))) 
}else{
  
  ## set up a grid of parameters to simulate with and collect all results:
  pReps = c(0.5, 0.75, 0.9, 0.95, 1)
  pSegs = c(0.5, 0.75, 0.9, 0.95, 1)
  param_grid <- expand_grid(pRep = pReps, pSeg = pSegs)
  
  extinction_3epi <- param_grid %>% 
    pmap(extinction, nTrials = ntrials, n_epi = n_epi)
  
  extinction_3epi_df <- 1:nrow(param_grid) %>% 
    map_df(function(i) cbind(param_grid[i,], extinction_3epi[[i]]$Totals)) 
  
  extinction_3epi_df_long <- pivot_extinction(extinction_3epi_df) %>% 
    group_by(time, trial, Pr = pRep, Ps = pSeg) %>% 
    mutate(total_cells = sum(number_of_cells)) %>% 
    ungroup %>% mutate(frac = number_of_cells/total_cells) %>% 
    arrange(desc(Pr), Ps) %>% 
    mutate(Pr = paste0(Pr*100, "%"),
           Ps = paste0(Ps*100, "%"),
           Pr = fct_inorder(factor(Pr)),
           Ps = fct_inorder(factor(Ps))) 
  
  
  # Plot distribution of episomes per cell over time
  p1 <- extinction_3epi_df_long %>% 
    ggplot(aes(time, frac, color = episomes_per_cell, group = interaction(trial, episomes_per_cell))) + 
    geom_line(alpha = 0.25) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) + 
    ggh4x::facet_grid2(Pr ~ Ps, labeller = "label_both", scales = "free_x", independent = "x") + 
    labs(x = "time (generations)", y = "Fraction of\npopulation", color = "episomes per cell")
  
  ggsave(here("results", "simulations", "extinction_3epi_lines.png"), p1, width = 8, height = 7)
  
  # Plot average number of episomes per cell over time
  averages <- extinction_3epi_df_long %>% 
    mutate(episomes_per_cell = as.numeric(as.character(episomes_per_cell))) %>% 
    group_by(time = round(time,1), Pr, Ps, trial) %>% 
    summarise(avg = sum(number_of_cells*episomes_per_cell)/sum(number_of_cells)) %>% 
    ungroup()
  
  # fill out Pr = Ps = 100%
  averages <- rbind(averages, 
                    data.frame(Pr = factor("100%", levels = levels(averages$Pr)), 
                               Ps = factor("100%", levels = levels(averages$Ps)), 
                               avg = 3, expand.grid(time = 41:(max(averages$time) + 30), trial = 1:100))
  )
  
  
  p2.1 <- averages %>% 
    ggplot(aes(time, avg, group = trial)) + 
    geom_line(alpha = 0.25, color = "gray") +
    ggh4x::facet_grid2(Pr ~ Ps, labeller = "label_both", scales = "free_x", independent = "x") + 
    labs(x = "time (generations)", y = "average number of episomes per cell") 
  
  ggsave(here("results", "simulations", "extinction_3epi_average_epi_freex.png"), p2.1, width = 8, height = 7)
  
  p2.2 <- averages %>% 
    ggplot(aes(time, avg, group = trial)) + 
    geom_line(alpha = 0.25, color = "gray") +
    facet_grid(Pr ~ Ps, labeller = "label_both") +
    scale_x_log10() +
    labs(x = "time (generations)", y = "average number of episomes per cell") 
  
  ggsave(here("results", "simulations", "extinction_3epi_average_epi_logscale.png"), p2.2, width = 8, height = 7)
  
  # Plot histograms of time to extinction
  time_to_extinction <- 1:nrow(param_grid) %>% 
    map_df(function(i) cbind(param_grid[i,], extinction_3epi[[i]]$ExtinctionTime))
  
  p3.1 <- time_to_extinction %>% 
    rename(Pr = pRep, Ps = pSeg) %>% 
    filter(Ps*Pr != 1) %>%
    arrange(desc(Pr), Ps) %>% 
    mutate(Pr = paste0(Pr*100, "%"),
           Ps = paste0(Ps*100, "%"),
           Pr = fct_inorder(factor(Pr)),
           Ps = fct_inorder(factor(Ps))) %>% 
    ggplot(aes(ExtinctionTime)) + 
    geom_histogram(bins = 30, aes(y = after_stat(density))) + 
    facet_grid(Pr ~ Ps, labeller = "label_both", scales = "free_y") +
    labs(x = "time to episome extinction (generations)")
  
  ggsave(here("results", "simulations", "time_to_exinction_3epi.png"), p3.1, width = 7, height = 7)
  
  p3.2 <- time_to_extinction %>% 
    rename(Pr = pRep, Ps = pSeg) %>% 
    filter(Ps*Pr != 1) %>%
    arrange(desc(Pr), Ps) %>% 
    mutate(Pr = paste0(Pr*100, "%"),
           Ps = paste0(Ps*100, "%"),
           Pr = fct_inorder(factor(Pr)),
           Ps = fct_inorder(factor(Ps))) %>% 
    ggplot(aes(ExtinctionTime)) + 
    geom_histogram(bins = 30, aes(y = after_stat(density))) + 
    ggh4x::facet_grid2(Pr ~ Ps, labeller = "label_both", scales = "free", independent = "x") +
    labs(x = "time to episome extinction (generations)")
  
  ggsave(here("results", "simulations", "time_to_exinction_3epi_freex.png"), p3.2, width = 7, height = 7)
  
  save(extinction_3epi, param_grid, 
       averages, time_to_extinction,
       file = here("results", "simulations", "constant_pop_3epi.RData"))
}


## Plot average episome over time for case with moderate defect
simple_plot <-  averages %>% 
  filter(Pr == "100%", Ps == "100%") %>% 
  ggplot(aes(time, avg, group = trial)) + 
  geom_line(alpha = 0.25, color= "gray") +
  labs(title = "100% Replication, 100% Segregation") +
  coord_cartesian(c(0,700)) +
  
  averages %>% 
  filter(Pr == "100%", Ps == "90%") %>% 
  ggplot(aes(time, avg, group = trial)) + 
  geom_line(alpha = 0.25, color= "gray") + 
  labs(title = "100% Replication, 90% Segregation") +
  coord_cartesian(c(0,700)) +
  
  averages %>% 
  filter(Pr == "90%", Ps == "100%") %>% 
  ggplot(aes(time, avg, group = trial)) + 
  geom_line(alpha = 0.25, color= "gray") + 
  labs(title = "90% Replication, 100% Segregation") +
  coord_cartesian(c(0,120)) +
  
  averages %>% 
  filter(Pr == "90%", Ps == "90%") %>% 
  ggplot(aes(time, avg, group = trial)) + 
  geom_line(alpha = 0.25, color= "gray") + 
  labs(title = "90% Replication, 90% Segregation") +
  coord_cartesian(c(0,120)) +
  
  plot_layout(ncol = 2, byrow = T) &
  # geom_rug(data = . %>% group_by(trial) %>% filter(time == max(time))) &
  ylim(0,5) & 
  theme(plot.margin = margin(15,15,15,15)) &
  labs(x = "Time (generations)", y = "Average number of episomes per cell")


ggsave(here("results", "simulations", "simple_averages.png"), simple_plot, width = 8, height = 8)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")


trajectory_plot <- function(Pr, Ps, xlim){
  PR <- Pr
  PS <- Ps
  
  plot <-  averages %>% 
    filter(Pr == PR, Ps == PS) %>% 
    ggplot(aes(time, avg, group = trial)) + 
    geom_line(alpha = 0.25, color= "gray") 
  
  if(!(PR == "100%" & PS == "100%")){
    temp <- time_to_extinction %>% 
      mutate(Pr = paste0(pRep*100, "%"),
             Ps = paste0(pSeg*100, "%")) %>% 
      filter(Pr == PR, Ps == PS) %>% 
      pull(ExtinctionTime) 
    
    mean_time <- round(mean(temp))
    sd_time <- round(sd(temp))
    
    label <- str_interp("T[E]: ${mean_time} %+-% ${sd_time}")
    
    plot <- plot + 
      geom_text(data = data.frame(NA), aes(mean_time, 0, 
                                           label = label),
                inherit.aes = F, parse = T, vjust = -0.3) +
      geom_point(data = data.frame(NA), aes(mean_time, 0),
                 inherit.aes = F) 
  }
  
  plot +
    labs(title = str_interp("${PR} Replication, ${PS} Segregation")) +
    labs(x = "Time (generations)",y = "Average number of episomes per cell") +
    ylim(c(0,5)) +
    coord_cartesian(c(0,xlim)) 
}

distribution_plot1 <- function(Pr, Ps, xlim){
  PR <- Pr
  PS <- Ps
  
  df <- extinction_3epi_df_long %>% 
    filter(Pr == PR, Ps == PS) %>% 
    group_by(time = round(time), episomes_per_cell) %>% 
    summarise(number_of_cells = mean(number_of_cells)) 
  
  if(PR == "100%" & PS == "100%"){
    df <- df %>% rbind(data.frame(time = 1:750, episomes_per_cell = factor("3", levels = levels(extinction_3epi_df_long$episomes_per_cell)), number_of_cells = 1000))
  } 
  
  df %>% 
    ggplot(aes(time, number_of_cells, fill = episomes_per_cell, color = episomes_per_cell)) + 
    geom_bar(position = position_fill(reverse = TRUE), stat = "identity") + 
    scale_fill_manual(values = safe_colorblind_palette) + 
    scale_color_manual(values = safe_colorblind_palette) + 
    coord_cartesian(c(0,xlim)) +
    scale_y_continuous(labels = scales::percent, breaks = c(0, 0.5, 1))  +
    labs(x = "Time (generations)", fill = "Episomes per cell", y = "Percent of\ncell population") +
    guides(color = "none") + 
    theme(panel.grid = element_blank())
}

new_plot <- trajectory_plot(Pr = "100%", Ps = "100%", 700) +
  
  trajectory_plot(Pr = "100%", Ps = "90%", 700) +
  
  distribution_plot1(Pr = "100%", Ps = "100%", 700) + 
  
  distribution_plot1(Pr = "100%", Ps = "90%", 700) +
  
  trajectory_plot(Pr = "90%", Ps = "100%", 120) +
  
  trajectory_plot(Pr = "90%", Ps = "90%", 120) +
  
  distribution_plot1(Pr = "90%", Ps = "100%", 120) +
  
  distribution_plot1(Pr = "90%", Ps = "90%", 120) +
  
  plot_layout(ncol = 2, byrow = T, heights = c(3,1, 3, 1), guides = "collect") &
  # theme(plot.margin = margin(15,15,15,15)) &
  labs(x = "Time (generations)", fill = "Episomes per cell") &
  guides(color = "none")

ggsave(here("results", "simulations", "new_simulations_fig.png"), new_plot, width = 10, height = 10)

#####
# Plot histograms with distribution of episomes in population at varying time points
#####
distribution_plot2 <- function(generation, PR, PS){
  df <- extinction_3epi_df_long %>% 
    filter(pRep == PR, pSeg == PS) %>% 
    group_by(trial, pRep, pSeg) %>% 
    filter((abs(time - generation)) == min(abs(time - generation))) %>% 
    mutate(episomes_per_cell = as.numeric(as.character(episomes_per_cell))) 
  
  # average <- df %>% group_by(trial) %>% summarise(avg = sum(number_of_cells*episomes_per_cell)/sum(number_of_cells)) 
  
  df %>% 
    group_by(episomes_per_cell, Pr, Ps) %>% 
    summarise(percent_of_cells = mean(frac)*100) %>% 
    ggplot(aes(episomes_per_cell, percent_of_cells)) + 
    geom_bar(stat = "identity") + 
    labs(x = "Number of episomes per cell", y = "Percent of cells with a\ngiven number of episomes",
         title = paste("Division", generation)) + 
    ylim(c(0, 100)) + 
    scale_x_continuous(breaks = seq(0,10, by = 2), limits = c(-1,10))
}

plot_epi_distribution <- function(Pr, Ps){
  plot <- distribution_plot2(1, Pr,Ps) + theme(axis.title.y = element_text(size = 10)) + 
    distribution_plot2(10, Pr,Ps) + 
    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + 
    distribution_plot2(20, Pr,Ps) + 
    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + 
    distribution_plot2(100, Pr,Ps) + 
    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + 
    plot_layout(nrow = 1) &
    # plot_annotation(title = str_interp("${Pr*100}% Replication, ${Ps*100}% Segregation")) &
    theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
          axis.title.x = element_blank(), plot.title = element_text(size = 8)) &
    scale_y_continuous(breaks = c(0, 50, 100), limits = c(0,100)) & 
    labs(y = "Percent of\ncell population")
  
  # # Use the tag label as an x-axis label
  # wrap_elements(panel = plot) +
  #   labs(tag = "Number of episomes per cell") +
  #   theme(
  #     plot.tag = element_text(size = rel(1)),
  #     plot.tag.position = "bottom"
  #   )
  plot
}

plot_epi_distribution(1,0.9)
plot_epi_distribution(0.9,1)
plot_epi_distribution(0.9,0.9)

x_label <- ggdraw() + 
  draw_label(
    "Number of episomes per cell",
    x = 0.5,
    hjust = 0.5,
    size = 10,
    y = 0.8
    # vjust = 1
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(-10, 0, 0, 0),
    plot.background = element_rect(color = "white")
  )
## add epi distributions to "simple plot"
new_plot2 <- plot_grid(
  trajectory_plot(Pr = "100%", Ps = "100%", 700),
  trajectory_plot(Pr = "100%", Ps = "90%", 700),
  plot_epi_distribution(1,1) ,#+ theme(plot.margin = margin(-10,5,-10,5)),
  plot_epi_distribution(1,0.9) ,#+ theme(plot.margin = margin(-10,5,-10,5)),
  x_label, x_label,
  trajectory_plot(Pr = "90%", Ps = "100%", 120),
  trajectory_plot(Pr = "90%", Ps = "90%", 120),
  plot_epi_distribution(0.9,1) ,#+ theme(plot.margin = margin(0,1,0,1)),
  plot_epi_distribution(0.9,0.9) ,#+ theme(plot.margin = margin(0,1,0,1)),
  x_label, x_label,
  
  byrow = T, ncol = 2, rel_heights = c(2,1, 0.1, 2,1, 0.1),
  axis = "lr", align = "v"
)

ggsave(here("results", "simulations", "new_simulations_fig_v2.png"), new_plot2, width = 10, height = 10, bg = "white")


#####
# Supplemental figures of simulation summary
#####

time_to_extinction_summary <- time_to_extinction %>% 
  mutate(Pr = paste0(pRep*100, "%"),
         Ps = paste0(pSeg*100, "%")) %>% 
  filter(!(pRep == 1 & pSeg == 1)) %>% 
  group_by(Pr, Ps, pRep, pSeg) %>% 
  summarise(sd = sd(ExtinctionTime), ExtinctionTime = mean(ExtinctionTime))

sup1 <- time_to_extinction_summary %>% 
  ggplot(aes(pSeg*100, ExtinctionTime, color = as.factor(pRep*100), group = pRep)) + 
  geom_point(size = 2) + geom_line(lwd = 1) + 
  geom_errorbar(aes(ymin = ExtinctionTime - sd, ymax = ExtinctionTime + sd), width = 0.5, lwd = 1) + 
  # add infity for Pr = Ps = 1
  geom_text(data = data.frame(pSeg = 1, pRep = 1, ExtinctionTime = Inf), vjust = 1, aes(label = "+"), size = 7, show.legend = F) +
  scale_color_manual(values = sort(scico_palette_data("acton", T))[c(1, 31, 51, 71, 81)]) + 
  labs(x = "Ps (%)", y = "Time to episome extinction (generations)", color = "Pr (%)") +
  theme(legend.position = c(0,1), legend.justification = c(-0.1,1.1))

sup2 <- time_to_extinction_summary %>% 
  mutate(pRep = ifelse(pSeg == 0.5, pRep - 0.01, 
                       ifelse(pSeg == 0.75, pRep - 0.005, 
                              ifelse(pSeg == 0.95, pRep + 0.005, 
                                     ifelse(pSeg == 1, pRep + 0.01, pRep))))) %>% 
  ggplot(aes(pRep*100, ExtinctionTime, color = as.factor(pSeg*100), group = pSeg)) + 
  geom_line(lwd = 1) +
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin = ExtinctionTime - sd, ymax = ExtinctionTime + sd), width = 0.5, lwd = 1) + 
  # add infity for Pr = Ps = 1
  geom_text(data = data.frame(pSeg = 1, pRep = 1.01, ExtinctionTime = Inf), vjust = 1, aes(label = "+"), size = 7, show.legend = F) +
  # scale_color_manual(values = sort(scico_palette_data("oslo", T))[c(1, 41, 61, 81)]) + 
  scale_color_manual(values = sort(scico_palette_data("oslo", T))[c(1, 31, 51, 71, 91)]) + 
  labs(x = "Pr (%)", y = "Time to episome extinction (generations)", color = "Ps (%)") +
  theme(legend.position = c(0,1), legend.justification = c(-0.1,1.1))


sup3 <- averages %>% 
  group_by(Pr, Ps, time = round(time, 2)) %>%
  summarise(avg = mean(avg)) %>%
  mutate(Ps = as.numeric(gsub("%", "", as.character(Ps)))) %>%
  ggplot(aes(time, avg, color = as.factor(Ps), group = Ps)) + 
  geom_line(size = 1) +
  facet_wrap(~Pr, ncol = 1, labeller = "label_both") +
  scale_color_manual(values = sort(scico_palette_data("oslo", T))[c(1, 31, 51, 71, 91)]) +
  labs(x = "Time (generations)", y = "Average number of episomes per cell", color = "Ps (%)") + 
  # add time to extinction annotation
  geom_point(data = time_to_extinction_summary %>%
               mutate(Ps = as.numeric(gsub("%", "", as.character(Ps))), 
                      Pr = factor(Pr, levels = levels(averages$Pr))), 
             aes(ExtinctionTime, 0), shape = 4, show.legend = F) +
  geom_rect(aes(xmin = 1850, xmax = 2250, ymin = 0.5, ymax = 3.2), fill = "white", alpha = 0.05, color = NA) +
  geom_text(data = time_to_extinction_summary %>%
              mutate(Ps = as.numeric(gsub("%", "", as.character(Ps))), 
                     Pr = factor(Pr, levels = levels(averages$Pr)),
                     y = ifelse(Ps == 100, 3, ifelse(Ps == 95, 2.5, ifelse(Ps == 90, 2, ifelse(Ps == 75, 1.5, 1))))), 
            aes(1850, y, label = paste("T[E]:", round(ExtinctionTime), "%+-%", round(sd))), 
            parse = T, hjust = 0, vjust = 1, 
            show.legend = F) 

sup <- ((plot_spacer() / sup1 / plot_spacer() / sup2 / plot_spacer() + plot_layout(heights = c(0.5,4,0.5,4,0.5))) | 
  plot_spacer() | 
  sup3) + plot_layout(widths = c(1,0.1,1))

ggsave(here("results", "simulations", "supplemental_figure.png"), sup , width = 10.5, height = 10)


sup3.2 <- averages %>% 
  group_by(Pr, Ps, time = round(time, 2)) %>%
  summarise(avg = mean(avg)) %>%
  mutate(Ps = as.numeric(gsub("%", "", as.character(Ps))),
         Pr = fct_rev(Pr)) %>%
  ggplot(aes(time, avg, color = as.factor(Ps), group = Ps)) + 
  geom_line(size = 1) +
  facet_wrap(~Pr, nrow = 1, labeller = "label_both") +
  scale_color_manual(values = sort(scico_palette_data("oslo", T))[c(1, 31, 51, 71, 91)]) +
  labs(x = "Time (generations)", y = "Average number of episomes per cell", color = "Ps (%)") + 
  # add time to extinction annotation
  geom_point(data = time_to_extinction_summary %>%
               mutate(Ps = as.numeric(gsub("%", "", as.character(Ps))), 
                      Pr = factor(Pr, levels = levels(averages$Pr))), 
             aes(ExtinctionTime, 0), shape = 4, show.legend = F) +
  geom_rect(aes(xmin = 1850, xmax = 2250, ymin = 0.5, ymax = 3.2), fill = "white", alpha = 0.05, color = NA) +
  geom_text(data = time_to_extinction_summary %>%
              mutate(Ps = as.numeric(gsub("%", "", as.character(Ps))), 
                     Pr = factor(Pr, levels = levels(averages$Pr)),
                     y = ifelse(Ps == 100, 3, ifelse(Ps == 95, 2.5, ifelse(Ps == 90, 2, ifelse(Ps == 75, 1.5, 1))))), 
            aes(1850, y, label = paste("T[E]:", round(ExtinctionTime), "%+-%", round(sd))), 
            parse = T, hjust = 0, vjust = 1, 
            show.legend = F)  + 
  theme(legend.position = "bottom")

sup.2 <-( (sup1 + plot_spacer() + sup2 + plot_layout(widths = c(1, 0.1, 1), nrow = 1)) / 
  plot_spacer() / sup3.2) + plot_layout(heights = c(1,0.1,1), ncol = 1)

ggsave(here("results", "simulations", "supplemental_figure_v2.png"), sup.2 , width = 15, height = 7)


## Calculate shannon index
shannon <- extinction_3epi_df_long %>% 
  filter(number_of_cells != 0) %>% 
  mutate(p = number_of_cells/total_cells) %>% 
  group_by(pRep, pSeg, trial, time) %>% 
  summarise(H = -sum(p*log(p)))

shannon %>% 
  mutate(Pr = pRep*100, Ps = pSeg*100) %>% 
  filter(Pr == 100, Ps == 95) %>% 
  group_by(Pr, Ps, time = round(time, 2)) %>% 
  summarise(avg = mean(H)) %>% head
  ggplot(aes(time, avg, color = as.factor(Ps), group = Ps)) + 
  geom_line() +
  # facet_wrap(~Ps, ncol = 1, labeller = "label_both") +
  scale_color_manual(values = sort(scico_palette_data("oslo", T))[c(1, 31, 51, 71, 91)]) +
  labs(x = "Time (generations)", y = "Average shannon index", color = "Ps (%)") 


