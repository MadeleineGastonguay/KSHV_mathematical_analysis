## This script runs and plots simulations of a cell population with constant size

library(tidyverse)
library(here)
library(patchwork)
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

## set up a grid of parameters to simulate with and collect all results:
pReps = c(0.5, 0.75, 0.9, 1)
pSegs = c(0.5, 0.75, 0.9, 1)
param_grid <- expand_grid(pRep = pReps, pSeg = pSegs)

ntrials <- 100
n_epi <- 3 

extinction_3epi <- param_grid %>% 
  pmap(extinction, nTrials = ntrials, n_epi = n_epi)

save(extinction_3epi, file = here("results", "simulations", "constant_pop_3epi.RData"))

extinction_3epi_df <- 1:nrow(param_grid) %>% 
  map_df(function(i) cbind(param_grid[i,], extinction_3epi[[i]]$Totals)) 

extinction_3epi_df_long <- pivot_extinction(extinction_3epi_df) %>% 
  group_by(time, trial, Pr = pRep, Ps = pSeg) %>% 
  mutate(total_cells = sum(number_of_cells)) %>% 
  ungroup %>% mutate(frac = number_of_cells/total_cells)


# Plot distribution of episomes per cell over time
p1 <- extinction_3epi_df_long %>% 
  arrange((desc(Pr))) %>% mutate(Pr = fct_inorder(factor(Pr))) %>% 
  ggplot(aes(time, frac, color = episomes_per_cell, group = interaction(trial, episomes_per_cell))) + 
  geom_line(alpha = 0.25) +  
  guides(color = guide_legend(override.aes = list(alpha = 1))) + 
  ggh4x::facet_grid2(Pr ~ Ps, labeller = "label_both", scales = "free_x", independent = "x") + 
  labs(x = "time (generations)", y = "fraction of cell population", color = "episomes per cell")

ggsave(here("results", "simulations", "extinction_3epi_lines.png"), p1, width = 8, height = 7)

# Plot average number of episomes per cell over time
averages <- extinction_3epi_df_long %>% 
  mutate(episomes_per_cell = as.numeric(as.character(episomes_per_cell))) %>% 
  group_by(time = round(time,1), Pr, Ps, trial) %>% 
  summarise(avg = sum(number_of_cells*episomes_per_cell)/sum(number_of_cells))

p2.1 <- averages %>% 
  ungroup %>% 
  arrange(desc(Pr)) %>% mutate(Pr = fct_inorder(factor(Pr))) %>% 
  ggplot(aes(time, avg, group = trial)) + 
  geom_line(alpha = 0.25, color = "gray") +
  ggh4x::facet_grid2(Pr ~ Ps, labeller = "label_both", scales = "free_x", independent = "x") + 
  labs(x = "time (generations)", y = "average number of episomes per cell") 

ggsave(here("results", "simulations", "extinction_3epi_average_epi_freex.png"), p2.1, width = 8, height = 7)

p2.2 <- averages %>% 
  ungroup %>% 
  arrange(desc(Pr)) %>% mutate(Pr = fct_inorder(factor(Pr))) %>% 
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
  arrange((desc(Pr))) %>% mutate(Pr = fct_inorder(factor(Pr))) %>% 
  ggplot(aes(ExtinctionTime)) + 
  geom_histogram(bins = 30, aes(y = after_stat(density))) + 
  facet_grid(Pr ~ Ps, labeller = "label_both", scales = "free_y") +
  labs(x = "time to episome extinction (generations)")

ggsave(here("results", "simulations", "time_to_exinction_3epi.png"), p3.1, width = 7, height = 7)

p3.2 <- time_to_extinction %>% 
  rename(Pr = pRep, Ps = pSeg) %>% 
  filter(Ps*Pr != 1) %>%
  arrange((desc(Pr))) %>% mutate(Pr = fct_inorder(factor(Pr))) %>% 
  ggplot(aes(ExtinctionTime)) + 
  geom_histogram(bins = 30, aes(y = after_stat(density))) + 
  ggh4x::facet_grid2(Pr ~ Ps, labeller = "label_both", scales = "free", independent = "x") +
  labs(x = "time to episome extinction (generations)")

ggsave(here("results", "simulations", "time_to_exinction_3epi_freex.png"), p3.2, width = 7, height = 7)


## Plot average episome over time for case with moderate defect
simple_plot <-  averages %>% 
  filter(Pr == 1, Ps == 1) %>% 
  ggplot(aes(time, avg, group = trial)) + 
  geom_line(alpha = 0.25, color= "gray") +
  labs(x = "time (generations)", y = "average number of episomes per cell",
       title = "100% Replication, 100% Segregation") +
  
  averages %>% 
  filter(Pr == 1, Ps == 0.9) %>% 
  ggplot(aes(time, avg, group = trial)) + 
  geom_line(alpha = 0.25, color= "gray") + 
  labs(x = "time (generations)", y = "average number of episomes per cell",
       title = "100% Replication, 90% Segregation") +
  
  
  averages %>% 
  filter(Pr == 0.9, Ps == 1) %>% 
  ggplot(aes(time, avg, group = trial)) + 
  geom_line(alpha = 0.25, color= "gray") + 
  labs(x = "time (generations)", y = "average number of episomes per cell",
       title = "90% Replication, 100% Segregation") +
  
  averages %>% 
  filter(Pr == 0.9, Ps == 0.9) %>% 
  ggplot(aes(time, avg, group = trial)) + 
  geom_line(alpha = 0.25, color= "gray") + 
  labs(x = "time (generations)", y = "average number of episomes per cell",
       title = "90% Replication, 90% Segregation") +

  plot_layout(ncol = 2, byrow = T) &
  # geom_rug(data = . %>% group_by(trial) %>% filter(time == max(time))) &
  ylim(0,5) &
  theme(plot.margin = margin(15,15,15,15))


ggsave(here("results", "simulations", "simple_averages.png"), simple_plot, width = 8, height = 8)

## Observed parameters

## Perfect for one and observed for the other

## Perfect for one and better than observed for  one


## Figures:

## time  verses total number of episomes (one faint line per  simulation so density adds up)




