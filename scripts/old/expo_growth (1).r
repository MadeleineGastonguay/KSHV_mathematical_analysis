require(ggplot2)
require(RColorBrewer)
require(Hmisc)
require(dplyr)
require(tidyverse)
require(patchwork)
theme_set(theme_bw())

#####
# Simulation functions
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

# Function to simulate exponential growth of a cell population with or without selection
expoSims <- function(pRep, pSeg, nIts, nRuns, selection, max_epi){
  # time points at which we want to record the number of each k-episome cell
  start = 1:100
  mid = 11:100*10
  end = 2:100*1000
  recordList = c(start, mid, end)
  
  # initialize matrix to collect results
  data <- matrix(NA,  ncol=5)
  colnames(data) = c("run", "time", "episomes", "frac", "total")
  
  for(run in 1:nRuns){ # iterate through runs
    print(run)
    deathVec = rep(0, 10) # define vector of death rates (all 0 in this case because cells don't die)
    # initial number of each k-episome cell
    cells = rep(0, max_epi + 1)
    cells[4] <- 1
    # number of time points to simulate 
    times = 0
    
    for(j in 1:length(cells)){
      # for each k-episome cell type, calculate the fraction of cells with in that category
      tmp = c(run, times[1], j-1, cells[j]/sum(cells), sum(cells))
      data <- rbind(data, tmp)
    }
    # calculate the average number of episomes per cell
    tmp = c(run, times[1], -1, sum((0:9)*cells)/sum(cells), sum(cells))
    data <- rbind(data, tmp)
    
    loop <- TRUE
    i = 1
    while(loop){
      # print(i)
      # for each time step, simulate a replication event
      birthVec = rep(1, 10)
      result = simStepFlex(pRep, pSeg, cells, birthVec, deathVec, selection, max_epi)
      # record the new number of each k-episome cell type
      cells = result[[2]]
      # record the new time
      times[i+1] = times[i] + result[[1]]
      # if(i%%1000 == 0){
      #   print(i)
      # }
      # if the number of cells is a point we want to keep track of:
      if(sum(cells) %in% recordList){
        #print(i)
        for(j in 1:length(cells)){
          # record the fraction of each k-episome cell type
          tmp = c(run, times[i+1], j-1, cells[j]/sum(cells), sum(cells)) #added sum cells
          data <- rbind(data, tmp)
        }
        # record average number of episomes per cell
        tmp = c(run, times[i+1], -1, sum((0:9)*cells)/sum(cells), sum(cells)) #totalEps[i+1])
        data <- rbind(data, tmp)
      }
      total = sum(cells)
      if(i > nIts) {
        loop <- FALSE 
      }else{
        loop <- total < 1e5
      }
      i <- i + 1
    }
    
  }
  data <- data[-1,]
  return(as_tibble(data))
}

#####
# plotting functions

plot_fraction <- function(sims, pRep = NA, pSeg = NA, selection = NA){
  colors = c("black", brewer.pal(n=10,name="Paired"), brewer.pal(n=10,name="Set3"))
  labs = c("average", 0:19)
  
  title <- ifelse(is.na(pRep), "", 
                  str_interp("pRep = ${pRep*100}%, pSeg =  ${pSeg*100}%, 3 episome exponential growth ${ifelse(selection, 'with', 'without')} selection"))
  
  ggplot(data = sims, aes(x=log10(total), y = log10(frac+1/100000), group=episomes, color=factor(episomes))) + 
    stat_summary(fun.data="mean_cl_normal", geom="smooth", se = T)+
    # stat_summary(fun.data="mean_cl_normal", geom="smooth", se = F)+
    theme(legend.title = element_text(colour="black", size=12, face="bold"))+
    xlab("Log10(Population Size)") + ylab("Log10(Fraction)")+
    scale_color_manual(values = colors, labels = labs, name = "k-episome cells, k = ...")+
    ggtitle(title)
}

plot_histogram <- function(sims, pRep = NA, pSeg = NA, selection = NA){
  
  hist_sizes = c(1, 10, 100, 1000, 10000, 100000)
  title <- ifelse(is.na(pRep), "", 
                  str_interp("pRep = ${pRep*100}%, pSeg =  ${pSeg*100}%, 3 episome exponential growth ${ifelse(selection, 'with', 'without')} selection"))
  browser()
  temp <- subset(sims, total %in% hist_sizes)
  temp2 <- subset(temp, episomes != -1)
  for(i in hist_sizes){
    
    temp2$frac[temp2$total == i] <- (temp2$frac[temp2$total==i]+1e-4)/sum(temp2$total == i)*10
  }
  
  cols2 = c(brewer.pal(n=10,name="Paired"))#, brewer.pal(n=10,name="Set3"))
  epis = 0:9
  kludge <- data.frame(NULL)
  for(i in epis){
    for(j in hist_sizes){
      t1 = temp2[temp2$total == j,]
      t2 = t1$frac[t1$episomes == i]
      avg = sum(t2)
      tmp = data.frame(episomes = i, total = paste("Population Size:", as.integer(j)), frac = avg, color = cols2[i+1])
      kludge <- rbind(kludge, tmp)
    }
  }
  
  ggplot(data = kludge, aes(x = episomes, y = frac, group = total, fill = factor(episomes)))+ geom_bar(stat="identity")+
    facet_wrap(~total) + scale_x_discrete(breaks = 0:9) + ylab("Fraction of population") + xlab("Number of episomes") +
    ggtitle(title) + 
    scale_fill_manual(values = cols2, labels = 0:9, name = "k-episome cells, k = ....")
  
}

plot_histogram2 <- function(sims, hist_size, pRep = NA, pSeg = NA, selection = NA){
  
  # hist_sizes = c(1, 10, 100, 1000, 10000, 100000)
  title <- ifelse(is.na(pRep), "", 
                  str_interp("pRep = ${pRep*100}%, pSeg =  ${pSeg*100}%, 3 episome exponential growth ${ifelse(selection, 'with', 'without')} selection"))
  cols2 = c(brewer.pal(n=10,name="Paired"))#, brewer.pal(n=10,name="Set3"))
  
  kludge <- sims %>% filter(total %in% hist_size, episomes != -1) %>% 
    group_by(Pr, Ps, total) %>% 
    mutate(frac = (frac + 1e-4)/n()*10) %>% 
    group_by(Pr, Ps, total, episomes) %>% 
    summarise(frac = sum(frac)) %>% 
    mutate(total = paste("Population Size:", as.integer(total)))
  
  ggplot(data = kludge, aes(x = episomes, y = frac, group = total, fill = factor(episomes)))+ geom_bar(stat="identity")+
    scale_x_discrete(breaks = 0:9) + ylab("Fraction of population") + xlab("Number of episomes") +
    ggtitle(title) + 
    scale_fill_manual(values = cols2, labels = 0:9, name = "k-episome cells, k = ....")
  
}

#####
# run simulations

nIts <- 150000

# start with MLE estimates with and without selection
pRep = 0.77
pSeg = 2*0.91 - 1 # convert segregation probability to probability of a mechanism for segregation

sims_with_MLE_parameters <- expoSims(pRep,pSeg,nIts, 200, TRUE, max_epi = 9)
sims_with_MLE_parameters_noSelection <- expoSims(pRep,pSeg, nIts, 200, FALSE, max_epi = 9)

save(sims_with_MLE_parameters, sims_with_MLE_parameters_noSelection, 
     file = "exponential_simulations.RData")

nIts <- 1e5
# keep MLE estimate of Pr and try upper and lower bounds of Ps
pRep = 0.77
pSeg = 2*1 - 1 # convert segregation probability to probability of a mechanism for segregation
sims2_noSelection <- expoSims(pRep,pSeg, nIts, 200, FALSE, max_epi = 9)

save(sims_with_MLE_parameters, sims_with_MLE_parameters_noSelection,
     sims2_noSelection,
      file = "exponential_simulations.RData")

# keep MLE estimate of Pr and try upper and lower bounds of Ps
pRep = 0.77
pSeg = 2*0.77-1 # convert segregation probability to probability of a mechanism for segregation

# sims3 <- expoSims(pRep,pSeg,nIts, 200, TRUE, max_epi = 9)
sims3_noSelection <- expoSims(pRep,pSeg,nIts, 200, FALSE, max_epi = 9)

save(sims_with_MLE_parameters, sims_with_MLE_parameters_noSelection,
     sims2_noSelection,
     sims3_noSelection,
     file = "exponential_simulations.RData")

# keep MLE estimate of Ps and try upper and lower bounds of Pr
pRep = 0.54
pSeg = 2*0.91-1 # convert segregation probability to probability of a mechanism for segregation

# sims4 <- expoSims(pRep,pSeg,nIts, 200, TRUE, max_epi = 9)
sims4_noSelection <- expoSims(pRep,pSeg, nIts, 200, FALSE, max_epi = 9)

save(sims_with_MLE_parameters, sims_with_MLE_parameters_noSelection,
     sims2_noSelection,
     sims3_noSelection,
     sims4_noSelection,
     file = "exponential_simulations.RData")

# keep MLE estimate of Ps and try upper and lower bounds of Pr
pRep = 0.89
pSeg = 2*0.91-1 # convert segregation probability to probability of a mechanism for segregation

# sims5 <- expoSims(pRep,pSeg, nIts, 200, TRUE, max_epi = 9)
sims5_noSelection <- expoSims(pRep,pSeg, nIts, 200, FALSE, max_epi = 9)

save(sims_with_MLE_parameters, sims_with_MLE_parameters_noSelection,
     sims2_noSelection,
     sims3_noSelection,
     sims4_noSelection,
     sims5_noSelection,
     file = "exponential_simulations.RData")

# upper bounds of both parameters
pRep = 0.89
pSeg = 2*1-1 # convert segregation probability to probability of a mechanism for segregation

# sims6 <- expoSims(pRep,pSeg, nIts, 200, TRUE, max_epi = 9)
sims6_noSelection <- expoSims(pRep,pSeg, nIts, 200, FALSE, max_epi = 9)

save(sims_with_MLE_parameters, sims_with_MLE_parameters_noSelection,
     sims2_noSelection,
     sims3_noSelection,
     sims4_noSelection,
     sims5_noSelection,
     sims6_noSelection,
     file = "exponential_simulations.RData")

# lower bounds of both parameters
pRep = 0.54
pSeg = 2*0.77-1 # convert segregation probability to probability of a mechanism for segregation

# sims7 <- expoSims(pRep,pSeg, nIts, 200, TRUE, max_epi = 9)
sims7_noSelection <- expoSims(pRep,pSeg, nIts,200, FALSE, max_epi = 9)

save(sims_with_MLE_parameters, sims_with_MLE_parameters_noSelection,
     sims2_noSelection,
     sims3_noSelection,
     sims4_noSelection,
     sims5_noSelection,
     sims6_noSelection,
     sims7_noSelection,
     file = "exponential_simulations.RData")

# lower bound of Pr and upper bound of Ps
pRep = 0.54
pSeg = 2*1-1 # convert segregation probability to probability of a mechanism for segregation

# sims8 <- expoSims(pRep,pSeg, nIts, 200, TRUE, max_epi = 9)
sims8_noSelection <- expoSims(pRep,pSeg,nIts, 200, FALSE, max_epi = 9)

save(sims_with_MLE_parameters, sims_with_MLE_parameters_noSelection,
     sims2_noSelection,
     sims3_noSelection,
     sims4_noSelection,
     sims5_noSelection,
     sims6_noSelection,
     sims7_noSelection,
     sims8_noSelection,
     file = "exponential_simulations.RData")

# lower bound of Ps and upper bound of Pr
pRep = 0.89
pSeg = 2*0.77-1 # convert segregation probability to probability of a mechanism for segregation

# sims9 <- expoSims(pRep,pSeg,nIts, 200, TRUE, max_epi = 9)
sims9_noSelection <- expoSims(pRep,pSeg,nIts, 200, FALSE, max_epi = 9)

save(sims_with_MLE_parameters, sims_with_MLE_parameters_noSelection,
     sims2_noSelection,
     sims3_noSelection,
     sims4_noSelection,
     sims5_noSelection,
     sims6_noSelection,
     sims7_noSelection,
     sims8_noSelection,
     sims9_noSelection,
     file = "exponential_simulations.RData")

exponential_sims_noSelection <- rbind(sims_with_MLE_parameters_noSelection %>% mutate(Pr = 0.77, Ps = 0.91),
      sims2_noSelection %>% mutate(Pr = 0.77, Ps = 1),
      sims3_noSelection %>% mutate(Pr = 0.77, Ps = 0.77),
      sims4_noSelection %>% mutate(Pr = 0.54, Ps = 0.91),
      sims5_noSelection %>% mutate(Pr = 0.89, Ps = 0.91),
      sims6_noSelection %>% mutate(Pr = 0.89, Ps = 1),
      sims7_noSelection %>% mutate(Pr = 0.54, Ps = 0.77),
      sims8_noSelection %>% mutate(Pr = 0.54, Ps = 1),
      sims9_noSelection %>% mutate(Pr = 0.89, Ps = 0.77)
)


#####
# figures
pdf("exponential_growth_simulations.pdf", width = 20, height = 10)
plot_fraction(exponential_sims_noSelection) + 
  facet_grid(Pr ~ Ps, labeller= "label_both") + 
  ggtitle("Exponential Growth With Varying Replication and Segregation Efficiencies")
dev.off()

pdf("exponential_growth_simulations_histogram.pdf", width = 20, height = 10)
plot_histogram2(exponential_sims_noSelection %>% filter(total == 1e5), hist_size = 1e5) + 
  facet_grid(Pr ~ Ps, labeller= "label_both") + 
  ggtitle("Long Term Population State for Exponential Growth With Varying Replication and Segregation Efficiencies")
dev.off()

