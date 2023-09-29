extinction <- function(pRep, pSeg, nTrials, n_epi){
  
  results = 1:nTrials*0
  for(z in 1:nTrials){
    print(z)
    time = 0
    deathVec = rep(1, 10)
    cells = rep(0, 10)
    cells[n_epi + 1] <- 1000
    total = sum((0:9)*cells)
    while(total >0){
      birthVec = rep(2*(1-sum(cells)/1000) + 1, 10)
      result = simStepFlex(pRep, pSeg, cells, birthVec, deathVec, FALSE, max_epi = 9)
      cells = result[[2]]
      time=time + result[[1]]
      #print(cells)
      total = sum((0:9)*cells)
    }
    results[z]= time
  }
  return(tibble(ExtinctionTime = results))
}


pRep = 0.77
pSeg = 2*0.91 - 1 # convert segregation probability to probability of a mechanism for segregation
extTimes3epi <- extinction(pRep, pSeg, 500, 3)
kludge3 = data.frame(stats = c(mean(extTimes3epi$ExtinctionTime), median(extTimes3epi$ExtinctionTime)), offset = c(1.2, -1.2))

extTimes1epi <- extinction(pRep, pSeg, 500, 1)
kludge1 = data.frame(stats = c(mean(extTimes1epi$ExtinctionTime), median(extTimes1epi$ExtinctionTime)), offset = c(1.2, -1.2))

pdf("extinctiontime.pdf", width=10, height=5)
ggplot(data=extTimes3epi) + geom_histogram(aes(x=ExtinctionTime, y=after_stat(density))) + 
  geom_vline(data=kludge3, aes(xintercept=stats, color=factor(stats)))+
  geom_text(data=kludge3, aes(x=stats+offset,y=.02, label=round(stats, digits = 1), color = factor(stats)), size=4, angle = 90, show.legend = FALSE)+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  ggtitle("3 Episome Extinction")+
  xlab("Time (generations)")+ 
  ylab("Density")+
  scale_color_manual(values = c("#A6CEE3", "#66CC99"), labels=c("Median","Mean"), name = "Statistics")+
  theme(text = element_text(size=20)) + 

ggplot(data=extTimes1epi) + geom_histogram(aes(x=ExtinctionTime, y=after_stat(density))) + 
  geom_vline(data=kludge1, aes(xintercept=stats, color=factor(stats)))+
  geom_text(data=kludge1, aes(x=stats+offset,y=.02, label=round(stats, digits = 1), color = factor(stats)), size=4, angle = 90, show.legend = FALSE)+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  ggtitle("1 Episome Extinction")+
  xlab("Time (generations)")+ 
  ylab("Density")+
  scale_color_manual(values = c("#A6CEE3", "#66CC99"), labels=c("Median","Mean"), name = "Statistics")+
  theme(text = element_text(size=20)) + 
  
  plot_layout(ncol = 2)

dev.off()


## Perfect replication and segregation

## Observed parameters

## Perfect for one and observed for the other

## Perfect for one and betterr than observed for  one


## Figures:

# as a function of generations, what is the total number of episomes
# distribution of number of episomes in cells over time


## time  verses total number of episomes (one faint line per  simulation so density adds up)

#interested in relative importance of replication and segregation (definite decrease in pop versus increased variance)



