require(ggplot2)
require(RColorBrewer)

makeChildren<- function(pRep, pSeg, numEpisomes){
	if(numEpisomes == 0){
		return(c(0, 0))
	}
	else{
		outcomes = c(1, 0, 0, 1, 1, 1, 2, 0, 0, 2)
		refMat = matrix(outcomes, nrow= 5, byrow=T) #possible episome fates
		# probVec = c(1/2*(1-pRep), 1/2*(1-pRep), pRep*(pSeg+1/2*(1-pSeg)), pRep*1/4*(1-pSeg), pRep*1/4*(1-pSeg)) # probabilities when pSeg = probability of tethering
		probVec = c(1/2*(1-pRep), 1/2*(1-pRep), pRep*pSeg, 1/2*pRep*(1-pSeg), 1/2*pRep*(1-pSeg)) # probabilities when pSeg = probability of segregation regardless of mechanism
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


expoDf <- function(pRep, pSeg, nIts, nRuns, max_epi){
	#try using matrix first and then convert to data frame
# rm(data)
dataCUT = 1
data <- matrix(NA, nrow=11*(nIts/dataCUT+1)*nRuns, ncol=4)
names(data) = c("run", "time", "episomes", "count")
z = 1
for(run in 1:nRuns){
	print(run)
	deathVec = rep(0, 10)
	cells = c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0)
	times = rep(0, nIts + 1)
	totalEps = append(c(sum((0:9)*cells)), rep(0, nIts))
	for(j in 1:length(cells)){
			tmp = c(run, times[1], j-1, cells[j])
			data[z,] <- tmp
			z=z+1
		}

		tmp = c(run, times[1], -1, sum(cells))
		data[z,] <- tmp
		z= z+1
	#data <- rbind(data, data.frame(run=run, time=times[1], episomes=-1, count = sum((0:9)*cells)))
	for(i in 1:nIts){

		birthVec = rep(1, 10)
		result = simStepFlex(pRep, pSeg, cells, birthVec, deathVec, max_epi = max_epi)
		cells = result[[2]]
		times[i+1] = times[i] + result[[1]]
		#print(cells)
		totalEps[i+1] = sum((0:9)*cells)
		if(i%%1000 == 0){
		#	print(i)
		}
		if(i%%dataCUT == 0){
			for(j in 1:length(cells)){
				tmp = c(run, times[i+1], j-1, cells[j]) #added sum cells
				data[z,] <- tmp
				z = z+1
			}
			tmp = c(run, times[i+1], -1, sum(cells))#totalEps[i+1])
			data[z,] <- tmp
			z= z+1
			#data <- rbind(data, data.frame(run=run, time=times[i+1], episomes=-1, count = sum((0:9)*cells)))
		}
	}
	
#plot(times, totalEps)
}
return(data)
}



data<- expoDf(1,.8,100, 2, 9)
test <- data.frame(data)

label_list <- list(
	"-1"="Total",
	"0"="0",
	"1"="1",
	"2"="2",
	"3"="3",
	"4"="4",
	"5"="5",
	"6"="6",
	"7"="7",
	"8"="8",
	"9"="9"
	)

labels <- function(variable, value){
	return(label_list[toString(value)])
}

colors = c("black", brewer.pal(n=10,name="Paired"))
labs = c("total", 0:9)
names(test) <- c("run", "time", "episomes", "count")
ggplot(data = test, aes(x = time, y = log10(count), group = episomes, color=factor(episomes)))+geom_point(size=1) + 
theme(legend.title = element_text(colour="black", size=12, face="bold"))+ facet_wrap(~episomes)+
xlab("Time (generations)") + 
  scale_color_manual(values = colors, labels = labs, name = "k-episome cells, k = ...")+ggtitle("pRep = 100%, pSeg =  80%, 3 episome exponential growth")



expoDf2 <- function(pRep, pSeg, nIts, nRuns, selection, max_epi){
	#try using matrix first and then convert to data frame
# rm(data)
dataCUT = 1
start = 1:100
mid = 11:100*10
end = 2:100*1000
recordList = c(start, mid, end)
data <- matrix(NA, nrow=2*(max_epi + 2)*(length(recordList)+1)*nRuns, ncol=5)
names(data) = c("run", "time", "episomes", "frac", "total")
z = 1
for(run in 1:nRuns){
	print(run)
	deathVec = rep(0, 10)
	cells = c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0)
	times = rep(0, nIts + 1)
	totalEps = append(c(sum((0:9)*cells)), rep(0, nIts))
	for(j in 1:length(cells)){
			tmp = c(run, times[1], j-1, cells[j]/sum(cells), sum(cells))
			data[z,] <- tmp
			z=z+1
		}

		tmp = c(run, times[1], -1, sum((0:9)*cells)/sum(cells), sum(cells))
		data[z,] <- tmp
		z= z+1
	#data <- rbind(data, data.frame(run=run, time=times[1], episomes=-1, count = sum((0:9)*cells)))
	for(i in 1:nIts){

		birthVec = rep(1, 10)
		result = simStepFlex(pRep, pSeg, cells, birthVec, deathVec, selection, max_epi)
		cells = result[[2]]
		times[i+1] = times[i] + result[[1]]
		#print(cells)
		totalEps[i+1] = sum((0:9)*cells)
		if(i%%1000 == 0){
		#	print(i)
		}
		if(sum(cells) %in% recordList){
			#print(i)
			for(j in 1:length(cells)){
				tmp = c(run, times[i+1], j-1, cells[j]/sum(cells), sum(cells)) #added sum cells
				data[z,] <- tmp
				z = z+1
			}
			tmp = c(run, times[i+1], -1, sum((0:9)*cells)/sum(cells), sum(cells)) #totalEps[i+1])
			data[z,] <- tmp
			z= z+1
			#data <- rbind(data, data.frame(run=run, time=times[i+1], episomes=-1, count = sum((0:9)*cells)))
		}
	}
	
#plot(times, totalEps)
}
return(data)
}


selection.against.zero=TRUE
max.epi = 9
# dataSelection<- expoDf2(.824,.882,120000, 200, selection.against.zero)
dataSelection<- exponential_growth(.824,.882,120000, 4, 1, 3, selection.against.zero, max.epi)
selection.against.zero=FALSE 
dataNoSelection <- expoDf2(.824, .882, 120000, 200, selection.against.zero)

dataPerfectReplication3 <- expoDf3(1, 0, 120000, 200, TRUE)
dataPerfectReplication4 <- expoDf3(.824, 0, 150000, 200, TRUE)
dataPerfectSegregation <- expoDf2(.824, 1, 120000, 200, TRUE)


#DO THIS FIX FIRST
test <- data.frame(dataSelection)


colors = c("black", brewer.pal(n=10,name="Paired"), brewer.pal(n=10,name="Set3"))
labs = c("average", 0:19)
names(test) <- c("run", "time", "episomes", "frac", "total")

ggplot(data = test, aes(x=total, y = log10(frac+1/10000), group=episomes, color=factor(episomes))) + geom_point(size=1)+
theme(legend.title = element_text(colour="black", size=12, face="bold"))+
xlab("Size") +
 scale_color_manual(values = colors, labels = labs, name = "k-episome cells, k = ...")+ggtitle("pRep = 100%, pSeg =  80%, 3 episome exponential growth") +
  scale_x_log10()
  # facet_wrap(~run)

require(Hmisc)

pdf("bad_segregation_bad_replication_final.pdf", width=12, height = 7.5)
ggplot(data = test, aes(x=log10(total), y = log10(frac+1/100000), group=episomes, color=factor(episomes))) + 
stat_summary(fun.data="mean_cl_normal", geom="smooth")+
theme(legend.title = element_text(colour="black", size=12, face="bold"))+
xlab("Log10(Population Size)") + ylab("Log10(Fraction)")+
  scale_color_manual(values = colors, labels = labs, name = "k-episome cells, k = ...")+ggtitle("pRep = 82.4%, pSeg =  50%, 3 episome exponential growth with selection")

dev.off()

#extract data and make histograms

hist_sizes = c(1, 10, 100, 1000, 10000, 100000)

temp <- subset(test, total %in% hist_sizes)
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
	#tmp = data.frame(episomes = 0, total = paste("Size:", as.integer(j)), frac = 0, color = cols2[1])
	#kludge <- rbind(kludge, tmp)
}

pdf("bad_segregation_bad_replication_histogram.pdf", width = 12, height = 7.5)
ggplot(data = kludge, aes(x = episomes, y = frac, group = total, fill = factor(episomes)))+ geom_bar(stat="identity")+
facet_wrap(~total) + scale_x_discrete(breaks = 0:9) + ylab("Fraction of population") + xlab("Number of episomes") +
ggtitle("pRep = 82.4%, pSeg =  50%, 3 episome exponential growth with selection")+ 
scale_fill_manual(values = cols2, labels = 0:9, name = "k-episome cells, k = ....")
dev.off()
