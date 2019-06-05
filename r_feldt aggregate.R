#setwd("C:\\Users\\James\\Dropbox\\Meta-analysis\\d-to-r conversion\\simulations and code")
#setwd("C:\\Users\\jep701\\Dropbox\\Meta-analysis\\d-to-r conversion\\simulations and code")
source("r_feldt simulation functions.R")
source("parameter values.R")
library(ggplot2)


# prepare array to store results

load(paste("d2r ",design[1]," fix=",fixed.cutoff[1]," p",1/p1[1]," n",n.T[1],".rData", sep=""))

allResults <- array(NA, dim = c(length(fixed.cutoff),length(design),length(p1), length(n.T), dim(results)), 
			dimnames = c(fixed=list(fixed.cutoff), design=list(design), p1=list(1/p1), n=list(n.T), dimnames(results)))


# read in results from separate files

for (f in 1:2)
  for (d in 1:2)
    for (p in 1:5)
      for (j in 1:5) {
        load(paste("d2r ",design[d]," fix=",fixed.cutoff[f]," p",1/p1[p]," n",n.T[j],".rData", sep=""))
        allResults[f,d,p,j,,,] <- results
	}

# turn array into data frame
allResults <- cast(melt(allResults), formula = ... ~ summary)

# fix factors
allResults$design <- ordered(allResults$design, levels=c("EG","Di"), 
                            labels=c("Extreme Group", "Dichotomization"))
allResults$stat <- ordered(allResults$stat, levels = stat.names())


# create smooth means and variances

var.sm <- melt(daply(allResults, .(fixed, design, p1, n, stat), 
                     function(x) predict(loess(var ~ rho, data = x))))
names(var.sm)[6:7] <- c("rho", "var.sm")
var.sm$rho <- (var.sm$rho - 1) / 100
allResults <- merge(var.sm, allResults)

mean.sm <- melt(daply(allResults, .(fixed, design, p1, n, stat), 
                      function(x) predict(loess(mean ~ rho, data = x))))
names(mean.sm)[6:7] <- c("rho", "mean.sm")
mean.sm$rho <- (mean.sm$rho - 1) / 100
allResults <- merge(mean.sm, allResults)


# save results to file

save(allResults, file="d2r results.rData")