#setwd("C:\\Users\\James\\Dropbox\\Meta-analysis\\d-to-r conversion\\simulations and code")
#setwd("C:\\Users\\jep701\\Dropbox\\Meta-analysis\\d-to-r conversion\\simulations and code")
source("r_feldt simulation functions.R")
source("parameter values.R")

f <- 2
d <- 1
p <- 5
j <- 4

rho.rand <- sample(rho)
system.time({
  
for (f in 1:2)
  for (d in 1:2)
    for (p in 1:5)
      for (j in 1:5) {

  results <- array(NA, dim=c(length(rho), length(stat.names()), 2), 
    dimnames = c(rho=list(rho.rand), stat=list(stat.names()), summary=list(c("mean","var"))))

	for (r in 1:length(rho)) {
    if (design[d] == "EG") {
      SMD <- simulate_SMD(rho.rand[r], p1[p], p1[p], n.T[j]/2, n.T[j]/2, iterations[j], fixed.cutoff[f])
      results[r,,] <- analyze(SMD, p1[p], p1[p], n.T[j]/2, n.T[j]/2, truncation=truncation, k=taylor.degree)
    } else {
      n1 <- round(p1[p] * n.T[j], 0)
      n2 <- n.T[j] - n1
      SMD <- simulate_SMD(rho.rand[r], p1[p], 1 - p1[p], n1, n2, iterations[j], fixed.cutoff[f])
      results[r,,] <- analyze(SMD, p1[p], 1 - p1[p], n1, n2, truncation = truncation, k=taylor.degree)
    }
		print(rho.rand[r])
		}
  results <- results[order(rho.rand),,]
  print(paste("d2r ",design[d]," fix=",fixed.cutoff[f]," p",1/p1[p]," n",n.T[j],".rData", sep=""))
	#save(results, file = paste("d2r ",design[d]," fix=",fixed.cutoff[f]," p",1/p1[p]," n",n.T[j],".rData", sep=""))

}
})