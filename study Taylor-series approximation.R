setwd("C:\\Documents and Settings\\James Pustejovsky\\My Documents\\My Dropbox\\Meta-analysis\\d-to-r conversion\\simulations and code")
source("r_feldt simulation functions.R")

##########################################
## look at a, b constants for varying p ##
##########################################

p <- seq(0, 0.5, 0.001)
c.EG <- sapply(p, function(x) ab.n(x, x, Inf))
c.Di <- sapply(p, function(x) ab.n(x, 1 - x, Inf))

plot(p, c.EG[1,], type = "l", ylab="a")
lines(p, c.Di[1,], lty="dashed")

plot(p, c.EG[2,], type = "l", ylim=c(1,2), ylab="b")
lines(p, c.Di[2,], lty="dashed")


########################################################
## compare truncation to Taylor series approximation  ##
## for Z transformation						##
########################################################

p1 <- c(1/2, 1/3, 1/4, 1/5, 1/8)
truncation <- 1 - 10^-4

b <- sapply(p1, function(x) ab.n(x, x, Inf)[2])
Z.trunc <- Vectorize(function(x, b) Z.Fisher(min(x*b, truncation)), "x")

curve(Z.trunc(x, b[1]), 0, 1)
curve(Z.Taylor(x, b[1]), add = TRUE, lty = "dashed")


