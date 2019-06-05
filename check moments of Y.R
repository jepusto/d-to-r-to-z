#setwd("C:\\Users\\James\\Dropbox\\Meta-analysis\\d-to-r conversion\\simulations and code")
#setwd("C:\\Users\\jep701\\Dropbox\\Meta-analysis\\d-to-r conversion\\simulations and code")
#setwd("C:\\Documents and Settings\\James Pustejovsky\\My Documents\\My Dropbox\\Meta-analysis\\d-to-r conversion\\simulations and code")
source("r_feldt simulation functions.R")
library(truncnorm)
library(ggplot2)

##########################################
## check cumulants from expression (11) ##
##########################################

n <- 10^7
rho <- 0.7
c <- 1 / 4
u <- qnorm(c)
v <- dnorm(u) / c

x <- rtruncnorm(n, b = beta)
e <- rnorm(n)
y <- rho * x + sqrt(1 - rho^2) * e

plot(density(y))

- rho * v
mean(y)

1 - rho^2 * v * (v + u)
mean((y + rho * v)^2)

- rho^3 * v * (2 * v^2 + 3 * u * v + u^2 - 1)
mean((y + rho * v)^3)

- rho^4 * v * (6 * v^3 + 12 * u * v^2 + 7 * u^2 * v - 4 * v - 3 * u + u^3)
mean((y + rho * v)^4) - 3 * mean((y + rho * v)^2)^2




#############################################
## Check asymptotic variance of y_i, s_i^2 ##
#############################################

rho <- 0.7
p1 <- 1/4
n1 <- 160
iterations <- 10^5


x <- matrix(rtruncnorm(n1 * iterations, b = qnorm(p1)), n1, iterations)

e <- matrix(rnorm(n1 * iterations), n1, iterations)
y <- rho * x + sqrt(1 - rho^2) * e

s.1.sq <- apply(y, 2, var)
y.bar.1 <- colMeans(y)

k1 <- cumulants(rho, p1)

mean(y.bar.1)
k1[1]
mean(s.1.sq)
k1[2]

n1 * var(cbind(y.bar.1, s.1.sq))
k1[2]
k1[3]
k1[4] + 2 * k1[2]^2


#####################################
## look at distribution of d, r, Z ##
#####################################


rho <- 0.8
p1 <- 1/4
iterations <- 10^4
n1 <- c(20,40,160,640)
colors <- c("orange","red","green","blue","brown","black")
delta <- r_to_d(rho, ab(p1,p1))
abC <- ab(p1, p1)

f <- function(n1) {
  d <- simulate_SMD(rho,p1,p1,n1,n1,iterations,TRUE)
  r <- d_to_r(d, abC)
  Z <- Z.Taylor(d_to_r(d, c(abC[1], 1)), abC[2])
  d.dens <- density(sqrt(2 * n1) * (d - delta))
  r.dens <- density(sqrt(2 * n1) * (r - rho))
  Z.dens <- density(sqrt(2 * n1) * (Z - Z.Fisher(rho)))
  return(list(d.dens,r.dens,Z.dens))
}
densities <- sapply(n1, f)


par(mfrow=c(1,3))

# d metric
curve(dnorm(x,mean=0,sd=sqrt(V.d.correct(rho,p1,p1))), 
      -8,8,ylim=c(0,.3),lty="dashed", ylab="", main = "d")
abline(v=0)
for (i in 1:length(n1)) lines(densities[1,i][[1]], col=colors[i])

# r metric
curve(dnorm(x,mean=0,sd=sqrt(V.r.Feldt(rho,p1,p1))), 
      -2*abC[2] - rho, 2*abC[2] - rho,ylim=c(0,.8),lty="dashed", ylab="", main = "r")
abline(v=0)
for (i in 1:length(n1)) lines(densities[2,i][[1]], col=colors[i])


# Z metric
curve(dnorm(x,mean=0,sd=sqrt(V.Z.Feldt(rho,p1,p1))), 
      -4, 4,ylim=c(0,.6),lty="dashed", ylab="", main="Z")
abline(v=0)
for (i in 1:length(n1)) lines(densities[3,i][[1]], col=colors[i])