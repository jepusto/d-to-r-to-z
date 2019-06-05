#setwd("C:\\Users\\James\\Dropbox\\Meta-analysis\\d-to-r conversion\\simulations and code")
#setwd("C:\\Users\\jep701\\Dropbox\\Meta-analysis\\d-to-r conversion\\simulations and code")


##########################
## simulation functions ##
##########################

Z.Fisher <- function(r) (log(1 + r) - log(1 - r)) / 2

simulation_run <- function(rho, n.T, iterations) {
  raw_data <- array(NA, dim=c(2,n.T, iterations / n.T))
  raw_data[1,,] <- rnorm(iterations)
  raw_data[2,,] <- raw_data[1,,] * rho + rnorm(iterations) * sqrt(1 - rho^2)

  r_p <- apply(raw_data, 3, function(x) cor(x[1,],x[2,]))
  Z_p <- Z.Fisher(r_p)
  V_r_p <- (1 - r_p^2)^2

  summary_stats <- matrix(NA, 2, 3)
  summary_stats[1,1] <- mean(r_p)
  summary_stats[2,1] <- var(r_p)
  summary_stats[1,2] <- mean(Z_p)
  summary_stats[2,2] <- var(Z_p)
  summary_stats[1,3] <- mean(V_r_p) / n.T
  summary_stats[2,3] <- var(V_r_p) / n.T^2

  return(summary_stats)
}


######################
## parameter values ##
######################

rho <- seq(0.0, 0.95, 0.01)
n.T <- c(20, 40, 80, 160, 1000)
iterations <- 10^7

#####################
## run simulations ##
#####################

sim_results <- array(NA, dim = c(length(rho), length(n.T), 2, 3), 
                            dimnames=list(rho=c(rho), n.T=c(n.T), summary=c("mean","var"), stat=c("r_p","Z_p","V_r_p")))

for (r in 1:length(rho)) 
   for (i in 1:length(n.T)) {
     print(c(rho[r], n.T[i]))
     sim_results[r,i,,] <- simulation_run(rho[r], n.T[i], iterations)
   }


library(reshape)
allResults <- cast(melt(sim_results), formula = ... ~ summary)

# create smooth means and variances

var.sm <- melt(daply(allResults, .(n.T, stat), function(x) predict(loess(var ~ rho, data = x))))
names(var.sm)[3:4] <- c("rho", "var.sm")
var.sm$rho <- (var.sm$rho - 1) / 100
allResults <- merge(var.sm, allResults)

mean.sm <- melt(daply(allResults, .(n.T, stat), function(x) predict(loess(mean ~ rho, data = x))))
names(mean.sm)[3:4] <- c("rho", "mean.sm")
mean.sm$rho <- (mean.sm$rho - 1) / 100
allResults <- merge(mean.sm, allResults)


# save results to file

save(allResults, file = "r_Pearson simulation results.rData")


#########################
## analyze simulations ##
#########################
load("r_Pearson simulation results.rData")
library(ggplot2)


# bias of r_p

r_p <- droplevels(subset(allResults, stat == "r_p" & n.T != 1000))
r_p$bias <- with(r_p, mean - rho)

ggplot(r_p, aes(rho, bias, linetype=factor(n.T))) +
         geom_smooth(method = "loess", se=FALSE, color = "black") + 
         labs(linetype = "n") + theme_bw() +
         scale_y_continuous(name=expression(Bias(r[p]))) + 
         scale_x_continuous(name=expression(rho))
ggsave(filename="r_p_bias.pdf", width=4, height=2.5)

summary(subset(r_p, n.T == 20)$bias)

# bias of Z_p

Z_p <- droplevels(subset(allResults, stat == "Z_p" & n.T != 1000))
Z_p$bias <- with(Z_p, mean - Z.Fisher(rho))

ggplot(Z_p, aes(rho, bias, linetype=factor(n.T))) +
  geom_smooth(method = "loess", se=FALSE, color = "black") + 
  labs(linetype = "n") + theme_bw() +
  scale_y_continuous(name=expression(Bias(z[p]))) + 
  scale_x_continuous(name=expression(rho))
ggsave(filename="Z_p_bias.pdf", width=4, height=2.5)

summary(subset(Z_p, n.T == 20)$bias)

# relative bias of V_r_p

Var_r_p <- subset(allResults, stat == "r_p", select = c("n.T","rho","var"))
V_r_p <- merge(Var_r_p, subset(allResults, stat == "V_r_p", select = c("n.T","rho","mean")))
V_r_p$RB <- with(V_r_p, mean / var - 1)

ggplot(subset(V_r_p, n.T != 1000),
  aes(rho, RB, linetype=factor(n.T))) +
  geom_smooth(method = "loess", se=FALSE, color = "black") + 
  labs(linetype = "n") + theme_bw() +
  scale_y_continuous(name=expression(Relative_Bias(V[p]^r)), limits = c(-0.15, 0)) + 
  scale_x_continuous(name=expression(rho))
ggsave(filename="V_r_p_relbias.pdf", width=4, height=2.5)


# relative bias of V_Z_p

V_Z_p <- subset(allResults, stat == "Z_p", select = c("n.T","rho","var"))
V_Z_p$RB <- with(V_Z_p, 1 / (var * (n.T - 3)) - 1)

ggplot(subset(V_Z_p, n.T != 1000),
       aes(rho, RB, linetype=factor(n.T))) +
         geom_smooth(method = "loess", se=FALSE, color = "black") + 
         labs(linetype = "n") + theme_bw() +
         scale_y_continuous(name=expression(Relative_Bias(V[p]^z))) + 
         scale_x_continuous(name=expression(rho))
ggsave(filename="V_Z_p_relbias.pdf", width=4, height=2.5)
