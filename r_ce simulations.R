#setwd("C:\\Users\\James\\Dropbox\\Meta-analysis\\d-to-r conversion\\simulations and code")
#setwd("C:\\Users\\jep701\\Dropbox\\Meta-analysis\\d-to-r conversion\\simulations and code")


##########################
## simulation functions ##
##########################

Z.Fisher <- function(r) (log(1 + r) - log(1 - r)) / 2

J <- function(x) 1 - 3 / (4 * x - 1)

simulation_run <- function(rho, lambda.diff, n.T, iterations) {
  d_ce <- rt(n = iterations / n.T, 
             df = n.T - 2, 
             ncp = rho * lambda.diff * sqrt(n.T) / (2 * sqrt(1 - rho^2))) * J(n.T - 2) * 2 / sqrt(n.T)
  r_ce <- d_ce / sqrt(d_ce^2 + lambda.diff^2)
  Z_ce <- Z.Fisher(r_ce)
  V_d_ce <- 4 / n.T + d_ce^2 / (2 * n.T)
  V_r_ce <- V_d_ce * lambda.diff^4 / (d_ce^2 + lambda.diff^2)^3
  V_Z_ce <- V_d_ce / (d_ce^2 + lambda.diff^2)
  
  summary_stats <- matrix(NA, 2, 5)
  summary_stats[1,1] <- mean(d_ce)
  summary_stats[2,1] <- var(d_ce)
  summary_stats[1,2] <- mean(r_ce)
  summary_stats[2,2] <- var(r_ce)
  summary_stats[1,3] <- mean(Z_ce)
  summary_stats[2,3] <- var(Z_ce)
  summary_stats[1,4] <- mean(V_r_ce)
  summary_stats[2,4] <- var(V_r_ce)
  summary_stats[1,5] <- mean(V_Z_ce)
  summary_stats[2,5] <- var(V_Z_ce)
  
  return(summary_stats)
}


######################
## parameter values ##
######################

rho <- seq(0.0, 0.95, 0.01)
n.T <- c(20, 40, 80, 160, 1000, 10000)
lambda.diff <- c(0.5, 1, 2)
iterations <- 10^7

#####################
## run simulations ##
#####################

sim_results <- array(NA, dim = c(length(rho), length(n.T), length(lambda.diff), 2, 5), 
                     dimnames=list(rho=c(rho), n.T=c(n.T), diff=c(lambda.diff), 
                                   summary=c("mean","var"), stat=c("d_ce","r_ce","Z_ce","V_r_ce","V_Z_ce")))

for (r in 1:length(rho)) 
  for (i in 1:length(n.T)) 
    for (l in 1:length(lambda.diff)) {
    print(c(rho[r], n.T[i], lambda.diff[l]))
    sim_results[r,i,l,,] <- simulation_run(rho[r], lambda.diff[l], n.T[i], iterations)
  }


library(reshape)
allResults <- cast(melt(sim_results), formula = ... ~ summary)

# create smooth means and variances

var.sm <- melt(daply(allResults, .(n.T, diff, stat), function(x) predict(loess(var ~ rho, data = x))))
names(var.sm)[4:5] <- c("rho", "var.sm")
var.sm$rho <- (var.sm$rho - 1) / 100
allResults <- merge(var.sm, allResults)

mean.sm <- melt(daply(allResults, .(n.T, diff, stat), function(x) predict(loess(mean ~ rho, data = x))))
names(mean.sm)[4:5] <- c("rho", "mean.sm")
mean.sm$rho <- (mean.sm$rho - 1) / 100
allResults <- merge(mean.sm, allResults)


# save results to file

save(allResults, file = "r_ce simulation results.rData")


#########################
## analyze simulations ##
#########################
load("r_ce simulation results.rData")
library(ggplot2)


# check distribution of d

d_ce <- subset(allResults, stat == "d_ce")
d_ce$g <- with(d_ce, mean * J(n.T - 2))
d_ce$delta <- with(d_ce, rho * diff / sqrt(1 - rho^2))
plot(mean ~ delta, data = d_ce)

d_ce$var_exact <- with(d_ce, (n.T - 2) * 4 / (n.T * (n.T - 4)) + delta^2 * ((n.T - 2) / (n.T - 4) - 1 / J(n.T - 2)^2))
plot(var ~ var_exact, data = d_ce)


# bias of r

r_ce <- droplevels(subset(allResults, stat == "r_ce" & n.T < 1000))
r_ce$bias <- with(r_ce, mean - rho)
r_ce$diff <- paste("w =", r_ce$diff)

ggplot(r_ce, aes(rho, bias, linetype=factor(n.T))) +
  geom_smooth(method = "loess", se=FALSE, color = "black") + 
  labs(linetype = "n") + theme_bw() +
  facet_wrap(~diff) +
  scale_y_continuous(name=expression(Bias(r[ce]))) + 
  scale_x_continuous(name=expression(rho))
ggsave(filename="r_ce_bias.pdf", width=7, height=2.25)


# bias of Z_ce

Z_ce <- droplevels(subset(allResults, stat == "Z_ce" & n.T < 1000))
Z_ce$bias <- with(Z_ce, mean - Z.Fisher(rho))
Z_ce$diff <- paste("w =", Z_ce$diff)

ggplot(Z_ce, aes(rho, bias, linetype=factor(n.T))) +
  geom_smooth(method = "loess", se=FALSE, color = "black") + 
  labs(linetype = "n") + theme_bw() +
  facet_wrap(~diff) +
  scale_y_continuous(name=expression(Bias(z[ce]))) + 
  scale_x_continuous(name=expression(rho))
ggsave(filename="Z_ce_bias.pdf", width=7, height=2.25)

# relative bias of V_r_ce

qplot(rho, var.sm, linetype = factor(n.T), data = subset(allResults, stat=="r_ce"),
      geom = "line",facets = ~diff)

Var_r_ce <- subset(allResults, stat == "r_ce", select = c("n.T","diff","rho","var"))
V_r_ce <- merge(Var_r_ce, subset(allResults, stat == "V_r_ce", select = c("n.T","diff","rho","mean")))
V_r_ce$RB <- with(V_r_ce, mean / var - 1)
V_r_ce$diff <- paste("w =", V_r_ce$diff)

ggplot(subset(V_r_ce, n.T < 1000), aes(rho, RB, linetype=factor(n.T))) +
         geom_smooth(method = "loess", se=FALSE, color = "black") + 
         labs(linetype = "n") + theme_bw() +
         facet_wrap(~diff) +
         scale_y_continuous(name=expression(Relative_Bias(V[ce]^r))) + 
         scale_x_continuous(name=expression(rho))
ggsave(filename="V_r_ce_relbias.pdf", width=7, height=2.25)


# relative bias of V_Z_p

Var_Z_ce <- subset(allResults, stat == "Z_ce", select = c("n.T","diff","rho","var"))
V_Z_ce <- merge(Var_Z_ce, subset(allResults, stat == "V_Z_ce", select = c("n.T","diff","rho","mean")))
V_Z_ce$RB <- with(V_Z_ce, mean / var - 1)
V_Z_ce$diff <- paste("w =", V_Z_ce$diff)


ggplot(subset(V_Z_ce, n.T < 1000), aes(rho, RB, linetype=factor(n.T))) +
         geom_smooth(method = "loess", se=FALSE, color = "black") + 
         labs(linetype = "n") + theme_bw() +
         facet_wrap(~diff) + 
         scale_y_continuous(name=expression(Relative_Bias(V[ce]^z))) + 
         scale_x_continuous(name=expression(rho))
ggsave(filename="V_Z_ce_relbias.pdf", width=7, height=2.25)
