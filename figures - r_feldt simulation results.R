#setwd("C:\\Users\\James\\Dropbox\\Meta-analysis\\d-to-r conversion\\simulations and code")
#setwd("C:\\Documents and Settings\\James Pustejovsky\\My Documents\\My Dropbox\\Meta-analysis\\d-to-r conversion\\simulations and code")
#setwd("C:\\Users\\jep701\\Dropbox\\Meta-analysis\\d-to-r conversion\\simulations and code")

source("r_feldt simulation functions.R")
library(ggplot2)
library(reshape)

load("d2r results.rData")
allResults$n <- ordered(allResults$n)
allResults$p.inv <- allResults$p1
allResults$p1 <- ordered(allResults$p1, 
                labels = paste("p1 = 1/",unique(allResults$p1), sep=""))
allResults$fixed <- ordered(allResults$fixed, levels=c("TRUE","FALSE"), 
                labels = c("Fixed percentiles","Sample percentiles"))

#################
## bias of r_F ##
#################

r_F <- droplevels(subset(allResults, stat=="r.i"))
levels(r_F$fixed) <- c("Pop. cutoff","Sample cutoff")
r_F$bias <- r_F$mean - r_F$rho
r_F$bias.sm <- r_F$mean.sm - r_F$rho
r_F$rmse <- sqrt(r_F$bias^2 + r_F$var)

# Dichotomization design
ggplot(droplevels(subset(r_F, design=="Dichotomization" & p.inv %in% c(2,3,5) & n != "1000")), 
  aes(rho, bias, linetype=n)) +
  geom_smooth(method="loess", se=FALSE, color = "black") + 
  facet_grid(fixed ~ p1) + theme_bw() +
  labs(linetype = "n") +
  scale_y_continuous(name=expression(Bias(r[eg]))) + 
  scale_x_continuous(name=expression(rho))
ggsave(filename="r_F_bias_Di.pdf", width=8, height=4)


# Extreme group design
ggplot(droplevels(subset(r_F, design=="Extreme Group" & p.inv %in% c(3,5,8) & n != "1000")), 
  aes(rho, bias, linetype = n)) +
  geom_smooth(method="loess", se=FALSE, color = "black") + 
  facet_grid(fixed ~ p1) + theme_bw() +
  labs(linetype = "n") +
  scale_y_continuous(name=expression(Bias(r[eg]))) + 
  scale_x_continuous(name=expression(rho))
ggsave(filename="r_F_bias_EG.pdf", width=8, height=4)


max.bias <- aggregate(abs(bias.sm) ~ n, 
                      data = subset(r_F, design=="Extreme Group" & p1!="p1 = 1/2"),
                      FUN = max)
max.bias[,2] * as.numeric(levels(max.bias$n))[max.bias$n]

# example: bias of r0 for 
# EG, sample percentiles, p1 = 1/5, n.T=160, rho=0.40
subset(r_F, design=="Extreme Group" & fixed=="Sample percentiles" & 
      p1=="p1 = 1/5" & n=="160" & rho==0.4)$mean.sm / ab(1/5,1/5)[2] - 0.4


#################################################
## Percentage of points falling outside (-1,1) ##
#################################################

InRange <- subset(allResults, stat=="InRange")

ggplot(droplevels(subset(InRange, fixed=="Sample percentiles" & 
  n!="1000" & p.inv %in% c(2,3,5) & rho >= 0.4 & design == "Dichotomization")),
       aes(rho, mean, linetype=n)) + 
         geom_line() + 
         facet_wrap( ~ p1) + theme_bw() +
         labs(linetype = "n") +
         scale_y_continuous(name=expression(Pr( group("|",r[eg],"|") < 1))) + 
         scale_x_continuous(name=expression(rho))
ggsave(filename="r_F_in_range_DI.pdf", width=8, height=3)

ggplot(droplevels(subset(InRange, fixed=="Sample percentiles" & 
  n!="1000" & p.inv %in% c(3,5,8) & rho >= 0.4 & design == "Extreme Group")),
       aes(rho, mean, linetype=n)) + 
         geom_line() + 
         facet_wrap( ~ p1) + theme_bw() +
         labs(linetype = "n") +
         scale_y_continuous(name=expression(Pr( group("|",r[eg],"|") < 1))) + 
         scale_x_continuous(name=expression(rho))
ggsave(filename="r_F_in_range_EG.pdf", width=8, height=3)

subset(InRange, fixed=="Sample percentiles" & design=="Dichotomization" & 
  n=="20" & p1=="p1 = 1/2" & rho >= 0.7, select=c("rho","mean.sm","mean"))

###############
## bias of Z ##
###############

Z <- droplevels(subset(allResults, substr(stat,1,3)=="z.i"))
levels(Z$stat) <- c("z.S","z.T")
levels(Z$n) <- paste("n =", levels(Z$n))
Z$bias <- Z$mean - Z.Fisher(Z$rho)
Z$rmse <- sqrt(Z$bias^2 + Z$var)

ggplot(droplevels(subset(Z, fixed=="Sample percentiles" & 
  design=="Dichotomization" & p.inv %in% c(2,3,5) & n != "n = 1000" &
  rho >= 0.5)),
  aes(rho, bias, linetype = stat)) +
  geom_smooth(method="loess", se=FALSE, color = "black") +
  facet_grid(n ~ p1, scale = "free_y") + theme_bw() +
  labs(linetype = "Estimator") +
  scale_y_continuous(name=expression(Bias(z[eg])), limits=c(-0.1, 0.6)) + 
  scale_x_continuous(name=expression(rho))
ggsave(filename="Z_bias_Di.pdf", width=8, height=5)


# Extreme group design

ggplot(droplevels(subset(Z, fixed=="Sample percentiles" & 
  design=="Extreme Group" & p.inv %in% c(3,5,8) & n != "n = 1000" &
  rho >= 0.5)),
  aes(rho, bias, linetype = stat)) +
  geom_smooth(method="loess", se=FALSE, color = "black") +
  facet_grid(n ~ p1, scale = "free_y") + theme_bw() +
  labs(linetype = "Estimator") +
  scale_y_continuous(name=expression(Bias(z[eg])), limits = c(-0.1, 0.3)) + 
  scale_x_continuous(name=expression(rho))
ggsave(filename="Z_bias_EG.pdf", width=8, height=5)


##########################
## Relative bias of V_r ##
##########################

Var.r <- subset(allResults, stat=="r.i", select=c(1:4,6,8,10))
V.r <- merge(Var.r, droplevels(subset(allResults, substr(stat,1,2)=="Vr", select=c(1:6,7,9))))
V.r$RB1 <- with(V.r, mean.sm / var.sm - 1)
V.r$RB2 <- with(V.r, mean / var.sm - 1)
V.r2 <- droplevels(subset(V.r, stat %in% c("Vr.asy","Vr.exp") & n != "1000" & p1 != "p1 = 1/4"))
levels(V.r2$n) <- paste("n =", levels(V.r2$n))
levels(V.r2$stat) <- c("V.eg","V.ce")
levels(V.r2$fixed) <- c("Pop. cutoff","Sample cutoff")

# Dichotomization, limit to sample percentiles, p1 = 1/2

ggplot(droplevels(subset(V.r2, design=="Dichotomization" & fixed == "Sample cutoff" 
                         & p1 == "p1 = 1/2" & n != "1000")),
  aes(rho, RB2, linetype = stat)) +
  #geom_line() +
  geom_smooth(method="loess", se=FALSE, color = "black") +
  facet_wrap(~ n,ncol = 4) + theme_bw() +
  labs(linetype = "Estimator") +
  scale_y_continuous(name=expression(Relative_Bias(V^r)), limits = c(-0.1, 0.3)) + 
  scale_x_continuous(name=expression(rho))
ggsave(filename="Vr_RB_Di.pdf", width=8, height=3)

# Extreme Groups, limit to p1 = 1/3

ggplot(droplevels(subset(V.r2, design=="Extreme Group"
                         & p1 == "p1 = 1/3" & n != "1000")),
       aes(rho, RB2, linetype = stat)) +
         geom_line() +
         facet_grid(fixed ~ n, scale = "free_y") + theme_bw() +
         labs(linetype = "Estimator") +
         scale_y_continuous(name=expression(Relative_Bias(V^r)), limits = c(-0.2, 0.3)) + 
         scale_x_continuous(name=expression(rho))
ggsave(filename="Vr_RB_EG.pdf", width=8, height=4)

##########################
## Relative bias of V.Z ##
##########################

Var.Z <- droplevels(subset(allResults, stat=="z.iT" & n != "1000" & p1 != "p1 = 1/4",
                select=c(1:4,6,8,10)))
V.Z <- merge(Var.Z, droplevels(subset(allResults, 
        substr(stat,1,2)=="VZ" & n != "1000" & p1 != "p1 = 1/4", 
        select=c(1:6,7,9))))
V.Z$RB <- with(V.Z, mean / var - 1)
VZ.sm <- melt(daply(V.Z, .(fixed, design, p1, n, stat), 
                     function(x) predict(loess(RB ~ rho, data = x,family="symmetric"))))
names(VZ.sm)[6:7] <- c("rho", "RB.sm")
VZ.sm$rho <- (VZ.sm$rho - 1) / 100
VZ.sm$n <- as.factor(VZ.sm$n)
levels(VZ.sm$n) <- paste("n =", levels(VZ.sm$n))
levels(VZ.sm$stat) <- c("V.eg","V.ce")
levels(VZ.sm$fixed) <- c("Pop. cutoff","Sample cutoff")

# Dichotomization, limit to sample percentiles, p1 = 1/2

ggplot(droplevels(subset(VZ.sm, design=="Dichotomization" & fixed == "Sample cutoff" &
  n != "n = 20" & p1 == "p1 = 1/2")),
  aes(rho, RB.sm, linetype = stat)) +
  geom_line() +
  facet_wrap( ~ n, ncol=3) + theme_bw() +
  labs(linetype = "Estimator") +
  scale_y_continuous(name=expression(Relative_Bias(V^z)), limits = c(-0.1, 0.4)) + 
  scale_x_continuous(name=expression(rho))
ggsave(filename="VZ_RB_Di.pdf", width=6, height=2.5)

ggplot(droplevels(subset(VZ.sm, design=="Dichotomization" & 
  p1 != "p1 = 1/8" & fixed == "Sample percentiles")),
       aes(rho, RB.sm, linetype = stat)) +
         geom_line() +
         facet_grid(n ~ p1, scale = "free_y") + theme_bw() +
         labs(linetype = "Estimator") +
         scale_y_continuous(name=expression(Relative_Bias(V^z))) + 
         scale_x_continuous(name=expression(rho))

ggplot(droplevels(subset(VZ.sm, design=="Dichotomization" & 
  p1 != "p1 = 1/8" & fixed == "Fixed percentiles" & n != "n = 20")),
       aes(rho, RB.sm, linetype = stat)) +
         geom_line() +
         facet_grid(n ~ p1, scale = "free_y") + theme_bw() +
         labs(linetype = "Estimator") +
         scale_y_continuous(name=expression(Relative_Bias(V^z))) + 
         scale_x_continuous(name=expression(rho))

# Extreme groups, limit to p1 = 1/3

ggplot(droplevels(subset(VZ.sm, design=="Extreme Group" & n != "n = 20" &
  p1 == "p1 = 1/3")),
       aes(rho, RB.sm, linetype = stat)) +
         geom_line() +
         facet_grid(fixed ~ n) + theme_bw() +
         labs(linetype = "Estimator") +
         scale_y_continuous(name=expression(Relative_Bias(V^z)), limits = c(-0.2, 0.3)) + 
         scale_x_continuous(name=expression(rho))
ggsave(filename="VZ_RB_EG.pdf", width=7, height=4)

