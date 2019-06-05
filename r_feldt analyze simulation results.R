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

#############################################################
## compare different conversion constants on natural scale ##
#############################################################

nat <- droplevels(subset(allResults, stat %in% c("r.i","r.n")))
nat$bias <- nat$mean.sm - nat$rho
nat$rmse <- sqrt((nat$mean - nat$rho)^2 + nat$var)


# bias
ggplot(subset(nat, fixed==TRUE & design=="Extreme Group"), 
  aes(rho, bias, color = stat)) + 
	geom_line() + 
  facet_grid(n ~ p1, scale="free_y") + theme_bw()
ggplot(subset(nat, fixed==FALSE & design=="Extreme Group"), 
  aes(rho, bias, color = stat)) + 
  geom_line() + 
  facet_grid(n ~ p1, scale="free_y") + theme_bw()
ggplot(subset(nat, fixed==TRUE & design=="Dichotomization"), 
  aes(rho, bias, color = stat)) + 
  geom_line() + 
  facet_grid(n ~ p1, scale="free_y") + theme_bw()
ggplot(subset(nat, fixed==FALSE & design=="Dichotomization"), 
  aes(rho, bias, color = stat)) + 
  geom_line() + 
  facet_grid(n ~ p1, scale="free_y") + theme_bw()



# rmse
ggplot(subset(nat, fixed==TRUE & design=="Extreme Group"), 
  aes(rho, rmse, color = stat)) + 
  stat_smooth(se=FALSE) + 
  facet_grid(n ~ p1, scale="free_y") + theme_bw()
ggplot(subset(nat, fixed==FALSE & design=="Extreme Group"), 
  aes(rho, rmse, color = stat)) + 
  stat_smooth(se=FALSE) + 
  facet_grid(n ~ p1, scale="free_y") + theme_bw()
ggplot(subset(nat, fixed==TRUE & design=="Dichotomization"), 
  aes(rho, rmse, color = stat)) + 
  stat_smooth(se=FALSE) + 
  facet_grid(n ~ p1, scale="free_y") + theme_bw()
ggplot(subset(nat, fixed==FALSE & design=="Dichotomization"), 
  aes(rho, rmse, color = stat)) + 
  stat_smooth(se=FALSE) + 
  facet_grid(n ~ p1, scale="free_y") + theme_bw()

rmse <- cast(subset(nat,select=c(names(nat)[1:6],"rmse")), 
             fixed + design + p1 + n ~ stat, mean, value="rmse")
rmse$diff <- rmse$r.n - rmse$r.i
rmse

#######################################################
## compare different conversion constants on Z scale ##
#######################################################

Z <- droplevels(subset(allResults, substr(stat,1,1)=="z"))
Z$bias <- Z$mean.sm - Z.Fisher(Z$rho)
Z$rmse <- sqrt((Z$mean - Z.Fisher(Z$rho))^2 + Z$var)

# bias
ggplot(subset(Z, fixed==TRUE & design=="Extreme Group"),
	aes(rho, bias, color = stat)) + theme_bw() + 
	geom_line() + 
	facet_grid(n ~ p1, scale="free_y")

bias <- cast(subset(Z,select=c(names(nat)[1:6],"bias")), 
             fixed + design + p1 + n ~ stat, mean, value="bias")
bias$diff.F <- bias$z.nF - bias$z.iF
bias$diff.T <- bias$z.nT - bias$z.iT
bias

# rmse
rmse <- cast(subset(Z,select=c(names(nat)[1:6],"rmse")), 
             fixed + design + p1 + n ~ stat, mean, value="rmse")
rmse$diff.F <- rmse$z.nF - rmse$z.iF
rmse$diff.T <- rmse$z.nT - rmse$z.iT
rmse


#################
## bias of r_F ##
#################

r_F <- droplevels(subset(allResults, stat=="r.i"))
levels(r_F$fixed) <- c("Pop. cutoff","Sample cutoff")
r_F$bias <- r_F$mean - r_F$rho
r_F$bias.sm <- r_F$mean.sm - r_F$rho
r_F$rmse <- sqrt(r_F$bias^2 + r_F$var)

ggplot(r_F, aes(rho, bias, color = n)) +
  stat_smooth(se=FALSE) + 
  facet_grid(design + fixed ~ p1, scale="free_y") + theme_bw()

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

ggplot(InRange, aes(rho, mean.sm, color=n)) +
    geom_line() + coord_cartesian(ylim=c(0.6, 1)) + 
    facet_grid(design + fixed ~ p1) + theme_bw()

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
levels(Z$stat) <- c("Z.S","Z.T")
levels(Z$n) <- paste("n =", levels(Z$n))
Z$bias <- Z$mean - Z.Fisher(Z$rho)
Z$rmse <- sqrt(Z$bias^2 + Z$var)


# Dichotomization design
ggplot(droplevels(subset(Z, design=="Dichotomization")),
  aes(rho, bias, color = n)) +
  geom_line() + 
  facet_grid(fixed + stat ~ p1, scale="free_y") + theme_bw()


ggplot(droplevels(subset(Z, fixed=="Sample percentiles" & 
  design=="Dichotomization" & p.inv %in% c(2,3,5) & n != "n = 1000" &
  rho >= 0.5)),
  aes(rho, bias, linetype = stat)) +
  geom_smooth(method="loess", se=FALSE, color = "black") +
  facet_grid(n ~ p1, scale = "free_y") + theme_bw() +
  labs(linetype = "Estimator") +
  scale_y_continuous(name=expression(Bias(Z[eg])), limits=c(-0.1, 0.6)) + 
  scale_x_continuous(name=expression(rho))
ggsave(filename="Z_bias_Di.pdf", width=8, height=5)


# Extreme group design
ggplot(droplevels(subset(Z, design=="Extreme Group")),
  aes(rho, bias, color = n)) +
  geom_line() +
  facet_grid(fixed + stat ~ p1, scale="free_y") + theme_bw()


ggplot(droplevels(subset(Z, fixed=="Sample percentiles" & 
  design=="Extreme Group" & p.inv %in% c(3,5,8) & n != "n = 1000" &
  rho >= 0.5)),
  aes(rho, bias, linetype = stat)) +
  geom_smooth(method="loess", se=FALSE, color = "black") +
  facet_grid(n ~ p1, scale = "free_y") + theme_bw() +
  labs(linetype = "Estimator") +
  scale_y_continuous(name=expression(Bias(Z[eg])), limits = c(-0.1, 0.3)) + 
  scale_x_continuous(name=expression(rho))
ggsave(filename="Z_bias_EG.pdf", width=8, height=5)



# look at rmse for each design & percentile

for (d in levels(Z$design))
for (f in levels(Z$fixed))
print(ggplot(droplevels(subset(Z, fixed== f & design==d & n != "1000")),
  aes(rho, rmse, linetype = stat)) +
  geom_line() + stat_smooth(se=FALSE) +
  facet_grid(n ~ p1, scale = "free_y") + theme_bw() +
  labs(linetype = "Estimator") + opts(title = paste(f, d)) +
  scale_y_continuous(name=expression(rmse(Z))) + 
  scale_x_continuous(name=expression(rho)))


###########################
## Compare V_d to Var(d) ##
###########################

Var.d <- subset(allResults, stat=="d", select=c(1:4,6,10,11))
Var.d$p2 <- ifelse(Var.d$design=="Dichotomization", 1 - 1 / Var.d$p.inv, 1 / Var.d$p.inv)
Var.d$V.correct <- mapply(V.d.correct, rho = Var.d$rho, 
                          p1 = 1/Var.d$p.inv, p2 = Var.d$p2) / 
                          as.numeric(levels(Var.d$n))[Var.d$n]
Var.d$V.ratio <- with(Var.d, var / V.correct)

ggplot(Var.d, aes(rho, V.ratio, color = n)) +
  #geom_line() +
  stat_smooth(se=FALSE) + theme_bw() +
  facet_grid(design + fixed ~ p1, scale="free_y")

#############################
## Compare V_r to Var(r_F) ##
#############################

Var.r <- subset(allResults, stat=="r.i", select=c(1:4,6,10,11))
Var.r$p2 <- ifelse(Var.r$design=="Dichotomization", 1 - 1 / Var.r$p.inv, 1 / Var.r$p.inv)
Var.r$V.correct <- mapply(V.r.Feldt, rho = Var.r$rho, p1 = 1/Var.r$p.inv, 
                    p2 = Var.r$p2) / as.numeric(levels(Var.r$n))[Var.r$n]
Var.r$V.exp <- mapply(V.r.Feldt, rho = Var.r$rho, p1 = 1/Var.r$p.inv, 
                          p2 = Var.r$p2, MoreArgs=list(V=V.d.normal)) / 
                          as.numeric(levels(Var.r$n))[Var.r$n]
Var.r$V.ratio.cor <- with(Var.r, var / V.correct)
Var.r$V.ratio.exp <- with(Var.r, var / V.exp)
Var.r$V.ratio.BVN <- with(Var.r, var * as.numeric(levels(n))[n]/ (1 - rho^2)^2)

ggplot(Var.r, aes(rho, V.ratio.cor, color = n)) +
  #geom_line() + 
  stat_smooth(se=FALSE) + theme_bw() +
  facet_grid(design + fixed ~ p1, scale="free_y")

ggplot(Var.r, aes(rho, V.ratio.exp, color = n)) +
  #geom_line() + 
  stat_smooth(se=FALSE) + theme_bw() +
  facet_grid(design + fixed ~ p1, scale="free_y")

ggplot(Var.r, aes(rho, V.ratio.BVN, color = n)) +
  geom_line() + 
  stat_smooth(se=FALSE) + theme_bw() +
  coord_cartesian(ylim=c(0,4)) +
  facet_grid(design + fixed ~ p1, scale="free_y")


###########################
## Compare V_Z to Var(Z) ##
###########################

Var.Z <- subset(allResults, stat == "z.iT", select=c(1:4,6,10,11))
Var.Z$p2 <- ifelse(Var.Z$design=="Dichotomization", 1 - 1 / Var.Z$p.inv, 1 / Var.Z$p.inv)
Var.Z$V.correct <- mapply(V.Z.Feldt, rho = Var.Z$rho, p1 = 1/Var.Z$p.inv, 
                          p2 = Var.Z$p2) / as.numeric(levels(Var.Z$n))[Var.Z$n]
Var.Z$V.exp <- mapply(V.Z.Feldt, rho = Var.Z$rho, p1 = 1/Var.Z$p.inv, 
                      p2 = Var.Z$p2, MoreArgs=list(V=V.d.normal)) / 
                        as.numeric(levels(Var.Z$n))[Var.Z$n]
Var.Z$V.ratio.cor <- with(Var.Z, var / V.correct)
Var.Z$V.ratio.exp <- with(Var.Z, var / V.exp)

# asymptotic variance
ggplot(Var.Z, aes(rho, V.correct, color = n)) +
  geom_line() + 
  theme_bw() + 
  facet_grid(design + fixed ~ p1, scale="free_y")

# actual variance
ggplot(Var.Z, aes(rho, var, color = n)) +
  geom_line() + stat_smooth(se=FALSE) + 
  theme_bw() + 
  coord_cartesian(ylim=c(0,2)) + 
  facet_grid(design + fixed ~ p1, scale="free_y")
ggplot(droplevels(subset(Var.Z,n!="20")), aes(rho, var, color = n)) +
  geom_line() + stat_smooth(se=FALSE) + 
  theme_bw() + 
  facet_grid(design + fixed ~ p1, scale="free_y")


ggplot(Var.Z, aes(rho, V.ratio.cor, color = n)) +
  #geom_line() + 
  stat_smooth(se=FALSE) + theme_bw() + 
  coord_cartesian(ylim=c(0,10)) + 
  facet_grid(design + fixed ~ p1, scale="free_y")
ggplot(droplevels(subset(Var.Z, n!="20")), aes(rho, V.ratio.cor, color = n)) +
  geom_line() + 
  stat_smooth(se=FALSE) + theme_bw() + 
  facet_grid(design + fixed ~ p1, scale="free_y")


ggplot(Var.Z, aes(rho, V.ratio.exp, color = n)) +
  #geom_line() + 
  stat_smooth(se=FALSE) + theme_bw() + 
  coord_cartesian(ylim=c(0,10)) + 
  facet_grid(design + fixed ~ p1, scale="free_y")
ggplot(droplevels(subset(Var.Z, n!="20")), aes(rho, V.ratio.exp, color = n)) +
  geom_line() + 
  stat_smooth(se=FALSE) + theme_bw() + 
  facet_grid(design + fixed ~ p1, scale="free_y")


#########################
## Look at bias of V_r ##
#########################

Var.r <- subset(allResults, stat=="r.i", select=c(1:4,6,8,10))
ggplot(Var.r, aes(rho, var, color = n)) +
  geom_line() +
  stat_smooth(se=FALSE) + theme_bw() + 
  facet_grid(design + fixed ~ p1, scale="free_y")

V.r <- merge(Var.r, droplevels(subset(allResults, substr(stat,1,2)=="Vr", select=c(1:6,7,9))))
V.r$RB1 <- with(V.r, mean.sm / var.sm - 1)
V.r$RB2 <- with(V.r, mean / var.sm - 1)
V.r2 <- droplevels(subset(V.r, stat %in% c("Vr.asy","Vr.exp") & n != "1000" & p1 != "p1 = 1/4"))
levels(V.r2$n) <- paste("n =", levels(V.r2$n))
levels(V.r2$stat) <- c("V.eg","V.ce")

# compare smoothing methods
V.r.long <- melt.data.frame(V.r, id.vars=c(1:5,8), measure.vars=c(11,12))
for (d in levels(V.r$design))
for (f in levels(V.r$fixed))
print(ggplot(subset(V.r.long, design==d & fixed==f),
       aes(rho, value, linetype=variable, color = stat)) +
         geom_line() + theme_bw() + opts(title=paste(d,f)) +
         facet_grid(p1 ~ n, scale="free_y"))


# limit to Vr.asy and Vr.exp
for (d in levels(V.r$design))
  for (f in levels(V.r$fixed))
    print(ggplot(subset(V.r2, design==d & fixed==f),
                 aes(rho, RB2, color = stat)) +
                   geom_line() + theme_bw() + 
                   stat_smooth() + opts(title=paste(d,f)) +
                   facet_grid(p1 ~ n, scale="free_y"))


# limit to sample percentiles, p1 = 1/3

ggplot(droplevels(subset(V.r2, design=="Extreme Group" & fixed == "Sample percentiles" 
                         & p1 == "p1 = 1/3" & n != "1000")),
  aes(rho, RB2, linetype = stat)) +
  geom_line() +
  facet_wrap(~ n, ncol = 4) + theme_bw() +
  labs(linetype = "Estimator") +
  scale_y_continuous(name=expression(Relative_Bias(V^r)), limits = c(-0.2, 0.3)) + 
  scale_x_continuous(name=expression(rho))
ggsave(filename="Vr_RB_EG.pdf", width=8, height=3)


ggplot(droplevels(subset(V.r2, design=="Dichotomization" & fixed == "Sample percentiles" 
                         & p1 == "p1 = 1/2" & n != "1000")),
  aes(rho, RB2, linetype = stat)) +
  #geom_line() +
  geom_smooth(method="loess", se=FALSE, color = "black") +
  facet_wrap(~ n,ncol = 4) + theme_bw() +
  labs(linetype = "Estimator") +
  scale_y_continuous(name=expression(Relative_Bias(V^r)), limits = c(-0.1, 0.3)) + 
  scale_x_continuous(name=expression(rho))
ggsave(filename="Vr_RB_Di.pdf", width=8, height=3)




#########################
## Look at bias of V.Z ##
#########################

Var.Z <- droplevels(subset(allResults, stat=="z.iT" & n != "1000" & p1 != "p1 = 1/4",
                select=c(1:4,6,8,10)))
ggplot(Var.Z, aes(rho, var, color = n)) +
  geom_line() +
  stat_smooth(se=FALSE, family="symmetric") + theme_bw() + 
  coord_cartesian(ylim=c(0,100)) +
  facet_grid(design + fixed ~ p1, scale="free_y")
ggplot(droplevels(subset(Var.Z, n!="20")), aes(rho, var, color = n)) +
  geom_line() +
  stat_smooth(se=FALSE, family="symmetric") + theme_bw() + 
  facet_grid(design + fixed ~ p1, scale="free_y")


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

for (d in levels(V.r$design))
  for (f in levels(V.r$fixed))
    print(ggplot(subset(VZ.sm, design==d & fixed==f),
                 aes(rho, RB.sm, color = stat)) +
                   geom_line() + theme_bw() + 
                   coord_cartesian(ylim=c(-1,2)) +
                   opts(title=paste(d,f)) +
                   facet_grid(p1 ~ n, scale="free_y"))


# limit to sample percentiles, p1 = 1/3

ggplot(droplevels(subset(VZ.sm, design=="Extreme Group" & fixed == "Sample percentiles" &
        p1 == "p1 = 1/3" & n != "n = 20")),
  aes(rho, RB.sm, linetype = stat)) +
  geom_line() +
  facet_wrap(~ n, ncol = 3) + theme_bw() +
  labs(linetype = "Estimator") +
  scale_y_continuous(name=expression(Relative_Bias(V^Z)), limits = c(-0.2, 0.3)) + 
  scale_x_continuous(name=expression(rho))
ggsave(filename="VZ_RB_EG.pdf", width=6, height=2.5)


# limit to sample percentiles, p1 = 1/2

ggplot(droplevels(subset(VZ.sm, design=="Dichotomization" & fixed == "Sample percentiles" &
  p1 == "p1 = 1/2")),
  aes(rho, RB.sm, linetype = stat)) +
  geom_line() +
  facet_wrap( ~ n, ncol=3) + theme_bw() +
  labs(linetype = "Estimator") +
  scale_y_continuous(name=expression(Relative_Bias(V^Z)), limits = c(-0.1, 0.4)) + 
  scale_x_continuous(name=expression(rho))
ggsave(filename="VZ_RB_Di.pdf", width=6, height=2.5)

ggplot(droplevels(subset(VZ.sm, design=="Dichotomization" & 
  p1 != "p1 = 1/8" & fixed == "Sample percentiles")),
       aes(rho, RB.sm, linetype = stat)) +
         geom_line() +
         facet_grid(n ~ p1, scale = "free_y") + theme_bw() +
         labs(linetype = "Estimator") +
         scale_y_continuous(name=expression(Relative_Bias(V^Z))) + 
         scale_x_continuous(name=expression(rho))

ggplot(droplevels(subset(VZ.sm, design=="Dichotomization" & 
  p1 != "p1 = 1/8" & fixed == "Fixed percentiles" & n != "n = 20")),
       aes(rho, RB.sm, linetype = stat)) +
         geom_line() +
         facet_grid(n ~ p1, scale = "free_y") + theme_bw() +
         labs(linetype = "Estimator") +
         scale_y_continuous(name=expression(Relative_Bias(V^Z))) + 
         scale_x_continuous(name=expression(rho))
