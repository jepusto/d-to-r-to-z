fixed.cutoff <- c(TRUE, FALSE)
design <- c("EG", "Di")
rho <- seq(0.0, 0.95, 0.01)
p1 <- c(1/2, 1/3, 1/4, 1/5, 1/8)
n.T <- c(20, 40, 80, 160, 1000)


iterations <- 10^7 / n.T

truncation <- 1 - 10^-12
taylor.degree <- 5


#for (p in 1:5)
#for (j in 1:5) {
#    n1 <- round(p1[p] * n.T[j],0)
#    print(paste("Di design: p1=",p1[p]," n1=",n1," n2=",n.T[j]-n1, " f.sample=",n1 / n.T[j]))
#}