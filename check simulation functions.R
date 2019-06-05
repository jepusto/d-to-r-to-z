
setwd("C:\\Users\\James\\Dropbox\\Meta-analysis\\d-to-r conversion\\simulations and code")
#setwd("C:\\Users\\jep701\\Dropbox\\Meta-analysis\\d-to-r conversion\\simulations and code")
source("r_feldt simulation functions.R")

rho <- 0.9
p1 <- 1/2
p2 <- 1 - p1
n1 <- 5
print(n2 <- p2 * n1 / p1)
iterations <- 10000
fixed <- TRUE

d <- simulate_SMD(rho, p1, p2, n1, n2, iterations, fixed)

results <- analyze(d, p1, p2, n1, n2)

results

results["d","var"]
V.d.correct(rho, p1, p2) / (n1 + n2)
V.d.correct(rho, p1, p2) / (n1 + n2) / results["d","var"]
V.d.normal(rho, p1, p2) / (n1 + n2)
V.d.normal(rho, p1, p2) / (n1 + n2) / results["d","var"]
results[c("Vd.asy","Vd.exp"),"mean"] / results["d","var"]

results["r.i","var"]
V.r.Feldt(rho, p1, p2) / (n1 + n2)
V.r.Feldt(rho, p1, p2) / (n1 + n2) / results["r.i","var"]
V.r.Feldt(rho, p1, p2, V=V.d.normal) / (n1 + n2)
V.r.Feldt(rho, p1, p2, V=V.d.normal) / (n1 + n2) / results["r.i","var"]
V.r.MVnormal(rho) / (n1 + n2)
V.r.MVnormal(rho) / (n1 + n2) / results["r.i","var"]
results[c("Vr.asy","Vr.exp","Vr.BVN"),"mean"] / results["r.i","var"]

results["z.iF","var"]
V.Z.Feldt(rho, p1, p2) / (n1 + n2)
V.Z.Feldt(rho, p1, p2) / (n1 + n2) / results["z.iF","var"]
V.Z.Feldt(rho, p1, p2, V=V.d.normal) / (n1 + n2)
V.Z.Feldt(rho, p1, p2, V=V.d.normal) / (n1 + n2) / results["z.iF","var"]
results[c("VZ.asy","VZ.exp"),"mean"] / results["z.iF","var"]