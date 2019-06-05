#setwd("C:\\Documents and Settings\\James Pustejovsky\\My Documents\\My Dropbox\\Meta-analysis\\d-to-r conversion\\simulations and code")
#setwd("C:\\Users\\James\\Dropbox\\Meta-analysis\\d-to-r conversion\\simulations and code")
#setwd("C:/Users/jep701/Dropbox/Meta-analysis/d-to-r conversion/simulations and code")
source("r_feldt simulation functions.R")

library(TeachingDemos)


########################################
## plot a, b constants for varying p1 ##
########################################
parfix <- par()
parfix$mar <- par()$mar + c(0,0,-3,0)


# a constant
pdf("aConst.pdf", width = 5, height = 5)
par(parfix)
curve( ab(x, x)[,1], 0, 0.5, 
       xlab=expression(p[1]), ylab = "a", ylim=c(0,35), bty="n")
curve( ab(x, 1 - x)[,1], add=TRUE, lty="dashed")
text(0.27, 10, expression(p[2]==p[1]))
text(0.1, 7, expression(p[2]==1 - p[1]))
dev.off()

# b constant
pdf("bConst.pdf", width = 5, height = 5)
par(parfix)
curve( ab(x, x)[,2], 0, 0.5, 
       xlab=expression(p[1]), ylab = "b",ylim=c(1, 2), bty="n")
curve( ab(x, 1 - x)[,2], add=TRUE, lty="dashed")
text(0.25, 1.10, expression(p[2]==p[1]))
text(0.27, 1.45, expression(p[2]==1 - p[1]))
dev.off()



##################################
## compare V.d.norm to true V.d ##
##################################

p <- c(2,3,5,8)
linetypes = c("solid","dashed", "dotted","dotdash")
parfix <- par()
parfix$mar <- par()$mar + c(0,1,-3,0)

D.ratio <- Vectorize(function(rho, p1, p2) V.d.normal(rho, p1, p2) / V.d.correct(rho, p1, p2), "rho")


# Extreme groups design
pdf("compareVd_EG.pdf", width = 5, height = 5)
par(parfix)
curve( D.ratio(x, 1 / p[1], 1 / p[1]), 0, 1,
	xlab = expression(rho), ylab = expression(V[ce]^" d" / V[eg]^" d"),
	ylim = c(0.5,2), bty = "n")
for (i in 2:4) curve( D.ratio(x, 1 / p[i], 1 / p[i]), add = TRUE, lty = linetypes[i])
legend("topleft", legend = parse(text=paste("p[1] == 1/",p, sep="")),
       lty=linetypes, bty="n", 
       horiz=FALSE, xjust = .5, yjust = .5, cex=1)
dev.off()

# Dichotomization design
pdf("compareVd_DI.pdf", width = 5, height = 5)
par(parfix)
curve( D.ratio(x, 1 / p[1], 1 - 1 / p[1]), 0, 1,
	xlab = expression(rho), ylab = expression(V[ce]^d / V[eg]^d),
	ylim = c(0.5,2), bty = "n")
for (i in 2:4) curve( D.ratio(x, 1 / p[i], 1 - 1 / p[i]), add = TRUE, lty = linetypes[i])
legend("topleft", legend = parse(text=paste("p[1] == 1/",p, sep="")),
       lty=linetypes, bty="n", 
       horiz=FALSE, xjust = .5, yjust = .5, cex=1)
dev.off()



#################################
## compare V.r.MVN to true V.r ##
#################################

R.ratio <- Vectorize(function(rho, p1, p2) V.r.MVnormal(rho) / V.r.Feldt(rho, p1, p2), "rho")



# Extreme groups design
pdf("compareVr_EG.pdf", width = 5, height = 5)
par(parfix)
curve( R.ratio(x, 1 / p[1], 1 / p[1]), 0, 1,
       xlab = expression(rho), ylab = expression(V[p]^" r" / V[eg]^" r"),
       ylim = c(0,3), bty = "n")
for (i in 2:4) curve( R.ratio(x, 1 / p[i], 1 / p[i]), add = TRUE, lty = linetypes[i])
legend("topright", legend = parse(text=paste("p[1] == 1/",p, sep="")),
       lty=linetypes, bty="n", 
       horiz=FALSE, xjust = .5, yjust = .5, cex=1)
dev.off()


# Dichotomization design
pdf("compareVr_DI.pdf", width = 5, height = 5)
par(parfix)
curve( R.ratio(x, 1 / p[1], 1 - 1 / p[1]), 0, 1,
       xlab = expression(rho), ylab = expression(V[p]^" r" / V[eg]^" r"),
       ylim = c(0,1), bty = "n")
for (i in 2:4) curve( R.ratio(x, 1 / p[i], 1 - 1 / p[i]), add = TRUE, lty = linetypes[i])
legend("topright", legend = parse(text=paste("p[1] == 1/",p, sep="")),
       lty=linetypes, bty="n", 
       horiz=FALSE, xjust = .5, yjust = .5, cex=1)
dev.off()



##################################
## Plot V.Z.correct and V.Z.exp ##
##################################

op <- par(mfrow=c(1,2), oma = c(5,0,0,0), mar = c(2,3,4,1) + 0.1)
  
# Extreme groups design
curve( V.Z.Feldt(x, 1 / p[1], 1 / p[1], V=V.d.correct), 0, 1,
       xlab = expression(rho), ylab = expression(V[Z]),
       ylim = c(0,10), bty = "n", main="Extreme group design")
for (i in 2:4) curve(V.Z.Feldt(x, 1 / p[i], 1 / p[i], V=V.d.correct), add = TRUE, lty = linetypes[i])
for (i in 1:4) curve(V.Z.Feldt(x, 1 / p[i], 1 / p[i], V=V.d.normal), add = TRUE, lty = linetypes[i], col="blue")

# Dichotomization design
curve( V.Z.Feldt(x, 1 / p[1], 1 - 1 / p[1], V=V.d.correct), 0, 1,
       xlab = expression(rho), ylab = expression(V[Z]),
       ylim = c(0,10), bty = "n", main="Dichotomization design")
for (i in 2:4) curve(V.Z.Feldt(x, 1 / p[i], 1 - 1 / p[i], V=V.d.correct), add = TRUE, lty = linetypes[i])
for (i in 1:4) curve(V.Z.Feldt(x, 1 / p[i], 1 - 1 / p[i], V=V.d.normal), add = TRUE, lty = linetypes[i], col="blue")

# common legend
par(xpd=NA)
leg.pos <- cnvrt.coords(0.5, 0.07, input='tdev' )$usr
legend(leg.pos$x,leg.pos$y, legend = parse(text=paste("p[1] == over(1,",p,")", sep="")),
       lty=linetypes, bty="n", 
       horiz=TRUE, xjust = .5, yjust = .5, cex=1.2)
par(op)
