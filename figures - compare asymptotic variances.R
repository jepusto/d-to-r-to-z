source("r_feldt simulation functions.R")

library(TeachingDemos)


########################################
## plot a, b constants for varying p1 ##
########################################

pdf("abConst.pdf", width = 10, height = 5)
op <- par(mfrow=c(1,2))

# a constant
curve( ab(x, x)[,1], 0, 0.5, xlab=expression(p[1]), ylab = "",
      bty="n", main = expression(a(p[1],p[2])))
curve( ab(x, 1 - x)[,1], add=TRUE, lty="dashed")
text(0.25, 10, expression(p[2]==p[1]))
text(0.1, 8, expression(p[2]==1 - p[1]))

# b constant
curve( ab(x, x)[,2], 0, 0.5, xlab=expression(p[1]), ylab = "",
       bty="n", main = expression(b(p[1],p[2])), ylim=c(1, 2))
curve( ab(x, 1 - x)[,2], add=TRUE, lty="dashed")
text(0.25, 1.10, expression(p[2]==p[1]))
text(0.25, 1.45, expression(p[2]==1 - p[1]))

par(op)
dev.off()



##################################
## compare V.d.norm to true V.d ##
##################################

p <- c(2,3,5,8)
linetypes = c("solid","dashed", "dotted","dotdash")

D.ratio <- Vectorize(function(rho, p1, p2) V.d.normal(rho, p1, p2) / V.d.correct(rho, p1, p2), "rho")

pdf("compareVd.pdf", width = 10, height = 5)
op <- par(mfrow=c(1,2), oma = c(2,0,0,6), mar = c(2,3,4,1) + 0.1)

# Extreme groups design
curve( D.ratio(x, 1 / p[1], 1 / p[1]), 0, 1,
	xlab = expression(rho), ylab = expression(V["d,normal"] / V[d]),
	ylim = c(0.5,2), bty = "n", main="Extreme group design")
for (i in 2:4) curve( D.ratio(x, 1 / p[i], 1 / p[i]), add = TRUE, lty = linetypes[i])

# Dichotomization design
curve( D.ratio(x, 1 / p[1], 1 - 1 / p[1]), 0, 1,
	xlab = expression(rho), ylab = "",
	ylim = c(0.5,2), bty = "n", main="Dichotomization design")
for (i in 2:4) curve( D.ratio(x, 1 / p[i], 1 - 1 / p[i]), add = TRUE, lty = linetypes[i])

# common legend
par(xpd=NA)
leg.pos <- cnvrt.coords(0.9, 0.50, input='tdev' )$usr
legend(leg.pos$x,leg.pos$y, legend = parse(text=paste("p[1] == 1/",p, sep="")),
	lty=linetypes, bty="n", 
	horiz=FALSE, xjust = .5, yjust = .5, cex=1.1)
par(op)
dev.off()



#################################
## compare V.r.MVN to true V.r ##
#################################
R.ratio <- Vectorize(function(rho, p1, p2) V.r.MVnormal(rho) / V.r.Feldt(rho, p1, p2), "rho")


pdf("compareVr.pdf", width = 10, height = 5)
op <- par(mfrow=c(1,2), oma = c(2,0,0,6), mar = c(2,3,4,1) + 0.1)

# Extreme groups design
curve( R.ratio(x, 1 / p[1], 1 / p[1]), 0, 1,
       xlab = expression(rho), ylab = expression(V["d,normal"] / V[d]),
       ylim = c(0,3), bty = "n", main="Extreme group design")
for (i in 2:4) curve( R.ratio(x, 1 / p[i], 1 / p[i]), add = TRUE, lty = linetypes[i])

# Dichotomization design
curve( R.ratio(x, 1 / p[1], 1 - 1 / p[1]), 0, 1,
       xlab = expression(rho), ylab = "",
       ylim = c(0,1), bty = "n", main="Dichotomization design")
for (i in 2:4) curve( R.ratio(x, 1 / p[i], 1 - 1 / p[i]), add = TRUE, lty = linetypes[i])

# common legend
par(xpd=NA)
leg.pos <- cnvrt.coords(0.9, 0.50, input='tdev' )$usr
legend(leg.pos$x,leg.pos$y, legend = parse(text=paste("p[1] == 1/",p, sep="")),
       lty=linetypes, bty="n", 
       horiz=FALSE, xjust = .5, yjust = .5, cex=1.1)
par(op)
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
