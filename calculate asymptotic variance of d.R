################
## simulate d ##
################

rho <- 0.7
p1 <- 1 / 4
p2 <- 3 / 4
n <- 1000
iterations <- 10000

x <- array(
	c(rtruncnorm(n * p1 * iterations, b = qnorm(p1)), - rtruncnorm(n * p2 * iterations, b = qnorm(1 - p2))),
	dim = c(n * p1, iterations, 2))
e <- array(rnorm(2 * n * p1 * iterations), dim = c(n * p1, iterations, 2))
y <- rho * x + sqrt(1 - rho^2) * e

y.mean <- apply(y, c(2,3), mean)
y.var <- apply(y, c(2,3), var)

delta <- function(y1, s1.sq, y2, s2.sq, f1) (y2 - y1) / sqrt(f1 * s1.sq + (1 - f1) * s2.sq)

d <- mapply(delta, y1 = y.mean[,1], s1.sq = y.var[,1], 
			 y2 = y.mean[,2], s2.sq = y.var[,2], 
			MoreArgs = list(f1 = p1 / (p1 + 1 - p2)))


###################################
## calculate cumulants of y1, y2 ##
###################################

cumulants <- function(rho, p) {
	u <- qnorm(p)
	v <- dnorm(u) / p
	k1 <- - rho * v
	k2 <- 1 - rho^2 * v * (u + v)
	k3 <- - rho^3  * v * (2 * v^2 + 3 * u * v + u^2 - 1)
	k4 <- - rho^4  * v * (6 * v^3 + 12 * u * v^2 + 7 * u^2 * v - 4 * v - 3 * u + u^3)
	return(c(k1, k2, k3, k4))
	}

cumulants(rho, 0.5)

k1 <- cumulants(rho, p1)
k2 <- cumulants(rho, 1 - p2)

###############################################
## matrix form of delta-method variance of d ##
###############################################

M <- matrix(0, 4, 4)
M[1:2, 1:2] <- c(k1[2], k1[3], k1[3], k1[4] + 2 * k1[2]^2) / p1
M[3:4, 3:4] <- c(k2[2], - k2[3], - k2[3], k2[4] + 2 * k2[2]^2) / (1 - p2)


# compare to simulated variances

y.sums <- cbind(y.mean, y.var)
colMeans(y.sums)
c(k1[1], -k2[1], k1[2], k2[2])
var(y.sums)[c(1,3,2,4),c(1,3,2,4)] * n
M


# gradient of d-function
delta.gradient <- function(y1, s1.sq, y2, s2.sq, f1) {
	g.y1 <- - 1 / sqrt(f1 * s1.sq + (1 - f1) * s2.sq)
	g.s1.sq <- -(f1 / 2) * (y2 - y1) / sqrt(f1 * s1.sq + (1 - f1) * s2.sq)^3
	g.y2 <-  1 / sqrt(f1 * s1.sq + (1 - f1) * s2.sq)
	g.s2.sq <- -((1 - f1) / 2) * (y2 - y1) / sqrt(f1 * s1.sq + (1 - f1) * s2.sq)^3
	return(c(g.y1, g.s1.sq, g.y2, g.s2.sq))
	}

g <- delta.gradient(k1[1], k1[2], -k2[1], k2[2], p1 / (p1 + 1 - p2))
t(g) %*% M %*% g


################################
## compare to formula for V.d ##
################################

sigma.d.sq <- p1 * k1[2] + (1 - p2) * k2[2]
V.d <- (p1 + 1 - p2) * (
	((1 - p2) * k1[2] + p1 * k2[2]) / (p1 * (1 - p2) * sigma.d.sq) - 
	(k1[1] + k2[1]) * (k1[3] + k2[3]) / sigma.d.sq^2 +
	(k1[1] + k2[1])^2 * (p1 * (2 * k1[2]^2 + k1[4]) + (1 - p2) * (2 * k2[2]^2 + k2[4])) / (4 * sigma.d.sq^3)
	)

V.d
var(d) * n

(4 * k1[2]^3 - 4 * k1[1] * k1[2] * k1[3] + k1[1]^2 * (k1[4] + 2 * k1[2]^2)) / (2 * p1 * k1[2]^3)

u <- qnorm(p1)
v <- dnorm(u) / p1
(4 - 2 * rho^2 * v * (6 * u + 5 * v) + 4 * rho^4 * v^2 * (2 * u^2 + 2 * u * v + 1) - rho^6 * u * v^3 * (u^2 + u * v + 1)) / (2 * p1 * (1 - rho^2 * v * (u + v))^3)


#################################
## check special case formulas ##
#################################

rho <- 0.9
p1 <- 1 / 10
p2 <- 1 / 10

k1 <- cumulants(rho, p1)
k2 <- cumulants(rho, 1 - p2)

M <- matrix(0, 4, 4)
M[1:2, 1:2] <- c(k1[2], k1[3], k1[3], k1[4] + 2 * k1[2]^2) / p1
M[3:4, 3:4] <- c(k2[2], - k2[3], - k2[3], k2[4] + 2 * k2[2]^2) / (1 - p2)
g <- delta.gradient(k1[1], k1[2], -k2[1], k2[2], p1 / (p1 + 1 - p2))
t(g) %*% M %*% g

sigma.d.sq <- p1 * k1[2] + (1 - p2) * k2[2]
V.d <- (p1 + 1 - p2) * (
	((1 - p2) * k1[2] + p1 * k2[2]) / (p1 * (1 - p2) * sigma.d.sq) - 
	(k1[1] + k2[1]) * (k1[3] + k2[3]) / sigma.d.sq^2 +
	(k1[1] + k2[1])^2 * (p1 * (2 * k1[2]^2 + k1[4]) + (1 - p2) * (2 * k2[2]^2 + k2[4])) / (4 * sigma.d.sq^3)
	)
V.d

u <- qnorm(p1)
v <- dnorm(u) / p1
(1 - rho^2 * v * (u * (1 - 2 * p1) / (1 - p1) + v * (1 - 3 * p1 + 3 * p1^2) / (1 - p1)^2)) / (p1 * (1 - p1 - rho^2 * v^2 * p1)) - 
rho^4  * v^2 * (2 * v^2 * (1 - 3 * p1 + 3 * p1^2) / (1 - p1)^2 + 3 * u * v * (1 - 2 * p1) / (1 - p1) + u^2 - 1) / (1 - p1 - rho^2 * v^2 * p1)^2 + 
rho^2 * v^2 * (1 - p1 - 2 * rho^2 * v^2 * p1 + rho^4 * v^2 * p1 * (u^2 + 2 * u * v * (1 - 2 * p1) / (1 - p1) + v^2 * (1 - 3 * p1 + 3 * p1^2) / (1 - p1)^2)) / (2 * (1 - p1 - rho^2 * v^2 * p1)^3) -
rho^6 * v^4 * p1 * (6 * v^2 * (1 - 3 * p1 + 3 * p1^2) / (1 - p1)^2 + 12 * u * v * (1 - 2 * p1) / (1 - p1) + 7 * u^2 - 4 ) / (4 * (1 - p1 - rho^2 * v^2 * p1)^3)


(4 * (1 - p1)^2 - rho^2 * v^2 * (6 * p1^2 - 6 * p1 + 4) - 
4 * rho^2 * u * v * (1 - p1) * (1 - 2 * p1) + 4 * rho^4 * v^2 * p1 * (1 - p1) - 
4 * rho^4 * u * v^3 * p1 * (1 - 2 * p1) - 4 * rho^4 * u^2 * v^2 * p1 * (1 - p1) - 
rho^6 * u^2 * v^4 * p1^2) /
(4 * p1 * (1 - p1 - rho^2 * v^2 * p1)^3)

######################################################################
## compare delta-method formula using V.d versus using V.r directly ##
######################################################################


ab.n <- function(p1, p2, n) {
	u1 <- qnorm(p1)
	v1 <- dnorm(u1) / p1
	u2 <- qnorm(1 - p2)
	v2 <- dnorm(u2) / (1 - p2)
	a.n <- ((p1 + 1 - p2) * n - 2) * ((v1 + v2)^2 - v1 * (u1 + v1) / (p1 * n) - v2 * (u2 + v2) / ((1 - p2) * n)) / 
		((p1 * n - 1) * v1 * (u1 + v1) + ((1 - p2) * n - 1) * v2 * (u2 + v2))
	b.n <- sqrt(a.n) / (v1 + v2)
	return(c(a.n, b.n))
}

r <- function(y1, s1.sq, y2, s2.sq, p1, p2, n) {
	ab <- ab.n(p1, p2, n)
	r <- ab[2] * (y2 - y1) / sqrt((y2 - y1)^2 + ab[1] * (p1 * s1.sq + (1 - p2) * s2.sq) / (p1 + 1 - p2))
	return(r)
	}

r.gradient <- function(y1, s1.sq, y2, s2.sq, p1, p2, n) {
	ab <- ab.n(p1, p2, n)
	f <- (n * p1 - 1) / ((p1 + 1 - p2) * n - 2)
	g.y1 <- - ab[2] / sqrt((y2 - y1)^2 + ab[1] * (f * s1.sq + (1 - f) * s2.sq)) + ab[2] * (y2 - y1)^2 / ((y2 - y1)^2 + ab[1] * (f * s1.sq + (1 - f) * s2.sq))^(3 / 2)
	g.s1.sq <- - (ab[1] * ab[2] * f / 2) * (y2 - y1) / ((y2 - y1)^2 + ab[1] * (f * s1.sq + (1 - f) * s2.sq))^(3 / 2)
	g.y2 <-  - g.y1
	g.s2.sq <- g.s1.sq * (1 - f) / f
	return(c(g.y1, g.s1.sq, g.y2, g.s2.sq))
	}

n <- 10^10
rho <- 0.9
p1 <- 1 / 4
p2 <- 3 / 4

k1 <- cumulants(rho, p1)
k2 <- cumulants(rho, 1 - p2)

M <- matrix(0, 4, 4)
M[1:2, 1:2] <- c(k1[2], k1[3], k1[3], k1[4] + 2 * k1[2]^2) / p1
M[3:4, 3:4] <- c(k2[2], - k2[3], - k2[3], k2[4] + 2 * k2[2]^2) / (1 - p2)
delta.g <- delta.gradient(k1[1], k1[2], -k2[1], k2[2], p1 / (p1 + 1 - p2))

(rho^2 * (p1 * k1[2] + (1 - p2) * k2[2])^3 / 
((k1[1] + k2[1])^2 * (p1 + 1 - p2)^3)) * 
t(delta.g) %*% M %*% delta.g

rho.g <- r.gradient(k1[1], k1[2], -k2[1], k2[2], p1, p2, n)
t(rho.g) %*% M %*% rho.g
