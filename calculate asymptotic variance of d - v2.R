setwd("C:\\Users\\James\\Dropbox\\Meta-analysis\\d-to-r conversion\\simulations and code")
#setwd("C:\\Users\\jep701\\Dropbox\\Meta-analysis\\d-to-r conversion\\simulations and code")
source("r_feldt simulation functions.R")

rho <- 0.7
p1 <- 1 / 5
p2 <- 1 / 5
n1 <- 10
n2 <- p2 * n1 / p1
# n2 <- 30
f <- n1 / (n1 + n2)
  
# calculate cumulants of y1, y2

k1 <- cumulants(rho, p1)
k2 <- cumulants(rho, p2)

###############################################
## matrix form of delta-method variance of d ##
###############################################

M <- matrix(0, 4, 4)
M[1:2, 1:2] <- c(k1[2], k1[3], k1[3], k1[4] + 2 * k1[2]^2) / f
M[3:4, 3:4] <- c(k2[2], - k2[3], - k2[3], k2[4] + 2 * k2[2]^2) / (1 - f)


# gradient of d-function
delta.gradient <- function(y1, s1.sq, y2, s2.sq, f) {
	g.y1 <- - 1 / sqrt(f * s1.sq + (1 - f) * s2.sq)
	g.s1.sq <- -(f / 2) * (y2 - y1) / sqrt(f * s1.sq + (1 - f) * s2.sq)^3
	g.y2 <-  1 / sqrt(f * s1.sq + (1 - f) * s2.sq)
	g.s2.sq <- -((1 - f) / 2) * (y2 - y1) / sqrt(f * s1.sq + (1 - f) * s2.sq)^3
	return(c(g.y1, g.s1.sq, g.y2, g.s2.sq))
	}

g <- delta.gradient(k1[1], k1[2], -k2[1], k2[2], f)
t(g) %*% M %*% g

V.d.correct(rho, p1, p2, f)

#################################
## check special case formulas ##
#################################

# equal p, equal n

rho <- 0.9
p <- 1 / 10
u <- qnorm(p)
v <- dnorm(u) / p

delta <- r_to_d(rho, ab(p,p))
print(V.d.delta <- 4 + delta^2 / 2 - 
	delta^4 * (2 * v^2 + 3 * u * v + u^2 - 1) / (4 * v^2) - 
	delta^6 * (6 * v^3 + 12 * u * v^2 + 7 * u^2 * v - 4 * v - 3 * u + u^3) / (64 * v^3))
print(V.d.rho <- (4 - 2 * rho^2 * v * (6 * u + 5 * v) + 
	4 * rho^4 * v^2 * (2 * u^2 + 2 * u * v + 1) -
	rho^6 * u * v^3 * (u^2 + u * v + 1)) / (1 - rho^2 * v * (u + v))^3)
V.d.correct(rho, p, p)

# dichotomization design

rho <- 0.9
p1 <- 1 / 5
p2 <- 1 - p1
u1 <- qnorm(p1)
v1 <- dnorm(u1) / p1
u2 <- qnorm(p2)
v2 <- dnorm(u2) / p2

V.d.correct(rho, p1, p2)
print(V.d.rho <- (4 * (1 - p1)^2 - rho^2 * v1^2 * (6 * p1^2 - 6 * p1 + 4) - 
	4 * rho^2 * u1 * v1 * (1 - 2 * p1) * (1 - p1 + rho^2 * v1^2 * p1) +
	4 * rho^4 * v1^2 * p1 * (1 - p1) * (1 - u1^2) -
	rho^6 * u1^2 * v1^4 * p1^2) / (4 * p1 * (1 - p1 - rho^2 * v1^2 * p1)^3))


# median split

rho <- 0.9
p <- 0.5
u <- qnorm(p)
v <- dnorm(u) / p

delta <- r_to_d(rho, ab(p,p))
print(V.d.delta <- 4 + delta^2 / 2 - 
	delta^4 * (4 - pi) / 8 + delta^6 * (pi - 3) / 32)
print(V.d.rho <- 4 * (pi^3 - rho^2 * pi^2 * (5 - 2 * rho^2)) / (pi - 2 * rho^2)^3)
V.d.correct(rho, p, p)

