library(truncnorm)

#########################
## determine constants ##
#########################

ab <- function(p1, p2, f = p1 / (p1 + p2), n) {
	u1 <- qnorm(p1)
	v1 <- dnorm(u1) / p1
	u2 <- qnorm(p2)
	v2 <- dnorm(u2) / p2
  if (missing(n)) {
    a <- (v1 + v2)^2 / (f * v1 * (u1 + v1) + (1 - f) * v2 * (u2 + v2))
    b <- sqrt(1 / (f * v1 * (u1 + v1) + (1 - f) * v2 * (u2 + v2)))
  } else {
    f <- (n[1] - 1) / (sum(n) - 2)
    a <- ((v1 + v2)^2 + v1 * (u1 + v1) / n[1] + v2 * (u2 + v2) / n[2])/ (f * v1 * (u1 + v1) + (1 - f) * v2 * (u2 + v2))
    b <- sqrt(a + 1/n[1] + 1/n[2]) / (v1 + v2)
  }
	return(cbind(a, b))
}


##########################
## conversion functions ##
##########################

d_to_r <- function(d, const) {
	r <- const[2] * d / sqrt(d^2 + const[1])
	return(r)
	}

r_to_d <- function(r, const) {
  d <- sqrt(const[1]) * r / sqrt(const[2]^2 - r^2)
  return(d)
}


##########################################
## calculate cumulants for given rho, p ##
##########################################

cumulants <- function(rho, p) {
  u <- qnorm(p)
  v <- dnorm(u) / p
  k1 <- - rho * v
  k2 <- 1 - rho^2 * v * (u + v)
  k3 <- - rho^3  * v * (2 * v^2 + 3 * u * v + u^2 - 1)
  k4 <- - rho^4  * v * (6 * v^3 + 12 * u * v^2 + 7 * u^2 * v - 4 * v - 3 * u + u^3)
  return(c(k1, k2, k3, k4))
}


######################################################
## calculate asymptotic variance and approximations ##
######################################################

V.d.correct <- function(rho, p1, p2, f = p1 / (p1 + p2)) {
  k1 <- cumulants(rho, p1)
  k2 <- cumulants(rho, p2)
  sigma.d.sq <- f * k1[2] + (1 - f) * k2[2]
  V.d <- (k1[2] / f + k2[2] / (1 - f)) / sigma.d.sq - 
    (k1[1] + k2[1]) * (k1[3] + k2[3]) / sigma.d.sq^2 +
    (k1[1] + k2[1])^2 * (2 * f * k1[2]^2 + f * k1[4] + 2 * (1 - f) * k2[2]^2 + (1 - f) * k2[4]) / (4 * sigma.d.sq^3)
return(V.d)
}
V.d.correct <- Vectorize(V.d.correct, "rho")

V.d.normal <- function(rho, p1, p2, f = p1 / (p1 + p2)) {
  const <- ab(p1, p2, f)
  delta <- rho * sqrt(const[1]) / sqrt(const[2]^2 - rho^2)
  V.d <- 1 / (f * (1 - f)) + delta^2 / 2 
  return(V.d)
}


V.r.Feldt <- function(rho, p1, p2, f = p1 / (p1 + p2), V=V.d.correct) {
  const <- ab(p1, p2, f)
  V.r <- V(rho, p1, p2, f) * (const[2]^2  - rho^2)^3 / (const[1] * const[2]^4)
  return(V.r)
}

V.r.MVnormal <- function(rho) (1 - rho^2)^2

V.Z.Feldt <- function(rho, p1, p2, f = p1 / (p1 + p2), V=V.d.correct) {
  const <- ab(p1, p2, f)
  V.Z <- V(rho, p1, p2, f) * (const[2]^2  - rho^2)^3 / (const[1] * const[2]^4 * (1 - rho^2)^2)
  return(V.Z)
}


#############################
## Fisher Z transformation ##
#############################

Z.Fisher <- function(r) (log(1 + r) - log(1 - r)) / 2

Z.Taylor <- function(r0, b, k=5) { 
	Z.der <- rbind(Z.Fisher(r0),
			(1 - r0^2)^-1,
			2 * r0 * (1 - r0^2)^-2,
			(2 + 6 * r0^2) * (1 - r0^2)^-3,
			(24 * r0 + 24 * r0^3) * (1 - r0^2)^-4,
			(24 + 240 * r0^2 + 120 * r0^4) * (1 - r0^2)^-5)
	r0.powers <- sapply(r0, function(r) r^(0:k))
	b.powers <- matrix(((b - 1)^(0:k)) / factorial(0:k), k+1, length(r0))
	Z.Taylor <- colSums(b.powers * r0.powers * Z.der[1:(k+1),])
	return(Z.Taylor)
	}


#########################
## generate random d's ##
#########################


simulate_SMD <- function(rho, p1, p2, n1, n2, iterations, fixed=TRUE) {
  if (fixed) {
    x <- matrix(NA, n1 + n2, iterations)
    x[1:n1,] <- matrix(rtruncnorm(n1 * iterations, b = qnorm(p1)), n1, iterations)
    x[n1 + 1:n2,] <- matrix(rtruncnorm(n2 * iterations, a = qnorm(1 - p2)), n2, iterations)
  } else {
    n <- round((n1 + n2)/ (p1 + p2), 0)
    x <- replicate(iterations, sort(rnorm(n))[c(1:n1,(n-n2+1):n)])
  }
	e <- matrix(rnorm((n1 + n2) * iterations), n1 + n2, iterations)
	y <- rho * x + sqrt(1 - rho^2) * e
	
  if (n1 > 1) {
    s.p.sq <- ((n1 - 1) * apply(y[1:n1,], 2, var) + 
      (n2 - 1) * apply(y[n1 + 1:n2,], 2, var)) / (n1 + n2 - 2)
  } else {
    s.p.sq <- apply(y[n1 + 1:n2,], 2, var)
  }
  
	d <- (colMeans(y[n1 + 1:n2,]) - colMeans(y[1:n1,])) / sqrt(s.p.sq)
	return(d)
	}


#############################
## calculate summary stats ##
#############################
#rho <- 0.7
#p1 <- 1/4
#p2 <- 1/4
#n1 <- 500
#n2 <- 500
#iterations <- 5000
#fixed <- TRUE
#d <- simulate_SMD(rho, p1, p2, n1, n2, iterations, fixed)
#f <- n1 / (n1 + n2)
#truncation <- 1 - 10^-4
#k <- 5

stat.names <- function() c("d", "r.i", "InRange", "r.n", "z.iF", "z.iT", "z.nF", "z.nT", 
                           "Vd.asy", "Vd.exp","Vr.asy", "Vr.exp", "Vr.BVN", "VZ.asy", "VZ.exp")

analyze <- function(d, p1, p2, n1, n2, f = n1 / (n1 + n2), truncation = 1 - 10^-4, k = 5) {
  
	iterations <- length(d)
  
	A <- matrix(NA, length(stat.names()), iterations, dimnames = list(stat=stat.names()))
  
	c.inf <- ab(p1, p2, f=f)
	c.n <- ab(p1, p2, n=c(n1,n2))
  
  A["d",] <- d
	A["r.i",] <- d_to_r(d, c.inf)
  r.i.trunc <- pmin(pmax(A["r.i",], -truncation), truncation)
  A["InRange",] <- abs(A["r.i",]) < 1
	A["r.n",] <- d_to_r(d, c.n)
	A["z.iF",] <- Z.Fisher(r.i.trunc)
	A["z.nF",] <- Z.Fisher(pmin(pmax(A["r.n",], -truncation), truncation))
	A["z.iT",] <- Z.Taylor(d_to_r(d, c(c.inf[1], 1)), c.inf[2])
	A["z.nT",] <- Z.Taylor(d_to_r(d, c(c.n[1], 1)), c.n[2])
	A["Vd.asy",] <- V.d.correct(r.i.trunc, p1, p2, f) / (n1 + n2)
	A["Vd.exp",] <- V.d.normal(r.i.trunc, p1, p2, f) / (n1 + n2)
	A["Vr.asy",] <- V.r.Feldt(r.i.trunc, p1, p2, f, V=V.d.correct) / (n1 + n2)
	A["Vr.exp",] <- V.r.Feldt(r.i.trunc, p1, p2, f, V=V.d.normal) / (n1 + n2)
  A["Vr.BVN",] <- V.r.MVnormal(r.i.trunc) / (n1 + n2)
  A[c("VZ.asy","VZ.exp"),] <- matrix(rep((1 - A["r.i",]^2)^-2, each=2), 2, iterations) * 
                                      A[c("Vr.asy","Vr.exp"),]
  
	summary <- matrix(NA, length(stat.names()), 2, dimnames = list(stat=stat.names(), summary=c("mean","var")))
	for (f in c("mean", "var")) summary[,f] <- apply(A, 1, f)
	
  return(summary)
	}