source("r_feldt simulation functions.R")

p <- 0.25
n <- 35
F.stat <- 18.2

print(d <- sqrt(F.stat * 2 / n))

print(u <- qnorm(p))
print(v <- dnorm(u) / p)

print(a <- 4 * v / (u + v))
print(b <- 1 / sqrt(v * (u + v)))

print(r.F <- b * d / sqrt(d^2 + a))
r.F / b
r.F - r.F / b

print(Z.F <- 0.5 * (log(1 + r.F) - log(1 - r.F)), digits = 12)
sapply(1:5, Z.Taylor, r0 = r.F / b, b = b)
print(Z.T <- Z.Taylor(r.F / b, b, k=2), digits = 12)
print(Z.F - Z.T)

print(V.d <- (4 + d^2 / 2 - 
	d^4 * (2 * v^2 + 3 * u * v + u^2 - 1) / (4 * v^2) -
	d^6 * (6 * v^3 + 12 * u * v^2 + 7 * u^2 * v - 4 * v - 3 * u + u^3) / (64 * v^3)) /
	(2 * n))

(4 - 2 * r.F^2 * v * (6 * u + 5 * v) + 
	4 * r.F^4 * v^2 * (2 * u^2 + 2 * u * v + 1) - 
	r.F^6 * u * v^3 * (u^2 + u * v + 1)) / 
	(2 * n * (1 - r.F^2 * v * (u + v))^3)

V.d.correct(r.F, p, p) / (2 * n)

print(V.d.norm <- 2 / n + d^2 / (4 * n))
V.d.normal(r.F, p, p) / (2 * n)


V.r.Feldt(r.F, p, p) / (2 * n)
V.Z.Feldt(r.F, p, p) / (2 * n)

sqrt(V.r.Feldt(r.F, p, p) / (2 * n))
sqrt(V.Z.Feldt(r.F, p, p) / (2 * n))

(1 - r.F^2)^2 / (2 * n)
sqrt((1 - r.F^2)^2 / (2 * n))
(1 - r.F^2)^2 / V.r.Feldt(r.F, p, p)
sqrt((1 - r.F^2)^2 / (2 * n)) / sqrt(V.r.Feldt(r.F, p, p) / (2 * n))

1 / (2 * n)
1 / V.Z.Feldt(r.F, p, p)
