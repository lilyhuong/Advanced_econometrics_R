
## Uniform distribution
# These functions provide information about the uniform distribution 
                                    #on the interval from min to max
# dunif: density
# punif: distribution function
# qunif: quantile function
# runif: generates random numbers  of deviates



## Normal distribution

# same 4 cases with Uniform distribution

X <- rnorm(100000, mean = 2, sd = 3) # nolint
hist(X, nclass = 100, probability = TRUE, col = "grey")
u <- seq(min(X), max(X), length.out = 100)
lines(u, dnorm(u, mean = 2, sd = 3), col = "red")


## Generalized Pareto Distribution (GPD)

dgpd <- function(x, k, s, a) {(a / s) * (1 + (x - k) / s) ^ (- a - 1)} # nolint
pgpd <- function(q, k, s, a) {1 - (1 + (q - k) / s)^(-a)}
qgpd <- function(p, k, s, a) {s * (1 - p)^(-1 / a) - s + k}
rgpd <- function(n, k, s, a) {qgpd(runif(n), k, s, a)}

X <- rgpd(100000, k = 1, s = 2, a = 2)
X <- X[X < 10]
hist(X, nclass = 100, probability = TRUE, col = "grey")
u <- seq(min(X), max(X), length.out = 100)
lines(u, dgpd(u, k = 1, s = 2, a = 2), col = "red")



# Monte Carlo experiments


# x ∼ 0.7N(−1, .4) + 0.3N(1, .4)
nobs = 100  #nb of observations, change to 100, 500 to see the difference when we increase the nb of observations
nrepet = 1e6 #cycle repetition, repeat something
m = numeric(nrepet)
for(i in 1:nrepet) {
  p = runif(nobs)
  X = (p<=.7) * rnorm(nobs,-1,.4)+(p>.7) * rnorm(nobs,+1,.4)
  m[i] = mean(X)
}
hist(m, nclass = 120, probability = TRUE, col = "grey")
u = seq(min(m), max(m), length.out = 100)
lines(u, dnorm(u, mean(m), sd(m)), col = "red")


## TCL: x from a Lognormal: MC
# Note : rlnorm/ The multivariate lognormal distribution
#Generates random amounts with a multivariate lognormal distribution,
# or gives the density of that distribution at a given point.


nobs = 500
nrepet = 1e6
m = numeric(nrepet)
for(i in 1:nrepet) {
  X = rlnorm(nobs,0,1)
  m[i] = mean(X)
}
hist(m,nclass = 120,probability = TRUE,col = "grey")
u = seq(min(m), max(m), length.out = 100)
lines(u, dnorm(u, mean(m), sd(m)), col = "red")

## MC experiment : t-stat in regression

#set.seed: Setting a seed in R means to initialize a pseudorandom number generator. 
#Most of the simulation methods in Statistics require the possibility to generate pseudorandom numbers 
#that mimic the properties of independent generations of a uniform distribution in the interval (0, 1)(0,1)


nobs = 100
nrepet = 50000
t = numeric(nrepet)
set.seed(1)
X = rlnorm(nobs)  # if dummy: X=(runif(nobs)<.7)*1
for (i in 1:nrepet){
  if (i %% 100 == 0) cat("i: ", i, "\r")
  eps = 1-(runif(nobs)<.5) * 2 
  y = 1 + 2 * X + eps
  a = summary(lm(y~X))
  t[i] = (a$coefficients[2] - 2) / a$coefficients[4]
}
hist(t,nclass = 120,probability = TRUE,col = "grey")
u = seq(min(t), max(t),length.out = 100)
lines(u, dnorm(u, 0, 1), col = "red")


## MC experiment : HCCME t-stat in regression

library(sandwich)
nobs = 100
nrepet = 50000
t = numeric(nrepet)
set.seed(1)
X = rlnorm(nobs)  
for (i in 1:nrepet){
  if (i%%100 == 0) cat("i: ", i, "\r")
  eps = 1-(runif(nobs)<.5) * 2 
  y = 1 + 2 * X + eps
  m = summary(lm(y~X))
  b.se = sqrt(vcovHC(lm(y~X))[4])
  t[i] = (m$coefficients[2]-2) / b.se
}
hist(t,nclass = 120,probability = TRUE,col = "grey")
u=seq(min(t),max(t),length.out = 100)
lines(u,dnorm(u, 0, 1), col = "red")


