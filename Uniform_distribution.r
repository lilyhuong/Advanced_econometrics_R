
## Uniform distribution
# These functions provide information about the uniform distribution 
                                    #on the interval from min to max
# dunif: density
# punif: distribution function
# qunif: quantile function
# runif: generates random numbers  of deviates

X <- runif(100000, min = 1, max = 2) # nolint #create 100 000 numbers from 1 to 2
hist(X, nclass = 100, probability = TRUE, col = "grey") #creat a histogram from X # nolint
u = seq(min(X), max(X), length.out = 100) # nolint
u <- seq(min(X), max(X), length.out = 100) # should use this form
# seq: create a sequence of elements in a Vector.
#       It takes the length and difference between values as optional argument.
lines(u, dunif(u, min = 1, max = 2), col = "red") #get the line of density


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
X = X[X < 10]
hist(X, nclass = 100, probability = TRUE, col = "grey")
u = seq(min(X), max(X), length.out = 100)
lines(u, dgpd(u, k = 1, s = 2, a = 2), col = "red")
