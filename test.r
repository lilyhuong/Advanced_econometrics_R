
dsinmad <- function(x, b, a, q) {(a * q / b) * ((x / b)^(a - 1)) * (1 + (x / b)^a )^(-q - 1)} # density
psinmad <- function(d, b, a, q)  {1 - (1 + (d / b)^a)^(-q)}  #distribution 
qsinmad <- function(p, b, a, q)  {b * ((1 - p)^( -1 / q) - 1)^(1 / a)} #quantile 
rsinmad <- function(n, b, a, q) {qsinmad(runif(n), b, a, q)}

X <- rsinmad(100000, b = 1, a = 2, q = 2)
X <- X[X < 5]
hist(X, nclass = 100, probability = TRUE, col = "grey")
u <- seq(min(X), max(X), length.out = 100)
lines(u, dsinmad(u, b = 1, a = 2, q = 2), col = "red")
