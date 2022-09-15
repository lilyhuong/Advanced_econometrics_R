
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



#15/09

## Bootstrap distribution: HCCME t-test

set.seed(1)
X = rlnorm(100) 
y = 1+2*X+1-(runif(100)<.5)*2 
library(sandwich) #this library we can use the fonction vovHC
t = numeric(50000)
mod=summary(lm(y~X))
res=mod$residuals
for (i in 1:50000){
  resb=sample(res,replace=TRUE) #vector of new bootstrap
  yb=mod$coefficients[1]+2*X+resb
  modb=summary(lm(yb~X))  #lm: fit the OLS
  b.se=sqrt(vcovHC(lm(yb~X))[4])
  t[i]=(modb$coefficients[2]-2)/b.se  #compute the t stat with true coeff 
}
hist(t,nclass=120,probability=TRUE,col="grey") #histogram of qtt of t
u=seq(min(t),max(t), by=.01) # asymptotic dist
lines(u,dnorm(u,0,1), lty=1, col="red")

b0.se=sqrt(vcovHC(lm(y~X))[4])
CI.inf=mod$coefficients[2]-quantile(t,0.975)*b0.se
CI.sup=mod$coefficients[2]-quantile(t,0.025)*b0.se

t0=(mod$coefficients[2]-2)/b0.se
p.value=sum(t>t0)/n    # unilateral p-value

n=NROW(t)              # bilateral p-value
p.value=sum(t>abs(t0))/n + sum(t<(-abs(t0)))/n


## Permutation tests for difference in means

set.seed(1) ; n=30 ; m=30
X=rlnorm(n,0,1)
Y=rlnorm(m,0,1)
t0=mean(X)-mean(Y)   #diff in mean
st0=(mean(X)-mean(Y))/sqrt(var(X)/n+var(Y)/m)  #standard deviation
p.value= 2*pnorm(-abs(st0))
nperm=5000
t=numeric(nperm)
for (i in 1:nperm){
  Zp=sample(c(X,Y),replace=FALSE)
  t[i]=mean(Zp[1:n])-mean(Zp[(n+1):(n+m)])  #diff in couple of sample 
}
hist(t,nclass=120,probability=TRUE,col="grey")
abline(v=t0, col="red")
p.value=sum(t>abs(t0))/nperm+sum(t<(-abs(t0)))/nperm