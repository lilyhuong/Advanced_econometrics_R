

#Resampling


############################### Pseudo random generator

## Uniform distribution

X <- runif(100000, min=1, max=2)
hist(X,nclass=100,probability=TRUE,col="grey")
u=seq(min(X),max(X),length.out =100)
lines(u,dunif(u, min=1, max=2), col="red")

## Normal distribution

X <- rnorm(100000, mean=2, sd=3)
hist(X,nclass=100,probability=TRUE,col="grey")
u=seq(min(X),max(X),length.out =100)
lines(u,dnorm(u,mean=2,sd=3), col="red")

## Generalized Pareto Distribution (GPD)

dgpd <- function(x,k,s,a) {(a/s)*(1+(x-k)/s)^(-a-1)}
pgpd <- function(q,k,s,a) {1-(1+(q-k)/s)^(-a)}
qgpd <- function(p,k,s,a) {s*(1-p)^(-1/a)-s+k}
rgpd <- function(n,k,s,a) {qgpd(runif(n),k,s,a)}

X <- rgpd(100000, k=1, s=2, a=2)
X=X[X<10]
hist(X,nclass=100,probability=TRUE,col="grey")
u=seq(min(X),max(X),length.out =100)
lines(u,dgpd(u, k=1, s=2, a=2), col="red")


############################### Monte Carlo experiments


## TCL: x from a mixture of Normals

nobs=10
nrepet=1e6
m=numeric(nrepet)
for(i in 1:nrepet) {
  p=runif(nobs)
  X=(p<=.7)*rnorm(nobs,-1,.4)+(p>.7)*rnorm(nobs,+1,.4)
  m[i]=mean(X)
}
hist(m,nclass=120,probability=TRUE,col="grey")
u=seq(min(m),max(m),length.out =100)
lines(u,dnorm(u, mean(m), sd(m)), col="red")


## TCL: x from a Lognormal

nobs=10
nrepet=1e6
m=numeric(nrepet)
for(i in 1:nrepet) {
  X=rlnorm(nobs,0,1)
  m[i]=mean(X)
}
hist(m,nclass=120,probability=TRUE,col="grey")
u=seq(min(m),max(m),length.out =100)
lines(u,dnorm(u, mean(m), sd(m)), col="red")


## MC experiment : t-stat in regression

nobs=100
nrepet=50000
t=numeric(nrepet)
set.seed(1)
X=rlnorm(nobs)  # if dummy: X=(runif(nobs)<.7)*1
for (i in 1:nrepet){
  if (i%%100==0) cat("i: ", i, "\r")
  eps=1-(runif(nobs)<.5)*2 
  y=1+2*X+eps
  a=summary(lm(y~X))
  t[i]=(a$coefficients[2]-2)/a$coefficients[4]
}
hist(t,nclass=120,probability=TRUE,col="grey")
u=seq(min(t),max(t),length.out =100)
lines(u,dnorm(u, 0, 1), col="red")

## MC experiment : HCCME t-stat in regression

library(sandwich)
nobs=100
nrepet=50000
t=numeric(nrepet)
set.seed(1)
X=rlnorm(nobs)  
for (i in 1:nrepet){
  if (i%%100==0) cat("i: ", i, "\r")
  eps=1-(runif(nobs)<.5)*2 
  y=1+2*X+eps
  m=summary(lm(y~X))
  b.se=sqrt(vcovHC(lm(y~X))[4])
  t[i]=(m$coefficients[2]-2)/b.se
}
hist(t,nclass=120,probability=TRUE,col="grey")
u=seq(min(t),max(t),length.out =100)
lines(u,dnorm(u, 0, 1), col="red")


########################### Bootstrap and Permutation tests


## Empirical Distribution Function EDF

X=sort(rnorm(30,0,1))
plot(ecdf(X), verticals = TRUE, do.points = FALSE)
points(X,seq(1/30,1,by=1/30))
lines(X,pnorm(X,0,1), col="red")


## Bootstrap distribution for the mean

nobs=30
X=rlnorm(nobs) 
nboot=50000
m=numeric(nboot)
for (i in 1:nboot){
  Xb=sample(X, replace=TRUE)  
  m[i]=mean(Xb)  
}
m=m-mean(m)+mean(X)  # recenter the bootstrap dist.
hist(m,nclass=120,probability=TRUE,col="grey")

abline(v=mean(exp(.5)), col="red") # population mean
points(mean(m),0, col="blue") # sample mean
u=seq(min(m),max(m), by=.01) # asymptotic dist
lines(u,dnorm(u,mean(X),sd(X)/sqrt(nobs)), lty=2)


## Bootstrap distribution for the median

nboot=50000
m=numeric(nboot)
set.seed(1)
for (i in 1:nboot){
  Xb=sample(X, replace=TRUE)  
  m[i]=median(Xb)  
}
m=m-mean(m)+median(X)  # recenter the bootstrap dist
hist(m,nclass=120,probability=TRUE,col="grey")
points(median(m),0, col="blue") # sample median
abline(v=1, col="red") # population median

## Bootstrap distribution: HCCME t-test

set.seed(1)
X=rlnorm(100) 
y=1+2*X+1-(runif(100)<.5)*2 
library(sandwich)
t=numeric(5000)
mod=summary(lm(y~X))
res=mod$residuals
for (i in 1:5000){
  resb=sample(res,replace=TRUE)
  yb=mod$coefficients[1]+2*X+resb
  modb=summary(lm(yb~X))
  b.se=sqrt(vcovHC(lm(yb~X))[4])
  t[i]=(modb$coefficients[2]-2)/b.se
}
hist(t,nclass=120,probability=TRUE,col="grey")
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
t0=mean(X)-mean(Y)
st0=(mean(X)-mean(Y))/sqrt(var(X)/n+var(Y)/m)
p.value= 2*pnorm(-abs(st0))
nperm=5000
t=numeric(nperm)
for (i in 1:nperm){
  Zp=sample(c(X,Y),replace=FALSE)
  t[i]=mean(Zp[1:n])-mean(Zp[(n+1):(n+m)])
}
hist(t,nclass=120,probability=TRUE,col="grey")
abline(v=t0, col="red")
p.value=sum(t>abs(t0))/nperm+sum(t<(-abs(t0)))/nperm


## Are permutation tests exact?

set.seed(1) ; n=30 ; m=30 ; nrepet=5000 ; nperm=500
p1=numeric(nrepet); p2=numeric(nrepet);t=numeric(nperm)
for(j in 1:nrepet) {   # MC experiment
  cat("j: ", j, "\r")
  X=rlnorm(n,0,1);   
  Y=rlnorm(m,0,1)
  s0=(mean(X)-mean(Y))/sqrt(var(X)/n+var(Y)/m)
  p1[j]= 2*pnorm(-abs(s0))  # asymptotic p-value
  for (i in 1:nperm){
    Zp=sample(c(X,Y),replace=FALSE)
    t[i]=mean(Zp[1:n])-mean(Zp[(n+1):(n+m)])
  }
  t0=mean(X)-mean(Y)        # permutation p-value
  p2[j]=sum(t>abs(t0))/nperm+sum(t<(-abs(t0)))/nperm
}
RP1=sum(p1<.05)/NROW(p1) #rejection proba asymp test
RP2=sum(p2<.05)/NROW(p2) #rejection proba permut test
c(RP1,RP2)

