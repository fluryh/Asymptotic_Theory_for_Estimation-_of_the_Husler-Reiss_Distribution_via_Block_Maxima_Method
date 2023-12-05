library(MASS)
library(ggplot2)
library(evd)
library(beepr)
set.seed(5838)
bm608 <- 2.96765
bm37138<-4.05082
bm1770297 <- 4.875154944523863
p <- function(l,m){
  return(1-l^2/log(m))
}
p2 <- function(l,bm){
    return((1-l^2/bm^2)/(1+l^2/bm^2))
}
p3 <- function(l,bm){
  return(.9 + .1*p2(l/sqrt(.1),bm))
}
grad <- function(l,X,Y){
  sum <- 0
  d <- length(X)
  for(i in 1:d){
    x <- X[i]
    y <- Y[i]
    sum <- sum -2*exp(-x)*dnorm(l + (y-x)/(2*l)) + exp(-x)*dnorm(l + (y-x)/(2*l))*(-1/2 - 1/(2*l^2)+(y-x)^2/(8*l^4) + pnorm(l + (x-y)/(2*l))*exp(-y)*(1-(y-x)/(2*l^2)) + pnorm(l + (y-x)/(2*l))*exp(-x)*(1-(x-y)/(2*l^2)))/(exp(-x-y)*pnorm(l + (x-y)/(2*l))*pnorm(l + (y-x)/(2*l)) + exp(-x)*dnorm(l + (y-x)/(2*l))/(2*l))
  }
  return(sum)
}

log_liks <- rep(0,1000)
for(j in 1:1000){
  samps <- matrix(data = rep(0,2*269),ncol= 2) 
  for(k in 1:269){
    XY_scaled <- apply(mvrnorm(37138, mu = c(0,0), Sigma = matrix(c(1, p2(1.25,bm37138), p2(1.25,bm37138), 1), ncol=2))
                       ,2,max)
    samps[k,] <- XY_scaled
  }
  log_liks[j] <- grad(1.25, samps[,1], samps[,2])/sqrt(269)
}
mean(log_liks)
hist(log_liks)
var(log_liks)

loghr <- rbvevd(100000,1/1.25, model = "hr")
# loghr <- rmev(100000,d=2,sigma = matrix(data = c(0,1/1.25^2,1/1.25^2,0),nrow=2), model = "hr")
grad(1.25,loghr[,1],loghr[,2])/100000
beep()

sum <- matrix(data = rep(0,2*1000),ncol= 2)
for(k in 1:1000){
  XY_scaled <- apply(mvrnorm(1770297, mu = c(-bm1770297^2,-bm1770297^2), Sigma = bm1770297^2*matrix(c(1, p2(1.25,bm1770297), p2(1.25,bm1770297), 1), ncol=2))
                     ,2,max)
  sum[k,] <- XY_scaled
}
beep()