library(MASS)
library(ggplot2)
library(evd)
library(beepr)
set.seed(5838)
bm608 <- 2.96765
bm5938<-4.05082
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
    sum <- sum -2*exp(-x)*dnorm(l + (y-x)/(2*l)) + exp(-x)*dnorm(l + (x-y)/(2*l))*(-1/2 - 1/(2*l^2)+(y-x)^2/(8*l^4) + pnorm(l + (x-y)/(2*l))*exp(-y)*(1-(y-x)/(2*l^2)) + pnorm(l + (y-x)/(2*l))*exp(-x)*(1-(x-y)/(2*l)))/(exp(-x-y)*pnorm(l + (x-y)/(2*l))*pnorm(l + (y-x)/(2*l)) + exp(-x)*dnorm(l + (y-x)/(2*l))/(2*l))
  }
  return(sum)
}

log_liks <- rep(0,100)
for(j in 1:100){
  samps <- matrix(data = rep(0,2*565),ncol= 2) 
  for(k in 1:565){
    XY_scaled <- apply(mvrnorm(1770297, mu = c(-bm1770297^2,-bm1770297^2), Sigma = bm1770297^2*matrix(c(1, p2(1.25,bm1770297), p2(1.25,bm1770297), 1), ncol=2))
                       ,2,max)
    samps[k,] <- XY_scaled
  }
  print(apply(samps,2,var))
  print(apply(samps,2,mean))
  # samps <- samps 
  # log_liks[j] <- sum/sqrt(565)
}
mean(log_liks)
hist(log)
var(log_liks)

loghr <- rbvhr(10000,1/1.332265)
# XY <- exp(loghr)
grad(1.25,loghr[,1],loghr[,2])/10000
beep()

sum <- matrix(data = rep(0,2*1000),ncol= 2)
for(k in 1:1000){
  XY_scaled <- apply(mvrnorm(1770297, mu = c(-bm1770297^2,-bm1770297^2), Sigma = bm1770297^2*matrix(c(1, p2(1.25,bm1770297), p2(1.25,bm1770297), 1), ncol=2))
                     ,2,max)
  sum[k,] <- XY_scaled
}
beep()