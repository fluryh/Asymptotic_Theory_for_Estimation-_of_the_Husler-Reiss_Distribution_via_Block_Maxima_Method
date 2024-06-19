library(MASS)
library(ggplot2)
library(evd)
library(progress)
library(beepr)
library(latex2exp)

set.seed(5838)
bm608 <- 2.96765
bm10221 <- 3.7400967411876915
bm37138<-4.05082
bm976624 <- 4.756753707412248
bm1770297 <- 4.875154944523863
bm1044 <- 3.12796

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
grad2 <- function(l,X,Y){
  sum <- 0
  d <- length(X)
  for(i in 1:d){
    x <- X[i]
    y <- Y[i]
    t1 <- 2*exp(-x)*dnorm(l + (y-x)/(2*l))*(l- (y-x)^2/(4*l^3))
    denom <- exp(-x-y)*pnorm(l + (x-y)/(2*l))*pnorm(l + (y-x)/(2*l)) + 
      exp(-x)*pnorm(l + (y-x)/(2*l))/(2*l)
    t2 <- exp(-x)*dnorm(l + (y-x)/(2*l)) * (2*exp(-x)*dnorm(l + (y-x)/(2*l))*(1 - (x-y)^2/(2*l^2)) +
                                              exp(-x)*pnorm(l + (y-x)/(2*l))*(-l+(x-y)/(2*l) - (y-x)^2/(4*l^3) - (x-y)^3/(8*l^5) + (x-y)/(l^3)) +
                                              exp(-y)*pnorm(l + (x-y)/(2*l))*(-l+(y-x)/(2*l) - (x-y)^2/(4*l^3) - (y-x)^3/(8*l^5) + (y-x)/(l^3)))
    t2 <- t2/denom
    t3 <- exp(-x)*dnorm(l + (y-x)/(2*l)) *(-l/2 - 1/(2*l)- 1/(l^3) + (y-x)^2/(4*l^3) + .625 * (y-x)^2/(l^5) - (y-x)^4/(32*l^7))
    t3 <- t3/denom
    t4 <- exp(-x)*dnorm(l + (y-x)/(2*l)) *(exp(-x)*pnorm(l + (y-x)/(2*l))*(1-(x-y)/(2*l^2)) +
                                             exp(-y)*pnorm(l + (x-y)/(2*l))*(1-(y-x)/(2*l^2)) -
                                             1/(2*l^2) - .5 + (x-y)^2/(8*l^4))
    t4 <- (t4/denom)^2 
    sum = sum + t1 + t2 + t3 + t4
  }
  return(sum)
}
log_lik <- function(l,X,Y){
  sum <- 0
  d <- length(X)
  for(i in 1:d){
    x <- X[i]
    y <- Y[i]
    sum <- sum - exp(-x)*pnorm(l + (y-x)/(2*l)) - exp(-y)*pnorm(l + (x-y)/(2*l)) - x +
      log(exp(-y)*pnorm(l + (y-x)/(2*l))*pnorm(l+(x-y)/(2*l)) + dnorm(l + (y-x)/(2*l))/(2*l))
  }
  return(-sum)
}
log_lik2 <- function(l,X,Y){
  data = cbind(X,Y)
  return(sum(dbvhr(data,l,log=T)))
}

log_liks <- rep(0,10000)
l_hats <- rep(0,10000)
pb <- progress_bar$new(total = 10000)
numdiff <- rep(0,10000)
for(j in 1:10000){
  pb$tick()
  samps <- matrix(data = rep(0,2*269),ncol= 2) 
  for(k in 1:269){
    XY_scaled <- apply(mvrnorm(37138, mu = c(-bm37138^2,-bm37138^2), Sigma = bm37138^2*matrix(c(1, p2(1.25,bm37138), p2(.75,bm37138), 1), ncol=2))
                       ,2,max)
    samps[k,] <- XY_scaled
  }
  log_liks[j] <- grad(1.25, samps[,1], samps[,2])/sqrt(269)
  l_hats[j] <- optimize(log_lik, interval = c(0,10), X = samps[,1], Y = samps[,2], maximum = T)$maximum
}
mean(log_liks)
median(log_liks)
var(log_liks)

(mean(l_hats) - 1.25) * sqrt(269)
(median(l_hats) - 1.25) * sqrt(269)
var(l_hats) * 269





lhatplot <- ggplot() + 
  geom_histogram(aes(x = l_hats, y = after_stat(density)),
                 breaks = seq(.5, 2, by = .025), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(l_hats), sd = sd(l_hats)), linewidth= 1.25) +
  stat_function(fun = dnorm, args = list(mean = 1.25 + (0.0107919/0.222535)/sqrt(269), sd = sqrt(1/(0.222535*269))), linewidth= 1.25, linetype = 5) +
  labs(x= TeX(r"($\hat{\lambda}$)"), title = TeX(r"($m = 37138$, $k = 269$, $\rho=0.826$)")) + 
  theme_bw()

llplot <- ggplot() + 
  geom_histogram(aes(x = log_liks, y = after_stat(density)),
                 breaks = seq(-2, 2, by = .1), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(log_liks), sd = sd(log_liks)), linewidth = 1.25) +
  stat_function(fun = dnorm, args = list(mean = -0.0107919, sd = sqrt(0.222535)), linewidth = 1.25, linetype = 5) +
  labs(x= "Score Function Value", title = TeX(r"($m = 37138$, $k = 269$, $\rho=0.826$)")) + 
  theme_bw()

pdf("hrplots.pdf")
lhatplot
llplot
dev.off()

beep()
lhatplot
llplot

# pdf("hrqqplots.pdf")
# ggplot() +
#   stat_qq(aes(sample = l_hats)) +
#   stat_qq_line(aes(sample = l_hats)) +
#   labs(x = "",y = "", title = TeX(r"($\hat{\lambda}$ QQ Plot)")) +
#   theme_bw()
#   
# ggplot() +
#   stat_qq(aes(sample = log_liks)) +
#   stat_qq_line(aes(sample = log_liks)) +
#   labs(x = "",y = "", title = TeX(r"(Score QQ Plot)")) +
#   theme_bw()
# dev.off()



log_liks2 <- rep(0,10000)
l_hats2 <- rep(0,10000)
grads2 <- matrix(nrow = 10000, ncol = 5, data = rep(0,50000))
pb2 <- progress_bar$new(total = 10000)
for(j in 1:10000){
  pb2$tick()
  samps <- matrix(data = rep(0,2*987),ncol= 2) 
  for(k in 1:987){
    XY_scaled <- apply(mvrnorm(10221, mu = c(-bm10221^2,-bm10221^2), Sigma = bm10221^2*matrix(c(1, p2(1.25,bm10221), p2(1.25,bm10221), 1), ncol=2))
                       ,2,max)
    samps[k,] <- XY_scaled
  }
  log_liks2[j] <- grad(1.25, samps[,1], samps[,2])/sqrt(987)
  l_hats2[j] <- optimize(log_lik, interval = c(0,10), X = samps[,1], Y = samps[,2], maximum = T)$maximum
}
mean(log_liks2)
median(log_liks2)
var(log_liks2)

(mean(l_hats2) - 1.25) * sqrt(987)
(median(l_hats2) - 1.25) * sqrt(987)
var(l_hats2) * 987





lhatplot2 <- ggplot() + 
  geom_histogram(aes(x = l_hats2, y = after_stat(density)),
                 breaks = seq(.5, 2, by = .025), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = 1.25 + (5 * 0.0107919/0.222535)/sqrt(987), sd = sqrt(1/(0.222535*987))), linetype = 5, linewidth = 1.25) +
  stat_function(fun = dnorm, args = list(mean = mean(l_hats2), sd = sd(l_hats2)), linewidth = 1.25) +
  theme_bw()

llplot2 <- ggplot() + 
  geom_histogram(aes(x = log_liks2, y = after_stat(density)),
                 breaks = seq(-2, 2, by = .1), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = 5 * -0.0107919, sd = sqrt(0.222535)), linewidth = 1.25, linetype = 5) +
  stat_function(fun = dnorm, args = list(mean = mean(log_liks2), sd = sd(log_liks2)), linewidth = 1.25) +
  theme_bw()

pdf("hrplotsincreasedbias.pdf")
lhatplot2
llplot2
dev.off()
beep()

lhatplot2
llplot2

mean(log_liks2)
var(log_liks2)

(mean(l_hats2) - 1.25) * sqrt(987)
var(l_hats2) * 987

# pdf("hrqqplots2.pdf")
# ggplot() +
#   stat_qq(aes(sample = l_hats2)) +
#   stat_qq_line(aes(sample = l_hats2)) +
#   labs(x = "",y = "", title = TeX(r"($\hat{\lambda}$ QQ Plot)")) +
#   theme_bw()
# 
# ggplot() +
#   stat_qq(aes(sample = log_liks2)) +
#   stat_qq_line(aes(sample = log_liks2)) +
#   labs(x = "",y = "", title = TeX(r"(Score QQ Plot)")) +
#   theme_bw()
# dev.off()

log_liks3 <- rep(0,10000)
l_hats3 <- rep(0,10000)
numdiff <-rep(0,10000)
pb3 <- progress_bar$new(total = 10000)
for(j in 1:10000){
  pb3$tick()
  samps <- matrix(data = rep(0,2*269),ncol= 2) 
  for(k in 1:269){
    XY_scaled <- apply(mvrnorm(37138, mu = c(-bm37138^2,-bm37138^2), Sigma = bm37138^2*matrix(c(1, p2(2,bm37138), p2(2,bm37138), 1), ncol=2))
                       ,2,max)
    samps[k,] <- XY_scaled
  }
  log_liks3[j] <- grad(2, samps[,1], samps[,2])/sqrt(269)
  ltemp <- optimize(log_lik, interval = c(0,10), X = samps[,1], Y = samps[,2], maximum = T)$maximum
  l_hats3[j] <- ltemp
  for(i in 1:269){
      if((samps[i,1]+log(.5))*(samps[i,2]+log(.5)) < 0){
        numdiff[j] = numdiff[j] + abs(samps[i,1] + samps[i,2])
      }
  }
}
mean(log_liks3)
median(log_liks3)
var(log_liks3)

(mean(l_hats3) - 2) * sqrt(269)
(median(l_hats3) - 2) * sqrt(269)
var(l_hats3) * 269



# 
lhatplot3 <- ggplot() +
  geom_histogram(aes(x = l_hats3, y = after_stat(density)),
                 breaks = seq(0, 10, by = .1),
                 colour = "black",
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(l_hats3), sd = sd(l_hats3)), linewidth= 1.25) +
  stat_function(fun = dnorm, args = list(mean = 2  +(0.009/0.284)/sqrt(269), sd = sqrt(1/(0.284*269))), linewidth= 1.25, linetype = 5) +
  labs(x= TeX(r"($\hat{\lambda}$)"), title = TeX(r"($m = 37138$, $k = 269$, $\rho=0.608$)")) +
  theme_bw()

llplot3 <- ggplot() + 
  geom_histogram(aes(x = log_liks3, y = after_stat(density)),
                 breaks = seq(-2, 2, by = .1), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = -0.009, sd = sqrt(0.284)), linewidth = 1.25, linetype = 5, xlim = c(-2,2)) +
  stat_function(fun = dnorm, args = list(mean = mean(log_liks3), sd = sd(log_liks3)), linewidth = 1.25, xlim = c(-2,2)) +
  theme_bw()

# ggplot() +
#   stat_qq(aes(sample = l_hats3)) +
#   stat_qq_line(aes(sample = l_hats3)) +
#   theme_bw()
# 
# ggplot() +
#   stat_qq(aes(sample = log_liks3)) +
#   stat_qq_line(aes(sample = log_liks3)) +
#   theme_bw()

mean(log_liks3)
var(log_liks3)

pdf("hrplotslowrho.pdf")
lhatplot3
llplot3
dev.off()

save.image(file = "HRSimData.RData")
load("HRSimData.RData")

log_liks4 <- rep(0,10000)
l_hats4 <- rep(0,10000)
for(j in 1:10000){
  samps <- matrix(data = rep(0,2*95),ncol= 2) 
  for(k in 1:95){
    XY_scaled <- apply(mvrnorm(1044, mu = c(-bm1044^2,-bm1044^2), Sigma = bm1044^2*matrix(c(1, p2(1.25,bm1044), p2(1.25,bm1044), 1), ncol=2))
                       ,2,max)
    samps[k,] <- XY_scaled
  }
  log_liks4[j] <- grad(1.25, samps[,1], samps[,2])/sqrt(95)
  ltemp <- optimize(log_lik, interval = c(0,10), X = samps[,1], Y = samps[,2], maximum = T)$maximum
  l_hats4[j] <- ltemp
}
mean(log_liks4)
median(log_liks4)
var(log_liks4)

(mean(l_hats4) - 2) * sqrt(95)
(median(l_hats4) - 2) * sqrt(95)
var(l_hats4-2) * 95


lhatplot4 <- ggplot() + 
  geom_histogram(aes(x = l_hats4, y = after_stat(density)),
                 breaks = seq(.5, 2, by = .025), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = 1.25 + (0.0107919/0.222535)/sqrt(95), sd = sqrt(1/(0.222535*95))), linetype = 5, linewidth = 1.25) +
  stat_function(fun = dnorm, args = list(mean = mean(l_hats4), sd = sd(l_hats4)), linewidth = 1.25) +
  theme_bw()

llplot4 <- ggplot() + 
  geom_histogram(aes(x = log_liks4, y = after_stat(density)),
                 breaks = seq(-2, 2, by = .1), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = -0.0107919, sd = sqrt(0.222535)), linewidth = 1.25, linetype = 5) +
  stat_function(fun = dnorm, args = list(mean = mean(log_liks4), sd = sd(log_liks4)), linewidth = 1.25) +
  theme_bw()


################################

l_hats5 <- rep(0,1000)
for(j in 1:1000){
  samps <- rbvhr(100,1/5)
  l_hats5[j] <- optimize(log_lik, interval = c(0,10), X = samps[,1], Y = samps[,2], tol=.Machine$double.eps^.85)$minimum
}
hist(l_hats5)
beep()


l_hats5 <- rep(0,100)
for(j in 1:100){
  samps <- rbvhr(100,1/5)
  l_hats5[j] <- optim(par=7, fn = log_lik, X = samps[,1], Y = samps[,2], upper= 10, lower = 0, method = "L-BFGS-B")$par
}
hist(l_hats5)

l_hats5 <- rep(0,100)
for(j in 1:100){
  samps <- rbvhr(100,1/5)
  l_hats5[j] <- optim(par=7, fn = log_lik, X = samps[,1], Y = samps[,2], upper= 10, lower = 0, method = "Brent")$convergence
}



condDist <- function(y,q,l,x){
  abs(exp(exp(-x)*(1-pnorm(l+(y-x)/(2*l))) - exp(-y) * pnorm(l +(x-y)/(2*l)))*pnorm(l + (y-x)/(2*l)) - q)
}

rhr <- function(n,l){
  out <- matrix(data = runif(2*n), nrow = n, ncol = 2)
  out[,1] <- -log(-log(out[,1]))
  for(i in 1:n){
    out[i,2] <- optim(par = 1, fn = condDist, q = out[i,2], l = l, x = out[i,1], lower = -10, upper = 10, method = "Brent")$par
  }
  return(out)
}

trysamp <- rhr(1000, 1/5)
saveVals <- seq(.1,10, by = .01)
for(i in 1:length(saveVals)){
  saveVals[i] <- log_lik(saveVals[i],trysamp[,1],trysamp[,2])
}
seq(.1,10, by = .01)[which.min(saveVals)]
optim(par = 7, fn = log_lik, lower = .01, upper = 10, X = trysamp[,1], Y = trysamp[,2], method = "Brent")$par


simFI <- function(l, n){
  samp <- rbvhr(n,1/l)
  sum <- 0
  for(i in 1:n){
    x <- samp[i,1]
    y <- samp[i,2]
    t1 <- 2 * exp(-x) * dnorm(l + (y-x)/(2*l)) * (l - (y-x)^2/(4*l^3))
    denom <- dbvhr(c(x,y), 1/l) / pbvhr(c(x,y),1/l)
    t2 <- exp(-x) * dnorm(l + (y-x)/(2*l)) * (2 * exp(-x) * dnorm(l + (y-x)/(2*l))*(1 - (x-y)^2/(2*l^2)) + 
                                                exp(-x)*pnorm(l + (y-x)/(2*l))*(-l + (x-y)/(2*l) + (y-x)^2/(4*l^3) - (x-y)^3/(8*l^5) + (x-y) / (l^3)) +
                                                exp(-y)*pnorm(l + (x-y)/(2*l))*(-l + (y-x)/(2*l) + (x-y)^2/(4*l^3) - (y-x)^3/(8*l^5) + (y-x) / (l^3)))
    t3 <- exp(-x) * dnorm(l + (y-x)/(2*l)) * (-l/2 -1/(2*l) - 1/(l^3) + (y-x)^2/(4*l^3) + 5*(y-x)^2/(32*l^5) - (y-x)^4/(32 * l^7))
    t4 <- exp(-x) * dnorm(l + (y-x)/(2*l)) * (exp(-x) * pnorm(l + (y-x)/(2*l)) * (1-(x-y)/(2*l^2)) +
                                                exp(-y) * pnorm(l + (x-y)/(2*l)) * (1-(y-x)/(2*l^2)) +
                                                -.5 - 1/(2*l^2) + (y-x)^2/(8*l^4))
    sum <- sum + t1 + t2/denom - t3/denom - t4^2/denom^2
  }
  return(sum/n)
}
simFI(.5,1000000)