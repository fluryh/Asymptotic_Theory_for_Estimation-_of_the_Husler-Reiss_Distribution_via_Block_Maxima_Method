library(evd)
library(MASS)
library(ggplot2)
library(progress)
library(beepr)
library(latex2exp)
set.seed(2293)

bm5938 <- 3.60239
bm1044 <- 3.12796
bm1128 <- 3.15023
bm334 <- 2.78177
bm4994 <- 3.55755
bm5471 <- 3.58125
bm37138 <- 4.05082
bm34092 <- 4.0308565224475466
bm31300 <- 4.010842635336875
bm17732 <- 3.87545

p2 <- function(l,bm){
  return((1-l^2/bm^2)/(1+l^2/bm^2))
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
  return(sum)
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

##### Simulations with l = .5, with differing sample sizes
log_liks <- rep(0,10000)
l_hats <- rep(0,10000)
pb <- progress_bar$new(total = 10000)
numdiff <- rep(0,10000)
for(j in 1:10000){
  pb$tick()
  samps <- matrix(data = rep(0,2*95),ncol= 2) 
  for(k in 1:95){
    XY_scaled <- apply(mvrnorm(1044, mu = c(-bm1044^2,-bm1044^2), Sigma = bm1044^2*matrix(c(1, p2(.5,bm1044), p2(.5,bm1044), 1), ncol=2))
                       ,2,max)
    samps[k,] <- XY_scaled
  }
  log_liks[j] <- grad(.5, samps[,1], samps[,2])/sqrt(95)
  l_hats[j] <- optimize(log_lik, interval = c(0,5), X = samps[,1], Y = samps[,2], maximum = T)$maximum
}
mean(log_liks)
var(log_liks)

(mean(l_hats) - .5) * sqrt(95)
(mean(l_hats) - .5) * sqrt(95) - 1.96 * sd(l_hats) * sqrt(95) / 100
(mean(l_hats) - .5) * sqrt(95) + 1.96 * sd(l_hats) * sqrt(95) / 100
var(l_hats) * 95

lhatplot <- ggplot() + 
  geom_histogram(aes(x = l_hats, y = after_stat(density)),
                 breaks = seq(.1, 1, by = .01), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(l_hats), sd = sd(l_hats)), linewidth= 1.25) +
  stat_function(fun = dnorm, args = list(mean = (0.0116526)/sqrt(95)+.5, sd = sqrt(1/(4.5533*95))), linewidth= 1.25, linetype = 5) +
  labs(x= TeX(r"($\hat{\lambda}$)"), title = TeX(r"($m = 1044$, $k = 95$, $\rho=0.95$)")) + 
  xlim(c(.3,.7))+
  theme_bw()
lhatplot

loglikeplot <- ggplot() + 
  geom_histogram(aes(x = log_liks, y = after_stat(density)),
                 breaks = seq(-10, 10, by = .25), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(log_liks), sd = sd(log_liks)), linewidth= 1.25) +
  stat_function(fun = dnorm, args = list(mean = 0.0530579, sd = sqrt(4.5533)), linewidth= 1.25, linetype = 5) +
  labs(x= "Score Function", title = TeX(r"($m = 1044$, $k = 95$, $\rho=0.95$)")) + 
  xlim(c(-8,10))+
  theme_bw()
loglikeplot

# pdf("Sim1Plots.pdf")
ggsave("Sim1lhat.jpg", plot = lhatplot, dpi = "retina")
# lhatplot
# dev.off()
# ggsave("Sim1lhat.jpg")
ggsave("Sim1ll.jpg", plot = loglikeplot, dpi = "retina")
loglikeplot
dev.off()

ggsave()

log_liks2 <- rep(0,10000)
l_hats2 <- rep(0,10000)
pb <- progress_bar$new(total = 10000)
numdiff <- rep(0,10000)
for(j in 1:10000){
  pb$tick()
  samps <- matrix(data = rep(0,2*168),ncol= 2) 
  for(k in 1:168){
    XY_scaled <- apply(mvrnorm(5938, mu = c(-bm5938^2,-bm5938^2), Sigma = bm5938^2*matrix(c(1, p2(.5,bm5938), p2(.5,bm5938), 1), ncol=2))
                       ,2,max)
    samps[k,] <- XY_scaled
  }
  log_liks2[j] <- grad(.5, samps[,1], samps[,2])/sqrt(168)
  l_hats2[j] <- optimize(log_lik, interval = c(0,5), X = samps[,1], Y = samps[,2], maximum = T)$maximum
}
mean(log_liks2)
var(log_liks2)

(mean(l_hats2) - .5) * sqrt(168)
(mean(l_hats2) - .5) * sqrt(168) - 1.96 * sd(l_hats2) * sqrt(168) / 100
(mean(l_hats2) - .5) * sqrt(168) + 1.96 * sd(l_hats2) * sqrt(168) / 100

var(l_hats2) * 168

lhatplot2 <- ggplot() + 
  geom_histogram(aes(x = l_hats2, y = after_stat(density)),
                 breaks = seq(.25, .75, by = .01), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(l_hats2), sd = sd(l_hats2)), linewidth= 1.25) +
  stat_function(fun = dnorm, args = list(mean = (0.0116526)/sqrt(168)+.5, sd = sqrt(1/(4.5533*168))), linewidth= 1.25, linetype = 5) +
  labs(x= TeX(r"($\hat{\lambda}$)"), title = TeX(r"($m = 4994$, $k = 168$, $\rho=0.96$)")) + 
  xlim(c(.35,.65))+
  theme_bw()
lhatplot2

loglikeplot2 <- ggplot() + 
  geom_histogram(aes(x = log_liks2, y = after_stat(density)),
                 breaks = seq(-7.5, 7.5, by = .25), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(log_liks2), sd = sd(log_liks2)), linewidth= 1.25) +
  stat_function(fun = dnorm, args = list(mean = 0.0530579, sd = sqrt(4.5533)), linewidth= 1.25, linetype = 5) +
  labs(x= "Score Function", title = TeX(r"($m = 4994$, $k = 168$, $\rho=0.96$)")) + 
  xlim(c(-7.5,7.5))+
  theme_bw()
loglikeplot2
beep()

# pdf("Sim2Plots.pdf")
ggsave("Sim2lhat.jpg", plot = lhatplot2, dpi = "retina")
ggsave("Sim2ll.jpg", plot = loglikeplot2, dpi = "retina")

jpeg("Sim2lhat.jpg")
lhatplot2
dev.off()
jpeg("Sim2ll.jpg")
loglikeplot2
dev.off()

log_liks3 <- rep(0,10000)
l_hats3 <- rep(0,10000)
pb <- progress_bar$new(total = 10000)
numdiff <- rep(0,10000)
for(j in 1:10000){
  pb$tick()
  samps <- matrix(data = rep(0,2*269),ncol= 2) 
  for(k in 1:269){
    XY_scaled <- apply(mvrnorm(37138, mu = c(-bm37138^2,-bm37138^2), Sigma = bm37138^2*matrix(c(1, p2(.5,bm37138), p2(.5,bm5938), 1), ncol=2))
                       ,2,max)
    samps[k,] <- XY_scaled
  }
  log_liks3[j] <- grad(.5, samps[,1], samps[,2])/sqrt(269)
  l_hats3[j] <- optimize(log_lik, interval = c(0,5), X = samps[,1], Y = samps[,2], maximum = T)$maximum
}
mean(log_liks3)
var(log_liks3)

(mean(l_hats3) - .5) * sqrt(269)
(mean(l_hats3) - .5) * sqrt(269) - 1.96 * sd(l_hats3) * sqrt(269) / 100
(mean(l_hats3) - .5) * sqrt(269) + 1.96 * sd(l_hats3) * sqrt(269) / 100
var(l_hats3) * 269

lhatplot3 <- ggplot() + 
  geom_histogram(aes(x = l_hats3, y = after_stat(density)),
                 breaks = seq(.35, .65, by = .01), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(l_hats3), sd = sd(l_hats3)), linewidth= 1.25) +
  stat_function(fun = dnorm, args = list(mean = (0.0116526)/sqrt(269)+.5, sd = sqrt(1/(4.5533*269))), linewidth= 1.25, linetype = 5) +
  labs(x= TeX(r"($\hat{\lambda}$)"), title = TeX(r"($m = 4994$, $k = 269$, $\rho=0.97$)")) +
  xlim(c(.35,.65)) + 
  theme_bw()
lhatplot3

loglikeplot3 <- ggplot() + 
  geom_histogram(aes(x = log_liks3, y = after_stat(density)),
                 breaks = seq(-7.5, 7.5, by = .25), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(log_liks3), sd = sd(log_liks3)), linewidth= 1.25) +
  stat_function(fun = dnorm, args = list(mean = 0.0530579, sd = sqrt(4.5533)), linewidth= 1.25, linetype = 5) +
  labs(x= "Score Function", title = TeX(r"($m = 4994$, $k = 168$, $\rho=0.97$)")) + 
  xlim(c(-7.5,7.5)) +
  theme_bw()
loglikeplot3

# pdf("Sim3Plots.pdf")
jpeg("Sim3lhat.jpg")
lhatplot3
dev.off()
jpeg("Sim3ll.jpg")
loglikeplot3
dev.off()

ggsave("Sim3lhat.jpg", plot = lhatplot3, dpi = "retina")
ggsave("Sim3ll.jpg", plot = loglikeplot3, dpi = "retina")


save.image("HRCode.R")
load("HRCode.R")

log_liks4 <- rep(0,10000)
l_hats4 <- rep(0,10000)
pb <- progress_bar$new(total = 10000)
numdiff <- rep(0,10000)
for(j in 1:10000){
  pb$tick()
  samps <- matrix(data = rep(0,2*56),ncol= 2) 
  for(k in 1:56){
    XY_scaled <- apply(mvrnorm(17732, mu = c(-bm17732^2,-bm17732^2), Sigma = bm17732^2*matrix(c(1, p2(.5,bm17732), p2(.5,bm17732), 1), ncol=2))
                       ,2,max)
    samps[k,] <- XY_scaled
  }
  log_liks4[j] <- grad(.5, samps[,1], samps[,2])/sqrt(56)
  l_hats4[j] <- optimize(log_lik, interval = c(0,5), X = samps[,1], Y = samps[,2], maximum = T)$maximum
}
mean(log_liks4)
var(log_liks4)

(mean(l_hats4) - .5) * sqrt(56)
(mean(l_hats4) - .5) * sqrt(56) - 1.96 * sd(l_hats4) * sqrt(56) / 100
(mean(l_hats4) - .5) * sqrt(56) + 1.96 * sd(l_hats4) * sqrt(56) / 100
var(l_hats4) * 56

lhatplot4 <- ggplot() + 
  geom_histogram(aes(x = l_hats4, y = after_stat(density)),
                 breaks = seq(.25, .75, by = .01), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(l_hats4), sd = sd(l_hats4)), linewidth= 1.25) +
  stat_function(fun = dnorm, args = list(mean = (0.0116526/3)/sqrt(56)+.5, sd = sqrt(1/(4.5533*56))), linewidth= 1.25, linetype = 5) +
  labs(x= TeX(r"($\hat{\lambda}$)"), title = TeX(r"($m = 4994$, $k = 56$, $\rho=0.95$)")) + 
  theme_bw()
lhatplot4

loglikeplot4 <- ggplot() + 
  geom_histogram(aes(x = log_liks4, y = after_stat(density)),
                 breaks = seq(-10, 10, by = .25), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(log_liks4), sd = sd(log_liks4)), linewidth= 1.25) +
  stat_function(fun = dnorm, args = list(mean = 0.0530579/3, sd = sqrt(4.5533)), linewidth= 1.25, linetype = 5) +
  labs(x= TeX(r"($\hat{\lambda}$)"), title = TeX(r"($m = 4994$, $k = 56$, $\rho=0.95$)")) + 
  theme_bw()
loglikeplot4

pdf("Sim4Plots.pdf")
lhatplot4
loglikeplot4
dev.off()


log_liks5 <- rep(0,10000)
l_hats5 <- rep(0,10000)
pb <- progress_bar$new(total = 10000)
numdiff <- rep(0,10000)
for(j in 1:10000){
  pb$tick()
  samps <- matrix(data = rep(0,2*3),ncol= 2) 
  for(k in 1:3){
    XY_scaled <- apply(mvrnorm(31300, mu = c(-bm31300^2,-bm31300^2), Sigma = bm31300^2*matrix(c(1, p2(.5,bm31300), p2(.5,bm31300), 1), ncol=2))
                       ,2,max)
    samps[k,] <- XY_scaled
  }
  log_liks5[j] <- grad(.5, samps[,1], samps[,2])/sqrt(3)
  l_hats5[j] <- optimize(log_lik, interval = c(0,5), X = samps[,1], Y = samps[,2], maximum = T)$maximum
}
mean(log_liks5)
var(log_liks5)

(mean(l_hats5) - .5) * sqrt(3)
var(l_hats5) * 3

lhatplot5 <- ggplot() + 
  geom_histogram(aes(x = l_hats5, y = after_stat(density)),
                 breaks = seq(0, 5, by = .1), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(l_hats5), sd = sd(l_hats5)), linewidth= 1.25) +
  stat_function(fun = dnorm, args = list(mean = (0.0116526/9)/sqrt(3)+.5, sd = sqrt(1/(4.5533*3))), linewidth= 1.25, linetype = 5) +
  labs(x= TeX(r"($\hat{\lambda}$)"), title = TeX(r"($m = 4994$, $k = 3$, $\rho=0.95$)")) + 
  theme_bw()
lhatplot5

loglikeplot5 <- ggplot() + 
  geom_histogram(aes(x = log_liks5, y = after_stat(density)),
                 breaks = seq(-10, 10, by = .25), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(log_liks5), sd = sd(log_liks5)), linewidth= 1.25) +
  stat_function(fun = dnorm, args = list(mean = 0.0530579/9, sd = sqrt(4.5533)), linewidth= 1.25, linetype = 5) +
  labs(x= TeX(r"($\hat{\lambda}$)"), title = TeX(r"($m = 4994$, $k = 3$, $\rho=0.95$)")) + 
  theme_bw()
loglikeplot5











log_liks <- rep(0,10000)
l_hats <- rep(0,10000)
pb <- progress_bar$new(total = 10000)
numdiff <- rep(0,10000)
for(j in 1:10000){
  pb$tick()
  samps <- matrix(data = rep(0,2*18),ncol= 2) 
  for(k in 1:18){
    XY_scaled <- apply(mvrnorm(1044, mu = c(-bm5471^2,-bm5471^2), Sigma = bm5471^2*matrix(c(1, p2(.5,bm5471), p2(.5,bm5471), 1), ncol=2))
                       ,2,max)
    samps[k,] <- XY_scaled
  }
  log_liks[j] <- grad(.5, samps[,1], samps[,2])/sqrt(18)
  l_hats[j] <- optimize(log_lik, interval = c(0,18), X = samps[,1], Y = samps[,2], maximum = T)$maximum
}
mean(log_liks)
beep()
var(log_liks)

(mean(l_hats) - .5) * sqrt(18)
var(l_hats) * 18

lhatplot <- ggplot() + 
  geom_histogram(aes(x = l_hats, y = after_stat(density)),
                 breaks = seq(.1, 1, by = .01), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(l_hats), sd = sd(l_hats)), linewidth= 1.25) +
  stat_function(fun = dnorm, args = list(mean = (0.0530579)/sqrt(18)+.5, sd = sqrt(1/(4.5533*18))), linewidth= 1.25, linetype = 5) +
  labs(x= TeX(r"($\hat{\lambda}$)"), title = TeX(r"($m = 1044$, $k = 18$, $\rho=0.95$)")) + 
  xlim(c(0,1))+
  theme_bw()
lhatplot

loglikeplot <- ggplot() + 
  geom_histogram(aes(x = log_liks, y = after_stat(density)),
                 breaks = seq(-10, 10, by = .25), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(log_liks), sd = sd(log_liks)), linewidth= 1.25) +
  stat_function(fun = dnorm, args = list(mean = 0.0530579, sd = sqrt(4.5533)), linewidth= 1.25, linetype = 5) +
  labs(x= "Score Function", title = TeX(r"($m = 1044$, $k = 18$, $\rho=0.95$)")) + 
  theme_bw()
loglikeplot




log_liks3 <- rep(0,100000)
l_hats3 <- rep(0,100000)
pb <- progress_bar$new(total = 100000)
numdiff <- rep(0,100000)
for(j in 1:100000){
  pb$tick()
  samps <- matrix(data = rep(0,2*886),ncol= 2) 
  for(k in 1:886){
    XY_scaled <- apply(mvrnorm(1128, mu = c(-bm1128^2,-bm1128^2), Sigma = bm1128^2*matrix(c(1, p2(.5,bm1128), p2(.5,bm1128), 1), ncol=2))
                       ,2,max)
    samps[k,] <- XY_scaled
  }
  log_liks3[j] <- grad(.5, samps[,1], samps[,2])/sqrt(886)
  l_hats3[j] <- optimize(log_lik, interval = c(0,10), X = samps[,1], Y = samps[,2], maximum = T)$maximum
}
mean(log_liks3)
var(log_liks3)

(mean(l_hats3) - .5) * sqrt(886)
var(l_hats3) * 886

lhatplot3 <- ggplot() + 
  geom_histogram(aes(x = l_hats3, y = after_stat(density)),
                 breaks = seq(.25, .75, by = .025), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(l_hats3), sd = sd(l_hats3)), linewidth= 1.25) +
  stat_function(fun = dnorm, args = list(mean = (sqrt(5)*0.00629098)/sqrt(886)+.5, sd = sqrt(1/(4.5533*886))), linewidth= 1.25, linetype = 5) +
  labs(x= TeX(r"($\hat{\lambda}$)"), title = TeX(r"($m = 1128$, $k = 886$, $\rho=0.937$)")) + 
  theme_bw()
lhatplot3


lhatplot3 <- ggplot() + 
  geom_histogram(aes(x = (l_hats3 - .5) * sqrt(886), y = after_stat(density)),
                 breaks = seq(-2, 2, by = .05), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = (mean(l_hats3) - .5) * sqrt(886), sd = sd((l_hats3 - .5) * sqrt(886))), linewidth= 1.25) +
  stat_function(fun = dnorm, args = list(mean = sqrt(5)*0.0286447, sd = sqrt(1/4.5533)), linewidth= 1.25, linetype = 5) +
  labs(x= TeX(r"($\hat{\lambda}$)"), title = TeX(r"($m = 1128$, $k = 886$, $\rho=0.937$)")) + 
  theme_bw()
lhatplot3
beep()


log_liks4 <- rep(0,10000)
l_hats4 <- rep(0,10000)
pb <- progress_bar$new(total = 10000)
numdiff <- rep(0,10000)
for(j in 1:10000){
  pb$tick()
  samps <- matrix(data = rep(0,2*168),ncol= 2) 
  for(k in 1:168){
    XY_scaled <- apply(mvrnorm(5938, mu = c(-bm5938^2,-bm5938^2), Sigma = bm5938^2*matrix(c(1, p2(.25,bm5938), p2(.25,bm5938), 1), ncol=2))
                       ,2,max)
    samps[k,] <- XY_scaled
  }
  log_liks4[j] <- grad(.25, samps[,1], samps[,2])/sqrt(168)
  l_hats4[j] <- optimize(log_lik, interval = c(0,10), X = samps[,1], Y = samps[,2], maximum = T)$maximum
}
mean(log_liks4)
var(log_liks4)

(mean(l_hats4) - .25) * sqrt(168)
var(l_hats4) * 168

lhatplot4 <- ggplot() + 
  geom_histogram(aes(x = l_hats4, y = after_stat(density)),
                 breaks = seq(.1, .4, by = .01), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(l_hats4), sd = sd(l_hats4)), linewidth= 1.25) +
  stat_function(fun = dnorm, args = list(mean = (0.00389052)/sqrt(168)+.25, sd = sqrt(1/(21.9023*168))), linewidth= 1.25, linetype = 5) +
  labs(x= TeX(r"($\hat{\lambda}$)"), title = TeX(r"($m = 37138$, $k = 269$, $\rho=0.826$)")) + 
  theme_bw()
lhatplot4

loglikeplot4 <- ggplot() + 
  geom_histogram(aes(x = log_liks4, y = after_stat(density)),
                 breaks = seq(-20, 20, by = .5), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(log_liks4), sd = sd(log_liks4)), linewidth= 1.25) +
  stat_function(fun = dnorm, args = list(mean = 0.0852113, sd = sqrt(21.9023)), linewidth= 1.25, linetype = 5) +
  labs(x= TeX(r"($\hat{\lambda}$)"), title = TeX(r"($m = 37138$, $k = 269$, $\rho=0.826$)")) + 
  theme_bw()
loglikeplot4
beep()

log_liks5 <- rep(0,10000)
l_hats5 <- rep(0,10000)
pb <- progress_bar$new(total = 10000)
numdiff <- rep(0,10000)
for(j in 1:10000){
  pb$tick()
  samps <- matrix(data = rep(0,2*168),ncol= 2) 
  for(k in 1:168){
    XY_scaled <- apply(mvrnorm(5938, mu = c(-bm5938^2,-bm5938^2), Sigma = bm5938^2*matrix(c(1, p2(.75,bm5938), p2(.75,bm5938), 1), ncol=2))
                       ,2,max)
    samps[k,] <- XY_scaled
  }
  log_liks5[j] <- grad(.75, samps[,1], samps[,2])/sqrt(168)
  l_hats5[j] <- optimize(log_lik, interval = c(0,10), X = samps[,1], Y = samps[,2], maximum = T)$maximum
}
mean(log_liks4)
var(log_liks4)

(mean(l_hats5) - .75) * sqrt(168)
var(l_hats5) * 168

lhatplot5 <- ggplot() + 
  geom_histogram(aes(x = l_hats4, y = after_stat(density)),
                 breaks = seq(.5, 1, by = .01), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(l_hats4), sd = sd(l_hats4)), linewidth= 1.25) +
  stat_function(fun = dnorm, args = list(mean = (0.0630814)/sqrt(168)+.75, sd = sqrt(1/(1.79065*168))), linewidth= 1.25, linetype = 5) +
  labs(x= TeX(r"($\hat{\lambda}$)"), title = TeX(r"($m = 37138$, $k = 269$, $\rho=0.826$)")) + 
  theme_bw()
lhatplot5

loglikeplot4 <- ggplot() + 
  geom_histogram(aes(x = log_liks4, y = after_stat(density)),
                 breaks = seq(-5, 5, by = .5), 
                 colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(log_liks4), sd = sd(log_liks4)), linewidth= 1.25) +
  stat_function(fun = dnorm, args = list(mean = 0.0630814, sd = sqrt(1.79065)), linewidth= 1.25, linetype = 5) +
  labs(x= TeX(r"($\hat{\lambda}$)"), title = TeX(r"($m = 37138$, $k = 269$, $\rho=0.826$)")) + 
  theme_bw()
loglikeplot4
beep()