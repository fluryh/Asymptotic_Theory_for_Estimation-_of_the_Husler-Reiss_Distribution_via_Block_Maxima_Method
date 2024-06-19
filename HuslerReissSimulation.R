library(MASS)
library(ggplot2)
library(evd)
library(latex2exp)
library(beepr)
set.seed(5838)
bm608 <- 2.96765
bm24480 <- 3.95278
bm3703 <- 3.47887
bm172117 <- 4.39468
bm249509 <- 4.47436
bm208 <- 2.6296
p <- function(l,m){
  return(1-l^2/log(m))
}
p3 <- function(l,bm){
  return(.9 + .1*p2(l/sqrt(.1),bm))
}
p2 <- function(l,bm){
  return((1-l^2/bm^2)/(1+l^2/bm^2))
}
d2 <- function(x,y,l){
  t1 <- 2*exp(-x)*dnorm(l + (y-x)/(2*l))*(l- (y-x)^2/(4*l^3))
  denom <- exp(-x-y)*pnorm(l + (x-y)/(2*l))*pnorm(l + (y-x)/(2*l)) + 
    exp(-x)*pnorm(l + (y-x)/(2*l))/(2*l)
  t2 <- exp(-x)*dnorm(l + (y-x)/(2*l)) * (2*exp(-x)*dnorm(l + (y-x)/(2*l))*(1 - (x-y)^2/(4*l^4)) +
                                            exp(-x)*pnorm(l + (y-x)/(2*l))*(-l+(x-y)/(2*l) + (y-x)^2/(4*l^3) - (x-y)^3/(8*l^5) + (x-y)/(l^3)) +
                                            exp(-y)*pnorm(l + (x-y)/(2*l))*(-l+(y-x)/(2*l) + (x-y)^2/(4*l^3) - (y-x)^3/(8*l^5) + (y-x)/(l^3)))
  t2 <- t2/denom
  t3 <- exp(-x)*dnorm(l + (y-x)/(2*l)) *(-l/2 - 1/(2*l)- 1/(l^3) + (y-x)^2/(4*l^3) + .625 * (y-x)^2/(l^5) - (y-x)^4/(32*l^7))
  t3 <- t3/denom
  t4 <- exp(-x)*dnorm(l + (y-x)/(2*l)) *(exp(-x)*pnorm(l + (y-x)/(2*l))*(1-(x-y)/(2*l^2)) +
                                           exp(-y)*pnorm(l + (x-y)/(2*l))*(1-(y-x)/(2*l^2)) -
                                           1/(2*l^2) - .5 + (x-y)^2/(8*l^4))
  t4 <- (t4/denom)^2
  return(c(t1,t2,t3,t4))
}
log_lik <- function(l,X,Y){
  sum <- 0
  d <- length(X)
  for(i in 1:d){
    x <- X[i]
    y <- Y[i]
    sum <- sum - exp(-x)*pnorm(l + (y-x)/(2*l)) - exp(-y)*pnorm(l + (x-y)/(2*l)) + 
      log(exp(-x-y)*pnorm(l + (y-x)/(2*l))*pnorm(l + (x-y)/(2*l)) +
            exp(-x)*dnorm(l + (y-x)/(2*l))/(2*l))
  }
  return(sum)
}
grad <- function(l,X,Y){
  sum <- 0
  d <- length(X)
  for(i in 1:d){
    x <- X[i]
    y <- Y[i]
    sum <- sum -2*exp(-x)*dnorm(l + (y-x)/(2*l)) + exp(-x)*dnorm(l + (x-y)/(2*l))*(1/2 + 1/(2*l^2)+(y-x)^2/(8*l^4) + pnorm(l + (x-y)/(2*l))*exp(-y)*(1-(y-x)/(2*l^2)) + pnorm(l + (y-x)/(2*l))*exp(-x)*(1-(x-y)/(2*l)))/(exp(-x-y)*pnorm(l + (x-y)/(2*l))*pnorm(l + (y-x)/(2*l)) + exp(-x)*dnorm(l + (y-x)/(2*l))/(2*l))
  }
  return(sum)
}

#Here we use 608 and 164. Real values are 164.388 and 608.3169
#Or 24,480 and 408
l_hats <- rep(0,1000)
for(j in 1:1000){
  XY_scaled <- data.frame(matrix(nrow = 164, ncol = 2))
  for(k in 1:164){
    XY_scaled[k,] <- apply(mvrnorm(608, mu = c(-bm608^2,-bm608^2), Sigma = bm608^2*matrix(c(1, p(2,608), p(2,608), 1), ncol=2))
                ,2,max)
  }
  l_hats[j] <- optimize(log_lik, interval = c(0, 10), X = XY_scaled[,1], Y = XY_scaled[,2],
           maximum = T, tol = .Machine$double.eps)$maximum
  if(l_hats[j] == 10){
    print("flag")
    k <- 2
    while(l_hats[j] == k*10 && k < 10){
      l_hats[j] <- optimize(log_lik, interval = c(0, k*10), X = XY_scaled[,1], Y = XY_scaled[,2],
                            maximum = T, tol = .Machine$double.eps)$maximum
    }
  }
}
l_hats_clean <- l_hats[sqrt(165)*(2-l_hats) > -20]
ggplot() + 
  geom_density(aes(x = sqrt(165)*(2-l_hats)),fill = 'lightgray', col = 'black')+
  stat_function(fun = dnorm, args = list(mean = mean(sqrt(165)*(2-l_hats)), sd = sqrt(1/.283969)), color ="red") +
  theme_bw()
ggplot() + 
  geom_density(aes(x = sqrt(165)*(2-l_hats_clean)),fill = 'lightgray', col = 'black')+
  stat_function(fun = dnorm, args = list(mean = mean(sqrt(165)*(2-l_hats_clean)), sd = sd(sqrt(165)*(2-l_hats_clean)))) +
  stat_function(fun = dnorm, args = list(mean = mean(sqrt(165)*(2-l_hats_clean)), sd = sqrt(1/.283969)), color ="red") +
  theme_bw()
var(sqrt(165)*(2-l_hats))
var(sqrt(165)*(2-l_hats_clean))

XY_scaled <- data.frame(matrix(nrow = 164, ncol = 2))
for(k in 1:164){
  XY <- apply(mvrnorm(608, mu = c(-bm608^2,-bm608^2), Sigma = bm608^2*matrix(c(1, p(2,608), p(2,608), 1), ncol=2))
              ,2,max)
  XY_scaled[k,] <- XY
}
ggplot() + 
  geom_density(aes(x = XY_scaled[,1]),fill = 'lightgray', col = 'black')+
  stat_function(fun = dgumbel, args = list(loc = 0, scale = 1), color ="red") +
  # stat_function(fun = dgumbel, args = list(loc = 3.975, scale = .08843979), color ="purple") +
  theme_bw()
# guplot(XY_scaled[,1])


# l_hats3 <- rep(0,200)
# for(j in 1:200){
#   XY_scaled <- data.frame(matrix(nrow = 270, ncol = 2))
#   for(k in 1:270){
#     XY <- apply(mvrnorm(3703, mu = c(0,0), Sigma = matrix(c(1, p(2,3703), p(2,3703), 1), ncol=2))
#                 ,2,max)
#     XY_scaled[k,] <- XY/bm3703 + bm3703
#   }
#   l_hats3[j] <- optimize(log_lik, interval = c(0, 10), X = XY_scaled[,1], Y = XY_scaled[,2],
#                         maximum = T, tol = .Machine$double.eps)$maximum
#   if(l_hats3[j] == 10){
#     print("flag")
#     k <- 2
#     while(l_hats3[j] == k*10 && k < 10){
#       l_hats3[j] <- optimize(log_lik, interval = c(0, k*10), X = XY_scaled[,1], Y = XY_scaled[,2],
#                             maximum = T, tol = .Machine$double.eps)$maximum
#     }
#   }
# }
# ggplot() + 
#   geom_density(aes(x = sqrt(270)*(2-l_hats3)),fill = 'lightgray', col = 'black')+
#   stat_function(fun = dnorm, args = list(mean = mean(sqrt(270)*(2-l_hats3)), sd = sd(sqrt(270)*(2-l_hats3)))) +
#   stat_function(fun = dnorm, args = list(mean = mean(sqrt(270)*(2-l_hats3)), sd = sqrt(1/.283969)), color ="red") +
#   theme_bw()
# ggplot() + geom_density(aes(x = l_hats3)) + theme_bw()
# var(sqrt(270)*(2-l_hats3))

# l_hats2 <- rep(0,200)
# flags <- rep(0,200)
# for(j in 1:200){
#   XY_scaled <- data.frame(matrix(nrow = 408, ncol = 2))
#   for(k in 1:408){
#     XY_scaled[k,] <- apply(mvrnorm(24480, mu = c(-bm24480^2,-bm24480^2), Sigma = bm24480^2*matrix(c(1, p(2,24480), p(2,24480), 1), ncol=2))
#                            ,2,max)
#   }
#   for(i in 1:408){
#     if((XY_scaled[i,1] > -log(.5) & XY_scaled[i,2] < -log(.5)) || (XY_scaled[i,1] < -log(.5) & XY_scaled[i,2] > -log(.5))){
#       flags[j] <- flags[j] + abs(XY_scaled[i,1] - XY_scaled[i,2])
#     }
#   }
#   l_hats2[j] <- optimize(log_lik, interval = c(0, 10), X = XY_scaled[,1], Y = XY_scaled[,2],
#                         maximum = T, tol = .Machine$double.eps)$maximum
#   if(l_hats2[j] == 10){
#     print("flag")
#     k <- 2
#     while(l_hats2[j] == k*10 && k < 10){
#       l_hats2[j] <- optimize(log_lik, interval = c(0, k*10), X = XY_scaled[,1], Y = XY_scaled[,2],
#                             maximum = T, tol = .Machine$double.eps)$maximum
#     }
#   }
# }
# # pdf("./Explanation for Second Mode")
# ggplot() + 
#   geom_point(aes(x= flags, y = l_hats2)) +
#   geom_hline(yintercept = 2, color = "red")+ 
#   labs(x = "Sum of Difference of Points on Opposite Sides of -ln(.5)", y = "Lambda Hat", title = "Explanation for Second Mode") +
#   theme_bw()
# l_hats2_clean <- l_hats2[sqrt(408)*(2-l_hats2) > -20]
# ggplot() + 
#   geom_density(aes(x = sqrt(408)*(2-l_hats2)),fill = 'lightgray', col = 'black')+
#   stat_function(fun = dnorm, args = list(mean = mean(sqrt(408)*(2-l_hats2)), sd = sd(sqrt(408)*(2-l_hats2)))) +
#   stat_function(fun = dnorm, args = list(mean = mean(sqrt(408)*(2-l_hats2)), sd = sqrt(1/.283969)), color ="red") +
#   theme_bw()
# ggplot() + 
#   geom_density(aes(x = sqrt(408)*(2-l_hats2_clean)),fill = 'lightgray', col = 'black')+
#   stat_function(fun = dnorm, args = list(mean = mean(sqrt(408)*(2-l_hats2_clean)), sd = sd(sqrt(408)*(2-l_hats2_clean)))) +
#   stat_function(fun = dnorm, args = list(mean = mean(sqrt(408)*(2-l_hats2_clean)), sd = sqrt(1/.283969)), color ="red") +
#   theme_bw()
# var(sqrt(408)*(2-l_hats2_clean))
# 
# # pdf("./Simulation-Plot.pdf")
# ggplot() + 
#   geom_density(aes(x = l_hats2), fill = "grey") +
#   geom_vline(xintercept = 2, color = "red", linetype = 2) + 
#   labs(x = TeX(r"($\hat{\lambda})"), y = "Density", title = TeX(r"(Density of $\hat{\lambda}, m = 24480, k = 408)")) +
#   xlim(c(0,10)) +
#   theme_bw()
# dev.off()

l_hats3 <- rep(0,200)
flags <- rep(0,200)
for(j in 1:200){
  XY_scaled <- data.frame(matrix(nrow = 581, ncol = 2))
  for(k in 1:581){
    XY_scaled[k,] <- apply(mvrnorm(172117, mu = c(-bm172117^2,-bm172117^2), Sigma = bm172117^2*matrix(c(1, p(2,172117), p(2,172117), 1), ncol=2))
                           ,2,max)
  }
  for(i in 1:581){
    if((XY_scaled[i,1] > -log(.5) & XY_scaled[i,2] < -log(.5)) || (XY_scaled[i,1] < -log(.5) & XY_scaled[i,2] > -log(.5))){
      flags[j] <- flags[j] + abs(XY_scaled[i,1] - XY_scaled[i,2])
    }
  }
  l_hats3[j] <- optimize(log_lik, interval = c(0, 10), X = XY_scaled[,1], Y = XY_scaled[,2],
                         maximum = T, tol = .Machine$double.eps)$maximum
  if(l_hats3[j] == 10){
    print("flag")
    k <- 2
    while(l_hats3[j] == k*10 && k < 10){
      l_hats3[j] <- optimize(log_lik, interval = c(0, k*10), X = XY_scaled[,1], Y = XY_scaled[,2],
                             maximum = T, tol = .Machine$double.eps)$maximum
    }
  }
}
# pdf("./Greater-m-Sim.pdf", res = )
ggplot() +
  geom_density(aes(x = l_hats3), fill = "grey") +
  geom_vline(xintercept = 2, color = "red", linetype = 2) +
  labs(x = TeX(r"($\hat{\lambda})"), y = "Density", title = TeX(r"(Density of $\hat{\lambda}, m = 172117, k = 581)")) +
  xlim(c(0,10)) +
  theme_bw()
dev.off()



l_hats4 <- rep(0,200)
flags <- rep(0,200)
# for(j in 1:200){
  # XY_scaled <- data.frame(matrix(nrow = 408, ncol = 2))
  # for(k in 1:408){
    # XY_scaled[k,] <- apply(mvrnorm(24480, mu = c(-bm24480^2,-bm24480^2), Sigma = bm24480^2*matrix(c(1, p3(2,bm24480), p3(2,bm24480), 1), ncol=2))
                           # ,2,max)
  # }
  # for(i in 1:408){
    # if((XY_scaled[i,1] > -log(.5) & XY_scaled[i,2] < -log(.5)) || (XY_scaled[i,1] < -log(.5) & XY_scaled[i,2] > -log(.5))){
      # flags[j] <- flags[j] + abs(XY_scaled[i,1] - XY_scaled[i,2])
    # }
  # }
  # l_hats4[j] <- optimize(log_lik, interval = c(0, 10), X = XY_scaled[,1], Y = XY_scaled[,2],
                         # maximum = T, tol = .Machine$double.eps)$maximum
  # if(l_hats4[j] == 10){
    # print("flag")
    # k <- 2
    # while(l_hats4[j] == k*10 && k < 10){
      # l_hats4[j] <- optimize(log_lik, interval = c(0, k*10), X = XY_scaled[,1], Y = XY_scaled[,2],
                             # maximum = T, tol = .Machine$double.eps)$maximum
    # }
  # }
# }
# ggplot() + 
  # geom_density(aes(x = l_hats4), fill = "grey") +
  # geom_vline(xintercept = 2, color = "red", linetype = 2) + 
  # labs(x = TeX(r"($\hat{\lambda})"), y = "Density", title = TeX(r"(Density of $\hat{\lambda}, m = 24480, k = 408)")) +
  # xlim(c(0,10)) +
  # theme_bw()


l_hats5 <- rep(0,200)
flags <- rep(0,200)
for(j in 1:200){
  XY_scaled <- data.frame(matrix(nrow = 581, ncol = 2))
  for(k in 1:581){
    XY_scaled[k,] <- apply(mvrnorm(172117, mu = c(-bm172117^2,-bm172117^2), Sigma = bm172117^2*matrix(c(1, p3(2,bm172117), p3(2,bm172117), 1), ncol=2))
                           ,2,max)
  }
  for(i in 1:581){
    if((XY_scaled[i,1] > -log(.5) & XY_scaled[i,2] < -log(.5)) || (XY_scaled[i,1] < -log(.5) & XY_scaled[i,2] > -log(.5))){
      flags[j] <- flags[j] + abs(XY_scaled[i,1] - XY_scaled[i,2])
    }
  }
  l_hats5[j] <- optimize(log_lik, interval = c(0, 10), X = XY_scaled[,1], Y = XY_scaled[,2],
                         maximum = T, tol = .Machine$double.eps)$maximum
  if(l_hats5[j] == 10){
    print("flag")
    k <- 2
    while(l_hats5[j] == k*10 && k < 10){
      l_hats5[j] <- optimize(log_lik, interval = c(0, k*10), X = XY_scaled[,1], Y = XY_scaled[,2],
                             maximum = T, tol = .Machine$double.eps)$maximum
    }
  }
}
# pdf("./Larger-rho.pdf")
ggplot() + 
  geom_density(aes(x = l_hats5), fill = "grey") +
  geom_vline(xintercept = 2, color = "red", linetype = 2) + 
  labs(x = TeX(r"($\hat{\lambda})"), y = "Density", title = TeX(r"(Density of $\hat{\lambda}, m = 172117, k = 581)")) +
  xlim(c(0,10)) +
  theme_bw()
# dev.off()
var(l_hats5)

pdf("./Larger-rho-qq.pdf")
qqnorm(l_hats5)
qqline(l_hats5)
dev.off()

#---------------------------------------------------------------

####test
samples <- matrix(data = rep(0,20000),ncol = 2)
for(j in 1:10000){
  XY_scaled <- apply(mvrnorm(172117, mu = c(-bm172117^2,-bm172117^2), Sigma = bm172117^2*matrix(c(1, p3(2,bm172117), p3(2,bm172117), 1), ncol=2))
                           ,2,max)
  x = XY_scaled[1]
  y = XY_scaled[2]
  samples[j,] = exp(-x)*dnorm(2 + (y-x)/(2*2))*(-2 + pnorm(2 + (x-y)/(2*2))*exp(-y)*(1-(y-x)/(2*2^2)) + pnorm(2 + (y-x)/(2*2))*exp(-x)*(1-(x-y)/(2*2^2)) + 1/(2*2^2)+1/2-(y-x)^2/(8*2^4))/(exp(-x-y)*pnorm(2 + (x-y)/(2*2))*pnorm(2 + (y-x)/(2*2)) + exp(-x)*dnorm(2+(y-x)/(2*2))/(2*2))
}
mean(samples)



l_hats6 <- rep(0,200)
flags <- rep(0,200)
for(j in 1:200){
  XY_scaled <- data.frame(matrix(nrow = 400, ncol = 2))
  for(k in 1:400){
    XY_scaled[k,] <- apply(mvrnorm(249509, mu = c(-bm249509^2,-bm249509^2), Sigma = bm249509^2*matrix(c(1, p2(1.25,bm249509), p3(1.25,bm249509), 1), ncol=2))
                           ,2,max)
  }
  for(i in 1:400){
    if((XY_scaled[i,1] > -log(.5) & XY_scaled[i,2] < -log(.5)) || (XY_scaled[i,1] < -log(.5) & XY_scaled[i,2] > -log(.5))){
      flags[j] <- flags[j] + abs(XY_scaled[i,1] - XY_scaled[i,2])
    }
  }
  l_hats6[j] <- optimize(log_lik, interval = c(0, 10), X = XY_scaled[,1], Y = XY_scaled[,2],
                         maximum = T, tol = .Machine$double.eps)$maximum
  if(l_hats6[j] == 10){
    print("flag")
    k <- 2
    while(l_hats6[j] == k*10 && k < 10){
      l_hats6[j] <- optimize(log_lik, interval = c(0, k*10), X = XY_scaled[,1], Y = XY_scaled[,2],
                             maximum = T, tol = .Machine$double.eps)$maximum
    }
  }
}
# pdf("./Larger-rho.pdf")
ggplot() + 
  geom_density(aes(x = l_hats6), fill = "grey") +
  geom_vline(xintercept = 1.25, color = "red", linetype = 2) + 
  labs(x = TeX(r"($\hat{\lambda})"), y = "Density", title = TeX(r"(Density of $\hat{\lambda}, m = 249509, k = 400)")) +
  xlim(c(0,5)) +
  theme_bw()
beep()

qqnorm(l_hats6)
qqline(l_hats6)



l_hats7 <- rep(0,200)
for(j in 1:200){
  XY_scaled <- data.frame(matrix(nrow = 400, ncol = 2))
  for(k in 1:400){
    XY_scaled[k,] <- apply(mvrnorm(249509, mu = c(-bm249509^2,-bm249509^2), Sigma = bm249509^2*matrix(c(1, p3(2,bm249509), p3(2,bm249509), 1), ncol=2))
                           ,2,max)
  }
  l_hats7[j] <- optimize(log_lik, interval = c(0, 10), X = XY_scaled[,1], Y = XY_scaled[,2],
                         maximum = T, tol = .Machine$double.eps)$maximum
}
ggplot() + 
  geom_density(aes(x = l_hats7), fill = "grey") +
  geom_vline(xintercept = 2, color = "red", linetype = 2) + 
  labs(x = TeX(r"($\hat{\lambda})"), y = "Density", title = TeX(r"(Density of $\hat{\lambda}, m = 249509, k = 400)")) +
  xlim(c(0,10)) +
  theme_bw()
var(l_hats7)

qqnorm(l_hats7)
qqline(l_hats7)
beep()

sim1 <- rep(0,400)
sim2 <- rep(0,400)
for(j in 1:40000){
  XY_scaled <- apply(mvrnorm(5938, mu = c(-bm5938^2,-bm5938^2), Sigma = bm5938^2*matrix(c(1, p(1.25,bm5938), p(1.25,bm5938), 1), ncol=2))
                           ,2,max)
  x = XY_scaled[1]
  y = XY_scaled[2]
  sim1[j] <- dnorm(1.25 + (y-x)/(2*1.25))*(-2*exp(-x) + (pnorm(1.25 + (x-y)/(2*1.25))*exp(-y)*(1-(y-x)/(2*1.25^2)) + pnorm(1.25 + (x-y)/(2*1.25))*exp(-x)*(1-(y-x)/(2*1.25^2)) + 1/(2* 1.25^2) + .5 - (y-x)^2/(8*1.25^4))/(exp(-y)*pnorm(1.25 + (y-x)/(2*1.25)) * pnorm(1.25 + (x-y)/(2*1.25)) + dnorm(1.25 + (y-x)/(2*1.25))/(2*1.25)))
  # sim2[j] <- dnorm(2 + (y-x)/(2*2))*(-2*exp(-x) + (pnorm(2 + (x-y)/(2*2))*exp(-y)*(1-(y-x)/(2*2^2)) + pnorm(2 + (x-y)/(2*2))*exp(-y)*(1-(y-x)/(2*2^2)) + 1/(2* 2^2) + .5 - (y-x)^2/(8*2^4))/(exp(-y)*pnorm(2 + (y-x)/(2*2)) * pnorm(2 + (x-y)/(2*2)) + dnorm(2 + (y-x)/(2*2))/(2*2)))
}
mean(sim1)
mean(sim2)
beep()



l_hats8 <- rep(0,200)
for(j in 1:200){
  XY_scaled <- data.frame(matrix(nrow = 400, ncol = 2))
  for(k in 1:400){
    XY_scaled[k,] <- apply(mvrnorm(249509, mu = c(-bm249509^2,-bm249509^2), Sigma = bm249509^2*matrix(c(1, p(1.25,249509), p(1.25,249509), 1), ncol=2))
                           ,2,max)
  }
  l_hats8[j] <- optimize(log_lik, interval = c(0, 10), X = XY_scaled[,1], Y = XY_scaled[,2],
                         maximum = T, tol = .Machine$double.eps)$maximum
}
ggplot() + 
  geom_density(aes(x = l_hats8), fill = "grey") +
  geom_vline(xintercept = 2, color = "red", linetype = 2) + 
  labs(x = TeX(r"($\hat{\lambda})"), y = "Density", title = TeX(r"(Density of $\hat{\lambda}, m = 249509, k = 400)")) +
  xlim(c(0,10)) +
  theme_bw()
var(l_hats8)
qqnorm(l_hats8)
qqline(l_hats8)
beep()

bm5938 <- 3.60239
l_hats9 <- rep(0,200)
for(j in 1:200){
  XY_scaled <- data.frame(matrix(nrow = 168, ncol = 2))
  for(k in 1:168){
    XY_scaled[k,] <- apply(mvrnorm(5938, mu = c(-bm5938^2,-bm5938^2), Sigma = bm5938^2*matrix(c(1, p(1.25,5938), p(1.25,5938), 1), ncol=2))
                           ,2,max)
  }
  l_hats9[j] <- optimize(log_lik, interval = c(0, 10), X = XY_scaled[,1], Y = XY_scaled[,2],
                         maximum = T, tol = .Machine$double.eps)$maximum
}
ggplot() + 
  geom_density(aes(x = l_hats9), fill = "grey") +
  geom_vline(xintercept = 2, color = "red", linetype = 1.25) + 
  labs(x = TeX(r"($\hat{\lambda})"), y = "Density", title = TeX(r"(Density of $\hat{\lambda}, m = 5938, k = 168)")) +
  xlim(c(0,10)) +
  theme_bw()
var(l_hats9)
qqnorm(l_hats9)
qqline(l_hats9)
