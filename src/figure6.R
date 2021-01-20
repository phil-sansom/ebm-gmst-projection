################################################
## Figure 6 - Equilibrium Climate Sensitivity ##
################################################

## Load libraries
library(coda)
library(MASS)

## Load data
load("output/samples-rcp45-10k.RData")
thetas = as.matrix(samples$theta)
phis   = as.matrix(samples$phi  )
zetas  = as.matrix(samples$zeta )
rm(samples)
gc()

sample.phis = function(x) {
  mu = x[1:13]
  D = diag(x[14:26])
  C = matrix(0,13,13)
  C[lower.tri(C)] = x[27:104]
  C = C + t(C)
  diag(C) = 1
  Sigma = D %*% C %*% D
  mvrnorm(1, mu, Sigma)
}

new.phis = t(apply(phis, 1, sample.phis))
ecs.models = exp(new.phis[,11] - new.phis[,5])
round(mean(ecs.models), 1)
round(median(ecs.models), 1)
round(quantile(ecs.models, c(0.05,0.95)), 1)
density.models = density(ecs.models)
round(density.models$x[which.max(density.models$y)], 1)

fco2.world = exp(zetas[,11])
k1.world   = exp(zetas[,5])
ecs.world = fco2.world/k1.world
round(mean(ecs.world), 1)
round(median(ecs.world), 1)
round(quantile(ecs.world, c(0.05,0.95)), 1)
density.world = density(ecs.world)
round(density.world$x[which.max(density.world$y)],1)

fco2s = exp(thetas[,seq(11,169,13)])
k1s   = exp(thetas[,seq( 5,169,13)])
ecss  = fco2s/k1s
ecss  = apply(ecss, 2, mean)

## Plot ECS distributions
width  = (210/25.4 - 2)/2
height = width/sqrt(2)
pdf("./fig/figure6.pdf", width = width, height = height, pointsize = 8 )
par(las = 1, mar = c(3,3,1,1)+0.1, mgp = c(2,0.5,0), tcl = -1/3,
    xaxs = "i", yaxs = "i")
plot (density.models, col = "black", lwd = 2, xlim = c(0,8), ylim = c(0,0.6), 
      main = "", xaxt = "n", #xlab = "Equilibrium Climate Sensitivity (K)")
      xlab = expression(paste("Equilibrium Climate Sensitivity (",
                              phantom()^degree ,"C)")))
lines(density.world , col = 2, lwd = 2)
rug(ecss)
axis(1, 0:8)
legend("topright", legend = c("CMIP5 ensemble","Real world"), 
       col = c(1,2), lwd = c(2,2), lty = c(1,1), bty = "n")
dev.off()
