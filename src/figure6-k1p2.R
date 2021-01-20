################################################
## Figure 6 - Equilibrium Climate Sensitivity ##
################################################

## Load libraries
library(coda)
library(MASS)

## Load data
load("./output/samples-rcp45-10k.RData")
zetas  = as.matrix(samples$zeta )
rm(samples)
gc()

fco2.world = exp(zetas[,11])
k1.world   = exp(zetas[,5])
ecs.world = fco2.world/k1.world
round(mean(ecs.world), 1)
round(median(ecs.world), 1)
round(quantile(ecs.world, c(0.05,0.95)), 1)
density.world = density(ecs.world)
round(density.world$x[which.max(density.world$y)],1)

## Load data
load("./output/samples-rcp45-k1p2-10k.RData")
zetas  = as.matrix(samples$zeta )
rm(samples)
gc()

fco2.k1p2 = exp(zetas[,11])
k1.k1p2   = exp(zetas[,5])
ecs.k1p2 = fco2.k1p2/k1.k1p2
round(mean(ecs.k1p2), 1)
round(median(ecs.k1p2), 1)
round(quantile(ecs.k1p2, c(0.05,0.95)), 1)
density.k1p2 = density(ecs.k1p2)
round(density.k1p2$x[which.max(density.k1p2$y)],1)


## Plot ECS distributions
width  = (210/25.4 - 2)/2
height = width/sqrt(2)
pdf("./fig/figure6-k1p2.pdf", width = width, height = height, pointsize = 8 )
par(las = 1, mar = c(3,3,1,1)+0.1, mgp = c(2,0.5,0), tcl = -1/3,
    xaxs = "i", yaxs = "i")
plot (density.world, col = 2, lwd = 2, 
      xlim = c(0,8), ylim = c(0,0.6), main = "", xaxt = "n", 
      #xlab = "Equilibrium Climate Sensitivity (K)")
      xlab = expression(paste("Equilibrium Climate Sensitivity (",
                              phantom()^degree ,"C)")))
lines(density.k1p2, col = 4, lwd = 2)
axis(1, 0:8)
legend("topright", legend = c(expression(paste(kappa," = 1.0")),
                              expression(paste(kappa," = 1.2"))), 
       col = c(2,4), lwd = c(2,2), lty = c(1,1), bty = "n")
dev.off()
