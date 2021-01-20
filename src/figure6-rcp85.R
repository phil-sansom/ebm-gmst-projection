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
load("./output/samples-rcp85-10k.RData")
zetas  = as.matrix(samples$zeta )
rm(samples)
gc()

fco2.rcp85 = exp(zetas[,11])
k1.rcp85   = exp(zetas[,5])
ecs.rcp85 = fco2.rcp85/k1.rcp85
round(mean(ecs.rcp85), 1)
round(median(ecs.rcp85), 1)
round(quantile(ecs.rcp85, c(0.05,0.95)), 1)
density.rcp85 = density(ecs.rcp85)
round(density.rcp85$x[which.max(density.rcp85$y)],1)


## Plot ECS distributions
width  = (210/25.4 - 2)/2
height = width/sqrt(2)
pdf("./fig/figure6-rcp85.pdf", width = width, height = height, pointsize = 8 )
par(las = 1, mar = c(3,3,1,1)+0.1, mgp = c(2,0.5,0), tcl = -1/3,
    xaxs = "i", yaxs = "i")
plot (density.world, col = 2, lwd = 2, 
      xlim = c(0,8), ylim = c(0,0.6), main = "", xaxt = "n", 
      #xlab = "Equilibrium Climate Sensitivity (K)")
      xlab = expression(paste("Equilibrium Climate Sensitivity (",
                              phantom()^degree ,"C)")))
lines(density.rcp85, col = 4, lwd = 2)
axis(1, 0:8)
legend("topright", legend = c("CMIP5 ensemble","Real world"), 
       col = c(2,4), lwd = c(2,2), lty = c(1,1), bty = "n")
dev.off()
