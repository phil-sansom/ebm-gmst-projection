####################################
## Figure 4 - Shared forcing bias ##
####################################

## Load libraries
library(coda)

## Load RCP4.5 data
load("output/samples-rcp45-10k.RData")
nus = as.matrix(samples$nus)

mus = nus
for (i in 1:nrow(mus))
  mus[i,] = cumsum(mus[i,])

fit = apply(mus, 2, mean)
lwr = apply(mus, 2, quantile, probs = 0.05)
upr = apply(mus, 2, quantile, probs = 0.95)

## Plot nus
width  = (210/25.4 - 2)/2
height = width/sqrt(2)
pdf("./fig/figure4.pdf", width = width, height = height, pointsize = 8 )
par(las = 1, mar = c(2.5,2.5,1,1)+0.1, mgp = c(1.5,0.5,0), 
    tcl = -1/3, xaxs = "i", yaxs = "i")

## RCP 4.5
plot (1850:2100, fit, type = "n", ylim = range(-0.3,+1.0), main = "", 
      xlab = "Time", ylab = expression(paste("Shared forcing discrepancy (", 
                                             mu[t],")")))
axis(4, seq(-0.2,+1.0,+0.2), labels = FALSE)
abline(h = 0, lty = "dashed")
polygon(c(rev(1850:2100),1850:2100), c(rev(upr),lwr), border = NA, 
        col = "grey")
lines(1850:2100, fit, lwd = 2)

dev.off()
