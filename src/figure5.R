#################################
## Figure 5 - Cross-validation ##
#################################

## Load libraries
library(coda)
library(parallel)

## Load source
source("./src/dlm.R")


##########
## Data ##
##########

## Load data
load("./data/cmip5.RData")
models = c("bcc-csm1-1","CanESM2","CCSM4","CNRM-CM5","FGOALS-s2","GFDL-ESM2G",
           "GISS-E2-R","HadGEM2-ES","IPSL-CM5A-MR","MIROC5","MPI-ESM-LR",
           "MRI-CGCM3","NorESM1-M")
data = lapply(data, function(x)
  list(tas = x$tas[models], toa = x$toa[models]))
models = names(data$abrupt4xCO2$tas)
n.models = length(models)
obs = read.table("./data/HadCRUT.4.6.0.0.annual_ns_avg.txt")
obs = obs[obs[,1] < 2020,]
sigmaz = (obs[,12] - obs[,11])/(qnorm(0.975) - qnorm(0.025))
obs = obs[,2] - mean(obs[1:50,2])

## abrupt4xCO2
y0 = list()
for (model in models)
  y0[[model]] = rbind(data$abrupt4xCO2$tas[[model]],
                      data$abrupt4xCO2$toa[[model]])
X0 = forcings$abrupt4xCO2$co2
X0 = matrix(X0, 1, length(X0))

## rcp45
y1 = list()
for (model in models)
  y1[[model]] = rbind(c(data$historical$tas[[model]],data$rcp45$tas[[model]]),
                      c(data$historical$toa[[model]],data$rcp45$toa[[model]]))
X1 = rbind( c(forcings$historical$co2,forcings$rcp45$co2),
            -c(forcings$historical$vol,forcings$rcp45$vol))
rm(model)

## Observations
z = matrix(obs, 1, length(obs))
Vz = array(sigmaz^2, c(1,1,length(sigmaz)))
kappa = 1.0


###############
## Load sims ##
###############

buffer = list()
for (model in models) {
  
  load(paste0("./output/samples-rcp45-",model,"-sims.RData"))
  buffer[[model]] = sims
  
}
sims = buffer
rm(buffer)

std.t = list()
std.n = list()
for (model in models) {
  t = sims[[model]][1,,]
  n = sims[[model]][2,,]
  std.t[[model]] = (y1[[model]][1,171:251] - apply(t, 1, mean, na.rm = TRUE)) / 
    apply(t, 1, sd, na.rm = TRUE)
  std.n[[model]] = (y1[[model]][2,171:251] - apply(n, 1, mean, na.rm = TRUE)) / 
    apply(n, 1, sd, na.rm = TRUE)
} ## m

std.t = simplify2array(std.t)
std.n = simplify2array(std.n)

mm = which(models == "HadGEM2-ES")
fit.t = apply(sims[[mm]][1,,], 1, mean, na.rm = TRUE)
lwr.t = apply(sims[[mm]][1,,], 1, quantile, probs = 0.05, na.rm = TRUE)
upr.t = apply(sims[[mm]][1,,], 1, quantile, probs = 0.95, na.rm = TRUE)
fit.n = apply(sims[[mm]][2,,], 1, mean, na.rm = TRUE)
lwr.n = apply(sims[[mm]][2,,], 1, quantile, probs = 0.05, na.rm = TRUE)
upr.n = apply(sims[[mm]][2,,], 1, quantile, probs = 0.95, na.rm = TRUE)



############
## Figure ##
############

width  = 210/25.4 - 2
height = width/sqrt(2)
pdf("./fig/figure5.pdf", width = width, height = height, pointsize = 8 )

par(las = 1, mar = c(2.5,2.5,1,1)+0.1, mgp = c(1.5,0.5,0), tcl = -1/3,
    xaxs = "i", yaxs = "i", mfrow = c(2,2))

## Temperature
plot(1850:2100, y1[[mm]][1,], type = "n", ylim = c(-0.5,+3.5), 
     xlab = "Time", #ylab = "GMST (K)",
     ylab = expression(paste("GMST anomaly (", phantom()^degree ,"C)")),
     yaxt = "n")

mtext("(a)", side = 3, line = 0, at = 1830)
axis(2, seq(-0.5,+3.5,+0.5))
axis(4, seq(-0.5,+3.5,+0.5), labels = FALSE)

lines(2020:2100, fit.t, col = "red", lwd = 2)
lines(2020:2100, lwr.t, col = "red")
lines(2020:2100, upr.t, col = "red")
lines(1850:2100, y1[[mm]][1,])

## Radiation
plot(1850:2100, y1[[mm]][2,], type = "n", ylim = c(-1.5,+2.0), 
     xlab = "Time", ylab = "Top-of-atmosphere radiation balance", 
     yaxt = "n")

mtext("(b)", side = 3, line = 0, at = 1830)
axis(2, seq(-1.5,+2.0,+0.5))
axis(4, seq(-1.5,+2.0,+0.5), labels = FALSE)

lines(2020:2100, fit.n, col = "red", lwd = 2)
lines(2020:2100, lwr.n, col = "red")
lines(2020:2100, upr.n, col = "red")
lines(1850:2100, y1[[mm]][2,])


## Temperature
plot(2020:2100, std.t[,1], type = "n", ylim = c(-4,+4), 
     xlab = "Time", ylab = "Standardised error", yaxt = "n")

mtext("(c)", side = 3, line = 0, at = 2014)
axis(2, seq(-4,+4,+1))
axis(4, seq(-4,+4,+1), labels = FALSE)

for (m in 1:n.models)
  lines(2020:2100, std.t[,m], col = "darkgrey")
lines(2020:2100, apply(std.t, 1, mean, na.rm = TRUE), lwd = 2)

abline(h =  0, lty = "dashed")
abline(h = -2, lty = "dotted")
abline(h = +2, lty = "dotted")

## Radiation
plot(2020:2100, std.n[,1], type = "n", ylim = c(-4,+4), 
     xlab = "Time", ylab = "Standardised error", yaxt = "n")

mtext("(d)", side = 3, line = 0, at = 2014)
axis(2, seq(-4,+4,+1))
axis(4, seq(-4,+4,+1), labels = FALSE)

for (m in 1:n.models)
  lines(2020:2100, std.n[,m], col = "darkgrey")
lines(2020:2100, apply(std.n, 1, mean, na.rm = TRUE), lwd = 2)

abline(h =  0, lty = "dashed")
abline(h = -2, lty = "dotted")
abline(h = +2, lty = "dotted")

dev.off()


pdf("./fig/figure5s.pdf", width = 210/25.4, height = 148/25.4, pointsize = 8)

for (m in 1:length(models)) {

  par(las = 1, mar = c(2.5,2.5,1,1)+0.1, mgp = c(1.5,0.5,0), tcl = -1/3,
      xaxs = "i", yaxs = "i")
  
  ## Temperature
  plot(1850:2100, y1[[m]][1,], type = "n", ylim = c(-1.0,+4.0), 
       ylab = expression(paste("GMST anomaly (", phantom()^degree ,"C)")),
       yaxt = "n", main = models[m])
  
  axis(2, seq(-1.0,+4.0,+1.0))
  axis(4, seq(-1.0,+4.0,+1.0), labels = FALSE)
  
  lines(2020:2100, apply(sims[[m]][1,,], 1, mean, na.rm = TRUE), 
        col = "red", lwd = 2)
  lines(2020:2100, apply(sims[[m]][1,,], 1, quantile, probs = 0.05, na.rm = TRUE), 
        col = "red")
  lines(2020:2100, apply(sims[[m]][1,,], 1, quantile, probs = 0.95, na.rm = TRUE), 
        col = "red")
  lines(1850:2100, y1[[m]][1,])
  
} ## m

for (m in 1:length(models)) {
  
  par(las = 1, mar = c(2.5,2.5,1,1)+0.1, mgp = c(1.5,0.5,0), tcl = -1/3,
      xaxs = "i", yaxs = "i")
  
  ## Temperature
  plot(1850:2100, y1[[m]][2,], type = "n", ylim = c(-3.0,+5.0), 
       ylab = "Top-of-atmosphere radiation balance",
       yaxt = "n", main = models[m])
  
  axis(2, seq(-3.0,+5.0,+1.0))
  axis(4, seq(-3.0,+5.0,+1.0), labels = FALSE)
  
  lines(2020:2100, apply(sims[[m]][2,,], 1, mean, na.rm = TRUE), 
        col = "red", lwd = 2)
  lines(2020:2100, apply(sims[[m]][2,,], 1, quantile, probs = 0.05, na.rm = TRUE),
        col = "red")
  lines(2020:2100, apply(sims[[m]][2,,], 1, quantile, probs = 0.95, na.rm = TRUE), 
        col = "red")
  lines(1850:2100, y1[[m]][2,])
  
} ## m

dev.off()
