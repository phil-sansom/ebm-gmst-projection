####################################
## Figure 3 - Standardized errors ##
####################################

## Load libraries
library(parallel)
library(MASS)

##########
## Data ##
##########

## Load data
load("./data/cmip5.RData")
models = c("bcc-csm1-1","BNU-ESM","CanESM2","CCSM4","CNRM-CM5","CSIRO-Mk3-6-0",
           "FGOALS-s2","GFDL-ESM2M","GISS-E2-R","HadGEM2-ES","inmcm4",
           "IPSL-CM5A-LR","MIROC5","MPI-ESM-LR","MRI-CGCM3","NorESM1-M")
data = lapply(data, function(x)
  list(tas = x$tas[models], toa = x$toa[models]))
n.models = length(models)

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

#########################
## Cummins et al. fits ##
#########################

library(EBM)

inits = list(
  gamma = 2.0,
  C = c(5.0,20.0,100.0),
  kappa = c(1.0,2.0,1.0),
  epsilon = 1.0,
  sigma_eta = 0.5,
  sigma_xi = 0.5,
  F_4xCO2 = 5.0
)

cummins = mclapply(y0, function(x)
  FitKalman(inits, x[1,], x[2,]), mc.cores = 4)

# data("CMIP5")
# CMIP5 = CMIP5[c("BCC","BNU","CCCMA","NCAR","CNRM","CSIRO","IAP","GFDL","GISS",
#                 "MOHC","INM","IPSL","MIROC","MPIM","MRI","NCC")]
# 
# cummins = mclapply(CMIP5, function(x)
#   FitKalman(inits, x$temp, x$flux), mc.cores = 4)



################
## Simulation ##
################

## Load source
source("./src/dlm.R")
source("./src/functions.R")

## Joint distribution
post.sim = function(ltheta, X, nu, V, m0 = NULL, C0 = NULL) {
  
  if (is.null(nu)) {
    XX = X
  } else{
    XX = rbind(X, nu)
  }
  
  ## Make matrices
  matrices = make.matrices(ltheta)
  F = matrices$F
  G = matrices$G
  W = matrices$W
  B = matrices$B
  
  ## Extract dimensions
  rr = nrow(F)
  tt = ncol(X)
  nn = ncol(F)
  pp = nrow(XX)
  
  ## Intercept terms
  c = matrix(0, rr, tt)
  d = B %*% XX
  
  if (rr < nrow(F))
    F = F[1,,drop=FALSE]
  
  ## Forecast uncertainty
  if (is.null(V))
    V = diag(rep(1e-12, rr))
  
  ## Initial conditions
  if (is.null(m0))
    m0 = rep(0, nrow(G))
  if (is.null(C0)) {
    C = solve(diag(4^2) - 
                kronecker(G[1:4,1:4],G[1:4,1:4]), as.vector(W[1:4,1:4]))
    C = matrix(C, 4, 4)
    if (nn > 4) {
      C0 = matrix(0, 5, 5)
      C0[1:4,1:4] = C
      C0[5,5] = 1e-12
    } else {
      C0 = C
    }
  }
  
  ## Transform matrices for use in fkf
  F  = array(F, c(rr,nn,tt))
  G  = array(G, c(nn,nn,tt))
  V  = array(V, c(rr,rr,tt))
  W  = array(W, c(nn,nn,tt))
  
  sim0.dlm(F, G, V, W, c, d, m0, C0)
  
}

nu = NULL
V  = diag(rep(1e-12), 2)

fn = function(i, m) {
  buffer = NULL
  while (!is.numeric(buffer)) {
    ltheta = mvrnorm(1, cummins[[m]]$mle, cummins[[m]]$vcov)
    ltheta[11] = log(0.5) + ltheta[11]
    buffer = try(post.sim(ltheta, X1[1,,drop=FALSE], nu, V)$f, silent = TRUE)
  }
  return(buffer)
}

sims = list()
for (m in 1:n.models) {
  print(m)
  buffer = mclapply(1:1000, fn, m = m, mc.cores = 4)
  sims[[m]] = simplify2array(buffer)
} ## m

std.t = list()
std.n = list()
for (m in 1:n.models) {
  t = sims[[m]][1,,]
  n = sims[[m]][2,,]
  std.t[[m]] = (y1[[m]][1,] - apply(t, 1, mean, na.rm = TRUE)) / 
    apply(t, 1, sd, na.rm = TRUE)
  std.n[[m]] = (y1[[m]][2,] - apply(n, 1, mean, na.rm = TRUE)) /
    apply(n, 1, sd, na.rm = TRUE)
} ## m

std.t = simplify2array(std.t)
std.n = simplify2array(std.n)

mm = which(models == "CanESM2")
fit.t = apply(sims[[mm]][1,,], 1, mean, na.rm = TRUE)
lwr.t = apply(sims[[mm]][1,,], 1, quantile, probs = 0.05, na.rm = TRUE)
upr.t = apply(sims[[mm]][1,,], 1, quantile, probs = 0.95, na.rm = TRUE)
fit.n = apply(sims[[mm]][2,,], 1, mean, na.rm = TRUE)
lwr.n = apply(sims[[mm]][2,,], 1, quantile, probs = 0.05, na.rm = TRUE)
upr.n = apply(sims[[mm]][2,,], 1, quantile, probs = 0.95, na.rm = TRUE)


width  = 210/25.4 - 2
height = width/sqrt(2)
pdf("./fig/figure3.pdf", width = width, height = height, pointsize = 8 )

par(las = 1, mar = c(2.5,2.5,1,1)+0.1, mgp = c(1.5,0.5,0), tcl = -1/3,
    xaxs = "i", yaxs = "i", mfrow = c(2,2))

## Temperature
plot(1850:2100, y1[[mm]][1,], type = "n", ylim = c(-0.5,+3.5), 
     xlab = "Time", ylab = "GMST (K)",
     # ylab = expression(paste("GMST (", phantom()^degree ,"C)")), 
     yaxt = "n")

mtext("(a)", side = 3, line = 0, at = 1830)
axis(2, seq(-0.5,+3.5,+0.5))
axis(4, seq(-0.5,+3.5,+0.5), labels = FALSE)

lines(1850:2100, fit.t, col = "red", lwd = 2)
lines(1850:2100, lwr.t, col = "red")
lines(1850:2100, upr.t, col = "red")
lines(1850:2100, y1[[mm]][1,])

## Radiation
plot(1850:2100, y1[[mm]][2,], type = "n", ylim = c(-1.5,+2.0), 
     xlab = "Time", ylab = "Top-of-atmosphere radiation balance", 
     yaxt = "n")

mtext("(b)", side = 3, line = 0, at = 1830)
axis(2, seq(-1.5,+2.0,+0.5))
axis(4, seq(-1.5,+2.0,+0.5), labels = FALSE)

lines(1850:2100, fit.n, col = "red", lwd = 2)
lines(1850:2100, lwr.n, col = "red")
lines(1850:2100, upr.n, col = "red")
lines(1850:2100, y1[[mm]][2,])


## Temperature
plot(1850:2100, std.t[,1], type = "n", ylim = c(-9,+9), 
     xlab = "Time", ylab = "Standardised error", yaxt = "n")

mtext("(c)", side = 3, line = 0, at = 1830)
axis(2, seq(-8,+8,+2))
axis(4, seq(-8,+8,+2), labels = FALSE)

for (m in 1:n.models)
  lines(1850:2100, std.t[,m], col = "darkgrey")
lines(1850:2100, apply(std.t, 1, mean, na.rm = TRUE), lwd = 2)

abline(h =  0, lty = "dashed")
abline(h = -2, lty = "dotted")
abline(h = +2, lty = "dotted")

## Radiation
plot(1850:2100, std.n[,1], type = "n", ylim = c(-13,+6), 
     xlab = "Time", ylab = "Standardised error", yaxt = "n")

mtext("(d)", side = 3, line = 0, at = 1830)
axis(2, seq(-12,+6,+2))
axis(4, seq(-12,+6,+2), labels = FALSE)

for (m in 1:n.models)
  lines(1850:2100, std.n[,m], col = "darkgrey")
lines(1850:2100, apply(std.n, 1, mean, na.rm = TRUE), lwd = 2)

abline(h =  0, lty = "dashed")
abline(h = -2, lty = "dotted")
abline(h = +2, lty = "dotted")

dev.off()


pdf("./fig/figure3s.pdf", width = 210/25.4, height = 148/25.4, pointsize = 8)

for (m in 1:length(models)) {
  
  par(las = 1, mar = c(2.5,2.5,1,1)+0.1, mgp = c(1.5,0.5,0), tcl = -1/3,
      xaxs = "i", yaxs = "i")
  
  ## Temperature
  plot(1850:2100, y1[[m]][1,], type = "n", ylim = c(-1.0,+6.0), 
       ylab = expression(paste("GMST anomaly (", phantom()^degree ,"C)")),
       yaxt = "n", main = models[m])
  
  axis(2, seq(-1.0,+6.0,+1.0))
  axis(4, seq(-1.0,+6.0,+1.0), labels = FALSE)
  
  lines(1850:2100, apply(sims[[m]][1,,], 1, mean), col = "red", lwd = 2)
  lines(1850:2100, apply(sims[[m]][1,,], 1, quantile, probs = 0.05), col = "red")
  lines(1850:2100, apply(sims[[m]][1,,], 1, quantile, probs = 0.95), col = "red")
  lines(1850:2100, y1[[m]][1,])
  
} ## m

dev.off()