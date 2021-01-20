############################
## Figure 7 - Projections ##
############################

## Load libraries
library(coda)
library(expm)
library(MASS)
library(parallel)

## Load source
source("./src/dlm.R")

## Load basic data
load("output/samples-rcp45-00.RData")
rm(burnin)
gc()

## Load RCP4.5 data
load("./output/samples-rcp45-10k.RData")
nus   = as.matrix(samples$nus )
phis  = as.matrix(samples$phi )
zetas = as.matrix(samples$zeta)
rm(samples)
gc()

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

V0 = diag(rep(1e-12,2))
V0 = array(1e-12, c(1, 1, 251))

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

fn = function(i) {
  ltheta = sample.phis(phis[i,])
  post.sim(ltheta, X1, nus[i,], V0)$f[1,]
}

gn = function(i) {
  mask = 1:ncol(z)
  lzeta = zetas[i,]
  fit = fit.model(lzeta, z, X1[,mask,drop=FALSE], Vz, nus[i,mask], rep(0,5))
  test = post.sim(lzeta, X1[,171:251], nus[i,171:251], array(1e-12, c(1,1,81)),
                  m0 = fit$m[,170], C0 = fit$C[,,170])$f[1,]
}

# model.sims = mclapply(1:nrow(phis), fn, mc.cores = 4)
# model.sims = simplify2array(model.sims)
# 
world.sims = mclapply(1:nrow(zetas), gn, mc.cores = 4)
world.sims = simplify2array(world.sims)

save(world.sims, file = "output/figure7.RData")

#load("output/figure7.RData")

##########
## IPCC ##
##########

## Load data
load("./data/cmip5-all-rcps.RData")
data = sapply(1:length(data$historical$tas), 
              function(m) c(data$historical$tas[[m]],data$rcp45$tas[[m]]))
X = data - 273.15
t = 1850:2100
T = nrow(X)
M = ncol(X)

## IPCC Reference period
tau = 1850:1900

## Compute anomalies relative to reference period
Xp = X
for (m in 1:M)
  Xp[,m] = X[,m] - mean(X[t %in% tau,m], na.rm = TRUE)
fit = apply(Xp, 1, mean, na.rm = TRUE)
sds = apply(Xp, 1, sd  , na.rm = TRUE)
lwr = fit + qnorm(0.05)*sds
upr = fit + qnorm(0.95)*sds


############
## Figure ##
############

width  = (210/25.4 - 2)/2
height = width/sqrt(2)
pdf("./fig/figure7.pdf", width = width, height = height, pointsize = 8)

par(las = 1, mar = c(2.5,2.5,1,1)+0.1, mgp = c(1.5,0.5,0), tcl = -1/3,
    xaxs = "i", yaxs = "i")

## RCP4.5
plot(1850:2100, y1[[1]][1,], type = "n", ylim = c(-0.5,+3.5), 
     xlab = "Time", #ylab = "GMST anomaly (K)",
     ylab = expression(paste("GMST anomaly(", phantom()^degree ,"C)")), 
     yaxt = "n")

axis(2, seq(0,3,1))
axis(4, seq(0,3,1), labels = FALSE)

# for (m in 1:n.models)
#   lines(1850:2100, y1[[m]][1,], col = "darkgrey")

lines(1850:2100, fit, lwd = 2)
lines(1850:2100, lwr)
lines(1850:2100, upr)

lines(1850:2019, z[1,], col = 2, lwd = 2)

lines(2020:2100, apply(world.sims, 1, mean), lwd = 2, col = 2)
lines(2020:2100, apply(world.sims, 1, quantile, probs = 0.05), col = 2)
lines(2020:2100, apply(world.sims, 1, quantile, probs = 0.95), col = 2)

legend("topleft", legend = c("IPCC method","EBM fit"), 
       col = c(1,2), lwd = c(2,2), lty = c(1,1), bty = "n")

dev.off()
