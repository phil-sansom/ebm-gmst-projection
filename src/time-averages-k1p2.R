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
load("./output/samples-rcp45-k1p2-0.RData")
rm(burnin)
gc()

## Load RCP4.5 data
load("./output/samples-rcp45-k1p2-10k.RData")
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

# fn = function(i) {
#   ltheta = sample.phis(phis[i,])
#   post.sim(ltheta, X1, nus[i,], V0)$f[1,]
# }

gn = function(i) {
  mask = 1:ncol(z)
  lzeta = zetas[i,]
  fit = fit.model(lzeta, z, X1[,mask,drop=FALSE], Vz, nus[i,mask], rep(0,5))
  test = post.sim(lzeta, X1[,171:251], nus[i,171:251], array(1e-12, c(1,1,81)),
                  m0 = fit$m[,170], C0 = fit$C[,,170])$f[1,]
}

# # model.sims = mclapply(1:nrow(phis), fn, mc.cores = 4)
# # model.sims = simplify2array(model.sims)

world.sims = mclapply(1:nrow(zetas), gn, mc.cores = 4)
world.sims = simplify2array(world.sims)

save(world.sims, file = "./output/figure7-k1p2.RData")

# load("./output/figure7-k1p2.RData")
k1p2.sims = world.sims

load("./output/figure7.RData")

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
tau = 1986:2005

## Compute anomalies relative to reference period
Xp = X
for (m in 1:M)
  Xp[,m] = X[,m] - mean(X[t %in% tau,m], na.rm = TRUE)
fit = apply(Xp, 1, mean, na.rm = TRUE)
sds = apply(Xp, 1, sd  , na.rm = TRUE)
lwr = fit + qnorm(0.05)*sds
upr = fit + qnorm(0.95)*sds

Xs = apply(Xp[t %in% 2081:2100,], 2, mean)
fit = mean(Xs)
lwr = fit + qnorm(0.05)*sd(Xs)
upr = fit + qnorm(0.95)*sd(Xs)



###################
## Time averages ##
###################

t.sims = 2020:2100
proj = apply(world.sims[t.sims %in% 2081:2100,], 2, mean)
proj.k1p2 = apply(k1p2.sims[t.sims %in% 2081:2100,], 2, mean)

t.obs = 1850:2019
ref = mean(z[1,t.obs %in% 1986:2005])

mean(proj - ref)
quantile(proj - ref, c(0.05,0.50,0.95))

mean(proj.k1p2 - ref)
quantile(proj.k1p2 - ref, c(0.05,0.50,0.95))
