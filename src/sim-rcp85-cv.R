############################
## Figure 7 - Projections ##
############################

## Load libraries
library(coda)
library(expm)
library(MASS)
library(parallel)
mc.cores = 4
options(mc.cores = mc.cores)

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

gn = function(i) {
  mask = 1:ncol(z)
  lzeta = zetas[i,]
  fit = fit.model(lzeta, z, X1[,mask,drop=FALSE], Vz, nus[i,mask], rep(0,5))
  sim = try(post.sim(lzeta, X1[,171:251], nus[i,171:251], 
                     array(1e-12, c(1,1,81)),
                     m0 = fit$m[,170], C0 = fit$C[,,170])$f, silent = TRUE)
  if(class(sim) == "try-error")
    sim = matrix(NA, 2, 81)
  return(sim)
}

## Load data
load("./data/cmip5.RData")
models = c("bcc-csm1-1","CanESM2","CCSM4","CNRM-CM5","FGOALS-s2","GFDL-ESM2G",
           "GISS-E2-R","HadGEM2-ES","IPSL-CM5A-MR","MIROC5","MPI-ESM-LR",
           "MRI-CGCM3","NorESM1-M")
X1 = rbind( c(forcings$historical$co2,forcings$rcp85$co2),
            -c(forcings$historical$vol,forcings$rcp85$vol))

for (model in models) {
  
  load(paste0("./output/samples-rcp85-",model,"-10k.RData"))
  zetas = as.matrix(samples$zeta)
  nus   = as.matrix(samples$nus)
  
  z  = c(data$historical$tas[[model]],data$rcp45$tas[[model]][1:15])
  z  = matrix(z, 1, length(z))
  Vz = array(1e-12, c(1,1,length(z)))
  
  sims = mclapply(1:nrow(zetas), gn, mc.cores = 4)
  sims = simplify2array(sims)
  
  save(sims, file = paste0("./output/samples-rcp85-",model,"-sims.RData"))
  
  rm(sims, samples, zetas, nus, z, Vz)
  gc()

}  
