###################################
## Table 2 - Parameter estimates ##
###################################

## Load libraries
library(coda)
library(MASS)

## Load source
source("./src/dlm.R")

## Load RCP4.5 data
load("./output/samples-rcp45-flat-10k.RData")
thetas = as.matrix(samples$theta)
phis   = as.matrix(samples$phi  )
zetas  = as.matrix(samples$zeta )
nus    = as.matrix(samples$nus  )
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

## Reshape
new.phis = t(apply(phis, 1, sample.phis))
phis   = new.phis
thetas = exp(thetas)
zetas  = exp(zetas)
phis   = exp(phis)
thetas = apply(thetas, 2, mean)
zetas  = apply(zetas , 2, mean)
phis   = apply(phis  , 2, mean)
rm(new.phis)

pars = data.frame(
  gamma   = thetas[seq( 1,length(thetas),13)],
  C1      = thetas[seq( 2,length(thetas),13)],
  C2      = thetas[seq( 3,length(thetas),13)],
  C3      = thetas[seq( 4,length(thetas),13)],
  k1      = thetas[seq( 5,length(thetas),13)],
  k2      = thetas[seq( 6,length(thetas),13)],
  k3      = thetas[seq( 7,length(thetas),13)],
  epsilon = thetas[seq( 8,length(thetas),13)],
  sigma_f = thetas[seq( 9,length(thetas),13)],
  sigma_t = thetas[seq(10,length(thetas),13)],
  F_CO2   = thetas[seq(11,length(thetas),13)],
  F_vol   = thetas[seq(12,length(thetas),13)],
  sigma_d = thetas[seq(13,length(thetas),13)]
)
pars = rbind(pars,phis,zetas)

models = c("BCC-CSM1.1","CanESM2","CCSM4","CNRM-CM5","FGOALS-s2","GFDL-ESM2G",
           "GISS-E2-R","HadGEM2-ES","IPSL-CM5A-MR","MIROC5","MPI-ESM-LR",
           "MRI-CGCM3","NorESM1-M")
rownames(pars) = c(models,"Ensemble","Observations")

pars$gamma   = formatC(pars$gamma  , digits = 2, format = "f")
pars$C1      = formatC(pars$C1     , digits = 2, format = "f")
pars$C2      = formatC(pars$C2     , digits = 1, format = "f")
pars$C3      = formatC(pars$C3     , digits = 0, format = "f")
pars$k1      = formatC(pars$k1     , digits = 2, format = "f")
pars$k2      = formatC(pars$k2     , digits = 2, format = "f")
pars$k3      = formatC(pars$k3     , digits = 2, format = "f")
pars$epsilon = formatC(pars$epsilon, digits = 2, format = "f")
pars$sigma_f = formatC(pars$sigma_f, digits = 2, format = "f")
pars$sigma_t = formatC(pars$sigma_t, digits = 2, format = "f")
pars$F_CO2   = formatC(pars$F_CO2  , digits = 2, format = "f")
pars$F_vol   = formatC(pars$F_vol  , digits = 1, format = "f")
pars$sigma_d = formatC(pars$sigma_d, digits = 3, format = "f")
rownames(pars) = formatC(rownames(pars), digits = max(nchar(rownames(pars))), 
                         format = "s", flag = "-")

write.table(pars, "./tab/table2-flat.txt", quote = FALSE, sep = " & ", 
            col.names = FALSE)
