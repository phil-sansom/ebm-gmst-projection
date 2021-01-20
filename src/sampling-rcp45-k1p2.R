## Load libraries
library(parallel)
mc.cores = 4
options(mc.cores = mc.cores)

## Load source
source("./src/dlm.R")
source("./src/functions.R")

## Sampling parameters
n.chains  = 4
n.burnin  = 25000
n.samples = 25000
n.blocks  = 4
n.thin    = 1

## Output parameters
output.dir = "output"
base.name  = "samples-rcp45-k1p2"

## Checks
print(paste("Chain samples =", n.samples*n.blocks))
print(paste("Total samples =", n.samples*n.chains*n.blocks))
print(paste("Output samples =", n.samples*n.chains*n.blocks/n.thin))

## Transform sampling parameters
a.burnin  = ceiling(n.burnin  / n.thin)
a.samples = ceiling(n.samples / n.thin)


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
kappa = 1.2


######################
## Hyper-parameters ##
######################

## Parameter names
theta.names = c("gamma","c1","c2","c3","k1","k2","k3","epsilon",
                "sigmaf","sigmat","fco2","fvol","sigmad")
phi.names   = c("gammam","c1m","c2m","c3m","k1m","k2m","k3m","epsilonm",
                "sigmafm","sigmatm","fco2m","fvolm","sigmadm",
                "gammap","c1p","c2p","c3p","k1p","k2p","k3p","epsilonp",
                "sigmafp","sigmatp","fco2p","fvolp","sigmadp",
                'rho1.2','rho1.3','rho1.4','rho1.5','rho1.6','rho1.7',
                'rho1.8','rho1.9','rho1.10','rho1.11','rho1.12','rho1.13',
                'rho2.3','rho2.4','rho2.5','rho2.6','rho2.7',
                'rho2.8','rho2.9','rho2.10','rho2.11','rho2.12','rho2.13',
                'rho3.4','rho3.5','rho3.6','rho3.7','rho3.8',
                'rho3.9','rho3.10','rho3.11','rho3.12','rho3.13',
                'rho4.5','rho4.6','rho4.7','rho4.8','rho4.9','rho4.10',
                'rho4.11','rho4.12','rho4.13','rho5.6','rho5.7',
                'rho5.8','rho5.9','rho5.10','rho5.11','rho5.12','rho5.13',
                'rho6.7','rho6.8','rho6.9','rho6.10','rho6.11',
                'rho6.12','rho6.13','rho7.8','rho7.9','rho7.10',
                'rho7.11','rho7.12','rho7.13','rho8.9','rho8.10',
                'rho8.11','rho8.12','rho8.13','rho9.10','rho9.11',
                'rho9.12','rho9.13','rho10.11','rho10.12','rho10.13',
                'rho11.12','rho11.13','rho12.13')

## Parameters dimensions
n.theta  = length(theta.names)
n.phi    = length(phi.names)
n.zeta   = length(theta.names)

## Hyper-parameters
hyper.parameters = list(
  mum = log(c(2,5,20,100,1,2,1,1,0.5,0.5,3,20,0.05)),
  mup = diag(rep(1/(log(10)/3)^2, n.theta)),
  Sigman  = n.theta  , Sigmas  = diag(rep(1e-3, n.theta)),
  sigmanm = log(0.10), sigmanp = 1/(log(10)/3)^2
)


##############
## Sampling ##
##############

## Initialization
print("Initialization...")
source("./src/initialization.R")

## Burn-in
print("Burn-in...")
burnin = mclapply(1:n.chains, wrapper, a.burnin, n.thin,
                  y0, X0, y1, X1, z, Vz, kappa,
                  n.theta, n.phi, n.zeta, initial.values, hyper.parameters, 
                  random.seeds, Sigma, adaptive = TRUE, jj = 0, 
                  theta.names, phi.names)

## Save burnin
save.image(paste0(output.dir,"/",base.name,"-",
                  formatC(0, width = nchar(n.blocks), flag = 0), ".RData"))

## Extract initial conditions for next block
random.seeds   = lapply(burnin, function(x) x$random.seed)
initial.values = lapply(burnin, function(x) x$initial.values)
Sigma          = lapply(burnin, function(x) x$Sigma)

## Garbage collection
rm(burnin)
gc()

## Sampling
print("Sampling...")
for (b in 1:n.blocks) {
  
  print(paste("Block", b))

  samples = mclapply(1:n.chains, wrapper, a.samples, n.thin, 
                     y0, X0, y1, X1, z, Vz, kappa,
                     n.theta, n.phi, n.zeta, initial.values, hyper.parameters, 
                     random.seeds, Sigma, adaptive = FALSE, jj = 0, 
                     theta.names, phi.names)

  file.name = paste0(output.dir, "/", base.name, "-",
                     formatC(b, width = nchar(n.blocks), flag = 0), ".RData")
  save(samples, file = file.name)
  
  ## Extract initial conditions for next block
  random.seeds   = lapply(samples, function(x) x$random.seed)
  initial.values = lapply(samples, function(x) x$initial.values)
  Sigma          = lapply(samples, function(x) x$Sigma)
  
  ## Garbage collection
  rm(samples)
  gc()

}
