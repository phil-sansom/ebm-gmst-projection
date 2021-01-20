## Load

## Load source
source("./src/load-wishart.R")

## Load samples
samples = load.samples("samples-rcp45", "output", 10, 1)

samples.theta = effectiveSize(samples$theta )
samples.phi   = effectiveSize(samples$phi   )
samples.zeta  = effectiveSize(samples$zeta  )
samples.nu    = effectiveSize(samples$nus   )
samples.sigma = effectiveSize(samples$sigman)

diag.theta  = gelman.diag(samples$theta, autoburnin = FALSE,
                          multivariate = FALSE)$psrf[,1]
diag.phi    = gelman.diag(samples$phi  , autoburnin = FALSE,
                          multivariate = FALSE)$psrf[,1]
diag.zeta   = gelman.diag(samples$zeta , autoburnin = FALSE,
                          multivariate = FALSE)$psrf[,1]
diag.nu     = gelman.diag(samples$nus  , autoburnin = FALSE,
                          multivariate = FALSE)$psrf[,1]
diag.sigma  = gelman.diag(samples$sigma, autoburnin = FALSE,
                          multivariate = FALSE)$psrf[,1]

acc.theta = samples$acceptance.theta
acc.zeta  = samples$acceptance.zeta
acc.sigma = samples$acceptance.sigman

save(samples.theta,samples.phi,samples.zeta,samples.nu,samples.sigma,
     diag.theta,diag.phi,diag.zeta,diag.nu,diag.sigma,
     acc.theta,acc.zeta,acc.sigma,
     file = "output/samples-rcp45-diag.RData")
