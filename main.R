####################################
## RCP4.5 analysis from main text ##
####################################

## Sampler
source("./src/sampling-rcp45.R")

## MCMC diagnostics
source("./src/check-rcp45.R")

## Thinning
source("./src/thin-rcp45.R")

## Figures
source("./src/figure2.R")
source("./src/figure3.R")
source("./src/figure4.R")
source("./src/figure6.R")
source("./src/figure7.R")
source("./src/figure8.R")

## Table 2
source("./src/table2.R")

## Time averages
source("./src/time-averages.R")

## Cross validation - huge run-time written to be split over cluster machine.
# source("./src/cross-validation-rcp45.R")
# source("./src/figure5.R")


#################################################
## RCP8.5 analysis from supplementary material ##
#################################################

## Sampler
source("./src/sampling-rcp85.R")

## MCMC diagnostics
source("./src/check-rcp85.R")

## Thinning
source("./src/thin-rcp85.R")

## Figures
source("./src/figure4-rcp85.R")
source("./src/figure6-rcp85.R")
source("./src/figure7-rcp85.R")
source("./src/figure8-rcp85.R")

## Table 2
source("./src/table2-rcp85.R")

## Time averages
source("./src/time-averages-rcp85.R")

## Cross validation - huge run-time written to be split over cluster machine.
# source("./src/cross-validation-rcp85.R")
# source("./src/figure5-rcp85.R")


######################################################
## Sensitivity analysis from supplementary material ##
######################################################

## All models

## Sampler
source("./src/sampling-rcp45-all.R")

## Thinning
source("./src/thin-rcp45-all.R")

## Figures
source("./src/figure6-all.R")
source("./src/figure7-all.R")
source("./src/figure8-all.R")

## Table 2
source("./src/table2-all.R")

## Time averages
source("./src/time-averages-all.R")


## Flat priors

## Sampler
source("./src/sampling-rcp45-flat.R")

## Thinning
source("./src/thin-rcp45-flat.R")

## Figures
source("./src/figure6-flat.R")
source("./src/figure7-flat.R")
source("./src/figure8-flat.R")

## Table 2
source("./src/table2-flat.R")

## Time averages
source("./src/time-averages-flat.R")


## kappa = 1.2

## Sampler
source("./src/sampling-rcp45-k1p2.R")

## Thinning
source("./src/thin-rcp45-k1p2.R")

## Figures
source("./src/figure6-k1p2.R")
source("./src/figure7-k1p2.R")
source("./src/figure8-k1p2.R")

## Table 2
source("./src/table2-k1p2.R")

## Time averages
source("./src/time-averages-k1p2.R")

