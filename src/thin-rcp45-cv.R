## Load source
source("./src/load-wishart.R")

models = c("bcc-csm1-1","CanESM2","CCSM4","CNRM-CM5","FGOALS-s2","GFDL-ESM2G",
           "GISS-E2-R","HadGEM2-ES","IPSL-CM5A-MR","MIROC5","MPI-ESM-LR",
           "MRI-CGCM3","NorESM1-M")

for (model in models) {
  
  samples = load.samples(paste0("samples-rcp45-",model), 
                         "./ebm-paper/output", 4, 40)
  save(samples, file = paste0("./output/samples-rcp45-",model,"-10k.RData"))
  
}
