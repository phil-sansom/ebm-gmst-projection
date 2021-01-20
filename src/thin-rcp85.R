## Load source
source("./src/load-wishart.R")

samples = load.samples("samples-rcp85", "output", 10, 200)

save(samples, file = "./output/samples-rcp85-10k.RData")
