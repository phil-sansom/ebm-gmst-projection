## Load source
source("./src/load-wishart.R")

samples = load.samples("samples-rcp45-flat", "output", 10, 200)

save(samples, file = "./output/samples-rcp45-flat-10k.RData")
