## Load source
source("./src/load-wishart.R")

samples = load.samples("samples-rcp45-all", "output", 10, 200)

save(samples, file = "./output/samples-rcp45-all-10k.RData")
