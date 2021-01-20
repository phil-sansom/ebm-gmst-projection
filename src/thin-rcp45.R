## Load source
source("./src/load-wishart.R")

samples = load.samples("samples-rcp45", "output", 10, 200)

save(samples, file = "./output/samples-rcp45-10k.RData")
