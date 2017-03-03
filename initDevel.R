source("geneMapping.R")
source("training.r")
source("prediction.R")
source("simulateData.R")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
