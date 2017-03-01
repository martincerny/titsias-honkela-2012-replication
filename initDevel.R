source("geneMapping.R")
source("training.r")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
