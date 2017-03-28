library(rstan)
options(mc.cores = parallel::detectCores())
options(warn=1)

data = readRDS("data.rds")

seed = 325845211
#this produces a large number of divergent transitions and very uncertain (sd > 50) values for interaction_weights[2,1]
testFit = stan('training.stan', data = data,control=list(adapt_delta = 0.99), seed = seed);
print(summary(testFit, pars = "interaction_weights")$summary)

#this does not produce divergent transitions and gives reasonable estimates for interaction_weights
testFitSimple = stan('training-simplified.stan', data = data,control=list(adapt_delta = 0.99), seed = seed);
print(summary(testFitSimple, pars = "interaction_weights")$summary)
