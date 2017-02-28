require(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

trainModel <- function(regulatorSpots, genesSpots, normalizedData, ...)
{
  regulatorIndices = rowIDsFromSpotNames(normalizedData, regulatorSpots);
  genesIndices = rowIDsFromSpotNames(normalizedData, genesSpots);
  
  interactionMatrix = array(1,c(length(regulatorIndices),length(genesIndices)));
  
  numRegulators = length(regulatorIndices);
  numGenes = length(genesIndices);
  
  modelData = list(num_time = length(drosophila_gpsim_normalized$times), 
                   num_integration_points = 10, 
                   num_regulators = numRegulators, 
                   regulator_profiles_observed = array(normalizedData$y[regulatorIndices, 1:12], c(numRegulators,12)), 
                   regulator_profiles_sigma = array(sqrt(normalizedData$yvar[regulatorIndices, 1:12]),c(numRegulators,12)), 
                   num_genes = numGenes, 
                   gene_profiles_observed = array(normalizedData$y[genesIndices, 1:12], c(numGenes,12)), 
                   gene_profiles_sigma = array(sqrt(normalizedData$yvar[genesIndices, 1:12]), c(numGenes,12)), 
                   interaction_matrix = interactionMatrix
                );
  return(stan('training.stan', data = modelData, ...));
}