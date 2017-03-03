require(rstan)

trainModel <- function(regulatorSpots, genesSpots, normalizedData, num_integration_points, ...)
{
  regulatorIndices = rowIDsFromSpotNames(normalizedData, regulatorSpots);
  genesIndices = rowIDsFromSpotNames(normalizedData, genesSpots);
  
  interactionMatrix = array(1,c(length(regulatorIndices),length(genesIndices)));
  
  numTime = length(normalizedData$times)
  
  numRegulators = length(regulatorIndices);
  numGenes = length(genesIndices);
  
  modelData = list(num_time = numTime, 
                   num_integration_points = num_integration_points, 
                   num_regulators = numRegulators, 
                    regulator_profiles_observed = array(normalizedData$y[regulatorIndices, 1:numTime], c(numRegulators,numTime)), 
                   regulator_profiles_sigma = array(sqrt(normalizedData$yvar[regulatorIndices, 1:numTime]),c(numRegulators,numTime)), 
                   num_genes = numGenes, 
                   gene_profiles_observed = array(normalizedData$y[genesIndices, 1:numTime], c(numGenes,numTime)), 
                   gene_profiles_sigma = array(sqrt(normalizedData$yvar[genesIndices, 1:numTime]), c(numGenes,numTime)), 
                   interaction_matrix = interactionMatrix
                );
  return(stan('training.stan', data = modelData, ...));
}