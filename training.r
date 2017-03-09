require(rstan)

trainModel <- function(regulatorSpots, genesSpots, normalizedData, num_integration_points, ...)
{
  regulatorIndices = rowIDsFromSpotNames(normalizedData, regulatorSpots);
  genesIndices = rowIDsFromSpotNames(normalizedData, genesSpots);
  
  interactionMatrix = array(1,c(length(regulatorIndices),length(genesIndices)));
  
  numTime = length(normalizedData$times)
  numReplicates = length(normalizedData$experiments)
  
  numRegulators = length(regulatorIndices);
  numGenes = length(genesIndices);
  
  modelData = list(num_time = numTime, 
                   num_integration_points = num_integration_points, 
                   num_regulators = numRegulators, 
                   num_replicates = numReplicates,
                   regulator_profiles_observed = array(normalizedData$y[,regulatorIndices,], c(numReplicates, numRegulators,numTime)), 
                   regulator_profiles_sigma = array(sqrt(normalizedData$yvar[,regulatorIndices,]),c(numReplicates,numRegulators,numTime)), 
                   num_genes = numGenes, 
                   gene_profiles_observed = array(normalizedData$y[, genesIndices ,], c(numReplicates, numGenes,  numTime)), 
                   gene_profiles_sigma = array(sqrt(normalizedData$yvar[,genesIndices, ]), c(numReplicates, numGenes,  numTime)), 
                   interaction_matrix = interactionMatrix
                );
  return(stan('training.stan', data = modelData, ...));
}

plotTrainingFit <- function(prediction, replicate, tfIndex, trueProtein, trueRNA, rnaSigma, numSamples = 20, title = "") {
  true_value = extract(prediction,pars="protein_profiles")$protein_profiles[,replicate,,tfIndex]
  
  numDetailedTime = length(trueProtein);
  numTime = length(trueRNA);
  
  detailedTime = ((1:numDetailedTime) - 1) * (numTime / (numDetailedTime + 1)) + 1;
  
  samplesToPlot = true_value[sample(1:(dim(true_value)[1]),numSamples),];
  
  matplot(detailedTime, t(samplesToPlot), type="l", main = title) 
  points(1:numTime - 0.1, trueRNA, pch=1, col = "lightblue");
  arrows(1:numTime - 0.1, trueRNA - rnaSigma, 1:numTime - 0.1,trueRNA + rnaSigma ,length=0.05, angle=90, code=3, col = "lightblue")
  points(detailedTime, trueProtein, pch=19);
}

plotTrainingTargetFit <- function(prediction, replicate, targetIndex, trueTarget, observedTarget, targetSigma, numSamples = 20, title = "")
{
  samples = extract(prediction,pars=c("protein_profiles","initial_conditions","basal_transcription","degradation","transcription_sensitivity","interaction_bias","interaction_weights"))
  sampleIndices = sample(1:(dim(samples$protein_profiles)[1]),numSamples);
  
  protein_profiles = samples$protein_profiles[,replicate,,]
  
  numDetailedTime = dim(protein_profiles)[2];
  numTime = length(trueTarget);
  
  detailedTime = ((1:numDetailedTime) - 1) * (numTime / (numDetailedTime + 1)) + 1;
  
  #solve the ODE to get the actual profiles of interest
  integrated_profile = array(-1, c(numSamples,numDetailedTime));
  
  for(sample in 1:numSamples) {
    sampleId = sampleIndices[sample];
    params = c(degradation = samples$degradation[sampleId,targetIndex], bias = samples$interaction_bias[sampleId,targetIndex], sensitivity = samples$transcription_sensitivity[sampleId,targetIndex], weight = samples$interaction_weights[sampleId,targetIndex,1], protein = approxfun(detailedTime, protein_profiles[sample,], rule=2));  
    
    integrated_profile[sample,] = ode( y = c(x = samples$initial_conditions[sampleId,replicate, targetIndex]), times = detailedTime, func = targetODE, parms = params, method = "ode45")[,"x"];
  }
  

  matplot(detailedTime, t(integrated_profile), type="l", main = title) 
  points(1:numTime - 0.1, observedTarget, pch=1, col = "lightblue");
  arrows(1:numTime - 0.1, observedTarget - targetSigma, 1:numTime - 0.1,observedTarget + targetSigma ,length=0.05, angle=90, code=3, col = "lightblue")
  points(1:numTime, trueTarget, pch=19);
  arrows(1:numTime, trueTarget - targetSigma, 1:numTime,trueTarget + targetSigma ,length=0.05, angle=90, code=3)
  
}

plotAllTargetFits <- function(prediction, simulatedData, simulatedDataIndices = 1:length(simulatedData$targetSpots))
{
  for(target in 1:length(simulatedDataIndices)){
    for(replicate in 1:length(simulatedData$experiments)){
      simulatedIndex = simulatedDataIndices[target]
      plotTrainingTargetFit(prediction,replicate,target,simulatedData$trueTargets[replicate,simulatedIndex,],simulatedData$y[replicate,simulatedIndex + 1,],simulatedData$yvar[replicate,simulatedIndex + 1,], title = paste0("Target ", target,"-",replicate))
    }
  }
    
}

plotRegulatorFit <- function(prediction, data, replicate = 1, tfIndex = 1, numSamples = 20, title = replicate) 
{
  true_value = extract(prediction,pars="regulator_profiles_true")$regulator_profiles_true[,replicate,tfIndex,]
  
  numDetailedTime = dim(true_value)[2];
  numTime = length(data$times);
  
  detailedTime = ((1:numDetailedTime) - 1) * (numTime / (numDetailedTime + 1)) + 1;
  
  samplesToPlot = true_value[sample(1:(dim(true_value)[1]),numSamples),];
  
  matplot(detailedTime, t(samplesToPlot), type="l", main = title) 
  values = data$y[replicate,tfIndex,];
  sigma = data$yvar[replicate,tfIndex,];
  points(1:numTime - 0.1, values, pch=19);
  arrows(1:numTime - 0.1, values - sigma, 1:numTime - 0.1,values + sigma ,length=0.05, angle=90, code=3)
  #points(detailedTime, trueProtein, pch=19);
}

testTraining <- function(simulatedData = NULL,numIntegrationPoints = 10, numTargets = 10, ...) {
  if(is.null(simulatedData))
  {
    simulatedData = simulateData(c(0.8,0.7,0.2,0.3,0.6,1.5,2.7,0.9,0.8,0.6,0.2,1.6), numIntegrationPoints, numTargets = numTargets)
  }
  
  trainResult = trainModel(simulatedData$regulatorSpots, simulatedData$targetSpots, simulatedData, numIntegrationPoints, ...);
  
  tryCatch({
  for(replicate in 1:length(simulatedData$experiments)){
    plotTrainingFit(trainResult,replicate,1,simulatedData$trueProtein[replicate,],simulatedData$y[replicate,1,], simulatedData$yvar[replicate,1,], title = replicate);
  }}, error = function(e) {});
  
  return(list(fit = trainResult, data = simulatedData));
}