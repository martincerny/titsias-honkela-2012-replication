require(rstan)

trainModel <- function(regulatorSpots, genesSpots, normalizedData, num_integration_points, modelFile = "training.stan", control = NULL,...)
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
  if(is.null(control))
  {
    control = list(adapt_delta = 0.95)
  }
  
  return(stan(modelFile, data = modelData, control = control, ...));
}

plotTrainingResult <- function(prediction, replicate, tfIndex, numSamples = 20, samples = NULL, title = "") {
  if(is.null(samples))
  {
    samples = extract(prediction,pars=c("log_tf_profiles"));
  }
  
  true_value = samples$log_tf_profiles[,replicate,,tfIndex]
    
  numDetailedTime = dim(true_value)[2];
  numTime = prediction@par_dims$gene_profiles_true[3];

  detailedTime = ((1:numDetailedTime) - 1) * (numTime / (numDetailedTime + 1)) + 1;
  
  sampleIndices = sample(1:(dim(true_value)[1]),numSamples);
  
  samplesToPlot = exp(true_value[sampleIndices,]);

  ylim = c(0, min(10,max(samplesToPlot)))
  matplot(detailedTime, t(samplesToPlot), type="l", main = title, ylim = ylim) 
}

plotAllTrainingResults <- function(prediction, numSamples = 20) {
  numReplicates = prediction@par_dims$regulator_profiles_true[1] 
  numTFs = prediction@par_dims$regulator_profiles_true[2]
  for(tf in 1:numTFs)
  {
    for(replicate in 1:numReplicates)
    {
      plotTrainingResult(prediction, replicate, tf, title = paste0(tf," - ", replicate))
    }
  }
}  

plotTrainingTargetResult <- function(prediction, replicate, targetIndex, observedProfile =NULL, observedSigma = NULL, numSamples = 20, title = "", samples = NULL) {
  if(is.null(samples))
  {
    samples = extract(prediction,pars=c("gene_profiles_true"));
  }
  true_value = samples$gene_profiles_true[,replicate,targetIndex,]
  
  numTime = prediction@par_dims$gene_profiles_true[3];
  
  sampleIndices = sample(1:(dim(true_value)[1]),numSamples);
  
  samplesToPlot = true_value[sampleIndices,];
  
  ylim = c(0, min(10,max(samplesToPlot)))
  matplot(1:numTime, t(samplesToPlot), type="l", main = title, ylim = ylim) 
  if(!is.null(observedProfile))
  {
    points(1:numTime, observedProfile, pch=19);
    if(!is.null(observedSigma)){
      arrows(1:numTime, observedProfile - observedSigma - 0.001, 1:numTime,observedProfile + observedSigma ,length=0.05, angle=90, code=3)
    }
  }
}

plotAllTrainingTargetResults <- function(prediction, data, geneSpots,  numSamples = 20) {
  geneIndices = rowIDsFromSpotNames(data, geneSpots);
  
  samples = extract(prediction,pars=c("gene_profiles_true"));
 
  numReplicates = prediction@par_dims$gene_profiles_true[1] 
  numTargets = prediction@par_dims$gene_profiles_true[2]
  for(target in 1:numTargets)
  {
    for(replicate in 1:numReplicates)
    {
      plotTrainingTargetResult(prediction, replicate, target, title = paste0(geneSpots[target]," - ", replicate), observedProfile = data$y[replicate, geneIndices[target],], observedSigma = data$yvar[replicate, geneIndices[target],])
    }
  }
  
}

plotTrainingFitGraphics <- function(samplesToPlot, trueProtein, trueRNA, rnaSigma,numTime, sampleTime, title){
  ylim = c(0, max(c(max(samplesToPlot), max(trueProtein), max(trueRNA))))
  matplot(sampleTime, t(samplesToPlot), type="l", main = title, ylim = ylim) 
  points(1:numTime - 0.1, trueRNA, pch=1, col = "lightblue");
  arrows(1:numTime - 0.1, trueRNA - rnaSigma, 1:numTime - 0.1,trueRNA + rnaSigma ,length=0.05, angle=90, code=3, col = "lightblue")
  points(sampleTime, trueProtein, pch=19);
  
}

plotTrainingFit <- function(prediction, data, replicate, tfIndex, numSamples = 20, title = "", useODE = FALSE) {
  trueProtein = data$trueProtein[replicate,];
  trueRNA = data$y[replicate,1,];
  rnaSigma = data$yvar[replicate,1,];

  if(useODE) 
  {
    samples = extract(prediction,pars=c("log_tf_profiles","regulator_profiles_true", "protein_initial_level","protein_degradation"));
  }
  else {
    samples = extract(prediction,pars=c("log_tf_profiles"));
  }
  true_value = samples$log_tf_profiles[,replicate,,tfIndex]
  
  numDetailedTime = length(trueProtein);
  numTime = length(trueRNA);
  
  detailedTime = ((1:numDetailedTime) - 1) * (numTime / (numDetailedTime + 1)) + 1;
  
  sampleIndices = sample(1:(dim(true_value)[1]),numSamples);
  
  
  samplesToPlot = exp(true_value[sampleIndices,]);
  
  
  # for(s in 1:numSamples){
  #   samplesToPlot[s,] = samplesToPlot[s,] - trueProtein;
  # }
  
  if(numDetailedTime == dim(samplesToPlot)[2])
  {
    sampleTime = detailedTime
  }
  else
  {
    sampleTime = 1:numTime
  }

  if(useODE)
  {
    odeResults = array(0, c(numSamples, numDetailedTime));
    for(sample in 1:length(sampleIndices)){
      sampleIndex = sampleIndices[sample];
      rnaProfile = samples$regulator_profiles_true[sampleIndex,replicate, tfIndex,];
      proteinODEParams = c(degradation = samples$protein_degradation[sampleIndex, tfIndex], regulator = approxfun(detailedTime, rnaProfile, rule=2));  
      
      odeResults[sample,] = ode( y = c(x = samples$protein_initial_level[sampleIndex,replicate,tfIndex]), times = detailedTime, func = proteinODE, parms = proteinODEParams, method = "ode45")[,"x"];
      
    }
    plotTrainingFitGraphics(odeResults, trueProtein, trueRNA, rnaSigma,numTime, sampleTime, paste0(title," - ODE"));  
  }
    
  plotTrainingFitGraphics(samplesToPlot, trueProtein, trueRNA, rnaSigma,numTime, sampleTime, title)  
}

plotTrainingTargetFit <- function(prediction, data, targetIndex, replicate, numSamples = 20, title = "", useODE = TRUE)
{
  observedTarget = data$y[replicate,targetIndex + 1,]  
  targetSigma = data$yvar[replicate,targetIndex + 1,]  
  trueTarget = data$trueTargets[replicate,targetIndex,]  
  
  samples = extract(prediction,pars=c("gene_profiles_true","log_tf_profiles","initial_condition","basal_transcription","degradation","transcription_sensitivity","interaction_bias","interaction_weights"))
  sampleIndices = sample(1:(dim(samples$log_tf_profiles)[1]),numSamples);
  
  protein_profiles = exp(samples$log_tf_profiles[sampleIndices,replicate,,])
  
  numDetailedTime = dim(protein_profiles)[2];
  numTime = length(trueTarget);
  
  detailedTime = ((1:numDetailedTime) - 1) * (numTime / (numDetailedTime + 1)) + 1;
  
  if(useODE) {
    #solve the ODE to get the actual profiles of interest
    integrated_profile = array(-1, c(numSamples,numTime));   
    
    for(sample in 1:numSamples) {
      sampleId = sampleIndices[sample];
      params = c(degradation = samples$degradation[sampleId,targetIndex], bias = samples$interaction_bias[sampleId,targetIndex], sensitivity = samples$transcription_sensitivity[sampleId,targetIndex], basalTranscription = samples$basal_transcription[sampleId,targetIndex], weight = samples$interaction_weights[sampleId,targetIndex,1],  protein = approxfun(detailedTime, protein_profiles[sample,], rule=2));

      integrated_profile[sample,] = ode( y = c(x = samples$initial_condition[sampleId,replicate, targetIndex]), times = 1:numTime, func = targetODE, parms = params, method = "ode45")[,"x"];
    }
    plotPredictFitGraphics(integrated_profile, observedTarget, trueTarget, targetSigma, numTime, paste0(title, " - R ODE"))
    
    # for(sample in 1:numSamples) {
    #    sampleId = sampleIndices[sample];
    #    integrated_profile[sample,] = numericalIntegration(
    #      initialCondition = samples$initial_condition[sampleId,replicate, targetIndex],
    #      degradation = samples$degradation[sampleId, targetIndex],
    #      bias = samples$interaction_bias[sampleId, targetIndex],
    #      sensitivity = samples$transcription_sensitivity[sampleId, targetIndex],
    #      basalTranscription = samples$basal_transcription[sampleId, targetIndex],
    #      weight = samples$interaction_weights[sampleId, targetIndex, 1],
    #      protein = protein_profiles[sample,],
    #      numTime = numTime,
    #      numIntegrationPoints = data$numIntegrationPoints)
    #  }
    #  plotPredictFitGraphics(integrated_profile, observedTarget, trueTarget, targetSigma, numTime, paste0(title, " - R Numeric"))
  }
  
  samplesToPlot = samples$gene_profiles_true[sampleIndices,replicate,targetIndex,];
  plotPredictFitGraphics(samplesToPlot, observedTarget, trueTarget, targetSigma, numTime, title)

}

plotAllTargetFits <- function(prediction, simulatedData, simulatedDataIndices = 1:length(simulatedData$targetSpots), useODE = FALSE)
{
  for(target in 1:length(simulatedDataIndices)){
    for(replicate in 1:length(simulatedData$experiments)){
      simulatedIndex = simulatedDataIndices[target]
      plotTrainingTargetFit(prediction,simulatedData, simulatedIndex, replicate, title = paste0("Target ", simulatedIndex,"-",replicate), useODE = useODE)
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
  trueValues = data$trueRegulator[replicate,];
  
  points(1:numTime - 0.1, values, pch=1, col = "lightblue");
  arrows(1:numTime - 0.1, values - sigma, 1:numTime - 0.1,values + sigma ,length=0.05, angle=90, code=3, col = "lightblue")
  points(1:numTime, trueValues, pch=19);
  arrows(1:numTime, trueValues - sigma, 1:numTime,trueValues + sigma ,length=0.05, angle=90, code=3)
  #points(detailedTime, trueProtein, pch=19);
}

testTraining <- function(simulatedData = NULL,numIntegrationPoints = 10, numTargets = 10, ...) {
  if(is.null(simulatedData))
  {
    simulatedData = simulateData(c(0.8,0.7,0.2,0.3,0.6,1.5,2.7,0.9,0.8,0.6,0.2,1.6), numIntegrationPoints, numTargets = numTargets)
  }
  
  trainResult = trainModel(simulatedData$regulatorSpots, simulatedData$targetSpots, simulatedData, numIntegrationPoints, ...);
  
  tryCatch({
    plotAllTargetFits(trainResult, simulatedData);
    for(replicate in 1:length(simulatedData$experiments)){
      plotTrainingFit(trainResult, simulatedData, replicate,1, title = replicate, useODE = FALSE);
    }
  }, error = function(e) {
    print(e);
  });
  
  
  return(list(fit = trainResult, data = simulatedData));
}