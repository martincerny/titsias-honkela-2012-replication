require(rstan)

predictModel <- function(tfProfiles, geneSpot, normalizedData, numIntegrationPoints, ...)
{
  geneIndex = rowIDsFromSpotNames(normalizedData, geneSpot);
  
  numTime = length(normalizedData$times)
  numReplicates = length(normalizedData$experiments)
  
  numDetailedTime = getNumDetailedTime(numTime, numIntegrationPoints)
  
  numRegulators = 1;

  modelData = list(num_time = numTime, 
                   num_integration_points = numIntegrationPoints, 
                   num_regulators = numRegulators, 
                   num_replicates = numReplicates,
                   tf_profiles = array(tfProfiles, c(numReplicates,numRegulators,numDetailedTime)), 
                   gene_profile_observed = array(normalizedData$y[,geneIndex, 1:numTime], c(numReplicates, numTime)), 
                   gene_profile_sigma = array(sqrt(normalizedData$yvar[,geneIndex, 1:numTime]), c(numReplicates, numTime))
  );
  return(stan('prediction.stan', data = modelData, ...));
}


plotPredictFit <- function(prediction, data, targetIndex, replicate, numSamples = 20, title = "", useODE = FALSE) {
  samples = extract(prediction,pars=c("gene_profile_true","initial_condition","basal_transcription","degradation","transcription_sensitivity","interaction_bias","interaction_weights"))
  sampleIndices = sample(1:(dim(samples$initial_condition)[1]),numSamples);
  
  
  numTime = length(data$times);
  
  
  if(useODE){
    #solve the ODE to get the actual profiles of interest
    samplesToPlot = array(-1, c(numSamples,numTime));   
    
    trueProtein = data$trueProtein[replicate,];
    numDetailedTime = length(trueProtein);
    detailedTime = ((1:numDetailedTime) - 1) * (numTime / (numDetailedTime + 1)) + 1;
    
    proteinFun = approxfun(detailedTime, trueProtein, rule=2)
    
    for(sample in 1:numSamples) {
      sampleId = sampleIndices[sample];
      params = c(degradation = samples$degradation[sampleId], bias = samples$interaction_bias[sampleId], sensitivity = samples$transcription_sensitivity[sampleId], basalTranscription = samples$basal_transcription[sampleId], weight = samples$interaction_weights[sampleId,1],  protein = proteinFun);  
      
      samplesToPlot[sample,] = ode( y = c(x = samples$initial_condition[sampleId,replicate]), times = 1:numTime, func = targetODE, parms = params, method = "ode45")[,"x"];
    }
  }
  else {
    samplesToPlot = samples$gene_profile_true[sampleIndices,replicate, ];
  }
  
  #Add the raw profile
  observedProfile = data$y[replicate,targetIndex + 1,]  
  sigma = data$yvar[replicate,targetIndex + 1,]  
  targetProfile = data$trueTargets[replicate,targetIndex,]  
  
  lineWidths = rep.int(1, numSamples);
  matplot(t(samplesToPlot), lwd = lineWidths,  type="l", main = title) 
  points(1:numTime - 0.1, observedProfile, pch=1, col = "lightblue");
  arrows(1:numTime - 0.1, observedProfile - sigma, 1:numTime - 0.1, observedProfile + sigma, length=0.05, angle=90, code=3, col = "lightblue")
  points(1:numTime, targetProfile, pch=19);
  arrows(1:numTime, targetProfile - sigma, 1:numTime, targetProfile + sigma,length=0.05, angle=90, code=3)
}

testSinglePrediction <- function(simulatedData, targetIndex, numIntegrationPoints, ...) {
  i = targetIndex;
  
  paramsOfInterest = c("initial_condition","basal_transcription","degradation","transcription_sensitivity","interaction_bias", "interaction_weights"); #, "model_mismatch_sigma"
  
  prediction = predictModel(simulatedData$trueProtein, paste0("t",i), simulatedData, numIntegrationPoints,  control = list(adapt_delta = 0.95), ...);
  resultSummary = summary(prediction, pars = paramsOfInterest)$summary[,c("mean","2.5%","25%","50%","75%","97.5%")];
  resultSummary = cbind(array(0,c(dim(resultSummary)[1],1)), resultSummary)
  
  colnames(resultSummary)[1] <- "true";
  for(replicate in 1:length(simulatedData$experiments)) {
    resultSummary[paste0("initial_condition[", replicate,"]"), "true"] = simulatedData$params$initialConditions[replicate,i];                                  
  }
  resultSummary["basal_transcription", "true"] = simulatedData$params$basalTranscription[i];                                  
  resultSummary["degradation", "true"] = simulatedData$params$degradation[i];                                  
  resultSummary["transcription_sensitivity", "true"] = simulatedData$params$sensitivity[i];                                  
  resultSummary["interaction_bias", "true"] = simulatedData$params$bias[i];                                  
  resultSummary["interaction_weights[1]", "true"] = simulatedData$params$weights[i];                                  
  
  title = paste0("Result ",i);
  cat(paste0("\n",title,"\n"));
  print(resultSummary);    
  
  for(replicate in 1:length(simulatedData$experiments))
  {
    plotPredictFit(prediction, simulatedData, i, replicate, numSamples = 20, title = paste0(title,"-",replicate, " Modelled"))
    plotPredictFit(prediction, simulatedData, i, replicate, numSamples = 20, title = paste0(title,"-",replicate, " ODE"), useODE = TRUE)
  }
  
  return(prediction);
}

testPrediction <- function(numTargets = 5, numIntegrationPoints = 10, ...) {
  ow <- options("warn")
  options(warn = 1)
  
  fits = vector("list", numTargets);
  
  simulatedData = simulateData(c(0.8,0.7,0.2,0.3,0.6,1.5,2.7,0.9,0.8,0.6,0.2,1.6), numIntegrationPoints, numTargets = numTargets)
  for(i in 1:numTargets){
        
    fits[[i]] = testSinglePrediction(simulatedData, i, numIntegrationPoints, ...);
  }
  options(ow);
  return(list(fits = fits, data = simulatedData));
}

averageSamplingTime <- function(fits)
{
  timeList = lapply(testResults$fits, get_elapsed_time)
  allTimes = Reduce(rbind,timeList, array(0,c(0,2)))
  warmupTimes = allTimes[,"warmup"]
  sampleTimes = allTimes[,"sample"]
  return(list(total = mean(warmupTimes + sampleTimes), sample = mean(sampleTimes)))
}