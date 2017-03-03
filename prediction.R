require(rstan)

predictModel <- function(tfProfiles, geneSpot, normalizedData, num_integration_points, ...)
{
  geneIndex = rowIDsFromSpotNames(normalizedData, geneSpot);
  
  numTime = length(normalizedData$times)
  
  numRegulators = 1;

  modelData = list(num_time = numTime, 
                   num_integration_points = num_integration_points, 
                   num_regulators = numRegulators, 
                   tf_profiles = array(tfProfiles, c(numRegulators,numTime * num_integration_points)), 
                   gene_profile_observed = array(normalizedData$y[geneIndex, 1:numTime], c(numTime)), 
                   gene_profile_sigma = array(sqrt(normalizedData$yvar[geneIndex, 1:numTime]), c(numTime))
  );
  return(stan('prediction.stan', data = modelData, ...));
}


plotPredictFit <- function(prediction, targetProfile, targetSigma, numSamples = 10, title = "") {
  true_value = extract(prediction,pars="gene_profile_true")$gene_profile_true
  
  samplesToPlot = true_value[sample(1:(dim(true_value)[1]),numSamples),];
  
  #Add the raw profile
  
  #plotData = t(rbind(targetProfile, samplesToPlot));
  #defaultWidth = 1;
  lineWidths = rep.int(1, numSamples);
  #lineWidths[1] = defaultWidth * 2;
  matplot(t(samplesToPlot), lwd = lineWidths,  type="l", main = title) 
  points(1:length(targetProfile), targetProfile, pch=19);
  arrows(1:12,targetProfile - targetSigma, 1:12, targetProfile + targetSigma,length=0.05, angle=90, code=3)
}

testSinglePrediction <- function(simulated_data, targetIndex, numIntegrationPoints, ...) {
  i = targetIndex;
  
  paramsOfInterest = c("initial_condition","basal_transcription","degradation","transcription_sensitivity","interaction_bias", "interaction_weights"); #, "model_mismatch_sigma"
  
  prediction = predictModel(simulated_data$trueProtein, paste0("t",i), simulated_data, numIntegrationPoints,  control = list(adapt_delta = 0.95), ...);
  resultSummary = summary(prediction, pars = paramsOfInterest)$summary[,c("mean","2.5%","25%","50%","75%","97.5%")];
  resultSummary = cbind(array(0,c(dim(resultSummary)[1],1)), resultSummary)
  
  colnames(resultSummary)[1] <- "true";
  resultSummary["initial_condition", "true"] = simulated_data$params$initialConditions[i];                                  
  resultSummary["basal_transcription", "true"] = simulated_data$params$basalTranscription[i];                                  
  resultSummary["degradation", "true"] = simulated_data$params$degradation[i];                                  
  resultSummary["transcription_sensitivity", "true"] = simulated_data$params$sensitivity[i];                                  
  resultSummary["interaction_bias", "true"] = simulated_data$params$bias[i];                                  
  resultSummary["interaction_weights[1]", "true"] = simulated_data$params$weights[i];                                  
  
  title = paste0("Result ",i);
  cat(paste0("\n",title,"\n"));
  print(resultSummary);    
  
  plotPredictFit(prediction, simulated_data$trueTargets[i,],simulated_data$targetSigma[i,], numSamples = 20, title = title)
  
  return(prediction);
}

testPrediction <- function(numTargets = 5, ...) {
  ow <- options("warn")
  options(warn = 1)
  
  fits = vector("list", numTargets);
  
  numIntegrationPoints = 5
  simulated_data = simulateData(c(0.8,0.7,0.2,0.3,0.6,1.5,2.7,0.9,0.8,0.6,0.2,1.6), numIntegrationPoints, numTargets = numTargets)
  for(i in 1:numTargets){
        
    fits[[i]] = testSinglePrediction(simulated_data, i, numIntegrationPoints, ...);
  }
  options(ow);
  return(list(fits = fits, data = simulated_data));
}
