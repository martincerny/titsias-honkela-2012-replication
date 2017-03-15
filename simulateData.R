require(deSolve);
require(abind)

proteinODE <- function(t, state, parameters)
{
  with(as.list(c(state, parameters)), {
    dX = regulator(t) - degradation * x;
    
    list(dX) 
  })
}

targetODE <- function(t, state, parameters)
{
  with(as.list(c(state, parameters)), {
    regulatoryInput = bias + weight * log(protein(t));
    dX = basalTranscription + (sensitivity/(1 + exp(-regulatoryInput))) - degradation * x;
    
    list(dX)
  })
}

bindReplicates <- function(a,b) {
  dimsA = dim(a);
  if(length(dimsA) == 2) {
    a = array(a,c(dimsA[1],1,dimsA[2]))
  }
  abind(a,b, along = 2);
}

simulateData <- function(regulatorProfile, numIntegrationPoints = 10,numTargets = 4, numReplicates = 3)
{
  time = 1:length(regulatorProfile);
  
  regulatorReplicates = array(0, c(numReplicates, length(time)));
  regulatorReplicates[1,] = regulatorProfile;
  for(i in 2:numReplicates)
  {
    regulatorReplicates[i,] = regulatorProfile + rnorm(length(time), 0,1);
    regulatorReplicates[regulatorReplicates < 0.05] = 0.05;
  }
  
  step = 1/numIntegrationPoints;
  integrationTime = seq(from = 1, to = length(regulatorProfile) + 0.00001, by = step);
  
  proteinDegradation = 0.4;#exp(rnorm(1, -0.5,0.3));
  proteinInitialLevel = exp(rnorm(1, -0.5,0.3));
  
  initialConditions = array(exp(rnorm(numTargets * numReplicates, -0.5,1)), c(numReplicates, numTargets));
  basalTranscription = exp(rnorm(numTargets, -0.8,0.2));
  degradation = exp(rnorm(numTargets, 2,1));
  sensitivity = exp(rnorm(numTargets, 2,1));
  
  bias = rnorm(numTargets, 0, 1);
  
  interactionWeights = rnorm(numTargets, 0, 2);
  
  sigmaGenerator <- function(len) {abs(rcauchy(len,0,0.1)) + 0.00001};
  
  regulatorSigma = array(0, c(numReplicates, length(time)));
  regulatorObserved = array(0, c(numReplicates, length(time)));
  for(i in 1:numReplicates) 
  {
    regulatorSigma[i,] =sigmaGenerator(length(time));
    regulatorObserved[i,] = regulatorReplicates[i,] + (rnorm(length(time)) * regulatorSigma[i,]);
  }

  
  regulatorProtein = array(0, c(numReplicates, length(integrationTime)));
  for(i in 1:numReplicates) 
  {
    proteinODEParams = c(degradation = proteinDegradation, regulator = approxfun(time, regulatorReplicates[i,], rule=2));  
    regulatorProtein[i,] = ode( y = c(x = proteinInitialLevel), times = integrationTime, func = proteinODE, parms = proteinODEParams, method = "ode45")[,"x"];
  
  }
  minProtein = min(regulatorProtein)
  maxProtein = max(regulatorProtein)
  
  regulatorProtein[regulatorProtein < 0.05] = 0.05;
  
  spots = character(numTargets + 1);
  spots[1] = "reg";

  targetValues = array(0, c(numReplicates, numTargets, length(time)));
  targetObserved = array(0, c(numReplicates, numTargets, length(time)));
  targetSigma = array(0, c(numReplicates, numTargets, length(time)));
  for(i in 1:numTargets)
  {
    #Rejection sampling to get full regulation potential
    while (sign(minProtein * interactionWeights[i] + bias[i]) == sign(maxProtein * interactionWeights[i] + bias[i])
           || abs(minProtein * interactionWeights[i] - maxProtein * interactionWeights[i]) < 2)
    {
      bias[i] = rnorm(1, 0,1);
      interactionWeights[i] = rnorm(1, -0,2);

    }
    while((sensitivity[i]) / degradation[i] < 0.5) {
      degradation[i] = exp(rnorm(1, 2,1));
      sensitivity[i] = exp(rnorm(1, 2,1));
    }
    for(j in 1:numReplicates)
    {
      params = c(degradation = degradation[i], bias = bias[i], sensitivity = sensitivity[i], weight = interactionWeights[i], basalTranscription = basalTranscription[i], protein = approxfun(integrationTime, regulatorProtein[j,], rule=2));  
      
      targetValues[j,i,] = ode( y = c(x = initialConditions[j,i]), times = time, func = targetODE, parms = params, method = "ode45")[,"x"];
      
      targetSigma[j,i,] = sigmaGenerator(length(time));
      
      targetObserved[j,i,] = targetValues[j,i,] + rnorm(length(time)) * targetSigma[j,i,]
      
    }    
    
    spots[i + 1] = paste0("t",i);
  }
  
  observed = bindReplicates(regulatorObserved, targetObserved);
  observed[observed <= 0.05] = 0.05;
  
  data = list(
    y = observed,
    yvar = bindReplicates(regulatorSigma, targetSigma),
    genes = spots,
    times = time,
    detailedTime = integrationTime,
    numIntegrationPoints = numIntegrationPoints,
    experiments = 1:numReplicates,
    trueProtein = regulatorProtein,
    trueTargets = targetValues,
    targetSigma = targetSigma,
    proteinDegradation = proteinDegradation,
    params = list(
      initialConditions = initialConditions,
      degradation = degradation,
      basalTranscription = basalTranscription,
      sensitivity = sensitivity,
      weights = interactionWeights,
      bias = bias
    )
    , regulatorSpots = spots[1]
    , targetSpots = spots[2:(numTargets + 1)]
  )
}

numericalIntegration <- function(basalTranscription, degradation, initialCondition, sensitivity, weight, bias, protein, numTime, numIntegrationPoints){
  
  numericalResult = numeric(numTime);
  numericalResult[1] = initialCondition
  integrationTimeStep = 1 / numIntegrationPoints
  
  basalOverDegradation = basalTranscription / degradation;
  
  regulationInput = bias + weight * log(protein);
  synthesis = integrationTimeStep * ( 1 / (1 + exp(-regulationInput)))
  
  residual = -0.5 * synthesis[1];
  degradationPerStep = exp(-degradation * integrationTimeStep)
  
  # for(detailedTimeIndex in 2:length(data$detailedTime)){
  #   residual = (residual + synthesis[detailedTimeIndex - 1])  * degradationPerStep;
  #   if((detailedTimeIndex - 1) %% numIntegrationPoints == 0)
  #   {
  #     integral = residual + 0.5 * synthesis[detailedTimeIndex];
  #     time = ((detailedTimeIndex - 1) %/% numIntegrationPoints) + 1
  #     numericalResult[time] = basalOverDegradation + (initialCondition - basalOverDegradation) * exp(-degradation * (time - 1)) + sensitivity * integral
  #   }
  # }
  
  
  for(time in 2:numTime){
    integral = 0;
    previousDetailed = (time - 1) * numIntegrationPoints + 1;
    for(previousIndex in 1:previousDetailed)
    {
      if(previousIndex == 1 || previousIndex == previousDetailed)
      {
        h = 0.5;
      }
      else
      {
        h = 1;
      }
      degradationCoeff = exp(-degradation * (previousDetailed - previousIndex) * integrationTimeStep)
      integral = integral + h * synthesis[previousIndex] * degradationCoeff;
    }
    numericalResult[time] = basalOverDegradation + (initialCondition - basalOverDegradation) * exp(-degradation * (time - 1)) + sensitivity * integral
  }
  return(numericalResult)
}

testNumericalIntegrationSingle <- function(data, targetIndex = 1, replicate = 1, doPlot = TRUE) {

  #The actual numerical integration
  basalTranscription = data$params$basalTranscription[targetIndex];
  degradation = data$params$degradation[targetIndex];
  initialCondition = data$params$initialConditions[replicate, targetIndex];
  sensitivity = data$params$sensitivity[targetIndex];
  weight = data$params$weights[targetIndex];
  bias = data$params$bias[targetIndex];
  
  numIntegrationPoints = data$numIntegrationPoints
  
  numericalResult = numericalIntegration(basalTranscription, degradation, initialCondition, sensitivity, weight, bias, data$trueProtein[replicate,], length(data$times), numIntegrationPoints)
  
  numericalData = data.frame(time = data$times, value = numericalResult);
  numericalData$type = "numeric";
  
  odeData = data.frame(time = data$times, value = data$trueTargets[replicate, targetIndex,]);
  odeData$type = "ode";
  
  allData = rbind(numericalData, odeData);
  
  if(doPlot){
    print(ggplot(allData, aes(x = time, y = value, color = type)) + geom_line())
  }
  errors = numericalResult - data$trueTargets[replicate, targetIndex,];
  return(sqrt(mean(errors ^ 2)))
}

testNumericalIntegration <- function(numTargets = 5, numIntegrationPoints = 10, doPlot = TRUE) {
  simulatedData = simulateData(c(0.8,0.7,0.2,0.3,0.6,1.5,2.7,0.9,0.8,0.6,0.2,1.6), numIntegrationPoints, numTargets = numTargets)
  totalError = 0
  for(target in 1:numTargets){
    for(replicate in 1:length(simulatedData$experiments))
    {
      totalError = totalError + testNumericalIntegrationSingle(simulatedData, target, replicate, doPlot) ^ 2; 
    }
  }
  return (sqrt(totalError / (numTargets * length(simulatedData$experiments))))
}