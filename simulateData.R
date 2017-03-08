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
    dX = sensitivity/(1 + exp(-regulatoryInput)) - degradation * x;
    
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
    while (sign(minProtein * interactionWeights[i] + bias[i]) == sign(maxProtein * interactionWeights[i] + bias[i])) 
    {
      bias[i] = rnorm(1, 0,1);
      interactionWeights[i] = rnorm(1, -0,2);
      
    }
    while((sensitivity[i] + basalTranscription[i]) / degradation[i] < 0.2) {
      degradation[i] = exp(rnorm(1, 2,1));
      sensitivity[i] = exp(rnorm(1, 2,1));
    }
    for(j in 1:numReplicates)
    {
      params = c(degradation = degradation[i], bias = bias[i], sensitivity = sensitivity[i], weight = interactionWeights[i], protein = approxfun(integrationTime, regulatorProtein[j,], rule=2));  
      
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