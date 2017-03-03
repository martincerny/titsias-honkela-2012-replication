require(deSolve);

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

simulateData <- function(regulatorProfile, numIntegrationPoints = 4,numTargets = 4)
{
  time = 1:length(regulatorProfile);
  step = 1/numIntegrationPoints;
  integrationTime = seq(from = 1, to = length(regulatorProfile) + step, by = step);
  
  proteinDegradation = 0.4;#exp(rnorm(1, -0.5,2));
  proteinInitialLevel = 0.8;#exp(rnorm(1, -0.5,2));
  
  initialConditions = exp(rnorm(numTargets, -0.5,2));
  basalTranscription = exp(rnorm(numTargets, -0.5,2));
  degradation = exp(rnorm(numTargets, -0.5,2));
  sensitivity = exp(rnorm(numTargets, -0.5,2));
  
  bias = rnorm(numTargets, 0, 1);
  
  interactionWeights = rnorm(numTargets, 0, 2);

  proteinODEParams = c(degradation = proteinDegradation, regulator = approxfun(time, regulatorProfile, rule=2));  
  regulatorProtein = ode( y = c(x = proteinInitialLevel), times = integrationTime, func = proteinODE, parms = proteinODEParams, method = "ode45")[,"x"];
  
  sigmaGenerator <- function(len) {abs(rcauchy(len,0,0.1)) + 0.00001};
  
  regulatorSigma =sigmaGenerator(length(time));
  regulatorObserved = regulatorProfile + (rnorm(length(time)) * regulatorSigma);
  
  spots = character(numTargets + 1);
  spots[1] = "reg";

  targetValues = array(0, c(numTargets, length(time)));
  targetObserved = array(0, c(numTargets, length(time)));
  targetSigma = array(0, c(numTargets, length(time)));
  for(i in 1:numTargets)
  {
    params = c(degradation = degradation[i], bias = bias[i], sensitivity = sensitivity[i], weight = interactionWeights[i], protein = approxfun(integrationTime, regulatorProtein, rule=2));  
    targetValues[i,] = ode( y = c(x = initialConditions[i]), times = time, func = targetODE, parms = params, method = "ode45")[,"x"];
    
    targetSigma[i,] = sigmaGenerator(length(time));
    
    targetObserved[i,] = targetValues[i,] + rnorm(length(time)) * targetSigma[i,]
    
    spots[i + 1] = paste0("t",i);
  }
  
  observed = rbind(regulatorObserved, targetObserved);
  observed[observed <= 0.05] = 0.05;
  
  data = list(
    y = observed,
    yvar = rbind(regulatorSigma, targetSigma),
    genes = spots,
    times = time,
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