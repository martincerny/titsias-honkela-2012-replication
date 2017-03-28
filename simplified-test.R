reps = 3
targets = 1:4

regulatorRep1 = trainResult1$data$trueProtein[1:reps,]

targetRep1 = trainResult1$data$y[1:reps,targets,]
targetRep1Var = trainResult1$data$yvar[1:reps,targets,]

repeatedData <- function(n) {
  modelData = list(num_time = 12, 
                   num_integration_points = 10, 
                   num_replicates = reps,
                   log_tf_profiles = log(array(regulatorRep1,c(reps, 111))),
                   num_genes = n * length(targets), 
                   gene_profiles_observed = array(targetRep1, c(reps, n * length(targets), 12)), 
                   gene_profiles_sigma = array(sqrt(targetRep1Var), c(reps, n * length(targets),  12))
  );
  return(modelData)
}

testData <- function(n) {
  modelData = list(num_time = 12, 
                   num_integration_points = 10, 
                   num_replicates = reps,
                   log_tf_profiles = log(array(regulatorRep1,c(reps, 111))),
                   log_tf_profiles_source = log(array(regulatorRep1,c(reps, 111))),
                   num_regulators = 1,
                   num_genes = length(n), 
                   gene_profiles_observed = array(trainResult1$data$y[1:reps,n + 1,], c(reps, length(n), 12)), 
                   gene_profiles_sigma = array(sqrt(trainResult1$data$yvar[1:reps,n + 1,]), c(reps, length(n),  12))
  );
  return(modelData)
}


testDataPred <- function(n) {
  modelData = list(num_time = 12, 
                   num_integration_points = 10, 
                   num_replicates = reps,
                   num_regulators = 1,
                   tf_profiles = array(regulatorRep1,c(reps, 1, 111)),
                   gene_profile_observed = array(trainResult1$data$y[1:reps,n + 1,], c(reps, 12)), 
                   gene_profile_sigma = array(sqrt(trainResult1$data$yvar[1:reps,n + 1,]), c(reps, 12))
  );
  return(modelData)
}

testFit = stan('prediction.stan', data = testDataPred(1));  


options(warn = 1)
for(i in 1:10)
{
  cat(paste0(i," Training-simple\n"))
  testFit = stan('prediction.stan', data = testDataPred(i));  
}
# testFit = stan('training-simplified.stan', data = repeatedData(1), control=list(adapt_delta = 0.95));
# 
options(warn = 1)
for(i in 1:10)
{
  cat(paste0(i," Training-simple\n"))
  testFit = stan('training-simplified.stan', data = testData(i),control=list(adapt_delta = 0.99));  
}

options(warn = 1)
for(i in 1:10)
{
  cat(paste0(i," Prediction\n"))
  testSinglePrediction(trainResult1$data, i, 10, silent=TRUE)
}
# 
# for(i in 1:10)
# {
#   cat(paste0(i," Training\n"))
#   testTraining(simulatedData = trainResult1$data,numIntegrationPoints = 10,targetIndices = i)
#   print(warnings())
#   cat(paste0(i," Prediction\n"))
#   testSinglePrediction(trainResult1$data, i, 10)
#   print(warnings())
# }