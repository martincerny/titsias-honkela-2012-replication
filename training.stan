functions {
  void transcription_params_prior_lp(vector v, int num_elements)
  {
    vector[num_elements] log_v = log(v);
    target += normal_lpdf(log_v | -0.5,2);
    target += -log_v;
  }
}

data {
  int num_time;
  
  int num_integration_points;
  
  int num_regulators; 
  vector<lower = 0>[num_time] regulator_profiles_observed[num_regulators];
  vector<lower = 0>[num_time] regulator_profiles_sigma[num_regulators];
  
  int num_targets;
  vector<lower = 0>[num_time] target_profiles_observed[num_targets];
  vector<lower = 0>[num_time] target_profiles_sigma[num_targets];
  
  //0 regulation impossible, 1 regulation possible
  int interaction_matrix[num_regulators, num_targets];
}

transformed data {
  int num_detailed_time = num_time * num_integration_points;
  real integration_step = 1.0 / num_integration_points;
  
  vector[num_detailed_time] zero_mean;

  for (i in 1:num_detailed_time)  
  {
    zero_mean[i] = 0;
  }
}

parameters {
  vector[num_detailed_time] regulator_profiles_true_gp[num_regulators];

  vector<lower = 0>[num_targets] model_mismatch_variance;

  vector<lower = 0>[num_targets] initial_conditions;
  vector<lower = 0>[num_targets] basal_transcription;
  vector<lower = 0>[num_targets] degradation;
  vector<lower = 0>[num_targets] transcription_sensitivity;
  vector[num_targets] interaction_bias;
  
  real interaction_weights[num_regulators, num_targets];

  vector<lower = 0>[num_regulators] protein_initial_level;
  vector<lower = 0>[num_regulators] protein_degradation;
  
  vector<lower=0>[num_regulators] gp_variance;
  vector<lower=0>[num_regulators] gp_length;

//  matrix<lower=0>[num_time * num_intermediate_points, num_regulators] protein_level;
//  matrix<lower=0>[num_time * num_intermediate_points, num_regulators] protein_level_sigma;
}

model {
  vector[num_detailed_time] regulator_profiles_true[num_regulators];
  vector[num_detailed_time] protein_profiles[num_regulators];
  vector[num_time] target_profiles_true[num_targets];
  

  //----------First compute transformed data------------
  
  //Transformation to ensure non-negativity
  for(regulator in 1:num_regulators)
  {
    regulator_profiles_true[regulator] = log(1.0 + regulator_profiles_true_gp[regulator]);
  }
  
  //TF protein synthesis
  for (regulator in 1:num_regulators)
  {
    for (time in 1:num_detailed_time) {
      real initial = protein_initial_level[regulator] * exp(-protein_degradation[regulator] * time * integration_step);
      
      real accumulator = 0;
      for(previous in 1:time){
        real decay_exponent = -protein_degradation[regulator] * (time - previous) * integration_step;
        accumulator = accumulator + regulator_profiles_true[regulator, previous] * integration_step * exp(decay_exponent);
      }
      
      protein_profiles[regulator, time] = initial + accumulator;
    }
  }
  
  //Target RNA synthesis
  for (target in 1:num_targets) {
    for (time in 1:num_time) {
      real basal_over_degradation = basal_transcription[target] / degradation[target];
      real initial = basal_over_degradation + (initial_conditions[target] - basal_over_degradation) * exp(-degradation[target] * time);
      
      real synthesis_accumulator = 0; 
      int detailed_time_index = (time - 1) * num_integration_points + 1;
      
      for (previous in 1:detailed_time_index) {
        real sigmoid_value;
        real decay_exponent = -degradation[target] * (detailed_time_index - previous) * integration_step;
        
        real regulation_input = interaction_bias[target];
        for (regulator in 1:num_regulators) 
        {
          regulation_input = regulation_input + interaction_matrix[regulator, target] * interaction_weights[regulator, target] * log(protein_profiles[regulator, previous]);
        }
        
        sigmoid_value = integration_step / (1.0 + exp(-regulation_input));
        synthesis_accumulator = synthesis_accumulator + sigmoid_value * exp(decay_exponent);
      }
      target_profiles_true[target, time] = initial + transcription_sensitivity[target] * synthesis_accumulator;
    }
  }  
  
  
  //-----------Now the actual model-------------------------------------
  
  //Observation model
  for (regulator in 1:num_regulators) {
    for (time in 1:num_time) {
      int true_time_index = (time - 1) * num_integration_points + 1;
      regulator_profiles_observed[regulator,time] 
        ~ normal(regulator_profiles_true[regulator, true_time_index], regulator_profiles_sigma[regulator, time]);
    }
  }
  
  for (target in 1:num_targets) {
    target_profiles_observed[target] ~ normal(target_profiles_true[target], target_profiles_sigma[target] + model_mismatch_variance[target]);
  }
  
  //GP prior on regulators
  for (regulator in 1: num_regulators) {
    matrix[num_detailed_time, num_detailed_time] regulator_gp_covariance;

    //Calculate covariance matrix    
    for(detailed_time_index in 1:num_detailed_time) {
      regulator_gp_covariance[detailed_time_index, detailed_time_index] = gp_variance[regulator];
      for(other_time_index in (detailed_time_index + 1):num_detailed_time) {
        real value = gp_variance[regulator] * exp(-0.5 * (1 / gp_length[regulator]) * square((detailed_time_index - other_time_index) * integration_step));
        regulator_gp_covariance[detailed_time_index, other_time_index] = value;
        regulator_gp_covariance[other_time_index, detailed_time_index] = value;
      }
    }
    
    regulator_profiles_true_gp[regulator] ~ multi_normal(zero_mean, regulator_gp_covariance);
  }
  
  //Other priors
  transcription_params_prior_lp(initial_conditions, num_targets);
  transcription_params_prior_lp(basal_transcription, num_targets);
  transcription_params_prior_lp(degradation, num_targets);
  transcription_params_prior_lp(transcription_sensitivity, num_targets);
  
  transcription_params_prior_lp(protein_degradation, num_regulators);
  
  interaction_bias ~ normal(0,2);

  //The following differs from the paper for simplicity (originally a conjugate inverse gamma)
  model_mismatch_variance ~ cauchy(0, 2); 

  for(regulator in 1:num_regulators) {
    interaction_weights[regulator] ~ normal(0,2);
  }
  
  gp_length ~ uniform(1, 144);
  
  gp_variance ~ cauchy(0, 2); //in the paper this prior is not specified
}
