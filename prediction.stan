functions {
  void transcription_params_prior_lp(real v)
  {
    real log_v = log(v);
    target += normal_lpdf(log_v | -0.5, sqrt2());
    target += -log_v;
  }
  
  int count_num_detailed_time(int num_time, int num_integration_points)
  {
    return (num_time - 1) * num_integration_points + 1;
  }
}

//TODO: replicates
data {
  int num_time;
  int num_replicates;
  
  int num_integration_points;
  
  int num_regulators; 
  vector<lower = 0>[count_num_detailed_time(num_time, num_integration_points)] tf_profiles[num_replicates, num_regulators];

  vector<lower = 0>[num_time] gene_profile_observed[num_replicates];
  vector<lower = 0>[num_time] gene_profile_sigma[num_replicates];
  
}

transformed data {
  int num_detailed_time = count_num_detailed_time(num_time, num_integration_points);
  real integration_step = 1.0 / num_integration_points;
  matrix[num_detailed_time, num_regulators] log_tf_profiles[num_replicates];
  
  for(replicate in 1:num_replicates) {
     for(regulator in 1:num_regulators){
       for(detailed_time in 1:num_detailed_time) {
        log_tf_profiles[replicate, detailed_time, regulator] =log(tf_profiles[replicate, regulator, detailed_time]);
       }
     }
  }
}

parameters {

  real<lower = 0> initial_condition[num_replicates];
  real<lower = 0> basal_transcription;
  real<lower = 0> degradation;
  real<lower = 0> transcription_sensitivity;
  real<lower = 0> interaction_bias;
  
  //real<lower = 0> model_mismatch_sigma;
  
  vector[num_regulators] interaction_weights;
}

transformed parameters {
  vector[num_time] gene_profile_true[num_replicates];
  
  //gene RNA synthesis
  { //new scope to make the variables local
    //real degradation_per_unit_time = exp(-degradation);
    real basal_over_degradation = basal_transcription / degradation;
    //real degradation_per_step = exp(-degradation * integration_step);

//!Init -> time 1!
    for (replicate in 1:num_replicates)
    {
      //real residual;
      vector[num_detailed_time] regulation_input;
      vector[num_detailed_time] synthesis;

      regulation_input = rep_vector(interaction_bias, num_detailed_time) + log_tf_profiles[replicate] * interaction_weights;

      for (i in 1:num_detailed_time)
      {
        synthesis[i] = integration_step / (1 + exp(-regulation_input[i]));
      }
      //synthesis = integration_step * inv(1 + exp(-regulation_input));
      
      gene_profile_true[replicate, 1] = initial_condition[replicate];

      for (time in 2:num_time)
      {
        real integral = 0;
        int detailed_time = (time - 1) * num_integration_points + 1;
        for(previous_detailed in 1:detailed_time){
          real degradation_coeff = exp(-degradation * (detailed_time - previous_detailed) * integration_step);
          real trapezoid_synthesis;
          if(previous_detailed == 1 || previous_detailed == detailed_time){
            trapezoid_synthesis = 0.5 * synthesis[previous_detailed];
          }
          else
          {
            trapezoid_synthesis = synthesis[previous_detailed];
          }
          integral = integral + trapezoid_synthesis * degradation_coeff;
          
        }
        gene_profile_true[replicate, time] = basal_over_degradation + (initial_condition[replicate] - basal_over_degradation) * exp(-degradation * (time - 1)) + transcription_sensitivity * integral;
        
      }
      /*
      //Calculating the integral by trapezoid rule in a single pass for all values
      residual = -0.5 * synthesis[1];
      for (detailed_time in 2:num_detailed_time)
      { 
        residual = (residual + synthesis[detailed_time - 1]) * degradation_per_step;
        
        if((detailed_time - 1) % num_integration_points == 0) {
          real integral_value = residual + 0.5 * synthesis[detailed_time];
          int time = ((detailed_time -1) / num_integration_points) + 1;
          gene_profile_true[replicate, time] = basal_over_degradation + (initial_condition[replicate] - basal_over_degradation) * exp(-degradation * (time - 1)) + transcription_sensitivity * integral_value;
        }
      }*/
    }
  }
}

model {
  

  
  //-----------Now the actual model-------------------------------------
  for(replicate in 1:num_replicates)
  {
    gene_profile_observed[replicate] ~ normal(gene_profile_true[replicate], gene_profile_sigma[replicate]); //+ model_mismatch_sigma);
    
    //Here we diverge from the original model. At least for simulated data, no prior on initial_condition works better
    //transcription_params_prior_lp(initial_condition[replicate]);
  }


  //Other priors
  transcription_params_prior_lp(basal_transcription);
  transcription_params_prior_lp(degradation);
  transcription_params_prior_lp(transcription_sensitivity);
  
  interaction_bias ~ normal(0,2);
  
  //The following differs from the paper for simplicity (originally a conjugate inverse gamma)
  //model_mismatch_sigma ~ cauchy(0, 2); 
  
  interaction_weights ~ normal(0,2);

}
