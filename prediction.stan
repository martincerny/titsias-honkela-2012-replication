functions {
  void transcription_params_prior_lp(real v)
  {
    real log_v = log(v);
    target += normal_lpdf(log_v | -0.5,10);
    target += -log_v;
  }
}

//TODO: replicates
data {
  int num_time;
  
  int num_integration_points;
  
  int num_regulators; 
  vector<lower = 0>[num_time * num_integration_points] tf_profiles[num_regulators];

  vector<lower = 0>[num_time] gene_profile_observed;
  vector<lower = 0>[num_time] gene_profile_sigma;
  
}

transformed data {
  int num_detailed_time = num_time * num_integration_points;
  real integration_step = 1.0 / num_integration_points;
  vector[num_detailed_time] zero_mean;
  real detailed_time[num_detailed_time];
  
  //print("RegO:", regulator_profiles_observed);
  //print("RegSig:", regulator_profiles_sigma);
  
  for (i in 1:num_detailed_time)  
  {
    detailed_time[i] = (i - 1) * integration_step;
    zero_mean[i] = 0;
  }
}

parameters {

  real<lower = 0> initial_condition;
  real<lower = 0> basal_transcription;
  real<lower = 0> degradation;
  real<lower = 0> transcription_sensitivity;
  real<lower = 0> interaction_bias;
  
  //real<lower = 0> model_mismatch_sigma;
  
  real interaction_weights[num_regulators];
}

transformed parameters {
  vector[num_time] gene_profile_true;
  
  //gene RNA synthesis
  { //new scope to make the variables local
    real degradation_per_unit_time = exp(-degradation);
    real basal_over_degradation = basal_transcription / degradation;
    real initial_residual = (initial_condition - basal_over_degradation);
    real previously_synthetized_residual = 0;

//!Init -> time 1!
    gene_profile_true[1] = initial_condition;
    for (time in 2:num_time) {
    
      real initial = basal_over_degradation + (initial_condition - basal_over_degradation) * exp(-degradation * time);
      
      real synthesis_accumulator = 0; 

      int previous_detailed_time = (time - 2) * num_integration_points + 1;
      
      for (mid_step in 1:num_integration_points) {
        real sigmoid_value;
        real decay_exponent = -degradation * (num_integration_points - mid_step) * integration_step;
        
        //TODO maybe put log(tf_profiles) in transformed data
        real regulation_input = interaction_bias + dot_product(interaction_weights, log(tf_profiles[,previous_detailed_time + mid_step]));

        sigmoid_value = integration_step / (1.0 + exp(-regulation_input));
        synthesis_accumulator = synthesis_accumulator + sigmoid_value * exp(decay_exponent);
      }
      
      initial_residual = initial_residual * degradation_per_unit_time;
      previously_synthetized_residual = previously_synthetized_residual * degradation_per_unit_time + synthesis_accumulator; 
      
      gene_profile_true[time] = basal_over_degradation + initial_residual + transcription_sensitivity * previously_synthetized_residual;
      
      /*
        if(is_nan(gene_profiles_true[gene, time]) || is_inf(gene_profiles_true[gene, time]))
        {
          print("gene ", gene, " weird: ", gene_profiles_true[gene, time], " init: ", initial_conditions[gene], " sensitivit: ", transcription_sensitivity[gene], " decay: ", degradation[gene], " weights: ", interaction_weights);
          gene_profiles_true[gene, time] = 100;
        }
      */
    }
  }
}

model {
  

  
  //-----------Now the actual model-------------------------------------
    
  gene_profile_observed ~ normal(gene_profile_true, gene_profile_sigma); //+ model_mismatch_sigma);


  //Other priors
  transcription_params_prior_lp(initial_condition);
  transcription_params_prior_lp(basal_transcription);
  transcription_params_prior_lp(degradation);
  transcription_params_prior_lp(transcription_sensitivity);
  
  interaction_bias ~ normal(0,2);
  
  //The following differs from the paper for simplicity (originally a conjugate inverse gamma)
  //model_mismatch_sigma ~ cauchy(0, 2); 
  
  interaction_weights ~ normal(0,2);

}
