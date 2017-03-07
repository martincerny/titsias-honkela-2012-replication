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
  /*
  for(i in 1:num_replicates) {
      print("Rep",i,": ",gene_profile_observed[i]);
  }*/
}

parameters {

  real<lower = 0> initial_condition[num_replicates];
  real<lower = 0> basal_transcription;
  real<lower = 0> degradation;
  real<lower = 0> transcription_sensitivity;
  real<lower = 0> interaction_bias;
  
  //real<lower = 0> model_mismatch_sigma;
  
  real interaction_weights[num_regulators];
}

transformed parameters {
  vector[num_time] gene_profile_true[num_replicates];
  
  //gene RNA synthesis
  { //new scope to make the variables local
    real degradation_per_unit_time = exp(-degradation);
    real basal_over_degradation = basal_transcription / degradation;

//!Init -> time 1!
    for (replicate in 1:num_replicates)
    {
      real previously_synthetized_residual = 0;
      real initial_residual = (initial_condition[replicate] - basal_over_degradation);
      
      gene_profile_true[replicate, 1] = initial_condition[replicate];
      for (time in 2:num_time) {
      
        real initial = basal_over_degradation + (initial_condition[replicate] - basal_over_degradation) * exp(-degradation * time);
        
        real synthesis_accumulator = 0; 
  
        int previous_detailed_time = (time - 2) * num_integration_points + 1;
        
        for (mid_step in 1:num_integration_points) {
          real sigmoid_value;
          real decay_exponent = -degradation * (num_integration_points - mid_step) * integration_step;
          
          //TODO maybe put log(tf_profiles) in transformed data
          real regulation_input = interaction_bias + dot_product(interaction_weights, log(tf_profiles[replicate,,previous_detailed_time + mid_step]));
  
          sigmoid_value = integration_step / (1.0 + exp(-regulation_input));
          synthesis_accumulator = synthesis_accumulator + sigmoid_value * exp(decay_exponent);
        }
        
        initial_residual = initial_residual * degradation_per_unit_time;
        previously_synthetized_residual = previously_synthetized_residual * degradation_per_unit_time + synthesis_accumulator; 
        
        gene_profile_true[replicate, time] = basal_over_degradation + initial_residual + transcription_sensitivity * previously_synthetized_residual;
        
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
