functions {
  void transcription_params_prior_lp(vector v)
  {
    vector[num_elements(v)] log_v = log(v); 
    target += normal_lpdf(log_v | -0.5,sqrt2());     
    target += -log_v;
  }
  
  int count_num_detailed_time(int num_time, int num_integration_points)
  {
    return (num_time - 1) * num_integration_points + 1;
  }
}

data {
  int num_time;
  int num_replicates;
  
  int num_integration_points;
  
  vector[count_num_detailed_time(num_time, num_integration_points)] log_tf_profiles_source[num_replicates];
  
  int num_genes;
  vector<lower = 0>[num_time] gene_profiles_observed[num_replicates, num_genes];
  vector<lower = 0>[num_time] gene_profiles_sigma[num_replicates, num_genes];
}

transformed data {
  int num_detailed_time = count_num_detailed_time(num_time, num_integration_points);
  real integration_step = 1.0 / num_integration_points;
  vector[num_detailed_time] log_tf_profiles[num_replicates];
  
  for(replicate in 1:num_replicates) {
    for(time in 1:num_detailed_time) {
        log_tf_profiles[replicate, time] = log_tf_profiles_source[replicate, time];
    }
  }
  
}

parameters {
  vector<lower = 0>[num_genes] initial_condition[num_replicates];
  vector<lower = 0>[num_genes] basal_transcription;
  vector<lower = 0>[num_genes] degradation;
  vector<lower = 0>[num_genes] transcription_sensitivity;
  vector[num_genes] interaction_bias;
  
  real interaction_weights[num_genes];
}

transformed parameters {
  vector[num_time] gene_profiles_true[num_replicates, num_genes];



  //----------First compute transformed data------------
  


  
  //gene RNA synthesis
  for (gene in 1:num_genes) {
    real basal_over_degradation = basal_transcription[gene] / degradation[gene];
    real degradation_per_step = exp(-degradation[gene] * integration_step);

    for (replicate in 1:num_replicates)
    {
      real residual;
      vector[num_detailed_time] regulation_input;
      vector[num_detailed_time] synthesis;

      for(t in 1:num_detailed_time)
      {
        regulation_input[t] = interaction_bias[gene];
        regulation_input[t] = regulation_input[t] + log_tf_profiles[replicate, t] * interaction_weights[gene];
      }

      synthesis = integration_step * inv(1 + exp(-regulation_input));
      
      gene_profiles_true[replicate, gene, 1] = initial_condition[replicate, gene];

      //Calculating the integral by trapezoid rule in a single pass for all values
      residual = -0.5 * synthesis[1];
      for (time in 2:num_time) 
      {
        int detailed_time_base = (time - 2) * num_integration_points + 1; //detailed_time needs to start at 2 (so this starts at 1)
        for (detailed_time_relative in 1:num_integration_points)
        { 
          int detailed_time_index = detailed_time_base + detailed_time_relative;
          residual = (residual + synthesis[detailed_time_index - 1]) * degradation_per_step;
        }
        
        { //new block to let me define new vars
          int detailed_time_index = detailed_time_base + num_integration_points;
          real integral_value = residual + 0.5 * synthesis[detailed_time_index];
          gene_profiles_true[replicate, gene, time] = basal_over_degradation + (initial_condition[replicate, gene] - basal_over_degradation) * exp(-degradation[gene] * (time - 1)) + transcription_sensitivity[gene] * integral_value;
        }
        
      }
    }
  }  
  
    
}

model {
  
  
  //-----------Now the actual model-------------------------------------
  
  for (replicate in 1:num_replicates) {
    for (gene in 1:num_genes) {
      gene_profiles_observed[replicate, gene] ~ normal(gene_profiles_true[replicate, gene], gene_profiles_sigma[replicate, gene]);// + model_mismatch_sigma[gene]);
    }
  }
  

  transcription_params_prior_lp(basal_transcription);
  transcription_params_prior_lp(degradation);
  transcription_params_prior_lp(transcription_sensitivity);
  
  interaction_bias ~ normal(0,2);

  interaction_weights ~ normal(0,2);

}
