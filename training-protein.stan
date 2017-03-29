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
  
  int num_regulators; 

  int num_genes;
  vector<lower = 0>[num_time] gene_profiles_observed[num_replicates, num_genes];
  vector<lower = 0>[num_time] gene_profiles_sigma[num_replicates, num_genes];
  
  //0 regulation impossible, 1 regulation possible
  int interaction_matrix[num_regulators, num_genes];
}

transformed data {
  int num_detailed_time = count_num_detailed_time(num_time, num_integration_points);
  real integration_step = 1.0 / num_integration_points;
  vector[num_detailed_time] zero_mean;
  real detailed_time[num_detailed_time];
  
  //print("RegO:", regulator_profiles_observed);
  //print("RegSig:", regulator_profiles_sigma);
  
  for (i in 1:num_detailed_time)  
  {
    detailed_time[i] = (i - 1) * integration_step + 1;
    zero_mean[i] = 0;
  }
}

parameters {
  vector[num_detailed_time] protein_profiles_true_gp[num_replicates, num_regulators];
  
  //vector<lower = 0>[num_genes] model_mismatch_sigma;
  
  vector<lower = 0>[num_genes] initial_condition[num_replicates];
  vector<lower = 0>[num_genes] basal_transcription;
  vector<lower = 0>[num_genes] degradation;
  vector<lower = 0>[num_genes] transcription_sensitivity;
  vector[num_genes] interaction_bias;
  
  vector[num_regulators] interaction_weights[num_genes];
  
  vector<lower=1, upper=num_time>[num_regulators] gp_length;
}

transformed parameters {
  matrix[num_detailed_time, num_regulators] log_tf_profiles[num_replicates];
  vector[num_time] gene_profiles_true[num_replicates, num_genes];
  
  
  //GP prior on protein levels
  for (regulator in 1: num_regulators) {
    matrix[num_detailed_time, num_detailed_time] covariance = cov_exp_quad(detailed_time, 1, gp_length[regulator]);
    matrix[num_detailed_time, num_detailed_time] covariance_cholesky;
    
    for (I in 1:num_detailed_time)
    {
      covariance[I, I] = covariance[I, I] + 1e-12;
    }
    
    covariance_cholesky = cholesky_decompose(covariance);
    
    for (replicate in 1: num_replicates) {
      log_tf_profiles[replicate, , regulator] = log1p_exp(covariance_cholesky * protein_profiles_true_gp[replicate, regulator]);
    } 
  }
  
  //gene RNA synthesis
  for (gene in 1:num_genes) {
    real basal_over_degradation = basal_transcription[gene] / degradation[gene];
    real degradation_per_step = exp(-degradation[gene] * integration_step);
    
    for (replicate in 1:num_replicates)
    {
      real residual;
      vector[num_detailed_time] regulation_input;
      vector[num_detailed_time] synthesis;
      
      regulation_input = rep_vector(interaction_bias[gene], num_detailed_time) + log_tf_profiles[replicate] * interaction_weights[gene];
      
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
  
  
  for (replicate in 1:num_replicates) {
    for (gene in 1:num_genes) {
      gene_profiles_observed[replicate, gene] ~ normal(gene_profiles_true[replicate, gene], gene_profiles_sigma[replicate, gene]);// + model_mismatch_sigma[gene]);
    }
  }
  
  //GP prior on regulators
  for (replicate in 1:num_replicates) {
    for (regulator in 1: num_regulators) {
      protein_profiles_true_gp[replicate, regulator] ~ normal(0, 1);
    }
  }
  
  //Other priors
  for (replicate in 1:num_replicates) {
    transcription_params_prior_lp(initial_condition[replicate]);
  }
  
  transcription_params_prior_lp(basal_transcription);
  transcription_params_prior_lp(degradation);
  transcription_params_prior_lp(transcription_sensitivity);
  
  interaction_bias ~ normal(0,2);
  
  //The following differs from the paper for simplicity (originally a conjugate inverse gamma)
  //model_mismatch_sigma ~ cauchy(0, 2); 
  
  for(gene in 1:num_genes) {
    interaction_weights[gene] ~ normal(0,2);
  }
  
}
