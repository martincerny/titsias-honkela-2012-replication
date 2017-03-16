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
  vector<lower = 0>[num_time] regulator_profiles_observed[num_replicates, num_regulators];
  vector<lower = 0>[num_time] regulator_profiles_sigma[num_replicates, num_regulators];
  
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
  vector[num_detailed_time] regulator_profiles_true_gp[num_replicates, num_regulators];

  //vector<lower = 0>[num_genes] model_mismatch_sigma;

  vector<lower = 0>[num_genes] initial_condition[num_replicates];
  vector<lower = 0>[num_genes] basal_transcription;
  vector<lower = 0>[num_genes] degradation;
  vector<lower = 0>[num_genes] transcription_sensitivity;
  vector[num_genes] interaction_bias;
  
  vector[num_regulators] interaction_weights[num_genes];

  vector<lower = 0>[num_regulators] protein_initial_level[num_replicates];
  vector<lower = 0>[num_regulators] protein_degradation;
  
  vector<lower=0>[num_regulators] gp_variance;
  vector<lower=0>[num_regulators] gp_length;
}

transformed parameters {
  vector[num_detailed_time] regulator_profiles_true[num_replicates, num_regulators];
  matrix[num_detailed_time, num_regulators] log_tf_profiles[num_replicates];
  vector[num_time] gene_profiles_true[num_replicates, num_genes];


  //GP prior on regulators
  for (regulator in 1: num_regulators) {
    matrix[num_detailed_time, num_detailed_time] covariance = cov_exp_quad(detailed_time, gp_variance[regulator] + 0.00001, gp_length[regulator]);
    matrix[num_detailed_time, num_detailed_time] covariance_cholesky;
    
    for (I in 1:num_detailed_time)
    {
      covariance[I, I] = covariance[I, I] + 1e-12;
    }
    
    covariance_cholesky = cholesky_decompose(covariance);
    
    for (replicate in 1: num_replicates) {
      regulator_profiles_true[replicate, regulator] = log1p_exp(covariance_cholesky * regulator_profiles_true_gp[replicate, regulator]);
    } 
  }
  
  //TF protein synthesis
  for(replicate in 1:num_replicates) {
    for(regulator in 1:num_regulators){
      real degradation_per_integration_step = exp(-protein_degradation[regulator] * integration_step);
      real residual = -0.5 * integration_step * regulator_profiles_true[replicate,regulator,1];
      
      log_tf_profiles[replicate,1, regulator] = log(protein_initial_level[replicate, regulator]);
      for (time in 2:num_detailed_time) {
        real initial_residual = protein_initial_level[replicate, regulator] * exp(-protein_degradation[regulator] * (time - 1) * integration_step);
        
        residual = (residual + integration_step * regulator_profiles_true[replicate,regulator,time - 1]) * degradation_per_integration_step;
        log_tf_profiles[replicate, time, regulator] = log(initial_residual + residual + 0.5 * integration_step * regulator_profiles_true[replicate, regulator, time]);
      }
    }
  }  
  
  

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
  
  
  //-----------Now the actual model-------------------------------------
  
  //Observation model
  for (replicate in 1:num_replicates) {
    for (regulator in 1:num_regulators) {
      for (time in 1:num_time) {
        int true_time_index = (time - 1) * num_integration_points + 1;
        regulator_profiles_observed[replicate, regulator,time] 
          ~ normal(regulator_profiles_true[replicate, regulator, true_time_index], regulator_profiles_sigma[replicate, regulator, time]);
      }
    }
  }
  
  for (replicate in 1:num_replicates) {
    for (gene in 1:num_genes) {
      gene_profiles_observed[replicate, gene] ~ normal(gene_profiles_true[replicate, gene], gene_profiles_sigma[replicate, gene]);// + model_mismatch_sigma[gene]);
    }
  }
  
  //GP prior on regulators
  for (replicate in 1:num_replicates) {
    for (regulator in 1: num_regulators) {
      regulator_profiles_true_gp[replicate, regulator] ~ normal(0, 1);
    }
  }
  
  //in the paper, the prior is uniform, between the smallest difference a and squared length
  gp_length ~ cauchy(5, 5);

  gp_variance ~ cauchy(0, 2); //in the paper this prior is not specified
  
  
  //Other priors
  for (replicate in 1:num_replicates) {
    transcription_params_prior_lp(initial_condition[replicate]);
    transcription_params_prior_lp(protein_initial_level[replicate]);
  }
  
  transcription_params_prior_lp(basal_transcription);
  transcription_params_prior_lp(degradation);
  transcription_params_prior_lp(transcription_sensitivity);
  
  transcription_params_prior_lp(protein_degradation);
  
  interaction_bias ~ normal(0,2);

  //The following differs from the paper for simplicity (originally a conjugate inverse gamma)
  //model_mismatch_sigma ~ cauchy(0, 2); 

  for(regulator in 1:num_regulators) {
    interaction_weights[regulator] ~ normal(0,2);
  }

}
