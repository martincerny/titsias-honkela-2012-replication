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

  vector<lower = 0>[num_genes] initial_conditions[num_replicates];
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
  vector[num_regulators] protein_profiles[num_replicates, num_detailed_time];

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
    for (regulator in 1:num_regulators)
    {
      real synthesis_accumulator = protein_initial_level[replicate, regulator];
      real degradation_per_integration_step = exp(-protein_degradation[regulator] * integration_step);

      for (time in 1:num_detailed_time) {
        protein_profiles[replicate, time, regulator] = synthesis_accumulator;
        
        synthesis_accumulator = synthesis_accumulator * degradation_per_integration_step + regulator_profiles_true[replicate, regulator, time] * integration_step;      
        }
    }
  }  
}

model {
  
  vector[num_time] gene_profiles_true[num_replicates, num_genes];
  

  //----------First compute transformed data------------
  


  
  //gene RNA synthesis
  for (gene in 1:num_genes) {
    real degradation_per_integration_step = exp(-degradation[gene] * integration_step);
    real basal_over_degradation = basal_transcription[gene] / degradation[gene];
    
    for(replicate in 1:num_replicates) {
      real previously_synthetized_residual = 0;
      real initial = initial_conditions[replicate, gene] - basal_over_degradation;
      
      gene_profiles_true[replicate, gene, 1] = initial_conditions[replicate, gene];
      for (detailed_time_index in 2:num_detailed_time) {
        //TODO: replaced midpoint with trapezoid rule
          
        //TODO maybe put log(tf_profiles) in transformed data
        real regulation_input = interaction_bias[gene] + dot_product(interaction_weights[gene], log(protein_profiles[replicate,detailed_time_index]));
  
        real sigmoid_value = integration_step / (1.0 + exp(-regulation_input));

        previously_synthetized_residual = previously_synthetized_residual * degradation_per_integration_step + sigmoid_value;
        if ((detailed_time_index - 1) % num_integration_points == 0)
        {
          int time = ((detailed_time_index - 1) / num_integration_points) + 1;
          real initial_residual = initial * exp(-degradation[gene] * (time - 1));
          gene_profiles_true[replicate, gene, time] = basal_over_degradation + initial_residual + transcription_sensitivity[gene] * previously_synthetized_residual;
        }
      }
    }
  }  
  
  
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
    transcription_params_prior_lp(initial_conditions[replicate]);
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
