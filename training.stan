functions {
  void transcription_params_prior_lp(vector v)
  {
    vector[num_elements(v)] log_v = log(v);
    target += normal_lpdf(log_v | -0.5,2);
    target += -log_v;
  }
}

//TODO: replicates
data {
  int num_time;
  
  int num_integration_points;
  
  int num_regulators; 
  vector<lower = 0>[num_time] regulator_profiles_observed[num_regulators];
  vector<lower = 0>[num_time] regulator_profiles_sigma[num_regulators];
  
  int num_genes;
  vector<lower = 0>[num_time] gene_profiles_observed[num_genes];
  vector<lower = 0>[num_time] gene_profiles_sigma[num_genes];
  
  //0 regulation impossible, 1 regulation possible
  int interaction_matrix[num_regulators, num_genes];
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
  vector[num_detailed_time] regulator_profiles_true_gp[num_regulators];

  //vector<lower = 0>[num_genes] model_mismatch_sigma;

  vector<lower = 0>[num_genes] initial_conditions;
  vector<lower = 0>[num_genes] basal_transcription;
  vector<lower = 0>[num_genes] degradation;
  vector<lower = 0>[num_genes] transcription_sensitivity;
  vector[num_genes] interaction_bias;
  
  real interaction_weights[num_regulators, num_genes];

  vector<lower = 0>[num_regulators] protein_initial_level;
  vector<lower = 0>[num_regulators] protein_degradation;
  
  vector<lower=0>[num_regulators] gp_variance;
  vector<lower=1>[num_regulators] gp_length;

//  matrix<lower=0>[num_time * num_intermediate_points, num_regulators] protein_level;
//  matrix<lower=0>[num_time * num_intermediate_points, num_regulators] protein_level_sigma;
}

transformed parameters {
  vector[num_detailed_time] regulator_profiles_true[num_regulators];
 
  //GP prior on regulators
  for (regulator in 1: num_regulators) {
    matrix[num_detailed_time, num_detailed_time] covariance = cov_exp_quad(detailed_time, gp_variance[regulator], gp_length[regulator]);
    
    matrix[num_detailed_time, num_detailed_time] covariance_cholesky;
    
    for (I in 1:num_detailed_time)
    {
      covariance[I, I] = covariance[I, I] + 1e-12;
    }
    
    covariance_cholesky = cholesky_decompose(covariance);
    
    regulator_profiles_true[regulator] = log1p_exp(covariance_cholesky * regulator_profiles_true_gp[regulator]);
  } 
}

model {
  
  vector[num_detailed_time] protein_profiles[num_regulators];
  vector[num_time] gene_profiles_true[num_genes];
  

  //----------First compute transformed data------------
  

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
  
  //gene RNA synthesis
  for (gene in 1:num_genes) {
    for (time in 1:num_time) {
      real basal_over_degradation = basal_transcription[gene] / degradation[gene];
      real initial = basal_over_degradation + (initial_conditions[gene] - basal_over_degradation) * exp(-degradation[gene] * time);
      
      real synthesis_accumulator = 0; 
      int detailed_time_index = (time - 1) * num_integration_points + 1;
      
      for (previous in 1:detailed_time_index) {
        real sigmoid_value;
        real decay_exponent = -degradation[gene] * (detailed_time_index - previous) * integration_step;
        
        real regulation_input = interaction_bias[gene];
        for (regulator in 1:num_regulators) 
        {
          regulation_input = regulation_input + interaction_matrix[regulator, gene] * interaction_weights[regulator, gene] * log(protein_profiles[regulator, previous]);
        }
        
        sigmoid_value = integration_step / (1.0 + exp(-regulation_input));
        synthesis_accumulator = synthesis_accumulator + sigmoid_value * exp(decay_exponent);
      }
      gene_profiles_true[gene, time] = initial + transcription_sensitivity[gene] * synthesis_accumulator;
      
	  /*
      if(is_nan(gene_profiles_true[gene, time]) || is_inf(gene_profiles_true[gene, time]))
      {
        print("gene ", gene, " weird: ", gene_profiles_true[gene, time], " init: ", initial_conditions[gene], " sensitivit: ", transcription_sensitivity[gene], " decay: ", degradation[gene], " weights: ", interaction_weights);
        gene_profiles_true[gene, time] = 100;
      }
	  */
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
  
  for (gene in 1:num_genes) {
    gene_profiles_observed[gene] ~ normal(gene_profiles_true[gene], gene_profiles_sigma[gene]);// + model_mismatch_sigma[gene]);
  }
  
  //GP prior on regulators
  for (regulator in 1: num_regulators) {
    regulator_profiles_true_gp[regulator] ~ normal(0, 1);
  }
  
  //in the paper, the prior is uniform, between the smallest difference a and squared length
  gp_length ~ cauchy(1, 5);

  gp_variance ~ cauchy(0, 2); //in the paper this prior is not specified
  
  
  //Other priors
  transcription_params_prior_lp(initial_conditions);
  transcription_params_prior_lp(basal_transcription);
  transcription_params_prior_lp(degradation);
  transcription_params_prior_lp(transcription_sensitivity);
  
  transcription_params_prior_lp(protein_degradation);
  transcription_params_prior_lp(protein_initial_level);
  
  interaction_bias ~ normal(0,2);

  //The following differs from the paper for simplicity (originally a conjugate inverse gamma)
  //model_mismatch_sigma ~ cauchy(0, 2); 

  for(regulator in 1:num_regulators) {
    interaction_weights[regulator] ~ normal(0,2);
  }

}
