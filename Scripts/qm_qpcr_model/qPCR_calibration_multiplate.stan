data {
  
  int Nplates;
  int Nobs;
  int NSamples; //number of unique biol samples, overall
  int NstdSamples; //number of unique biol samples with known concentrations (standards)
  int plate_idx[Nobs];
  int std_idx[NstdSamples]; //index relative to NSamples; which ones are the standards?
  int unkn_idx[NSamples-NstdSamples];
  int plateSample_idx[Nobs]; //index of unique combinations of plate and biological sample
  
  
  vector[Nobs] y; //Ct observations
  vector[NstdSamples] known_concentration;
  
  real stdCurvePrior_intercept[2];
  real stdCurvePrior_slope[2];
  
  //vector[NSamples-NstdSamples] dilutionFactor;
  
}

transformed data {
  vector[NstdSamples] known_conc;
    known_conc = log10(known_concentration);
}

parameters {
  vector<lower=0>[Nplates] beta_std_curve_0;
  vector<upper=0>[Nplates] beta_std_curve_1;
  vector<lower=0>[Nplates]  gamma_0; //intercept to scale variance w the mean
  vector<lower=0>[Nplates]  gamma_1; //slope to scale variance w the mean
  vector[NSamples-NstdSamples] envir_concentration;
  
}

transformed parameters{
  vector[NSamples] Concentration;
  vector[NSamples] mu;
  vector<lower=0>[NSamples] sigma;

  Concentration[std_idx] = known_conc; //log10 scale
  Concentration[unkn_idx] = envir_concentration;
  
  for(i in 1:Nobs){
    mu[plateSample_idx[i]] = beta_std_curve_0[plate_idx[i]] + 
                              beta_std_curve_1[plate_idx[i]] * Concentration[plateSample_idx[i]];
                              
    sigma[plateSample_idx[i]] = exp(gamma_0[plate_idx[i]] + 
                                  gamma_1[plate_idx[i]]*mu[plateSample_idx[i]]);
  }
}


model {
   for(i in 1:Nobs){
      y[i] ~ normal(mu[plateSample_idx[i]], sigma[plateSample_idx[i]]);
    }

  beta_std_curve_0 ~ normal(stdCurvePrior_intercept[1], stdCurvePrior_intercept[2]);
  beta_std_curve_1 ~ normal(stdCurvePrior_slope[1], stdCurvePrior_slope[2]);
  envir_concentration ~ normal(0, 5); //log10 scale
  gamma_0 ~ normal(0,2);
  gamma_1 ~ normal(0,.1);
}

