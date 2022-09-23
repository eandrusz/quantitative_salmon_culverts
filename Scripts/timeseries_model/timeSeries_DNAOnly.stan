
data {
  int<lower=0> N;
  int<lower=0> Ntime;
  int time_idx[N];
  // int y_fish[N];
  vector[N] y_logeDNA;
}

parameters {
  // real log_lambda_0;
  real mu_0; //starting concentration for eDNA (log scale)
  // real beta_0;
  // real beta_1;
  real phi[Ntime-1]; //difference between mean at time = t-1 and mean at time = t
  // real eta[Ntime];
  real<lower=0> sigma_dna; //SD of DNA concentration (log scale)
  real<lower=0> sigma_phi; //SD of values of phi across the whole dataset
}

transformed parameters {
  // real log_lambda[Ntime];
  real mu[Ntime]; //mean eDNA (log scale) from which obs was drawn
  
  // log_lambda[1] = log_lambda_0;
   mu[1] = mu_0;
  for (t in 2:Ntime){
    // log_lambda[t] = log_lambda[t-1] + phi[t-1];
    mu[t] = mu[t-1] + phi[t-1];
  }
  
}


model {

  for (i in 1:N){
    // y_fish[i] ~ poisson(exp(log_lambda[time_idx[i]]));
    y_logeDNA[i] ~ normal(mu[time_idx[i]], sigma_dna);
  }
  
  phi ~ normal(0, sigma_phi);
  // log_lambda_0 ~ normal(0, sigma_fish);
  mu_0 ~ normal(0, 10);

  sigma_phi ~ gamma(1,2);
  sigma_dna ~ gamma(1,2);
  // eta ~ normal(0, 1);
  
  // beta_0 ~ normal(0, 5);
  // beta_1 ~ normal(0, 1);
  
}

