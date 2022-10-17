
data {
  int<lower=0> Ncreek;
  int<lower=0> Ntime;
  int<lower=0> Nobs; //total observations
  int<lower=0> Nstations; //probably 2, downstream and upstream
  int<lower=0> Nspecies; 
  int<lower=0> time_idx[Nobs];
  int<lower=0> creek_idx[Nobs];
  int<lower=0> species_idx[Nobs];
  int<lower=0> station_idx[Nobs]; //downstream = 1, upstream = 2
  vector[Nobs] y_logeDNA;
  
  int<lower=0>N_unobserved;
  int<lower=0>unobserved_time_idx[N_unobserved];
  int<lower=0>unobserved_creek_idx[N_unobserved];
  int<lower=0>unobserved_station_idx[N_unobserved];
  int<lower=0>unobserved_species_idx[N_unobserved];

}

parameters {
  
  vector[Ncreek] mu_0[Nstations, Nspecies]; //starting concentration for eDNA (log scale)
  vector[Ncreek] eta[Ntime-1, Nstations, Nspecies]; //
  vector<lower=-1, upper = 1>[Ntime-1] beta_1; //slope relating present concentration to past concentration
  // vector<lower=0>[Ncreek] sigma_dna; //SD of DNA concentration (log scale) within creek for a given time
  // vector<lower=0>[Ncreek] sigma_eta; //SD of values of eta within creek across time
   vector<lower=0>[Nspecies] sigma_dna; //SD of DNA concentration (log scale)
  // real<lower=0> sigma_dna; //SD of DNA concentration (log scale)
  // vector<lower=0>[Nspecies] sigma_eta; //SD of values of eta 
   real<lower=0> sigma_eta;
  vector[Ntime-1] phi[Ncreek, Nspecies]; //random effect for upstream vs. downstream within each creek 
  vector[N_unobserved] mu_unobserved; //vector of unobserved samples to be estimated as params
  vector[Nspecies] mu_phi;
  real<lower=0> sigma_phi;
}

transformed parameters {
  vector[Ncreek] mu[Ntime, Nstations, Nspecies]; //mean eDNA (log scale) from which obs was drawn // dimensions = c(time, station, creek) 
  
  // for (d in 1:Nstations){
  //   mu[1,d,,] = mu_0[d,,];  //slot in param vector m0 at the starting position (t == 1) in the mu array  
  // }
  mu[1,,,] = mu_0;  //slot in param vector m0 at the starting position (t == 1) in the mu array  

  //here, slot in unobserved array elements (params), for sites/samples we wish to estimate as true params
  for (i in 1:N_unobserved){
    mu[unobserved_time_idx[i], unobserved_station_idx[i], unobserved_species_idx[i], unobserved_creek_idx[i]] = mu_unobserved[i];
  }
  
  //apply time-series estimation to array
    for (t in 2:Ntime){
      for (d in 1:Nstations){
        for (j in 1:Nspecies){
          for (i in 1:Ncreek){
    
            mu[t,d,j,i] =  //mu at t+1 and creek i...
                // mu[t-1, d, j, i] +  //is mu at t and creek i...
                beta_1[t-1]*mu[t-1, d, j, i] + //plus some increment determined by slope beta_1
                phi[i, j, t-1] + //
                eta[t-1, d, j, i]; //plus some additional random increment, which varies by upstream/downstream position and time
            
          }
        }
      }
    }
    
}


model {

  for (i in 1:Nobs){
    y_logeDNA[i] ~ normal(mu[time_idx[i], station_idx[i], species_idx[i], creek_idx[i]], 
                          sigma_dna[species_idx[i]]); 
  }

for (d in 1:Nstations){
  for (i in 1:Ncreek){
    for (j in 1:Nspecies){
      eta[,d,j,i] ~ normal(0, sigma_eta); //sigma_eta[species_idx[j]]
      mu_0[d,j,i] ~ normal(0, 5);
    }}}
  
  sigma_eta ~ gamma(1,1);  //TODO: play around w sigma_eta; perhaps make hierarchical across time within creek
  
  // for (t in 1:Ntime){
    sigma_dna ~ gamma(1,1); //[t,]
  // }
    //TODO: do per-species, or else tie to mean; lower means should have higher variances

  // for (j in 1:Nspecies){
      // beta_1[j,] ~ normal(0, 1);  
      beta_1 ~ normal(0, 1);  
  // }
  
  
   for (i in 1:Ncreek){
    for (j in 1:Nspecies){
      phi[i,j] ~ normal(mu_phi[j], sigma_phi);  
    }
    }
  
  mu_phi ~ normal(0,5);
  sigma_phi ~ gamma(1,1);
  mu_unobserved ~ normal(0, 5);
  
}

generated quantities{
 vector[Nobs] log_lik;
   for (i in 1:Nobs){
    log_lik[i]  =  normal_lpdf(y_logeDNA[i] | mu[time_idx[i], station_idx[i], species_idx[i], creek_idx[i]], 
                          sigma_dna[species_idx[i]]);

  }

 
}
