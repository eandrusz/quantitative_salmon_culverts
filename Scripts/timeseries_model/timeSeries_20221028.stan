
data {
  int<lower=0> Ncreek;
  int<lower=0> Ntime;
  int<lower=0> Nobs; //total observations
  int<lower=0> Nstations; //probably 2, downstream and upstream
  int<lower=0> Nspecies;
  int<lower=0> Nconstruction; //N culvert states; removed (=2) vs. nonremoved (=1)
  int<lower=0> MinconstructionTimepoint; //first timepoint with construction
  int<lower=0> NconstructionTimepoints;
  int<lower=0> time_idx[Nobs];
  int<lower=0> creek_idx[Nobs];
  int<lower=0> species_idx[Nobs];
  int<lower=0> station_idx[Nobs]; //downstream = 1, upstream = 2
  int<lower=0> construction_idx[Nobs];
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
  vector<lower=-1, upper = 1>[Nspecies] beta_1; //slope relating present concentration to past concentration
  vector<lower=0>[Nspecies] sigma_dna; //SD of DNA concentration (log scale)
  real<lower=0> sigma_eta;
  vector[Ntime-1] alpha[Ncreek, Nspecies]; //random effect for upstream vs. downstream within each creek 
  real gamma_unk[Ntime-1, Nspecies, Nconstruction-1];
  vector[N_unobserved] mu_unobserved; //vector of unobserved samples to be estimated as params
  vector[Nspecies] mu_alpha;
  real<lower=0> sigma_alpha;
}

transformed parameters {
  vector[Ncreek] mu[Ntime, Nstations, Nspecies]; //mean eDNA (log scale) from which obs was drawn // dimensions = c(time, station, creek) 
  real gamma[Ntime-1, Nspecies, Nconstruction];
  
  mu[1,,,] = mu_0;  //slot in param vector m0 at the starting position (t == 1) in the mu array  
  //here, slot in unobserved array elements (params), for sites/samples we wish to estimate as true params
  for (i in 1:N_unobserved){
    mu[unobserved_time_idx[i], unobserved_station_idx[i], unobserved_species_idx[i], unobserved_creek_idx[i]] = mu_unobserved[i];
  }

  //from Ole by email Oct 18, 2022:
  //To be identifiable, you need to define one of the gammas equal to 0.  
  //I'd suggest  gamma_(t,d=1,j,r=not-removed) = 0 so the gamma term is interpretable as the effect of going upstream from the culvert before you removed the culvert.
      
      gamma[,,1] = rep_array(0.0, Ntime-1, Nspecies); //define gamma to be 0 for all species/time where construction is absent
      //then fill in the rest of the array w param values
       for (t in 2:Ntime){
          for (j in 1:Nspecies){
            // for (i in 1:Ncreek){
            //for (d in 1:Nstations){
              // for (r in 2:Nconstruction){
              gamma[t-1,j,2] = gamma_unk[t-1,j,1];
            }}

       //define gamma 
       
       
       // only estimate gamma at time points where construction has happened, otherwise define as 0
       //gamma[,,2,1:(MinconstructionTimepoint-2)] = rep_array(0.0, Nspecies, Nstations, (Ntime-1)-NconstructionTimepoints);
  
  //apply time-series estimation to array
    for (t in 2:Ntime){
      for (d in 1:Nstations){
        for (j in 1:Nspecies){
          for (i in 1:Ncreek){
             for (r in 1:Nconstruction){
    
            mu[t,d,j,i] =  
                beta_1[j]*mu[t-1, d, j, i] + //slope; autocorr
                alpha[i, j, t-1] + // intercept, shared across stations
                gamma[t-1, j, r] + // effect of culvert removal, shared across creeks
                eta[t-1, d, j, i]; //plus some additional random increment, which varies by upstream/downstream position and time
            
            }
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


// Priors
for (d in 1:Nstations){
  for (i in 1:Ncreek){
    for (j in 1:Nspecies){
      // for (r in 1:Nconstruction){
      eta[,d,j,i] ~ normal(0, sigma_eta); 
      mu_0[d,j,i] ~ normal(0, 5);
      gamma_unk[,j,1] ~ normal(0, 5);  
    }}}
    // }
  
  // for (i in 1:Ncreek){
  //   for (j in 1:Nspecies){
  //     for (d in 1:Nstations){
  //       
  //   }}}

  sigma_eta ~ gamma(1,1); 
  sigma_dna ~ gamma(1,1); 

  beta_1 ~ normal(0, 1);  

    for (i in 1:Ncreek){
    for (j in 1:Nspecies){
      alpha[i,j] ~ normal(mu_alpha[j], sigma_alpha);  
    }
    }
  
  mu_alpha ~ normal(0,5);
  sigma_alpha ~ gamma(1,1);
  mu_unobserved ~ normal(0, 5);
}

generated quantities{
 vector[Nobs] log_lik;
   for (i in 1:Nobs){
    log_lik[i]  =  normal_lpdf(y_logeDNA[i] | mu[time_idx[i], station_idx[i], species_idx[i], creek_idx[i]], 
                          sigma_dna[species_idx[i]]);

  }

 
}
