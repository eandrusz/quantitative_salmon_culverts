#Clean Multicreek, Multispecies time-series analysis

library(tidyverse)
library(here)
library(rstan)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


d <-readRDS(here("Output/salmonids_abs_abundance_bio_posterior_rpk.RDS")) %>% 
  filter(creek != "4Pad5") %>%  #omit upstream of Padden 5, for symmetry w other creeks
  filter(creek != "2Brn") %>%   #omit Barnes because we have few datapoints for that creek in terms of abs concentration
  #filter(species != "Oncorhynchus tshawytscha") %>%  #for now, filter out as rare
  mutate(station_idx = as.integer(station),
         bio = as.integer(bio)) %>% 
  dplyr::select(-station) %>% 
  arrange(creek) %>% 
  mutate(month = substr(time, 1,2),
         year = substr(time, 3,4)) %>% 
  arrange(year, month) %>% 
  dplyr::select(-c(month, year))
#index outside of tidy notation, so I know it works
d$creek_idx <- match(d$creek, unique(d$creek))
d$time_idx <- match(d$time, unique(d$time))
d$species_idx <- match(d$species, unique(d$species))

#insert dummy data to reflect non-observation at Barnes upstream for cutthroat
# d <- d %>% 
#   full_join(expand_grid(creek ="2Brn",
#                         creek_idx = 2,
#                         station_idx = 2,
#                         species_idx = 1,
#                         time_idx = 1:length(unique(d$time_idx)),
#                         meandnaconc = 0.01))



# create duplicate object, useful for filtering to subset, etc
f <- d %>% arrange(time_idx, creek_idx) %>% 
  ungroup()

#find elements to treat as missing data
missing <- f %>% 
  select(creek_idx, time_idx, station_idx, species_idx) %>% 
  distinct() %>% 
  mutate(present = 1) %>% 
  right_join(expand_grid(creek_idx = 1:length(unique(f$creek_idx)),
                         time_idx = 1:length(unique(f$time_idx)),
                         station_idx = 1:length(unique(f$station_idx)),
                         species_idx = 1:length(unique(f$species_idx)))) %>% 
  filter(is.na(present)) %>% 
  dplyr::select(-present)

# f %>% 
#   filter(creek == "2Brn")


#format for stan input
stan_data <- list(
  Nobs = nrow(f),
  Ntime = length(unique(f$time_idx)),
  Ncreek = length(unique(f$creek_idx)),
  Nstations = length(unique(f$station_idx)),
  Nspecies = length(unique(f$species_idx)),
  time_idx = f$time_idx,
  creek_idx = f$creek_idx,
  station_idx = f$station_idx,
  species_idx = f$species_idx,
  y_logeDNA = log(f$meandnaconc),
  N_unobserved = nrow(missing),
  unobserved_time_idx= as.array(missing$time_idx),
  unobserved_creek_idx= as.array(missing$creek_idx),
  unobserved_station_idx = as.array(missing$station_idx),
  unobserved_species_idx = as.array(missing$species_idx)
)

#######RUN Stan Models#######
#here, one model with multiple betas per time point (one for each creek) and one model with a single beta across creeks in a time point

# stanMod_multiBeta = stan(file = here("Scripts/timeseries_model/timeSeries_multispecies4.stan") ,data = stan_data,
#                          verbose = FALSE, chains = 3, thin = 1,
#                          warmup = 500, iter = 1000,
#                          control = list(adapt_init_buffer = 175,
#                                         max_treedepth=12,
#                                         stepsize=0.01,
#                                         adapt_delta=0.7,
#                                         metric="diag_e"),
#                          #pars = stan_pars,
#                          refresh = 10,
#                          boost_lib = NULL,
#                          init = 0,
#                          #sample_file = "temp/tmp.csv"
# )
stanMod_oneBeta = stan(file = here("Scripts/timeseries_model/timeSeries_oneBeta_phi3.stan") ,data = stan_data,
                       verbose = FALSE, chains = 3, thin = 1,
                       warmup = 500, iter = 1500,
                       control = list(adapt_init_buffer = 175,
                                      max_treedepth=12,
                                      stepsize=0.01,
                                      adapt_delta=0.7,
                                      metric="diag_e"),
                       #init = 0,
                       #pars = stan_pars,
                       refresh = 10,
                       boost_lib = NULL
                       #init = stan_init_f2(n.chain=N_CHAIN,N_species=N_species),
                       #sample_file = "temp/tmp.csv"
)
# stanMod_phi = stan(file = here("Scripts/timeseries_model/timeSeries_oneBeta_phi.stan") ,data = stan_data,
#                        verbose = FALSE, chains = 3, thin = 1,
#                        warmup = 500, iter = 1000,
#                        control = list(adapt_init_buffer = 175,
#                                       max_treedepth=12,
#                                       stepsize=0.01,
#                                       adapt_delta=0.7,
#                                       metric="diag_e"),
#                        init = 0,
#                        #pars = stan_pars,
#                        refresh = 10,
#                        boost_lib = NULL
#                        #init = stan_init_f2(n.chain=N_CHAIN,N_species=N_species),
#                        #sample_file = "temp/tmp.csv"
# )

#model selection
# loo::loo_compare(loo(stanMod_multiBeta),
#                  loo(stanMod_oneBeta))
#we prefer the simpler model, slightly

plot(stanMod_oneBeta, par = c("mu"))
plot(stanMod_oneBeta, par = c("sigma_eta", "sigma_dna"))
#plot(stanMod_oneBeta, par = c("eta"))
plot(stanMod_oneBeta, par = c("beta_1"))
plot(stanMod_oneBeta, par = c("phi"))

# 
# traceplot(stanMod_oneBeta, par = c("sigma_eta", "sigma_dna"))
# ############################

# Save, if desired
# saveRDS(list(stan_data,
#              stanMod_oneBeta,
#              f), file = "modelFit_20221011.RDS")



