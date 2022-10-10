#Test effect of construction

library(tidyverse)
library(here)
library(rstan)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


d <-readRDS(here("Output/salmonids_abs_abundance_bio_posterior_rpk.RDS")) %>% 
  filter(creek != "4Pad5") %>%  #omit upstream of Padden 5, for symmetry w other creeks
  filter(creek != "2Brn") %>%   #omit Barnes because we have few datapoints for that creek in terms of abs concentration
  filter(species != "Oncorhynchus tshawytscha") %>%  #for now, filter out as rare
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

construction <- which(d$creek_idx == 3 & d$time_idx %in% c(8:10))
f <- d %>%  #use this to filter out anything you don't want to model
  arrange(time_idx, creek_idx) %>% 
  ungroup()
f <- f[-construction,]

## Run new model without those construction observations
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



stan_data <- list(
  Nobs = nrow(f),
  Ntime = length(unique(f$time)),
  Ncreek = length(unique(f$creek)),
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

stanMod_noConstr = stan(file = here("Scripts/timeseries_model/timeSeries_multispecies5.stan") ,data = stan_data,
                        verbose = FALSE, chains = 3, thin = 1,
                        warmup = 500, iter = 1000,
                        control = list(adapt_init_buffer = 175,
                                       max_treedepth=12,
                                       stepsize=0.01,
                                       adapt_delta=0.7,
                                       metric="diag_e"),
                        #pars = stan_pars,
                        refresh = 10,
                        boost_lib = NULL,
                        init = 0,
                        #sample_file = "temp/tmp.csv"
)

# saveRDS(list(stan_data,
#              stanMod_noConstr,
#              f), file = "modelFit_20221009_noConstruction.RDS")


