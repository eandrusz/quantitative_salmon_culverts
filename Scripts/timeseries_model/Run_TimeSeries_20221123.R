#Clean Multicreek, Multispecies time-series analysis

library(tidyverse)
library(here)
library(rstan)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(shinystan)

d <- readRDS(here("Output","20221123_abundance_flowcorrected.RDS")) %>% 
  filter(creek != "Barnes") %>%   #omit Barnes because we have few datapoints for that creek in terms of abs concentration
  mutate(station_idx = ifelse(station == "Up", 2, 1),
         bio = as.integer(bio)) %>% 
  dplyr::select(-station) %>% 
  arrange(creek) %>% 
  arrange(newtime)
  

#index outside of tidy notation, so I know it works
d$creek_idx <- match(d$creek, unique(d$creek))
d$time_idx <- match(d$newtime, unique(d$newtime))
d$species_idx <- match(d$species, unique(d$species))

# insert dummy data to reflect non-observation at Barnes upstream for cutthroat
# d <- d %>%
#   full_join(expand_grid(creek ="2Brn",
#                         creek_idx = 2,
#                         station_idx = 2,
#                         species_idx = 1,
#                         time_idx = 1:length(unique(d$time_idx)),
#                         meandnaconc = 0.01))
# 


# create duplicate object, useful for filtering to subset, etc
f <- d %>% arrange(time_idx, creek_idx) %>% 
  ungroup() 
  # filter(station_idx == 1) %>%  #just look downstream
  # filter(log(meandnaconc) > 2.3) #impose minimum of 10 copies/L, because this is the approx limit of reliability

# add construction index
construction <- which(f$creek == "Padden" & f$time_idx %in% c(7:12))
  f$construction_idx <- 1
  f$construction_idx[construction] <- 2

#find elements to treat as missing data
missing <- f %>% 
  dplyr::select(creek_idx, time_idx, station_idx, species_idx) %>% 
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
  Nconstruction = length(unique(f$construction_idx)),
  MinconstructionTimepoint = f[f$construction_idx==2,]$time_idx %>% min(),
  NconstructionTimepoints = f[f$construction_idx==2,]$time_idx %>% unique() %>% length(),
  time_idx = f$time_idx,
  creek_idx = f$creek_idx,
  station_idx = f$station_idx,
  species_idx = f$species_idx,
  construction_idx = f$construction_idx,
  # construction_indicator = f$construction_idx-1,
  #y_logeDNA = log(f$meandnaconc),
  y_logeDNA = log(f$meandnaconcflow),
  N_unobserved = nrow(missing),
  unobserved_time_idx= as.array(missing$time_idx),
  unobserved_creek_idx= as.array(missing$creek_idx),
  unobserved_station_idx = as.array(missing$station_idx),
  unobserved_species_idx = as.array(missing$species_idx),
  unobserved_construction_idx = ifelse(missing$creek_idx == 2 & missing$time_idx > 6, 2, 1) #NOTE assumes Padden is creek_idx = 2
)

#######RUN Stan Model#######

stanMod = stan(file = here("Scripts/timeseries_model/timeSeries_20221029.stan") ,data = stan_data,
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

# plot(stanMod, par = c("mu"))
# plot(stanMod, par = c("sigma_eta", "sigma_dna"))
# #plot(stanMod, par = c("eta"))
# plot(stanMod, par = c("delta"))
# rstan::traceplot(stanMod, par = c("beta_1"))
# plot(stanMod, par = c("alpha"))
# plot(stanMod, par = c("gamma"))
# rstan::traceplot(stanMod, par = c("gamma[1,1,1,7]",
#                            "gamma[2,1,1,7]",
#                            "gamma[1,1,2,5]",
#                            "gamma[2,1,2,7]"
#                            ))
# 
#gamma[station, species, construction, time]
summary(stanMod, par = "gamma")$summary[,1] %>% sort() %>% head()
summary(stanMod, par = "gamma")$summary[,"Rhat"] %>% sort() %>% tail()
summary(stanMod)$summary[,"Rhat"] %>% sort() %>% tail()
# shinystan::launch_shinystan(stanMod)
#drop_parameters(as.shinystan(stanMod), pars = c("eta", "alpha", "mu", "mu_0")) %>% launch_shinystan()
# 
# traceplot(stanMod, par = c("sigma_dna"))
# ############################

# Save, if desired
# saveRDS(list(stan_data,
#              stanMod,
#              f), file = here("Scripts/timeseries_model/20221123_modelFit_flowcorrected.RDS"))
# 


