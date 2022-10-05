library(tidyverse)
library(rstan)
library(here)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# source(here("Scripts", "functions", "calibrate_qPCR.R"))
# 
# a <- read.csv(here("Output","qpcr","cut_final.csv"))
# qMod_out <- run_qPCR_model("../../Output/qpcr/cut_final.csv",
#                           here("Scripts", "qm_qpcr_model", "qPCR_calibration_enchilada.stan"))

qMod_out <- readRDS(here("Output","qpcr","cut_qMod_out.RDS"))

res <- qMod_out$results_qPCR %>% 
  dplyr::select(time, creek, station, mean_concentration_est) %>% 
  mutate(creek = case_when(creek == "4Pad" & station == "Dn" ~ "4Pad11",
                           creek == "4Pad" & station == "Up11" ~ "4Pad11",
                           TRUE ~ creek)) %>% 
  mutate(creek = case_when(creek == "4Pad" & station == "Up5" ~ "4Pad5",
                           TRUE ~ creek)) %>% 
  mutate(station = case_when(creek == "4Pad11" & station == "Up11" ~ "Up",
                             TRUE ~ station)) %>% 
  mutate(station = case_when(creek == "4Pad5" & station == "Up5" ~ "Up",
                             TRUE ~ station))


allcreek <- c("1Prt","2Brn","3Chk","4Pad11", "4Pad5","5Sqm")
allstation <- c("Dn","Up")
allmodeled <- list()
counter=1

for (i in 1:length(allcreek)) {
  for (j in 1:length(allstation)) {
    
    x <- res %>% 
      filter(creek == allcreek[i] & station == allstation[j]) %>% 
      mutate(month = substr(time, 1,2),
             year = substr(time, 3,4)) %>% 
      arrange(year, month)
    x$time_idx = match(x$time, unique(x$time))
    
    #plot(log(padden$mean_concentration_est) ~ padden$time)
    
    stan_data <- list(
      N = nrow(x),
      Ntime = length(unique(x$time)),
      time_idx = match(x$time, unique(x$time)),
      # y_fish = d$count,
      y_logeDNA = log(x$mean_concentration_est)
    )
    
    
    
    stanMod = stan(file = here("Scripts", "timeseries_model", "timeSeries_DNAOnly.stan"), data = stan_data,
                   verbose = FALSE, chains = 3, thin = 1,
                   warmup = 500, iter = 1000,
                   control = list(adapt_init_buffer = 175,
                                  max_treedepth=12,
                                  stepsize=0.01,
                                  adapt_delta=0.7,
                                  metric="diag_e"),
                   #pars = stan_pars,
                   refresh = 10,
                   boost_lib = NULL
                   #init = stan_init_f2(n.chain=N_CHAIN,N_species=N_species),
                   #sample_file = "temp/tmp.csv"
    )
    
    #plot(stanMod, par = c("mu", "sigma_phi", "sigma_dna"))
    
    allmodeled[[counter]] <- data.frame(
      time_idx = unique(stan_data$time_idx),
      mu_mean = summary(stanMod, par = "mu")$summary[,1],
      ci_25 = summary(stanMod, par = "mu")$summary[,5],
      ci_75 = summary(stanMod, par = "mu")$summary[,7]
    ) %>% 
      left_join(x) 
    
    counter = counter + 1
    
  }
  
}


allmodeled2 <- do.call(rbind, allmodeled)
  
allmodeled2 %>%   
  ggplot(aes(x = time_idx, y = log(mean_concentration_est))) +
    geom_point() +
    geom_point(aes(x = time_idx, y = mu_mean), color = "red") +
    geom_segment(aes(x = time_idx, xend = time_idx, y = ci_25, yend = ci_75), color = "red") + 
    facet_grid(~creek ~ station)

