library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("/Volumes/GoogleDrive/My Drive/Kelly_Lab/Functions/qPCR_calibration/calibrate_qPCR.R")

a <- read.csv("../../Output/qpcr/cut_final.csv")
qMod_out <- run_qPCR_model("../../Output/qpcr/cut_final.csv",
                          here("Scripts", "qm_qpcr_model", "qPCR_calibration_enchilada.stan"))

res <- qMod_out$results_qPCR %>% 
  dplyr::select(time, creek, station, mean_concentration_est)

padden <- res %>% 
  filter(creek == "4Pad" & station == "Up11") %>% 
  mutate(month = substr(time, 1,2),
         year = substr(time, 3,4)) %>% 
  arrange(year, month)
padden$time_idx = match(padden$time, unique(padden$time))

plot(log(padden$mean_concentration_est) ~ padden$time)

stan_data <- list(
  N = nrow(padden),
  Ntime = length(unique(padden$time)),
  time_idx = match(padden$time, unique(padden$time)),
  # y_fish = d$count,
  y_logeDNA = log(padden$mean_concentration_est)
)



stanMod = stan(file = "timeSeries_DNAOnly.stan" ,data = stan_data,
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

plot(stanMod, par = c("mu", "sigma_phi", "sigma_dna"))

data.frame(
  time_idx = unique(stan_data$time_idx),
  mu_mean = summary(stanMod, par = "mu")$summary[,1],
  ci_25 = summary(stanMod, par = "mu")$summary[,5],
  ci_75 = summary(stanMod, par = "mu")$summary[,7]
) %>% 
  left_join(padden) %>% 
  ggplot(aes(x = time_idx, y = log(mean_concentration_est))) +
    geom_point() +
    geom_point(aes(x = time_idx, y = mu_mean), color = "red") +
    geom_segment(aes(x = time_idx, xend = time_idx, y = ci_25, yend = ci_75), color = "red")

