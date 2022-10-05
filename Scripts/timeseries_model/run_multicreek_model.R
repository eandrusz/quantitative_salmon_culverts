library(tidyverse)
library(here)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source(here("Scripts", "functions","calibrate_qPCR.R"))

qMod_out <- readRDS(here("Output","qpcr","cut_qMod_out.RDS"))

res <- qMod_out$results_qPCR %>%
  dplyr::select(time, creek, station, mean_concentration_est, biorep)

d <- res %>% 
  #filter(station %in% c("Dn")) %>% 
  #filter(creek == "4Pad") %>% 
  filter(!station %in% c("Up5")) %>% 
  filter(creek != "2Brn") %>%  #Barnes has too much missing data for the moment
  mutate(station_idx = case_when(station == "Dn" ~ 1,
                                 station %in% c("Up", "Up11") ~ 2)) %>% 
  dplyr::select(-c(station)) %>% 
 # filter(station_idx == 1) %>% 
  arrange(creek) %>% 
  mutate(creek_idx = match(creek, unique(creek))) %>% 
  mutate(month = substr(time, 1,2),
         year = substr(time, 3,4)) %>% 
  arrange(year, month) %>% 
  mutate(time_idx = match(time, unique(time))) %>% 
  dplyr::select(-c(month, year)) %>% 
  arrange(time_idx, creek_idx)

f <- d  #use this to filter out anything you don't want to model

### USE THIS TO SIMULATE WHAT YOU WOULD FIND IN PADDEN WITHOUT CONSTRUCTION
# f_construct <- f %>%
#   mutate(keep = case_when(creek == "4Pad" & time_idx > 7 ~ FALSE,
#                           creek == "4Pad" & time_idx < 8 ~ TRUE,
#                           creek != "4Pad" ~ TRUE)) %>%
#   filter(keep == TRUE) %>%
#   select(-keep)
# 
# f <- f_construct

#3D find zeros
missing <- f %>% 
  select(creek_idx, time_idx, station_idx) %>% 
  distinct() %>% 
  mutate(present = 1) %>% 
  right_join(expand_grid(creek_idx = 1:length(unique(f$creek_idx)),
                         time_idx = 1:length(unique(f$time_idx)),
                         station_idx = 1:length(unique(f$station_idx)))) %>% 
  filter(is.na(present)) %>% 
  dplyr::select(-present)

stan_data <- list(
  Nobs = nrow(f),
  Ntime = length(unique(f$time)),
  Ncreek = length(unique(f$creek)),
  Nstations = length(unique(f$station_idx)),
  time_idx = f$time_idx,
  creek_idx = f$creek_idx,
  station_idx = f$station_idx,
  y_logeDNA = log(f$mean_concentration_est),
  N_unobserved = nrow(missing),
  unobserved_time_idx= as.array(missing$time_idx),
  unobserved_creek_idx= as.array(missing$creek_idx),
  unobserved_station_idx = as.array(missing$station_idx)
)

stanMod = stan(file = here("Scripts", "timeseries_model","timeSeries_multicreek_wCulvertEffect.stan") ,data = stan_data,
               verbose = FALSE, chains = 3, thin = 1,
               warmup = 500, iter = 1500,
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

plot(stanMod, par = c("mu"))
plot(stanMod, par = c("sigma_eta", "sigma_dna", "beta_1"))
plot(stanMod, par = c("eta"))
# plot(stanMod, par = c("phi"))

#posterior predictions for observed DNA concentrations
a <- extract(stanMod, par = "mu")
b <- extract(stanMod, par = "sigma_dna") %>% unlist()

pp <- expand_grid(creek_idx = 1:length(unique(f$creek_idx)),
                  time_idx = 1:length(unique(f$time_idx)),
                  station_idx = 1:length(unique(f$station_idx)))
pp <- pp %>%
  bind_cols(matrix(NA, nrow = nrow(pp), ncol = length(a$mu[,1,1,1]))) %>%
  nest(data = 4:ncol(.))

for (i in 1:nrow(pp)){
  pp$data[i][[1]] <- rnorm(length(a$mu[,1,1,1]),
                           a$mu[,pp$time_idx[i],pp$station_idx[i],pp$creek_idx[i]],
                           b)
}

pp <- pp %>%
  mutate(creekname = case_when(creek_idx == 1 ~ "Portage",
                               creek_idx == 2 ~ "Chuckanut",
                               creek_idx == 3 ~ "Padden",
                               creek_idx == 4 ~ "Squalicum"))

## TO PLOT DOWN AND UP ON SAME PLOT
pp_plot <- pp 
pp_plot$mean_conc <- 0
pp_plot$ci.25 <- 0
pp_plot$ci.75 <- 0

for (i in 1:nrow(pp_plot)){
  l.model <- lm(unlist(data[[i]])~1, pp)
  cis <- confint(l.model, level=0.50)
  pp_plot$mean_conc[i] <- as.numeric(l.model$coefficients)
  pp_plot$ci.25[i] <- cis[1]
  pp_plot$ci.75[i] <- cis[2]
}

pp_plot %>% 
  #filter(creekname == "Padden") %>% 
  ggplot(aes(x=time_idx, y=mean_conc, color=as.factor(station_idx))) + 
  geom_point() +
  geom_segment(aes(x = time_idx, xend = time_idx, y = ci.25, yend = ci.75), size = .7, alpha = .5) +
  facet_wrap(~creekname) +
  labs(color="Downstream (1) or Upstream (2)")

## TO CALCULATE DIFFERENCE BETWEEN DOWN AND UP AND THEN PLOT 
#to derive differences in up vs. down, something like this?
pp_diff <- pp %>% 
  unnest(cols = data) %>%
  pivot_wider(names_from = station_idx, values_from = data, values_fn = list) %>%
  rename(Down = '1') %>%
  rename(Up = '2') %>%
  mutate(diff = map2(.x = Down, .y = Up, .f = function(x,y) y-x)) %>% 
  select(-c(Down, Up))

pp_diff$mean_diff <- 0
pp_diff$ci.25 <- 0
pp_diff$ci.75 <- 0

for (i in 1:nrow(pp_diff)){
  l.model <- lm(unlist(diff[[i]])~1, pp_diff)
  cis <- confint(l.model, level=0.50)
  pp_diff$mean_diff[i] <- as.numeric(l.model$coefficients)
  pp_diff$ci.25[i] <- cis[1]
  pp_diff$ci.75[i] <- cis[2]
}


ggplot(pp_diff, aes(x=time_idx, y = mean_diff)) +
  geom_point() +
  theme_bw() +
  geom_segment(aes(x = time_idx, xend = time_idx, y = ci.25, yend = ci.75), size = .7, alpha = .5) +
  geom_hline(yintercept=0, linetype='dotted' ) +
  facet_wrap(~creekname) +
  labs(y="Log Up - Log Down")






################################################## GRAVEYARD ####################################################
#shinystan::launch_shinystan(stanMod)

#dimensions = c(time, station, creek) 

resOut <- expand_grid(time_idx = 1:length(unique(f$time_idx)),
                      station_idx = 1:length(unique(f$station_idx)),
                      creek_idx = 1:length(unique(f$creek_idx))) %>%
  mutate(mean_est = summary(stanMod, par = "mu")$summary[,1],
         ci25 = summary(stanMod, par = "mu")$summary[,5],
         ci75 = summary(stanMod, par = "mu")$summary[,7]) %>%
  left_join(f) %>%
  left_join(data.frame(sigma_dna = summary(stanMod, par = "sigma_dna")$summary[,1],
                       creek_idx = 1:1:length(unique(f$creek_idx))))

#plot
resOut %>%
  ggplot(aes(x = time_idx, y = log(mean_concentration_est))) +
    geom_point() +
    geom_point(aes(x = time_idx, y = mean_est), color = "red") +
    geom_segment(aes(x = time_idx, xend = time_idx, y = ci25, yend = ci75), color = "red") +
    facet_grid(station_idx~creek_idx) +
    ggtitle("Estimated Mean Concentrations w Interquartile Range")
# 
# 
# #posterior predictive check:
# #95% CI for mu, sigma; w observed data
# 
p2 <- f %>%
  # filter(creek_idx == focalCreek) %>%
  mutate(logY = log(mean_concentration_est)) %>%
  right_join(resOut) %>%
  mutate(month = ifelse(time_idx < 11, time_idx + 2, time_idx - 10),
         month = as.factor(month)) %>%
  mutate(Sigma_low = mean_est - (summary(stanMod, par = "sigma_dna")$summary[,1]),
         TwoSigma_low = mean_est - 2*(summary(stanMod, par = "sigma_dna")$summary[,1]),
         Sigma_high = mean_est + (summary(stanMod, par = "sigma_dna")$summary[,1]),
         TwoSigma_high = mean_est + 2*(summary(stanMod, par = "sigma_dna")$summary[,1])) %>%
  mutate(station = ifelse(station_idx == 1, "Downstream", "Upstream")) %>%
  mutate(creekname = case_when(creek_idx == 1 ~ "Portage",
                               creek_idx == 2 ~ "Chuckanut",
                               creek_idx == 3 ~ "Padden",
                               creek_idx == 4 ~ "Squalicum"))

#p2_observed <- p2
#p2_expected <- p2


  ggplot(p2, aes(x = month, y = log(mean_concentration_est))) +
    geom_point() +
    geom_point(aes(x = month, y = mean_est), color = "red", size = 1.5, alpha = .5) +
    geom_segment(aes(x = month, xend = month, y = Sigma_low, yend = Sigma_high), color = "red", size = .7, alpha = .5) +
    geom_segment(aes(x = month, xend = month, y = TwoSigma_low, yend = TwoSigma_high), color = "red", size = .3, alpha = .5) +
    facet_grid(station~creekname) +
    ggtitle("Posterior Predictive Check \n(predicted mean +/- 1 and 2SD)")
# #ggsave(p2, filename = "pp_check_cutthroat.pdf")
# #ggsave(p2, filename = "pp_check_cutthroat.jpeg")
# 
#   
# padden_expected <- p2_expected %>% 
#   filter(creekname == "Padden") %>% 
#   unite(sample, c("time_idx", "creek_idx", "station_idx"), remove=FALSE) %>% 
#   select(c("sample","mean_est", "month", "station")) %>% 
#   rename(expected_est = mean_est)
# 
# padden_observed <- p2_observed %>% 
#   filter(creekname == "Padden") %>% 
#   unite(sample, c("time_idx", "creek_idx", "station_idx")) %>% 
#   select(c("sample","mean_est")) 
# 
# padden_compare <- padden_expected %>% 
#   left_join(padden_observed, by = "sample")
# 
# ggplot(padden_compare, aes(x=month, y=log(mean_est))) +
#   geom_point() + 
#   geom_point(aes(x = month, y = log(expected_est)), color = "blue", size = 1.5, alpha = .5) +
#   facet_wrap(~station)
#   
#   
# 
# # to see effect of culvert at Squalicum, subtract downstream from upstream in the eta terms:
# #dimensions of eta are c(time, station, creek)
# unlist(extract(stanMod, pars = "eta[8, 2, 4]")) - unlist(extract(stanMod, pars = "eta[8, 1, 4]")) %>% hist()
# 
# s <- extract(stanMod, pars = "eta")
# data.frame(August = s$eta[,6,2,4] - s$eta[,6,1,4],
#            September = s$eta[,7,2,4] - s$eta[,7,1,4],
#            October = s$eta[,8,2,4] - s$eta[,8,1,4],
#            November = s$eta[,9,2,4] - s$eta[,9,1,4],
#            December = s$eta[,10,2,4] - s$eta[,10,1,4]) %>% 
#   mutate(idx = 1:n()) %>% 
#   pivot_longer(-idx, names_to = "Month", values_to = "CulvertDifference") %>% 
#   filter(Month == "August") %>% 
#   ggplot(aes(x = Month, y = CulvertDifference)) +
#     geom_boxplot() +
#     geom_point(alpha = .2) +
#     ylab("(log qPCR concentration \nUpstream - Downstream)")




