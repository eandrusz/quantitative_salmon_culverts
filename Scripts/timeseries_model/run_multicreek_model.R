library(tidyverse)
library(here)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("/Volumes/GoogleDrive/My Drive/Kelly_Lab/Functions/qPCR_calibration/calibrate_qPCR.R")

a <- read.csv("../../Output/qpcr/cut_final.csv")
qMod_out <- run_qPCR_model("../../Output/qpcr/cut_final.csv",
                          here("Scripts", "qm_qpcr_model", "qPCR_calibration_enchilada.stan"))

res <- qMod_out$results_qPCR %>%
  dplyr::select(time, creek, station, mean_concentration_est, biorep)



# b <- readRDS("/Users/rpk/Downloads/salmonids_abs_abundance_bio_posterior.RDS")
# head(a)
# b %>%
#   filter(species == "Oncorhynchus clarkii" & creek == "4Pad5" & station == 2) %>%
#   pivot_longer(-c(time, creek, station, bio, species)) %>% 
#   group_by(time, creek, station, bio, species) %>% 
#   summarize(mean_conc = mean(value)) %>% 
#   ggplot(aes(x = time, y = log(mean_conc))) +
#     geom_point() +
#     facet_grid(~creek)

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

# construction <- which(d$creek_idx == 4 & d$time_idx %in% c(8:10))
f <- d  #use this to filter out anything you don't want to model
  
  
#note unobserved datapoints
# datamat <- f %>% 
#   select(creek_idx, time_idx) %>% 
#   distinct() %>% 
#   mutate(present = 1) %>% 
#   pivot_wider(names_from = time_idx, values_from = present, values_fill = 0) %>% 
#   arrange(creek_idx) %>% 
#   select(-creek_idx) %>% 
#   as.matrix()
# 
# findZeros <- function(mat){
#   rowidx <- NA
#   colidx <- NA
#   
#   for(i in 1:nrow(mat)){
#     tmp <- which(mat[i,]==0)
#     if (length(tmp) > 0) {
#       rowidx <- c(rowidx, rep(i, length(tmp)))
#       colidx <- c(colidx, tmp)
#     }
#   }
#   return(data.frame(rowidx, colidx)[-1,])
# }
# missing <- findZeros(datamat)

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

stanMod = stan(file = "timeSeries_multicreek_wCulvertEffect.stan" ,data = stan_data,
               verbose = FALSE, chains = 3, thin = 1,
               warmup = 500, iter = 700,
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


#posterior predictive check:
#95% CI for mu, sigma; w observed data

(p2 <- f %>% 
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
                               creek_idx == 4 ~ "Squalicum")) %>% 
  ggplot(aes(x = month, y = log(mean_concentration_est))) +
    geom_point() +
    geom_point(aes(x = month, y = mean_est), color = "red", size = 1.5, alpha = .5) +
    geom_segment(aes(x = month, xend = month, y = Sigma_low, yend = Sigma_high), color = "red", size = .7, alpha = .5) +
    geom_segment(aes(x = month, xend = month, y = TwoSigma_low, yend = TwoSigma_high), color = "red", size = .3, alpha = .5) +
    facet_grid(station~creekname) +
    ggtitle("Posterior Predictive Check \n(predicted mean +/- 1 and 2SD)"))
#ggsave(p2, filename = "pp_check_cutthroat.pdf")
#ggsave(p2, filename = "pp_check_cutthroat.jpeg")


# to see effect of culvert at Squalicum, subtract downstream from upstream in the eta terms:
#dimensions of eta are c(time, station, creek)
unlist(extract(stanMod, pars = "eta[8, 2, 4]")) - unlist(extract(stanMod, pars = "eta[8, 1, 4]")) %>% hist()

s <- extract(stanMod, pars = "eta")
data.frame(August = s$eta[,6,2,4] - s$eta[,6,1,4],
           September = s$eta[,7,2,4] - s$eta[,7,1,4],
           October = s$eta[,8,2,4] - s$eta[,8,1,4],
           November = s$eta[,9,2,4] - s$eta[,9,1,4],
           December = s$eta[,10,2,4] - s$eta[,10,1,4]) %>% 
  mutate(idx = 1:n()) %>% 
  pivot_longer(-idx, names_to = "Month", values_to = "CulvertDifference") %>% 
  filter(Month == "August") %>% 
  ggplot(aes(x = Month, y = CulvertDifference)) +
    geom_boxplot() +
    geom_point(alpha = .2) +
    ylab("(log qPCR concentration \nUpstream - Downstream)")



