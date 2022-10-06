library(tidyverse)
library(here)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# source("/Volumes/GoogleDrive/My Drive/Kelly_Lab/Functions/qPCR_calibration/calibrate_qPCR.R")
# 
# #a <- read.csv("../../Output/qpcr/cut_final.csv")
# qMod_out <- run_qPCR_model("../../Output/qpcr/cut_final.csv",
#                           here("Scripts", "qm_qpcr_model", "qPCR_calibration_enchilada.stan"))
# qMod_out2 <- run_qPCR_model("../../Output/qpcr/coho_final.csv",
#                            here("Scripts", "qm_qpcr_model", "qPCR_calibration_enchilada.stan"))
# 

# 
# 
# res <- qMod_out$results_qPCR %>%
#   dplyr::select(time, creek, station, mean_concentration_est, biorep) %>% 
#   mutate(species_idx = 1) %>% 
#   bind_rows(qMod_out2$results_qPCR %>% dplyr::select(time, creek, station, mean_concentration_est, biorep) %>% 
#               mutate(species_idx = 2))

d <-readRDS("/Volumes/GoogleDrive/My Drive/RPKDesktop/github_repos/quantitative_salmon_culverts/Output/salmonids_abs_abundance_bio_posterior_rpk.RDS")
  #dplyr::select(-c(propreads, absconc))
  
d <- d %>% 
  filter(creek != "4Pad5") %>% 
  filter(creek != "2Brn") %>% 
  filter(species != "Oncorhynchus tshawytscha") %>%  #for now, filter out as rare
  mutate(station_idx = as.integer(station),
         bio = as.integer(bio)) %>% 
  dplyr::select(-station) %>% 
  #rename(station_idx = station) %>% 
  #filter(creek != "2Brn") %>%  #Barnes has too much missing data for the moment
  # mutate(station_idx = case_when(station == "Dn" ~ 1,
  #                                station %in% c("Up", "Up11") ~ 2)) %>% 
  # dplyr::select(-c(station)) %>% 
 # filter(station_idx == 1) %>% 
  arrange(creek) %>% 
  mutate(month = substr(time, 1,2),
         year = substr(time, 3,4)) %>% 
  arrange(year, month) %>% 
  dplyr::select(-c(month, year))


d$creek_idx <- match(d$creek, unique(d$creek))
d$time_idx <- match(d$time, unique(d$time))
d$species_idx <- match(d$species, unique(d$species))

  
d %>% 
  filter(species_idx == 3) %>% 
  ggplot(aes(x = time_idx, y = log(meandnaconc))) +
    geom_point() +
    facet_grid(creek_idx~station_idx)

# construction <- which(d$creek_idx == 4 & d$time_idx %in% c(8:10))
f <- d %>%  #use this to filter out anything you don't want to model
  arrange(time_idx, creek_idx) %>% 
  ungroup()
  
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

stanMod = stan(file = here("Scripts/timeseries_model/timeSeries_multispecies4.stan") ,data = stan_data,
               verbose = FALSE, chains = 3, thin = 1,
               warmup = 300, iter = 700,
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


#plot(stanMod, par = c("mu"))
plot(stanMod, par = c("sigma_eta", "sigma_dna"))
#plot(stanMod, par = c("eta"))
plot(stanMod, par = c("beta_1"))
traceplot(stanMod, par = c("sigma_eta", "sigma_dna"))

#shinystan::launch_shinystan(stanMod)

#dimensions = c(time, station, species, creek) 

resOut <- expand_grid(time_idx = 1:length(unique(f$time_idx)),
                      station_idx = 1:length(unique(f$station_idx)),
                      species_idx = 1:length(unique(f$species_idx)),
                      creek_idx = 1:length(unique(f$creek_idx))) %>% 
  mutate(mean_est = summary(stanMod, par = "mu")$summary[,1],
         ci25 = summary(stanMod, par = "mu")$summary[,5],
         ci75 = summary(stanMod, par = "mu")$summary[,7]) %>% 
  left_join(f) %>% 
  left_join(data.frame(sigma_dna = summary(stanMod, par = "sigma_dna")$summary[,1], 
                       species_idx = 1:1:length(unique(f$species_idx))))

# resOut %>% 
#   filter(creek_idx == 1  & species_idx == 2 & time_idx == 8) %>% 
#   as.data.frame()

#plot
resOut %>% 
  filter(species_idx == 3) %>% 
  ggplot(aes(x = time_idx, y = log(meandnaconc))) +
    geom_point() +
    geom_point(aes(x = time_idx, y = mean_est), color = "red") +
    geom_segment(aes(x = time_idx, xend = time_idx, y = ci25, yend = ci75), color = "red") +
    facet_grid(~station_idx ~creek_idx) +
    ggtitle("Estimated Mean Concentrations w Interquartile Range")


#posterior predictive check:
#95% CI for mu, sigma; w observed data
#dimensions = c(time, station, species, creek) 

getPostPred <- function(time_idx, station_idx, species_idx, creek_idx){
  n = length(unlist(extract(stanMod, par = "mu[1,1,1,1]")))
  f_mu = paste0("mu[",
                time_idx,",",
                station_idx,",",
                species_idx,",",
                creek_idx,"]")
  f_s = paste0("sigma_dna[",species_idx,"]")
  
  postsamples <- rnorm(n, 
        unlist(extract(stanMod, par = f_mu)),
        unlist(extract(stanMod, par = f_s))
  ) 
  
  return(quantile(postsamples, c(.025, .05, .25, .5, .75, .9, .975)))
  
}

ppOut <- expand_grid(time_idx = 1:length(unique(f$time_idx)),
                     station_idx = 1:length(unique(f$station_idx)),
                     species_idx = 1:length(unique(f$species_idx)),
                     creek_idx = 1:length(unique(f$creek_idx)))
ppOut[,5:11] <- 0.0
colnames(ppOut)[5:11] <- paste0("pp_", c("025", "05", 25, 50, 75, 90, 975))

for (i in 1:nrow(ppOut)){
  ppOut[i,5:11] <- as.list(getPostPred(ppOut$time_idx[i], 
                               ppOut$station_idx[i], 
                               ppOut$species_idx[i], 
                               ppOut$creek_idx[i]))
}


(p2 <- f %>%
  mutate(logY = log(meandnaconc)) %>%
  right_join(resOut) %>%
  mutate(month = ifelse(time_idx < 11, time_idx + 2, time_idx - 10),
         month = as.factor(month)) %>% 
  left_join(ppOut) %>% 
  mutate(station = ifelse(station_idx == 1, "Downstream", "Upstream")) %>%
  mutate(creekname = case_when(creek_idx == 1 ~ "Portage",
                               creek_idx == 2 ~ "Chuckanut",
                               creek_idx == 3 ~ "Padden",
                               creek_idx == 4 ~ "Squalicum")) %>%
  mutate(species = case_when(species_idx == 1 ~ "Oncorhynchus clarkii",
                             species_idx == 2 ~ "Oncorhynchus kisutch",
                             species_idx == 3 ~ "Oncorhynchus mykiss")) %>%  ##to fix NAs in unobserved samples
    filter(species_idx == 3) %>%
    ggplot(aes(x = month, y = log(meandnaconc))) +
    geom_point() +
    geom_point(aes(x = month, y = mean_est), color = "red", size = 1.5, alpha = .5) +
    geom_segment(aes(x = month, xend = month, y = pp_05, yend = pp_975), color = "red", size = .3, alpha = .5) +
    geom_segment(aes(x = month, xend = month, y = pp_25, yend = pp_75), color = "red", size = .7, alpha = .5) +
    facet_grid(~creekname ~station) +
    ggtitle("Posterior Predictive Check \n(predicted mean +/- 95 CI)"))
# #ggsave(p2, filename = here("Figures/pp_check.pdf"))
# #ggsave(p2, filename = "pp_check_cutthroat.jpeg")
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


#posterior predictions for observed DNA concentrations
#   a <- extract(stanMod, par = "mu")
#   #b <- extract(stanMod, par = "sigma_dna") %>% unlist()
#   
#   pp <- expand_grid(creek_idx = 1:length(unique(f$creek_idx)),
#               time_idx = 1:length(unique(f$time_idx)),
#               station_idx = 1:length(unique(f$station_idx)),
#               species_idx = 1:length(unique(f$species_idx))) 
#   pp <- pp %>% 
#     bind_cols(matrix(NA, nrow = nrow(pp), ncol = length(a$mu[,1,1,1,1]))) %>% 
#     nest(data = 5:ncol(.))
#   
#   for (i in 1:nrow(pp)){
#     pp$data[i][[1]] <- rnorm(length(a$mu[,1,1,1]), 
#                              a$mu[,pp$time_idx[i],pp$station_idx[i],pp$creek_idx[i]], 
#                              b)
#   }
# 
#   
# pp <- pp %>% unnest(cols = data) %>% 
#   pivot_wider(names_from = station_idx, values_from = data, values_fn = list) 
# 
# map2(pp[1,3], pp[1,4], `-`)
# 
#  unlist(pp[1,3]) - unlist(pp[1,4])
# 
# 
# 
