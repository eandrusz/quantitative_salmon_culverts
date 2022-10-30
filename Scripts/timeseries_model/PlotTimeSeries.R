library(tidyverse)
library(here)

####
stanMod <- readRDS(here("Scripts/timeseries_model/modelFit_20221028.RDS"))

f <- stanMod[3][[1]]  #data in
s <- stanMod[[2]] #fitted model itself

####


resOut <- expand_grid(time_idx = 1:length(unique(f$time_idx)),
                      station_idx = 1:length(unique(f$station_idx)),
                      species_idx = 1:length(unique(f$species_idx)),
                      creek_idx = 1:length(unique(f$creek_idx))) %>% 
  mutate(mean_est = rstan::summary(s, par = "mu")$summary[,1],
         ci25 = rstan::summary(s, par = "mu")$summary[,5],
         ci75 = rstan::summary(s, par = "mu")$summary[,7]) %>% 
  left_join(f) %>% 
  left_join(data.frame(sigma_dna = rstan::summary(s, par = "sigma_dna")$summary[,1], 
                       species_idx = 1:1:length(unique(f$species_idx))))


#plot means
resOut %>% 
  filter(species_idx == 8) %>% 
  ggplot(aes(x = time_idx, y = log(meandnaconc))) +
  geom_point() +
  geom_point(aes(x = time_idx, y = mean_est), color = "red") +
  geom_segment(aes(x = time_idx, xend = time_idx, y = ci25, yend = ci75), color = "red") +
  facet_grid(~station_idx ~creek_idx) +
  ggtitle("Estimated Mean Concentrations w Interquartile Range")


#posterior predictive check:
#95% CI for mu, sigma; w observed data
#dimensions = c(time, station, species, creek) 

getPostPred <- function(modname, time_idx, station_idx, species_idx, creek_idx){
  s = modname
  n = length(unlist(rstan::extract(s, par = "mu[1,1,1,1]")))
  f_mu = paste0("mu[",
                time_idx,",",
                station_idx,",",
                species_idx,",",
                creek_idx,"]")
  f_s = paste0("sigma_dna[",species_idx,"]")
  
  postsamples <- rnorm(n, 
                       unlist(rstan::extract(s, par = f_mu)),
                       unlist(rstan::extract(s, par = f_s))
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
  ppOut[i,5:11] <- as.list(getPostPred(s,
                                       ppOut$time_idx[i], 
                                       ppOut$station_idx[i], 
                                       ppOut$species_idx[i], 
                                       ppOut$creek_idx[i]))
}

for (i in 1:length(unique(f$species))){
  assign(paste0("p",i), value = 
    f %>%
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
       # mutate(species = case_when(species_idx == 1 ~ "Oncorhynchus clarkii",
                                  # species_idx == 2 ~ "Oncorhynchus kisutch",
                                  # species_idx == 3 ~ "Oncorhynchus mykiss")) %>%  ##to fix NAs in unobserved samples
       filter(species_idx == i) %>%
       ggplot(aes(x = month, y = log(meandnaconc))) +
       geom_point() +
       geom_point(aes(x = month, y = mean_est), color = "red", size = 1.5, alpha = .5) +
       geom_segment(aes(x = month, xend = month, y = pp_025, yend = pp_975), color = "red", size = .3, alpha = .5) +
       geom_segment(aes(x = month, xend = month, y = pp_25, yend = pp_75), color = "red", size = .7, alpha = .5) +
       facet_grid(~creekname ~station) +
       ggtitle(paste0(unique(f$species)[i], "\nPosterior Predictive Check \n(predicted mean +/- 95 CI)"))
  )
}
#visualize
p1
p2
p3
p4
p5
p6

f %>% 
  filter(species == "Oncorhynchus tshawytscha") %>% 
  filter(creek == "3Chk") %>% 
  as.data.frame()

##One way to plot effect of culverts:
f %>%
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
  filter(species_idx == 1) %>%
  ggplot(aes(x = as.numeric(month), y = log(meandnaconc))) +
  # geom_point() +
  geom_point(aes(x = as.numeric(month), y = mean_est, color = station), size = 1.5, alpha = .5) +
  geom_smooth(aes(x = as.numeric(month), y = mean_est, color = station), size = 1.5, alpha = .5, se = F) +
  #geom_segment(aes(x = as.numeric(month), xend = as.numeric(month), y = pp_05, yend = pp_975, color = station),  size = .3, alpha = .5) +
  geom_segment(aes(x = as.numeric(month), xend = as.numeric(month), y = pp_25, yend = pp_75, color = station),  size = .3, alpha = .5) +
  facet_grid(~creekname ) +
  ggtitle("Posterior Predictive Check \n(predicted mean +/- 95 CI)")





gamma_means <- expand_grid(time_idx = 2:(length(unique(f$time_idx))),
                           species_idx = 1:length(unique(f$species_idx)),
                           #creek_idx = 1:length(unique(f$creek)),
                           #station_idx = 1:length(unique(f$station_idx)),
                           construction_idx = 1:length(unique(f$construction_idx))
                           ) %>% 
  mutate(gamma_mean = NA,
         gamma_025 = NA,
         gamma_975 = NA)

for (i in 1:nrow(gamma_means)){
  a <- paste0("gamma[", 
              gamma_means$time_idx[i]-1, ",",
              gamma_means$species_idx[i], ",",
              #gamma_means$creek_idx[i], ",",
              #gamma_means$station_idx[i], ",",
              gamma_means$construction_idx[i], "]")
  
  gamma_means$gamma_mean[i] <- rstan::summary(s, par = a)$summary[,1]
  gamma_means$gamma_025[i] <- rstan::summary(s, par = a)$summary[,4]
  gamma_means$gamma_975[i] <- rstan::summary(s, par = a)$summary[,8]
}

# gamma_means %>% 
#   mutate(month = time_idx) %>% 
#   filter(construction_idx == 2 & time_idx > 0 & species_idx < 4) %>% 
#   ggplot(aes(x = as.factor(time_idx), 
#              y = gamma_mean)) +
#     geom_point() +
#     geom_segment(aes(x = time_idx, xend = time_idx, y = gamma_025, yend = gamma_975)) +
#     facet_grid(~species_idx) +
#     xlab("time_idx") + ylab("Gamma")
  
gamma_means <- f %>% select(time, time_idx) %>%
  distinct() %>%
  mutate(month = ifelse(time_idx < 11, time_idx + 2, time_idx-10)) %>% 
  left_join(gamma_means)


gamma_means %>% 
  filter(construction_idx == 2, species_idx < 6 & month %in% c(1:12)) %>% 
  ggplot(aes(x = as.factor(month), y = gamma_mean)) +
    geom_point() +
    geom_segment(aes(x = as.factor(month), xend = as.factor(month), y = gamma_025, yend = gamma_975)) +
    facet_grid(.~species_idx)


