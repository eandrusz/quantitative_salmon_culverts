library(tidyverse)
library(here)

####
stanMod <- readRDS(here("Scripts/timeseries_model/modelFit_20221009.RDS"))

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
  filter(species_idx == 1) %>% 
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
       mutate(species = case_when(species_idx == 1 ~ "Oncorhynchus clarkii",
                                  species_idx == 2 ~ "Oncorhynchus kisutch",
                                  species_idx == 3 ~ "Oncorhynchus mykiss")) %>%  ##to fix NAs in unobserved samples
       filter(species_idx == i) %>%
       ggplot(aes(x = month, y = log(meandnaconc))) +
       geom_point() +
       geom_point(aes(x = month, y = mean_est), color = "red", size = 1.5, alpha = .5) +
       geom_segment(aes(x = month, xend = month, y = pp_05, yend = pp_975), color = "red", size = .3, alpha = .5) +
       geom_segment(aes(x = month, xend = month, y = pp_25, yend = pp_75), color = "red", size = .7, alpha = .5) +
       facet_grid(~creekname ~station) +
       ggtitle(paste0(unique(f$species)[i], "\nPosterior Predictive Check \n(predicted mean +/- 95 CI)"))
  )
}
#visualize
p1
p2
p3


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
  filter(species_idx == 3) %>%
  ggplot(aes(x = as.numeric(month), y = log(meandnaconc))) +
  # geom_point() +
  geom_point(aes(x = as.numeric(month), y = mean_est, color = station), size = 1.5, alpha = .5) +
  geom_smooth(aes(x = as.numeric(month), y = mean_est, color = station), size = 1.5, alpha = .5, se = F) +
  #geom_segment(aes(x = as.numeric(month), xend = as.numeric(month), y = pp_05, yend = pp_975, color = station),  size = .3, alpha = .5) +
  geom_segment(aes(x = as.numeric(month), xend = as.numeric(month), y = pp_25, yend = pp_75, color = station),  size = .3, alpha = .5) +
  facet_grid(~creekname ) +
  ggtitle("Posterior Predictive Check \n(predicted mean +/- 95 CI)")


###CONSTRUCTION Causal Analysis Plot:



fit_nc <- readRDS(here("Scripts/timeseries_model/modelFit_20221009_noConstruction.RDS"))

f1 <- fit_nc[3][[1]]  #data in
s1 <- fit_nc[[2]] #fitted model itself

resOutNoConstr <- expand_grid(time_idx = 1:length(unique(f1$time_idx)),
                              station_idx = 1:length(unique(f1$station_idx)),
                              species_idx = 1:length(unique(f1$species_idx)),
                              creek_idx = 1:length(unique(f1$creek_idx))) %>% 
  mutate(mean_est = rstan::summary(s1, par = "mu")$summary[,1],
         ci25 = rstan::summary(s1, par = "mu")$summary[,5],
         ci75 = rstan::summary(s1, par = "mu")$summary[,7]) %>% 
  left_join(f1) %>% 
  left_join(data.frame(sigma_dna = rstan::summary(s1, par = "sigma_dna")$summary[,1], 
                       species_idx = 1:1:length(unique(f1$species_idx))))

ppOutNoConstr <- expand_grid(time_idx = 1:length(unique(f1$time_idx)),
                             station_idx = 1:length(unique(f1$station_idx)),
                             species_idx = 1:length(unique(f1$species_idx)),
                             creek_idx = 1:length(unique(f1$creek_idx)))
ppOutNoConstr[,5:11] <- 0.0
colnames(ppOutNoConstr)[5:11] <- paste0("pp_", c("025", "05", 25, 50, 75, 90, 975))

for (i in 1:nrow(ppOutNoConstr)){
  ppOutNoConstr[i,5:11] <- as.list(getPostPred(s1, 
                                               ppOutNoConstr$time_idx[i], 
                                               ppOutNoConstr$station_idx[i], 
                                               ppOutNoConstr$species_idx[i], 
                                               ppOutNoConstr$creek_idx[i]))
}


expected <- f1 %>%
  mutate(logY = log(meandnaconc)) %>%
  right_join(resOutNoConstr) %>%
  mutate(month = ifelse(time_idx < 11, time_idx + 2, time_idx - 10),
         month = as.factor(month)) %>% 
  left_join(ppOutNoConstr) %>% 
  mutate(station = ifelse(station_idx == 1, "Downstream", "Upstream")) %>%
  mutate(creekname = case_when(creek_idx == 1 ~ "Portage",
                               creek_idx == 2 ~ "Chuckanut",
                               creek_idx == 3 ~ "Padden",
                               creek_idx == 4 ~ "Squalicum")) %>%
  mutate(species = case_when(species_idx == 1 ~ "Oncorhynchus clarkii",
                             species_idx == 2 ~ "Oncorhynchus kisutch",
                             species_idx == 3 ~ "Oncorhynchus mykiss")) %>% 
  filter(creek_idx == 3 & month %in% c(8,9,10,11,12))

observed <- f %>% 
  arrange(time_idx, creek_idx) %>% 
  ungroup() %>% 
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
                             species_idx == 3 ~ "Oncorhynchus mykiss")) %>% 
  filter(creek_idx == 3 & month %in% c(8,9,10,11,12))

observed <- observed %>% 
  mutate(state = "observed")
expected <- expected %>% 
  mutate(state = "expected")

(ConstrEffect <- 
  observed %>% 
    bind_rows(expected) %>% 
    #filter(species_idx == 3) %>%
    ggplot(aes(x = month, y = log(meandnaconc), color = state)) +
    # geom_jitter(width = 0.05) +
    geom_point(aes(x = month, y = mean_est, color = state), size = 1.5, alpha = .5) +
    #geom_segment(aes(x = month, xend = month, y = pp_05, yend = pp_975, color = state), size = .3, alpha = .5) +
    geom_segment(aes(x = month, xend = month, y = pp_25, yend = pp_75, color = state), size = .7, alpha = .5) +
    facet_grid(~species ~station))

ggsave(ConstrEffect)
