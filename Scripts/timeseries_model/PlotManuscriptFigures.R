library(tidyverse)
library(here)

####
#read in time series model fit
stanMod <- readRDS(here("Scripts/timeseries_model/modelFit_20221102.RDS"))

f <- stanMod[3][[1]]  #data in
s <- stanMod[[2]] #fitted model itself

####

# summarize expected values at each time/station/species/creek
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
# resOut %>% 
#   filter(species_idx == 8) %>% 
#   ggplot(aes(x = time_idx, y = log(meandnaconc))) +
#   geom_point() +
#   geom_point(aes(x = time_idx, y = mean_est), color = "red") +
#   geom_segment(aes(x = time_idx, xend = time_idx, y = ci25, yend = ci75), color = "red") +
#   facet_grid(~station_idx ~creek_idx) +
#   ggtitle("Estimated Mean Concentrations w Interquartile Range")


#posterior predictive check:
#95% CI for mu, sigma; w observed data
#dimensions = c(time, station, species, creek) 

  #sample posterior for expected value +/- observation variance
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

    
    #Make combined-species plot of trends over time
    
    (multispeciesTrends <- f %>%
      mutate(logY = log(meandnaconc)) %>%
      right_join(resOut) %>%
      mutate(month = ifelse(time_idx < 11, time_idx + 2, time_idx - 10),
             # month = as.factor(month)
      ) %>% 
      left_join(ppOut) %>% 
      drop_na() %>% 
      mutate(station = ifelse(station_idx == 1, "Downstream", "Upstream")) %>%
      # mutate(creekname = case_when(creek_idx == 1 ~ "Portage",
      #                              creek_idx == 2 ~ "Chuckanut",
      #                              creek_idx == 3 ~ "Padden",
      #                              creek_idx == 4 ~ "Squalicum")) %>%
      # mutate(species = case_when(species_idx == 1 ~ "Oncorhynchus clarkii",
      # species_idx == 2 ~ "Oncorhynchus kisutch",
      # species_idx == 3 ~ "Oncorhynchus mykiss")) %>%  ##to fix NAs in unobserved samples
      ungroup() %>% 
      # filter(species_idx == i) %>%
      ggplot(aes(x = month, y = log(meandnaconc), color = station)) +
      geom_point(alpha = .2) +
      geom_point(aes(x = month, y = mean_est, color = station), size = 1.5, alpha = .5) +
      geom_smooth(aes(x = month, y = mean_est, color = station), se = F, size = 1, alpha = .5, span = .2) +
      geom_segment(aes(x = month, xend = month, y = pp_025, yend = pp_975, color = station), size = .3, alpha = .2) +
      geom_segment(aes(x = month, xend = month, y = pp_25, yend = pp_75, color = station), size = .7, alpha = .4) +
      facet_grid(species~creek) +
      xlab("Month") + ylab("Log DNA Concentration (copies/L)") +
        labs(color='') +
      # ggtitle(paste0(unique(f$species)[i], "\n(predicted mean +/- 95 CI)")) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            legend.position = "bottom",
            legend.box.spacing = unit(0, "mm")) +
      scale_x_continuous(breaks = 1:12,
                         labels = as.character(1:12)))
    
    
    
    
    #make species-specific plots and save as ggplot objects
    for (i in 1:length(unique(f$species))){
      assign(paste0("p",i), value = 
               f %>%
               mutate(logY = log(meandnaconc)) %>%
               right_join(resOut) %>%
               mutate(month = ifelse(time_idx < 11, time_idx + 2, time_idx - 10),
                      # month = as.factor(month)
               ) %>% 
               left_join(ppOut) %>% 
               mutate(station = ifelse(station_idx == 1, "Downstream", "Upstream")) %>%
               mutate(creekname = case_when(creek_idx == 1 ~ "Portage",
                                            creek_idx == 2 ~ "Chuckanut",
                                            creek_idx == 3 ~ "Padden",
                                            creek_idx == 4 ~ "Squalicum")) %>%
               # mutate(species = case_when(species_idx == 1 ~ "Oncorhynchus clarkii",
               # species_idx == 2 ~ "Oncorhynchus kisutch",
               # species_idx == 3 ~ "Oncorhynchus mykiss")) %>%  ##to fix NAs in unobserved samples
               ungroup() %>% 
               filter(species_idx == i) %>%
               ggplot(aes(x = month, y = log(meandnaconc), color = station)) +
               geom_point(alpha = .2) +
               geom_point(aes(x = month, y = mean_est, color = station), size = 1.5, alpha = .5) +
               geom_smooth(aes(x = month, y = mean_est, color = station), se = F, size = 1, alpha = .5, span = .2) +
               geom_segment(aes(x = month, xend = month, y = pp_025, yend = pp_975, color = station), size = .3, alpha = .2) +
               geom_segment(aes(x = month, xend = month, y = pp_25, yend = pp_75, color = station), size = .7, alpha = .4) +
               facet_grid(~creekname) +
               xlab("Month") + ylab("Log DNA Concentration (copies/L)") +
               ggtitle(paste0(unique(f$species)[i], "\n(predicted mean +/- 95 CI)")) +
               theme_bw() +
               theme(panel.grid.minor = element_blank(),
                     panel.grid.major.x = element_blank()) +
               scale_x_continuous(breaks = 1:12,
                                  labels = as.character(1:12))
             
      )
    }
#visualize
p1
p2
p3
# p4
# p5
# p6

ggsave(p1, file = here("Output/Figures/Oclarkii_timeseries.png"))
ggsave(p2, file = here("Output/Figures/Okisutch_timeseries.png"))
ggsave(p3, file = here("Output/Figures/Omykiss_timeseries.png"))

## Calculate and plot the effect of culvert-removal construction (Padden Creek)

#recover posterior estimates of gamma term
    gamma_est <- expand_grid(time_idx = 2:(length(unique(f$time_idx))),
                               species_idx = 1:length(unique(f$species_idx)),
                               #creek_idx = 1:length(unique(f$creek)),
                               #station_idx = 1:length(unique(f$station_idx)),
                               construction_idx = 1:length(unique(f$construction_idx))
                               ) %>% 
      mutate(gamma_mean = NA,
             gamma_025 = NA,
             gamma_975 = NA)
    
    for (i in 1:nrow(gamma_est)){
      a <- paste0("gamma[", 
                  gamma_est$time_idx[i]-1, ",",
                  gamma_est$species_idx[i], ",",
                  #gamma_est$creek_idx[i], ",",
                  #gamma_est$station_idx[i], ",",
                  gamma_est$construction_idx[i], "]")
      
      gamma_est$gamma_mean[i] <- rstan::summary(s, par = a)$summary[,1]
      gamma_est$gamma_025[i] <- rstan::summary(s, par = a)$summary[,4]
      gamma_est$gamma_975[i] <- rstan::summary(s, par = a)$summary[,8]
    }

  # link to information about sampling dates  
  gamma_est <- f %>% select(newtime, time_idx) %>%
    distinct() %>%
    mutate(month = ifelse(time_idx < 11, time_idx + 2, time_idx-10)) %>% 
    left_join(gamma_est)

  Nspp <- length(unique(f$species))
  

gamma_est %>% 
  left_join(tibble(species_idx = 1:Nspp, jitter = (1:Nspp - median(1:Nspp))*0.07)) %>% 
  left_join(f %>% dplyr::select(species_idx, species) %>% distinct()) %>% 
  mutate(jitter_time = time_idx + (jitter)) %>% 
  filter(gamma_mean != 0) %>% 
  # filter(construction_idx == 2, species_idx < 4 & month %in% c(1:12)) %>% 
  ggplot(aes(x = jitter_time, y = gamma_mean, color = species)) +
    geom_point() +
    geom_segment(aes(x = jitter_time, xend = jitter_time, y = gamma_025, yend = gamma_975)) +
    scale_x_continuous(breaks = c(7:12),
                     labels = c(9:12,1,2)) +
    xlab("Month") + ylab("Gamma")

colorPal <- c("#8b62ca",
              "#55b04a",
              "#c760a3",
              "#9cb24e",
              "#688bcd",
              "#c7a03d",
              "#4bb092",
              "#6d732f")

#percent change due to gamma, as expressed by:
# gamma_mean/(mean_est-gamma_mean)
(construction_effect <- gamma_est %>% 
  filter(gamma_mean!= 0) %>% 
  left_join(tibble(species_idx = 1:Nspp, jitter = (1:Nspp - median(1:Nspp))*0.07)) %>% 
  left_join(f %>% dplyr::select(species_idx, species) %>% distinct()) %>% 
  mutate(jitter_time = time_idx + (jitter)) %>% 
  left_join(resOut %>% filter(creek_idx == 3, station_idx == 1) %>% dplyr::select(time_idx, species_idx, mean_est)) %>% distinct() %>% 
  mutate(gamma_mean_pct = gamma_mean/(mean_est-gamma_mean),
         gamma_025_pct = gamma_025/(mean_est-gamma_mean),
         gamma_975_pct = gamma_975/(mean_est-gamma_mean)) %>% 
  ggplot(aes(x = jitter_time, y = gamma_mean_pct, color = species)) +
  geom_point() +
  geom_segment(aes(x = jitter_time, xend = jitter_time, y = gamma_025_pct, yend = gamma_975_pct)) +
  scale_x_continuous(breaks = c(7:12),
                     labels = c(9:12,1,2)) +
  xlab("Month") + ylab("Percent Change due to Culvert Removal") +
  scale_color_discrete(type = colorPal))
# ggsave(construction_effect, file = here("Output/Figures/Construction_effect.png"))


### Effect of culvert -- difference between eta terms above vs. below


eta <- expand_grid(time_idx = 2:(length(unique(f$time_idx))),
                         station_idx = 1:length(unique(f$station_idx)),
                         species_idx = 1:length(unique(f$species_idx)),
                         creek_idx = 1:length(unique(f$creek))
                   ) %>% 
  mutate(eta_mean = NA,
         eta_025 = NA,
         eta_975 = NA)

for (i in 1:nrow(eta)){
  a <- paste0("eta[", 
              eta$time_idx[i]-1, ",",
              eta$station_idx[i], ",",
              eta$species_idx[i], ",",
              eta$creek_idx[i],"]")
  
  eta$eta_mean[i] <- rstan::summary(s, par = a)$summary[,1]
  eta$eta_025[i] <- rstan::summary(s, par = a)$summary[,4]
  eta$eta_975[i] <- rstan::summary(s, par = a)$summary[,8]
}

  # link to information about sampling dates  
  eta <- f %>% select(newtime, time_idx) %>%
    distinct() %>%
    mutate(month = ifelse(time_idx < 11, time_idx + 2, time_idx-10)) %>% 
    left_join(eta) %>% 
    left_join(resOut %>% dplyr::select(time_idx, species_idx, creek_idx, station_idx, mean_est)) %>% distinct() %>% 
    mutate(eta_mean_pct = eta_mean/(mean_est-eta_mean),
           eta_025_pct = eta_025/(mean_est-eta_mean),
           eta_975_pct = eta_975/(mean_est-eta_mean))
  
  eta_downstream <- eta %>% 
    filter(station_idx==1)
  
  eta_upstream <- eta %>% 
    filter(station_idx==2)
  
  eta_diff <- eta_downstream %>% 
    dplyr::select(newtime, time_idx, month, species_idx, creek_idx) %>% 
    distinct() %>% 
    mutate(eta_mean_diff = eta_upstream$eta_mean - eta_downstream$eta_mean,
           eta_mean_diff_pct = eta_upstream$eta_mean_pct - eta_downstream$eta_mean_pct) %>% 
    left_join(f %>% dplyr::select(species_idx, species, creek_idx, creek))
  
  (culvert_effect <-  eta_diff %>% 
      ggplot(aes(x = as.factor(month), y = eta_mean_diff_pct)) +
      geom_boxplot()) 
  
  eta_diff %>% 
    group_by(month) %>% 
    summarise(mean(eta_mean_diff_pct))
      
# ggsave(culvert_effect, file = here("Output/Figures/culvert_effect.png"))

  (culvert_effect_creek_species <-  eta_diff %>% 
      ggplot(aes(x = as.factor(month), y = eta_mean_diff_pct)) +
      geom_point() +
      facet_grid(species~creek))
  
  # ggsave(culvert_effect_creek_species, file = here("Output/Figures/culvert_effect_creek_species.png"))
  
  