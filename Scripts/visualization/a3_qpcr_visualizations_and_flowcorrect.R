## NGN Visualizing qPCR results before / after Stan model
# Author: Eily Allan 
# Person running: Eily
# Last modified: 10/19/22 by Eily
# Date of run: 10/19/22 by Eily 

# Overview 
# This script 

# Inputs: 
# 1) 

# Outputs: 
# 1) 

####################################################################
# Set up
####################################################################

# Load packages
library(here)
library(tidyverse)

cut_data_for_stan <- read_csv(here("Output","qpcr","cut_data_for_stan.csv"))
cut_modeled_conc <- readRDS(here("Output","qpcr","cut_modeled_conc.RDS"))
monthlyflow <- read.csv(here("Output","qpcr","monthly_flow.csv"))
closestflow <- read.csv(here("Output","qpcr","closest_flow.csv"))

####################################################################
# qPCR: standard curves before Stan model
####################################################################

cutplot <- cut_data_for_stan %>% 
  filter(str_detect(Sample, "St")) %>% 
  mutate(Quantity = 0) %>% 
  mutate(Quantity = case_when(Sample == "St1" ~ 100000,
                                 Sample == "St2" ~ 10000,
                                 Sample == "St3" ~ 1000,
                                 Sample == "St4" ~ 100,
                                 Sample == "St5" ~ 10,
                                 Sample == "St6" ~ 5,
                                 Sample == "St7" ~ 3,
                                 Sample == "St8" ~ 1,
                              TRUE ~ Quantity)) 
  #filter(Plate == 2) %>% 
  #filter(z == 1) %>% 
  ggplot(cutplot, aes(x= log(Quantity), y = Ct, color = as.factor(Plate))) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()


intercepts <- cutplot %>% 
  filter(Quantity == 10) %>% 
  filter( !is.na(Ct))

### GET EFFICIENCIES TO REPORT RANGE

####################################################################
# qPCR: example of stan model output standard curve
####################################################################

# qMod <- readRDS("/Users/elizabethandruszkiewicz/Desktop/20221129_model_output/cut_qMod_out.RDS")
# 
# PLATE = 2
# pp_beta0 <- extract(qMod, "beta_std_curve_0")$beta_std_curve_0[,PLATE]
# pp_beta1 <- extract(qMod, "beta_std_curve_1")$beta_std_curve_1[,PLATE]
# # pp_sigma <- extract(qMod, "sigma_std_curve")$sigma_std_curve
# pp_gamma_0 <- extract(qMod, "gamma_0")$gamma_0 #[,PLATE]
# pp_gamma_1 <- extract(qMod, "gamma_1")$gamma_1[,PLATE]
# 
# conc <- runif(length(pp_beta0), 0, 6)
# mu <- pp_beta0 + pp_beta1*conc
# sigma <- exp(pp_gamma_0 + pp_gamma_1*conc)
# y<-NA
# for(i in 1:length(pp_beta0)){y[i] <- rnorm(1, mean = mu[i], sd = sigma[i])}
# 
# (plotstd <- data.frame(y, conc) %>%
#     ggplot(aes(x = conc, y = y)) +
#     geom_point(alpha = .1) +
#     ggtitle(paste("Plate ", PLATE)))
# 
# (p1 <- plotstd +
#     geom_point(data = qPCRdata %>% filter(Plate == PLATE & z == 1), aes(x = log10(conc), y = Ct), color = "red"))
# #geom_point(data = qPCRdata %>% filter(Plate == PLATE & z == 0), aes(x = log10(conc), y = 20), color = "blue") +

# ggsave(p1, file = here("Output/SupplementalFigures/qPCR_calibration_supplemental.png"))    

####################################################################
# qPCR: cutthroat modeled concentration 
####################################################################

cut_modeled_conc %>% 
  filter(station != "Up5") %>% 
  mutate(station = case_when(station == "Up11" ~ "Up",
                             TRUE ~ station)) %>% 
  separate(time, into = c("month","year"), sep = 2) %>% 
  unite(newtime, c(year,month), sep="-") %>% 
  mutate(creek = case_when(creek == "1Prt" ~ "Portage",
                           creek == "2Brn" ~ "Barnes",
                           creek == "3Chk" ~ "Chuckanut",
                           creek == "4Pad" ~ "Padden",
                           creek == "5Sqm" ~ "Squalicum",
                           TRUE ~ creek)) %>% 
  mutate(station = case_when(station == "Dn" ~ "Down",
                           TRUE ~ station)) %>% 
  ggplot(aes(x = newtime, y = log(mean_concentration_est))) +
  geom_point() +
  geom_segment(aes(x = newtime, xend = newtime, y = log(ci25_concentration_est), yend = log(ci75_concentration_est))) +
  facet_grid(~creek~station) +
  labs(x="Date (YY-MM)", y= "Log copies/L water", title = "Cutthroat Trout") + 
  theme_bw() + 
  scale_x_discrete(guide = guide_axis(angle = -45))

#ggsave(here("Output","Figures", "modeled_cut_qpcr.png"))


####################################################################
# qPCR: cutthroat modeled concentration -- up/down on same plot
####################################################################

cut_modeled_conc %>% 
  filter(station != "Up5") %>% 
  mutate(station = case_when(station == "Up11" ~ "Up",
                             TRUE ~ station)) %>% 
  separate(time, into = c("month","year"), sep = 2) %>% 
  unite(newtime, c(year,month), sep="-") %>% 
  mutate(creek = case_when(creek == "1Prt" ~ "Portage",
                           creek == "2Brn" ~ "Barnes",
                           creek == "3Chk" ~ "Chuckanut",
                           creek == "4Pad" ~ "Padden",
                           creek == "5Sqm" ~ "Squalicum",
                           TRUE ~ creek)) %>% 
  mutate(station = case_when(station == "Dn" ~ "Down",
                             TRUE ~ station)) %>% 
  mutate(facetorder = factor(creek, levels=c('Padden','Portage','Chuckanut','Squalicum', 'Barnes')))%>% 
  ggplot(aes(x = newtime, y = log(mean_concentration_est), color=station)) +
  geom_point() +
  geom_segment(aes(x = newtime, xend = newtime, y = log(ci25_concentration_est), yend = log(ci75_concentration_est))) +
  facet_grid(rows=vars(facetorder)) +
  labs(x="Date (YY-MM)", y= "Log copies/L water", title = "Cutthroat Trout", color="Station") + 
  theme_bw() + 
  scale_x_discrete(guide = guide_axis(angle = -45)) +
  scale_color_manual(values=c("darkgrey", "black"))

#ggsave(here("Output","Figures", "20221129_modeled_cut_qpcr_updown.png"))


# ####################################################################
# # qPCR: cutthroat modeled concentration MULTIPLIED BY FLOW
# ####################################################################
# 
# ## for now, use padden for all of them 
# closestflow2 <- closestflow %>% 
#   filter(str_detect(Sample, "Dn")) %>% 
#   filter(creek=="Padden") %>% 
#   select(c(yearplot,monthplot,dayplot,flow_m3s)) %>% 
#   distinct() %>% 
#   slice(rep(1:n(), 3)) %>% 
#   mutate(creek=rep(c("Padden","Chuckanut","Squalicum"), each=12))
# 
# flowmerge <- 
#   closestflow2 %>% 
#   #closestflow %>% 
#   dplyr::select(c(creek, yearplot, monthplot, flow_m3s)) %>%
#   #dplyr::select(c(creek, yearplot, monthplot, flow_m3s)) %>%
#   #monthlyflow %>% 
#   #dplyr::select(c(creek, yearplot, monthplot, meanmonthflow)) %>%
#   separate(yearplot, into = c("xx","year"), sep = 2) %>% 
#   dplyr::rename(month = monthplot) %>% 
#   mutate(year = as.numeric(year)) %>% 
#   mutate(month = as.numeric(month)) %>% 
#   dplyr::select(-xx) %>% 
#   unite(creekyrmo, c(creek,year,month))
#   
# 
# cut_flow_corrected <- cut_modeled_conc %>% 
#   filter(station != "Up5") %>% 
#   mutate(station = case_when(station == "Up11" ~ "Up",
#                              TRUE ~ station)) %>% 
#   separate(time, into = c("month","year"), sep = 2) %>% 
#   unite(newtime, c(year,month), sep="-", remove=FALSE) %>% 
#   mutate(year = as.numeric(year)) %>% 
#   mutate(month = as.numeric(month)) %>%
#   mutate(creek = case_when(creek == "1Prt" ~ "Portage",
#                            creek == "2Brn" ~ "Barnes",
#                            creek == "3Chk" ~ "Chuckanut",
#                            creek == "4Pad" ~ "Padden",
#                            creek == "5Sqm" ~ "Squalicum",
#                            TRUE ~ creek)) %>% 
#   mutate(station = case_when(station == "Dn" ~ "Down",
#                              TRUE ~ station)) %>% 
#   unite(creekyrmo, c(creek,year,month), remove=FALSE) %>% 
#   left_join(flowmerge, by="creekyrmo") %>% 
#   #filter(!is.na(meanmonthflow)) %>% 
#   #mutate(mean_conc_flow_correct = mean_concentration_est*meanmonthflow) %>% 
#   #mutate(ci25_flow_correct = ci25_concentration_est*meanmonthflow) %>% 
#   #mutate(ci75_conc_flow_correct = ci25_concentration_est*meanmonthflow)
#   filter(!is.na(flow_m3s)) %>%
#   mutate(mean_conc_flow_correct = mean_concentration_est*flow_m3s) %>% 
#   mutate(ci25_flow_correct = ci25_concentration_est*flow_m3s) %>% 
#   mutate(ci75_conc_flow_correct = ci25_concentration_est*flow_m3s)
#  
# #write_rds(cut_flow_corrected, here("Output","qPCR", "modeled_cut_qpcr_flowcorrected_monthlyavg.RDS"))
# #write_rds(cut_flow_corrected, here("Output","qPCR", "modeled_cut_qpcr_flowcorrected_allpadden.RDS"))
# 
# ggplot(cut_flow_corrected, aes(x=newtime)) +
#   geom_point(aes(y=log10(mean_concentration_est))) + 
#   geom_point(aes(y=log10(mean_conc_flow_correct)), color="blue") + 
#   scale_y_continuous(name = "Log10 copies/L water",
#     sec.axis = sec_axis(~.*1, name="Log10 copies/s (flow corrected)")) +
#   facet_grid(~creek~station) +
#   labs(x="Date (YY-MM)", y= "Log10 copies/L water", title = "Cutthroat Trout") + 
#   theme_bw() + 
#   theme(axis.title.y.right = element_text(color = "blue")) + 
#   scale_x_discrete(guide = guide_axis(angle = -45))
# 
# #ggsave(here("Output","Figures", "modeled_cut_qpcr_updown_flowcorrected_allpadden.png"))
# 
