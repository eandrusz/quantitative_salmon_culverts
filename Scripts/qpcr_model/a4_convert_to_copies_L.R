## NGN Running stan model for QPCR data
# Author: Eily Allan but really Ryan Kelly
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
library(tidyverse)
library(rstan)
library(here)
library(readxl)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Call functions for Stan model
source(here("Scripts","functions", "calibrate_qPCR.R")) 

# Find metadata and all data for cutthroat and coho qpcr
#### Metadata file is for both cutthroat and coho 
qPCRmeta <- here("Input","qpcr","qPCR_samples_ALL.xlsx")

#### Cutthroat and coho file directory 
cut_files <- here("Input","qpcr","Results","CUT")
coho_files <- here("Input","qpcr","Results","COHO")
n_cut <- length(list.files(path = cut_files, pattern = "*data", recursive = T, full.names = T))
n_coho <- length(list.files(path = coho_files, pattern = "*data", recursive = T, full.names = T))


####################################################################
# Write function to read in files and format for Stan model
####################################################################

## cutthroat 
cutcleandata <- list()
for (i in 1:n_cut) {
  plate_num <- i
  cutcleandata[[i]] <- input_qPCR_data(cut_files, qPCRmeta, plate_num)
}
cut_data_for_stan <- do.call(rbind, cutcleandata)  

check <- cut_data_for_stan %>% 
  separate(Sample, into=c("time","creek","stn","bio")) %>% 
  unite(creekstnbio, c(creek,stn,bio)) 

check2 <- check %>% 
  group_by(creekstnbio) %>% 
  summarize(n=n())

# write.csv(cut_data_for_stan, here("Output","qpcr","cut_data_for_stan.csv"), row.names=FALSE)


## coho
cohocleandata <- list()
for (i in 1:n_coho) {
  plate_num <- i
  cohocleandata[[i]] <- input_qPCR_data(coho_files, qPCRmeta, plate_num)
}
coho_data_for_stan <- do.call(rbind, cohocleandata)  

check <- coho_data_for_stan %>% 
  separate(Sample, into=c("time","creek","stn","bio")) %>% 
  unite(creekstnbio, c(creek,stn,bio)) 

check2 <- check %>% 
  group_by(creekstnbio) %>% 
  summarize(n=n())

# write.csv(coho_data_for_stan, here("Output","qpcr","coho_data_for_stan.csv"), row.names=FALSE)


####################################################################
# Run Stan model for cutthroat and coho 
####################################################################

#cut_data_for_stan <- read_csv(here("Output","qpcr","cut_data_for_stan.csv"))
#coho_data_for_stan <- read_csv(here("Output","qpcr","coho_data_for_stan.csv"))


### RUN MODEL FOR CUTTHROAT 
cut_qMod_out <- run_qPCR_model(here("Output","qpcr","cut_data_for_stan.csv"),
                               here("Scripts", "functions", "qPCR_calibration_enchilada.stan"))
#write_rds(cut_qMod_out, "/Users/elizabethandruszkiewicz/Desktop/20221019-ngn-model-output/cut_qMod_out.RDS")
cut_modeled_conc <- cut_qMod_out$results_qPCR
write_rds(cut_modeled_conc, here("Output","qpcr","cut_modeled_conc.RDS"))

### RUN MODEL FOR COHO 
coho_qMod_out <- run_qPCR_model(here("Output","qpcr","coho_data_for_stan.csv"),
                                here("Scripts", "functions", "qPCR_calibration_enchilada.stan"))
write_rds(coho_qMod_out, "/Users/elizabethandruszkiewicz/Desktop/20221019-ngn-model-output/coho_qMod_out.RDS")
coho_modeled_conc <- coho_qMod_out$results_qPCR
write_rds(coho_modeled_conc, here("Output","qpcr","coho_modeled_conc.RDS"))


####################################################################
# Check output to make sure it looks reasonable  
####################################################################

checkcut <- cut_modeled_conc %>% 
  unite(Sample, c(time,creek,station,biorep), sep =".") %>% 
  select(-c(Plate, dilution, Adj_Vol)) %>% 
  left_join(cut_data_for_stan, by = "Sample") %>% 
  separate(Sample, into=c("time","creek","station","biorep"))

ggplot(checkcut, aes(y=Ct, x=mean_concentration_est, color=Plate)) +
  geom_point() +
  scale_x_log10() +
  scale_color_continuous(type = "viridis") +
  facet_wrap(~dilution)

####################################################################
# Plot cutthroat environmental samples 
####################################################################

cut_modeled_conc %>% 
  ggplot(aes(x = time, y = log(mean_concentration_est), col = station)) +
  geom_point() +
  geom_segment(aes(x = time, xend = time, y = log(ci25_concentration_est), yend = log(ci75_concentration_est))) +
  facet_grid(rows=vars(creek)) +
  labs(y= "Log copies/L water", title = "Cutthroat Trout") + 
  theme_bw()

####################################################################
# Plot coho environmental samples 
####################################################################

coho_modeled_conc %>% 
  ggplot(aes(x = time, y = log(mean_concentration_est), col = station)) +
  geom_point() +
  geom_segment(aes(x = time, xend = time, y = log(ci25_concentration_est), yend = log(ci75_concentration_est))) +
  facet_grid(rows=vars(creek)) +
  labs(y= "Log copies/L water", title = "Coho Salmon") + 
  theme_bw()




