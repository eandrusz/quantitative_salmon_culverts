# Person running: Eily
# Last modified: 10/19/22 by Eily
# Date of run: 10/19/22 by Eily 

# Overview 
# This script takes 

# Inputs: 
# 1) 

# Outputs: 
# 1)

####################################################################
# Set up
####################################################################

# Load packages
library(tidyverse)
library(here)

# Read in files 
cut_output <- readRDS(here("Output","qpcr", "cut_modeled_conc.RDS"))
coho_output <- readRDS(here("Output","qpcr", "coho_modeled_conc.RDS"))
cut_input <- read.csv(here("Output", "qpcr", "cut_data_for_stan.csv"))
coho_input <- read.csv(here("Output", "qpcr", "coho_data_for_stan.csv"))

####################################################################
# Number of samples per assay
####################################################################

checkcutin <- cut_input %>% 
  filter(!str_detect(Sample, "Up5")) %>% 
  filter(!str_detect(Sample, "St")) %>% 
  group_by(Sample) %>% 
  summarize(n=n()) 
  # %>% 
  # ungroup() %>% 
  # separate(Sample, into=c("time","creek","station","biorep")) %>% 
  # group_by(time,creek) %>% 
  # summarize(timecreekn = n())

checkcutout <- cut_output %>%
  unite(Sample, c(time,creek,station,biorep)) %>% 
  filter(!str_detect(Sample, "Up5")) %>% 
  group_by(Sample) %>% 
  summarize(n=n())

checkcohoin <- coho_input %>% 
  filter(!str_detect(Sample, "Up5")) %>% 
  filter(!str_detect(Sample, "St")) %>% 
  group_by(Sample) %>% 
  summarize(n=n()) 

checkcohoout <- coho_output %>%
  unite(Sample, c(time,creek,station,biorep)) %>% 
  filter(!str_detect(Sample, "Up5")) %>% 
  group_by(Sample) %>% 
  summarize(n=n())


####################################################################
# Standard curves 
####################################################################

cutstds <- cut_input %>% 
  filter(str_detect(Sample, "St")) %>% 
  group_by(Plate, Sample) %>% 
  summarize(meanCt = mean(Ct))

cohostds <- coho_input %>% 
  filter(str_detect(Sample, "St")) %>% 
  group_by(Plate, Sample) %>% 
  summarize(meanCt = mean(Ct))


####################################################################
# Inhibition
####################################################################

inhibition <- cut_input %>%
  filter(!str_detect(Sample, "Up5")) %>% 
  filter(!str_detect(Sample, "St")) %>% 
  group_by(dilution) %>% 
  summarize(n=n())

perc_not_inhibited = inhibition$n[1]/sum(inhibition$n)
perc_of_inihibied_1to10 = inhibition$n[4]/sum(inhibition$n[2:8])

####################################################################
# Model output stats
####################################################################

cut_output %>% 
  summarize(min = min(mean_concentration_est), max = max(mean_concentration_est), mean = mean(mean_concentration_est))

coho_output %>% 
  summarize(min = min(mean_concentration_est), max = max(mean_concentration_est), mean = mean(mean_concentration_est))

