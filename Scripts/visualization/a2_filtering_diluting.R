## NGN Visualizing field / filtering data
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
library(ggplot2)

cut_data_for_stan <- read_csv(here("Output","qpcr","cut_data_for_stan.csv"))
coho_data_for_stan <- read_csv(here("Output","qpcr","coho_data_for_stan.csv"))


####################################################################
# volume filtered
####################################################################

cut_data_for_stan %>% 
  filter(! str_detect(Sample, "St")) %>% 
  select(c(Sample, dilution, Adj_Vol)) %>% 
  separate(Sample, into=c("time","creek","station","biol")) %>% 
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
  ggplot(aes(x= newtime, y = Adj_Vol)) +
  geom_point() +
  facet_grid(~creek ~station, scales = "free") +
  theme_bw() +
  labs(x="Date (YY-MM)", y= "Volume filtered (L)") + 
  scale_y_continuous(limits=c(0,3), breaks=c(0,1,2,3)) + 
  scale_x_discrete(guide = guide_axis(angle = -45))

ggsave(here("Output","Figures","volume_filtered_by_sample.png"))

####################################################################
# dilution factor during inhibiton testing 
####################################################################

cut_data_for_stan %>% 
  filter(! str_detect(Sample, "St")) %>% 
  select(c(Sample, dilution, Adj_Vol)) %>% 
  separate(Sample, into=c("time","creek","station","biol")) %>% 
  mutate(creek = case_when(creek == "1Prt" ~ "Portage",
                           creek == "2Brn" ~ "Barnes",
                           creek == "3Chk" ~ "Chuckanut",
                           creek == "4Pad" ~ "Padden",
                           creek == "5Sqm" ~ "Squalicum",
                           TRUE ~ creek)) %>% 
  filter(station != "Up5") %>% 
  mutate(station = case_when(station == "Up11" ~ "Up",
                             TRUE ~ station)) %>% 
  separate(time, into = c("month","year"), sep = 2) %>% 
  unite(newtime, c(year,month), sep="-") %>% 
  mutate(station = case_when(station == "Dn" ~ "Down",
                             TRUE ~ station)) %>% 
  ggplot(aes(x= newtime, y = log10(dilution))) +
  geom_point() +
  facet_grid(~creek ~station, scales = "free") +
  theme_bw() +
  labs(x="Date", y= "Log 10 Dilution Factor") + 
  lims(y=c(0,4)) + 
  scale_x_discrete(guide = guide_axis(angle = -45))

ggsave(here("Output","Figures","inhibition_dilution_factors.png"))
