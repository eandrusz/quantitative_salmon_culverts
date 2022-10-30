## NGN Filtration data
# Author: Eily Allan 
# Person running: Eily
# Last modified: 10/17/22 by Eily
# Date of run: 10/17/22 by Eily 

# Overview 
# This script takes the all of the individual backpack files generated in the field and cleans up to get volume filtered for each sample.

# Inputs: 
# 1) 161 individual csv files for each sample run on the backpack (note: some were single samples, some were sets of 3 with trident)
# 2) csv file with March 2021 files (because no backpack files) and information on the remainder of the year if single or trident samples

# Outputs: 
# 1) table of samples with volume of water per filter

####################################################################
# Set up
####################################################################

# Load packages
library(here)
library(tidyverse)

# Set file path where individual csv files live
backpackpath <- here("Input","qpcr","backpack")
files <- list.files(path = backpackpath, pattern = "*.csv", recursive = T, full.names = T)

####################################################################
# Read in all files, clean up, store revelant information
####################################################################

# create empty matrix, then for loop to keep just volume filtered
file_vol <- matrix(NA, nrow=length(files), ncol=2)
for (i in 1:length(files)) {
  filename = files[i]
  
  lognum = unlist(strsplit(filename, "/"))[11]
  lognum = str_sub(lognum, start=1, end=5)
  lognum = as.numeric(lognum)
  
  data = read_csv(filename, n_max = 10, col_names = FALSE)
  vol = as.numeric(data[5,2])
  
  file_vol[i,1]=lognum
  file_vol[i,2]=vol
}
file_vol <- as_tibble(file_vol) %>% 
  rename(LogID = V1, Volume = V2)
  
# march 2021 files are not on here - so use this spreadsheet to read in values for march and then add the rest 
backpackmeta <- read_csv(here("Input","qpcr","backpack","vol_filtered.csv"))

# add april 2021 to feb 2022
backpackmeta <- backpackmeta %>% 
  left_join(file_vol, by="LogID") %>% 
  mutate(Volume = case_when(is.na(Volume.x) ~ Volume.y,
                           !is.na(Volume.x) ~ Volume.x)) %>% 
  dplyr::select(-c(Volume.x, Volume.y)) %>% 
  mutate(Adj_Vol = case_when(Trident == 0 ~ Volume,
                             Trident == 1 ~ Volume/3)) %>% 
  dplyr::select(c(Sample,Adj_Vol))


####################################################################
# Write output
####################################################################

write_rds(backpackmeta, here("Output","qpcr","backpack","adj_vol_filtered.RDS"), row.names=FALSE)
