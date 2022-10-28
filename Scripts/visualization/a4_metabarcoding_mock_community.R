## NGN Visualizing metabarcoding results before / after Stan model
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



####################################################################
# Plot proportions of salmonids before QM correction
####################################################################


quants_to_plot <- mergedc %>% 
  filter(creek != "4Pad5") %>% 
  mutate(station = case_when(station == 1 ~ "Down",
                             station == 2 ~ "Up")) %>%
  separate(time, into = c("month","year"), sep = 2) %>% 
  unite(newtime, c(year,month), sep="-") %>% 
  mutate(creek = case_when(creek == "1Prt" ~ "Portage",
                           creek == "2Brn" ~ "Barnes",
                           creek == "3Chk" ~ "Chuckanut",
                           creek == "4Pad11" ~ "Padden",
                           creek == "5Sqm" ~ "Squalicum",
                           TRUE ~ creek))  

ggplot(quants_to_plot, aes(x=newtime, y=log10(meandnaconc), color=species)) +
  geom_point() +
  facet_grid(~creek ~station) +
  scale_color_manual(values = pal_okabe_ito) +
  guides(color= 'none') + 
  labs(x="Date (YY-MM)", y= "Log10 copies/L water") + 
  theme_bw() + 
  scale_x_discrete(guide = guide_axis(angle = -45))