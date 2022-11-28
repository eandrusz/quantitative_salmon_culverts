## COMPARE CUTTHROAT TOTAL DNA TO COHO TOTAL DNA 

library(tidyverse)
library(here)
library(ggplot2)

cut_total_est <- readRDS(here("Output","salmonids_abs_abundance_bio_posterior_EAA_20201020.RDS"))
coho_total_est <- readRDS(here("Output","salmonids_abs_abundance_fromCOHO.RDS"))

cut <- cut_total_est %>% 
  select(c(bottle, newtime, creek, station, bio, meantotdna, meanpropreads)) %>% 
  distinct() %>% 
  mutate(anchor="cutthroat")

coho <- coho_total_est %>% 
  select(c(bottle, meantotdna)) %>% 
  distinct() %>% 
  mutate(anchor="coho")

both <- coho %>% 
  left_join(cut, by="bottle")

both %>% 
  ggplot(aes(x=log10(meantotdna.x), y=log10(meantotdna.y), color=as.factor(meanpropreads))) +
  guides(color=FALSE) +
  #guides(fill = guide_colourbar(nbin = 10)) + 
  geom_point() + 
  geom_abline() + 
  labs(x="Coho as anchor", y="Cutthroat as anchor", title="Estimated total DNA (log10)")
