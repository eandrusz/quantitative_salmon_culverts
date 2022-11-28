

library(here)
library(tidyverse)
library(ggplot2)
library(lubridate)

stocking <- read.csv(here("Input","stocking_data.csv"))

# 2021 numbers
stocking2021 <- stocking %>% 
  filter(y==2021) %>% 
  group_by(species) %>% 
  summarize(totalstocked = sum(quant))

stock <- stocking %>% 
  mutate(d=replace_na(d,15)) %>% 
  mutate(dateplot = make_datetime(y, m, d)) %>% 
  mutate(location="Padden") 

ggplot(stock, aes(x=dateplot, y=quant)) + 
  geom_point() +
  theme_bw() +
  facet_wrap(~species) +
  labs(x="Date", y="Number of Fish Stocked", title = "Padden Stocking")
#ggsave(here("Output","SupplementalFigures","historical_stocking_by_species_ts.png"))

stock %>% 
  filter(dateplot > "2021-01-01") %>% 
  ggplot(aes(x=dateplot, y=quant, fill=factor(species))) + 
  geom_col() +
  theme_bw() + 
  #facet_wrap(~species) +
  labs(x="Date", y="Number of Fish Stocked", title = "2021 Padden Stocking", fill="Species")
#ggsave(here("Output","SupplementalFigures","stocking_2021only.png"))


stock2 <- stocking %>% 
  mutate(d=replace_na(d,15)) %>% 
  mutate(dateplot = make_datetime(2022, m, d)) %>% 
  mutate(location="Padden")  

ggplot(stock2, aes(x=dateplot, y=quant, color=factor(y))) + 
  geom_point() +
  facet_wrap(~species) +
  theme_bw() +
  scale_y_log10() +
  labs(x="Date", y="Number of Fish Stocked", title = "Padden Stocking", color="Year")
#ggsave(here("Output","SupplementalFigures","historical_stocking_by_species_year.png"))

