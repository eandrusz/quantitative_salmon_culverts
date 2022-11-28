

library(here)
library(tidyverse)
library(ggplot2)
library(lubridate)

spawnerdata <- read.csv(here("Input","spawner_data.csv"))

spawner <- spawnerdata %>% 
  filter(! Species=="") %>% 
  filter(!Species=="Unknown") %>% 
  filter(! is.na(Fish.Count)) %>% 
  mutate(timeplot = parse_date_time(Survey.Date, "ymd", tz = "US/Pacific")) %>% 
  mutate(timeplot = with_tz(timeplot, "UTC")) %>% 
  mutate(yearplot = year(timeplot)) %>% 
  mutate(monthplot = month(timeplot)) %>% 
  mutate(dayplot = day(timeplot)) %>% 
  mutate(datetime = lubridate::make_datetime(yearplot, monthplot, dayplot))

ggplot(spawner, aes(x=datetime, y=Fish.Count)) + 
  geom_point() +
  facet_wrap(~Species) +
  theme_bw() +
  guides(color='none') +
  #scale_x_date(date_labels = "%b") +
  labs(x="Date", y="Number of Fish Counted", title = "Padden Spawner Surveys")
#ggsave(here("Output","SupplementalFigures","historical_spawner_data.png"))
