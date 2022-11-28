## NGN Flow data
# Author: Eily Allan 
# Person running: Eily
# Last modified: 10/17/22 by Eily
# Date of run: 10/17/22 by Eily 

# Overview 
# This script takes the flow data from three gauges (obtained by request from the City of Bellingham) and cleans it up and plots it for the time of the study. 

# Inputs: 
# 1-3) flow data for Padden Creek, Squalicum Creek, Chuckanut Creek 

# Outputs: 
# 1-3) only relevant dates and converted to m3/s

####################################################################
# Set up
####################################################################

# Load packages
library(here)
library(tidyverse)
library(lubridate) # package that makes handling dates much easier 
library(ggplot2)
library(gridExtra)

# Read in files
padflow <- read.csv(here("Input","qpcr","flow","Padden_Creek_Computed_Discharge_(cfs).Time_Series_Data.2022102422231948.csv"), skip=29)[,1:2]
sqmflow <- read.csv(here("Input","qpcr","flow","Squalicum_Creek_Computed_Discharge_(cfs).Time_Series_Data.2022102422132690.csv"), skip=27)[,1:2]
chkflow <- read.csv(here("Input","qpcr","flow","Chuckanut_Creek_Computed_Discharge_(cfs).Time_Series_Data.2022102422194411.csv"), skip=30)[,1:2]

# read in when filtering occurred 
filtermeta <- readRDS(here("Output","qpcr","backpack","adj_vol_filtered.RDS"))

####################################################################
# Clean up data, trim to correct dates, convert units
####################################################################

# clean up flow gauge data
pad.flow <- padflow %>% 
  add_column(creek = "Padden") %>% 
  dplyr::rename(flow_cfs = Discharge.Fairhaven.Park..cfs.) %>% 
  mutate(flow_m3s = flow_cfs*0.028316847) %>% 
  mutate(timeplot = parse_date_time(TimeStamp, "mdy_HM", tz = "US/Pacific")) %>% 
  mutate(timeplot = with_tz(timeplot, "UTC")) %>% 
  mutate(daysampletoavg = floor_date(timeplot, unit="day")) %>% 
  mutate(yearplot = year(timeplot)) %>% 
  mutate(monthplot = month(timeplot)) %>% 
  mutate(dayplot = day(timeplot)) %>% 
  mutate(hourplot = hour(timeplot)) %>% 
  mutate(minplot = minute(timeplot)) %>% 
  filter(yearplot > 2020) %>% 
  filter(timeplot < "2022-04-01 00:00:00 UTC")
  
sqm.flow <- sqmflow %>% 
  dplyr::rename(flow_cfs = Discharge.West.Street..cfs.) %>% 
  add_row(TimeStamp="4/1/2022 0:00", flow_cfs =0) %>% 
  add_row(TimeStamp="11/16/2021 9:30", flow_cfs =0) %>% 
  add_column(creek = "Squalicum") %>% 
  mutate(flow_m3s = flow_cfs*0.028316847) %>% 
  mutate(timeplot = parse_date_time(TimeStamp, "mdy_HM", tz = "US/Pacific")) %>% 
  mutate(timeplot = with_tz(timeplot, "UTC")) %>% 
  mutate(daysampletoavg = floor_date(timeplot, unit="day")) %>% 
  mutate(yearplot = year(timeplot)) %>% 
  mutate(monthplot = month(timeplot)) %>% 
  mutate(dayplot = day(timeplot)) %>% 
  mutate(hourplot = hour(timeplot)) %>% 
  mutate(minplot = minute(timeplot)) %>%  
  filter(yearplot > 2020) 

chk.flow <- chkflow %>% 
  dplyr::rename(flow_cfs = Discharge.Arroyo.Park..cfs.) %>%
  add_row(TimeStamp="4/1/2022 0:00", flow_cfs =0) %>%
  add_row(TimeStamp="11/12/2021 00:00", flow_cfs =0) %>%
  add_column(creek = "Chuckanut") %>% 
  mutate(flow_m3s = flow_cfs*0.028316847) %>% 
  mutate(timeplot = parse_date_time(TimeStamp, "mdy_HM", tz = "US/Pacific")) %>% 
  mutate(timeplot = with_tz(timeplot, "UTC")) %>% 
  mutate(daysampletoavg = floor_date(timeplot, unit="day")) %>% 
  mutate(yearplot = year(timeplot)) %>% 
  mutate(monthplot = month(timeplot)) %>% 
  mutate(dayplot = day(timeplot)) %>% 
  mutate(hourplot = hour(timeplot)) %>% 
  mutate(minplot = minute(timeplot)) %>% 
  filter(yearplot > 2020)

# clean up sample collection metadata 

filtermetaplot <- filtermeta %>% 
  mutate(DateStamp = gsub("UTC","", DateStamp)) %>% 
  mutate(timeplot = parse_date_time(DateStamp, "ymd_HMS")) %>% 
  select(-DateStamp) %>% 
  filter(!str_detect("Up5", Sample)) %>% 
  separate(Sample, into=c("time","creek","station","biol"), remove=FALSE) %>% 
  mutate(creek = case_when(creek == "3Chk" ~ "Chuckanut",
                           creek == "4Pad" ~ "Padden",
                           creek == "5Sqm" ~ "Squalicum",
                           TRUE ~ creek)) %>% 
  select(-c(biol,station,time, Adj_Vol)) %>% 
  mutate(timeplot = round_date(timeplot, unit="15 minutes")) %>% 
  mutate(daysampletoavg = floor_date(timeplot, unit="day"))


#pad.flow.save <- pad.flow %>% select(c(timeplot, flow_m3s)) %>% ungroup() %>% mutate(creek="Padden")
#pad.flow.save <- pad.flow %>% select(c(timeplot, flow_m3s)) %>% ungroup()
#pad.flow.save <- pad.flow %>% select(c(timeplot, flow_m3s)) %>% ungroup()

#write.csv(pad.flow.save, here("Output", "qpcr", "pad_hourly_flow.csv"), row.names=FALSE)

####################################################################
# Visualize 
####################################################################

padplot <- ggplot(pad.flow, aes(x=timeplot, y=flow_m3s)) + 
  geom_line() +
  theme_bw() +
  labs(x="Date", y="Flow, m^3/s", title = "Padden Flow")

sqmplot <- ggplot(sqm.flow, aes(x=timeplot, y=flow_m3s)) + 
  geom_line() +
  theme_bw() +
  labs(x="Date", y="Flow, m^3/s", title = "Squalicum Flow")

chkplot <- ggplot(chk.flow, aes(x=timeplot, y=flow_m3s)) + 
  geom_line() +
  theme_bw() +
  labs(x="Date", y="Flow, m^3/s", title = "Chuckanut Flow")

grid.arrange(padplot, sqmplot, chkplot, ncol=1, nrow =3)

####################################################################
# Time average per month so we get one value per sample  
####################################################################

pad.monthly <- pad.flow %>% 
  group_by(monthplot) %>% 
  mutate(meanmonthflow = mean(flow_m3s)) %>% 
  select(c(creek, yearplot, monthplot, meanmonthflow)) %>% 
  distinct() %>% 
  add_column(dayplot=15) %>% 
  mutate(timestring = paste0(yearplot,"-",monthplot,"-",dayplot," 12:00")) %>% 
  mutate(timeplot = ymd_hm(timestring))

chk.monthly <- chk.flow %>% 
  group_by(monthplot) %>% 
  mutate(meanmonthflow = mean(flow_m3s)) %>% 
  select(c(creek, yearplot, monthplot, meanmonthflow)) %>% 
  distinct() %>% 
  add_column(dayplot=15) %>% 
  mutate(timestring = paste0(yearplot,"-",monthplot,"-",dayplot," 12:00")) %>% 
  mutate(timeplot = ymd_hm(timestring))

sqm.monthly <- sqm.flow %>% 
  group_by(monthplot) %>% 
  mutate(meanmonthflow = mean(flow_m3s)) %>% 
  select(c(creek, yearplot, monthplot, meanmonthflow)) %>% 
  distinct() %>% 
  add_column(dayplot=15) %>% 
  mutate(timestring = paste0(yearplot,"-",monthplot,"-",dayplot," 12:00")) %>% 
  mutate(timeplot = ymd_hm(timestring))

####################################################################
# Time average over each day of sampling (still one per month)
####################################################################

paddailyflow <- pad.flow %>% 
  group_by(daysampletoavg) %>% 
  mutate(meandayflow = mean(flow_m3s)) %>% 
  select(c(creek, daysampletoavg, meandayflow)) %>% 
  distinct() 

pad.daily <- filtermetaplot %>% 
  filter(creek == "Padden") %>% 
  left_join(paddailyflow, by = "daysampletoavg")

chkdailyflow <- chk.flow %>% 
  group_by(daysampletoavg) %>% 
  mutate(meandayflow = mean(flow_m3s)) %>% 
  select(c(creek, daysampletoavg, meandayflow)) %>% 
  distinct() 

chk.daily <- filtermetaplot %>% 
  filter(creek == "Chuckanut") %>% 
  left_join(chkdailyflow, by = "daysampletoavg")

sqmdailyflow <- sqm.flow %>% 
  group_by(daysampletoavg) %>% 
  mutate(meandayflow = mean(flow_m3s)) %>% 
  select(c(creek, daysampletoavg, meandayflow)) %>% 
  distinct() 

sqm.daily <- filtermetaplot %>% 
  filter(creek == "Squalicum") %>% 
  left_join(sqmdailyflow, by = "daysampletoavg")

####################################################################
# Get discrete measurement at time of sampling (still one per month)
####################################################################

pad.discrete <- filtermetaplot %>% 
  filter(creek=="Padden") %>% 
  left_join(pad.flow %>% select(c(timeplot, flow_m3s)), by = "timeplot") %>% 
  mutate(yearplot = year(timeplot)) %>% 
  mutate(monthplot = month(timeplot)) %>% 
  mutate(dayplot = day(timeplot)) 

chk.discrete <- filtermetaplot %>% 
  filter(creek=="Chuckanut") %>% 
  left_join(chk.flow %>% select(c(timeplot, flow_m3s)), by = "timeplot") %>% 
  mutate(yearplot = year(timeplot)) %>% 
  mutate(monthplot = month(timeplot)) %>% 
  mutate(dayplot = day(timeplot)) 

sqm.discrete <- filtermetaplot %>% 
  filter(creek=="Squalicum") %>% 
  left_join(sqm.flow %>% select(c(timeplot, flow_m3s)), by = "timeplot") %>% 
  mutate(yearplot = year(timeplot)) %>% 
  mutate(monthplot = month(timeplot)) %>% 
  mutate(dayplot = day(timeplot)) 

####################################################################
# Visualize now with discrete flow at time of sampling
####################################################################

padplot <- ggplot() + 
  geom_line(data=pad.flow, aes(x=timeplot, y=flow_m3s), color="darkgrey") +
  geom_point(data=pad.discrete, aes(x=timeplot, y=flow_m3s), size=3,  pch=21) +
  #geom_point(data=pad.daily, aes(x=timeplot, y=meandayflow), size=3, bg="cyan", pch=22) +
  #geom_point(data=pad.monthly, aes(x=timeplot, y=meanmonthflow), size=3, bg="aquamarine", pch=23) +
  labs(x="Date", y=bquote('Average Discharge '(m^3/s)), title = "Padden Discharge") +
  theme_bw() 

chkplot <- ggplot() + 
  geom_line(data=chk.flow, aes(x=timeplot, y=flow_m3s), color="darkgrey") +
  geom_point(data=chk.discrete, aes(x=timeplot, y=flow_m3s), size=3,  pch=21) +
  #geom_point(data=chk.daily, aes(x=timeplot, y=meandayflow), size=3, bg="cyan", pch=22) +
  #geom_point(data=chk.monthly, aes(x=timeplot, y=meanmonthflow), size=3, bg="aquamarine", pch=23) +
  labs(x="Date", y=bquote('Average Discharge '(m^3/s)), title = "Chuckanut Discharge") +
  theme_bw() 

sqmplot <- ggplot() + 
  geom_line(data=sqm.flow, aes(x=timeplot, y=flow_m3s), color="darkgrey") +
  geom_point(data=sqm.discrete, aes(x=timeplot, y=flow_m3s), size=3, pch=21) +
  #geom_point(data=sqm.daily, aes(x=timeplot, y=meandayflow), size=3, bg="cyan", pch=22) +
  #geom_point(data=sqm.monthly, aes(x=timeplot, y=meanmonthflow), size=3, bg="aquamarine", pch=23) +
  labs(x="Date", y=bquote('Average Discharge '(m^3/s)), title = "Squalicum Discharge") +
  theme_bw() 

grid.arrange(padplot, chkplot, sqmplot,ncol=1, nrow =3)

#save
g <- arrangeGrob(padplot, chkplot, sqmplot, nrow=3) #generates g
#ggsave(file=here("Output","Figures","flow_gauges.png"),g)

####################################################################
# Write output 
####################################################################

all_monthly_avg <- rbind(pad.monthly, chk.monthly, sqm.monthly)
write.csv(all_monthly_avg, here("Output", "qpcr", "monthly_flow.csv"), row.names=FALSE)

all_daily_avg <- rbind(pad.daily, chk.daily, sqm.daily)
write.csv(all_monthly_avg, here("Output", "qpcr", "daily_flow.csv"), row.names=FALSE)

all_closest <- rbind(pad.discrete, chk.discrete, sqm.discrete)
write.csv(all_closest, here("Output", "qpcr", "closest_flow.csv"), row.names=FALSE)


#write.csv(pad.discrete, here("Output", "qpcr", "pad_sampling_flow.csv"), row.names=FALSE)


####################################################################
# some stats....  
####################################################################


padflowstats <- pad.flow %>% 
  summarize(meanflow = mean(flow_m3s), maxflow = max(flow_m3s))

chkflowstats <- chk.flow %>% 
  summarize(meanflow = mean(flow_m3s), maxflow = max(flow_m3s))

sqmflowstats <- sqm.flow %>% 
  summarize(meanflow = mean(flow_m3s), maxflow = max(flow_m3s))

padlow <- pad.flow %>% 
  filter(flow_cfs < 0.11)
padperclow = nrow(padlow)/nrow(pad.flow)*100

chklow <- chk.flow %>% 
  filter(flow_cfs < 0.11)
chkperclow = nrow(chklow)/nrow(chk.flow)*100

sqmlow <- sqm.flow %>% 
  filter(flow_cfs < 0.21)
sqmperclow = nrow(sqmlow)/nrow(sqm.flow)*100

