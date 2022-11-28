


####################################################################
# Set up
####################################################################

# Load packages
library(here)
library(tidyverse)
library(lubridate) # package that makes handling dates much easier 
library(ggplot2)
library(gridExtra)

source(here("Scripts","functions","flow_functions.R"))

# Read in files
padflow1 <- read.csv(here("Input","qpcr","flow","old", "Padden_Creek_Discharge.Time_Series_Data.2021110822572368.csv"), skip=29)[,1:2]
sqmflow1 <- read.csv(here("Input","qpcr","flow","old","Squalicum_Creek_Discharge.Time_Series_Data.2021110822504744.csv"), skip=27)[,1:2]
chkflow1 <- read.csv(here("Input","qpcr","flow","old","Chuckanut_Creek.Time_Series_Data.2021110823201271.csv"), skip=30)[,1:2]

padflow2 <- read.csv(here("Input","qpcr","flow","Padden_Creek_Computed_Discharge_(cfs).Time_Series_Data.2022102422231948.csv"), skip=29)[,1:2]
sqmflow2 <- read.csv(here("Input","qpcr","flow","Squalicum_Creek_Computed_Discharge_(cfs).Time_Series_Data.2022102422132690.csv"), skip=27)[,1:2]
chkflow2 <- read.csv(here("Input","qpcr","flow","Chuckanut_Creek_Computed_Discharge_(cfs).Time_Series_Data.2022102422194411.csv"), skip=30)[,1:2]

# read in when filtering occurred 
filtermeta <- readRDS(here("Output","qpcr","backpack","adj_vol_filtered.RDS"))

####################################################################
# Clean up gauge data
####################################################################

padflow <- rbind(padflow1, padflow2)
chkflow <- rbind(chkflow1, chkflow2)
sqmflow <- rbind(sqmflow1, sqmflow2)

pad.flow <- padflow %>% 
  distinct() %>% 
  add_column(creek = "Padden") %>% 
  dplyr::rename(flow_cfs = Discharge.Fairhaven.Park..cfs.) 

chk.flow <- chkflow %>% 
  distinct() %>% 
  add_column(creek = "Chuckanut") %>% 
  dplyr::rename(flow_cfs = Discharge.Arroyo.Park..cfs.) 

sqm.flow <- sqmflow %>% 
  distinct() %>% 
  add_column(creek = "Squalicum") %>% 
  dplyr::rename(flow_cfs = Discharge.West.Street..cfs.) 

####################################################################
# Clean up sampling data
####################################################################

filtermetaplot <- filtermeta %>% 
  mutate(DateStamp = gsub("UTC","", DateStamp)) %>% 
  mutate(timeplot = parse_date_time(DateStamp, "ymd_HMS")) %>% 
  dplyr::select(-DateStamp) %>% 
  filter(!str_detect("Up5", Sample)) %>% 
  separate(Sample, into=c("time","creek","station","biol"), remove=FALSE) %>% 
  mutate(creek = case_when(creek == "3Chk" ~ "Chuckanut",
                           creek == "4Pad" ~ "Padden",
                           creek == "5Sqm" ~ "Squalicum",
                           TRUE ~ creek)) %>% 
  dplyr::select(-c(biol,station,time, Adj_Vol)) %>% 
  mutate(timeplot = round_date(timeplot, unit="15 minutes")) %>% 
  mutate(daysampletoavg = floor_date(timeplot, unit="day"))


####################################################################
# Use functions to manipulate and plot 
####################################################################

pad.past.flow <- format_flow(pad.flow)
pad.sampling.flow <- format_flow2(pad.flow)
pad.year.avg <- yearly_avg_flow(pad.past.flow)
pad.discrete <- flow_discrete(filtermetaplot, "Padden", pad.sampling.flow)
pad.ts.plot <- plot_ts(pad.past.flow, "Padden")
pad.ts.facet.plot <- plot_facet_year(pad.past.flow, "Padden")
pad.ts.year.plot <- plot_year_avg(pad.year.avg, "Padden")
pad.ts.year.discrete.plot <- plot_year_avg_discrete(pad.year.avg, pad.discrete, "Padden")
pad.monthly.avg <- monthly_avg_flow(pad.past.flow)
pad.ts.monthly.discrete.plot <- plot_year_avg_discrete(pad.monthly.avg, pad.discrete, "Padden")
pad.ts.sampling.discrete.plot <-  plot_sampling_discrete(pad.sampling.flow, pad.discrete, "Padden")


sqm.past.flow <- format_flow(sqm.flow)
sqm.sampling.flow <- format_flow2(sqm.flow)
sqm.year.avg <- yearly_avg_flow(sqm.past.flow)
sqm.discrete <- flow_discrete(filtermetaplot, "Squalicum", sqm.sampling.flow)
sqm.ts.plot <- plot_ts(sqm.past.flow, "Squalicum")
sqm.ts.facet.plot <- plot_facet_year(sqm.past.flow, "Squalicum")
sqm.ts.year.plot <- plot_year_avg(sqm.year.avg, "Squalicum")
sqm.ts.year.discrete.plot <- plot_year_avg_discrete(sqm.year.avg, sqm.discrete, "Squalicum")
sqm.monthly.avg <- monthly_avg_flow(sqm.past.flow)
sqm.ts.monthly.discrete.plot <- plot_year_avg_discrete(sqm.monthly.avg, sqm.discrete, "Squalicum")

chk.past.flow <- format_flow(chk.flow)
chk.sampling.flow <- format_flow2(chk.flow)
chk.year.avg <- yearly_avg_flow(chk.past.flow)
chk.discrete <- flow_discrete(filtermetaplot, "Chuckanut", chk.sampling.flow)
chk.ts.plot <- plot_ts(chk.past.flow, "Chuckanut")
chk.ts.facet.plot <- plot_facet_year(chk.past.flow, "Chuckanut")
chk.ts.year.plot <- plot_year_avg(chk.year.avg, "Chuckanut")
chk.ts.year.discrete.plot <- plot_year_avg_discrete(chk.year.avg, chk.discrete, "Chuckanut")
chk.monthly.avg <- monthly_avg_flow(chk.past.flow)
chk.ts.monthly.discrete.plot <- plot_year_avg_discrete(chk.monthly.avg, chk.discrete, "Chuckanut")

flow_CFs <- pad.monthly.avg %>%
  dplyr::rename(padflow = val) %>% 
  left_join(chk.monthly.avg, by="datetime") %>% 
  dplyr::rename(chkflow = val) %>% 
  left_join(sqm.monthly.avg, by="datetime") %>% 
  dplyr::rename(sqmflow = val) %>% 
  mutate(chk_cf = chkflow/padflow)%>% 
  mutate(sqm_cf = sqmflow/padflow)

use_padden_and_cfs <- pad.discrete %>% 
  filter(str_detect(Sample, "Dn")) %>% 
  mutate(monthplot = month(timeplot)) %>% 
  mutate(datetime = lubridate::make_datetime(2022, monthplot, 15)) %>%
  left_join(flow_CFs %>% dplyr::select(c(datetime, chk_cf, sqm_cf)), by= "datetime") %>% 
  mutate(flow_chk_by_padden = flow_m3s*chk_cf) %>% 
  mutate(flow_sqm_by_padden = flow_m3s*sqm_cf) 

compare_chk_sqm <- use_padden_and_cfs %>% 
  dplyr::select(c(daysampletoavg, flow_m3s,flow_chk_by_padden,flow_sqm_by_padden)) %>%
  dplyr::rename(pad_flowm3s = flow_m3s) %>% 
  left_join(chk.discrete %>% dplyr::select(c(daysampletoavg, flow_m3s)), by= "daysampletoavg") %>% 
  dplyr::rename(chk_flowm3s = flow_m3s) %>% 
  left_join(sqm.discrete %>% dplyr::select(c(daysampletoavg, flow_m3s)), by= "daysampletoavg") %>% 
  dplyr::rename(sqm_flowm3s = flow_m3s) %>% 
  distinct()

ggplot(newflow2, aes(x = timeplot, y = flow_m3s)) + 
  geom_line() + 
  #scale_x_date(date_labels = "%b") +
  #scale_x_datetime(breaks = pad.flow$datetime, labels = monthplot) +
  labs(y=bquote('Average Discharge '(m^3/s)), x="Month", title = paste0(creek, " Flow Rate (m3/s): March 2021-March 2022"), colour = "Year") + 
  theme_bw() +
  geom_point(data=flow.discrete, aes(x=timeplot, y=flow_m3s), size=3, bg="blue", pch=21)


flow_cf_plot <- ggplot(compare_chk_sqm) + 
  geom_point(aes(x=daysampletoavg, y=flow_chk_by_padden), color="red", size=3, pch=4) +
  geom_point(aes(x=daysampletoavg, y=chk_flowm3s), color="red", size=3, pch=1) + 
  geom_point(aes(x=daysampletoavg, y=flow_sqm_by_padden), color="blue", size=3, pch=4) +
  geom_point(aes(x=daysampletoavg, y=sqm_flowm3s), color="blue", size=3, pch=1) +
  geom_line(data=pad.sampling.flow, aes(x = timeplot, y = flow_m3s), color="black", alpha=.1) +
  #geom_line(data=chk.sampling.flow, aes(x = timeplot, y = flow_m3s), color="red", alpha=.1) +
  #geom_line(data=sqm.sampling.flow, aes(x = timeplot, y = flow_m3s), color="blue", alpha=.1) +
  #geom_point(data=pad.discrete, aes(x=timeplot, y=flow_m3s), size=3, bg="black", pch=1) +
  theme_bw() +
  labs(y=bquote('Average Discharge '(m^3/s)), x="Date of Sampling")

chk_cf_plot <- ggplot(compare_chk_sqm) + 
  geom_point(aes(x=daysampletoavg, y=flow_chk_by_padden), color="red", size=3, pch=4) +
  geom_point(aes(x=daysampletoavg, y=chk_flowm3s), color="red", size=3, pch=1) + 
  #geom_point(aes(x=daysampletoavg, y=flow_sqm_by_padden), color="blue", size=3, pch=4) +
  #geom_point(aes(x=daysampletoavg, y=sqm_flowm3s), color="blue", size=3, pch=1) +
  #geom_line(data=pad.sampling.flow, aes(x = timeplot, y = flow_m3s), color="black", alpha=.1) +
  geom_line(data=chk.sampling.flow, aes(x = timeplot, y = flow_m3s), color="red", alpha=.5) +
  #geom_line(data=sqm.sampling.flow, aes(x = timeplot, y = flow_m3s), color="blue", alpha=.1) +
  #geom_point(data=pad.discrete, aes(x=timeplot, y=flow_m3s), size=3, bg="black", pch=1) +
  theme_bw() +
  labs(y=bquote('Average Discharge '(m^3/s)), x="Date of Sampling", title="Chuckanut Predictions")

sqm_cf_plot <- ggplot(compare_chk_sqm) + 
  #geom_point(aes(x=daysampletoavg, y=flow_chk_by_padden), color="red", size=3, pch=4) +
  #geom_point(aes(x=daysampletoavg, y=chk_flowm3s), color="red", size=3, pch=1) + 
  geom_point(aes(x=daysampletoavg, y=flow_sqm_by_padden), color="blue", size=3, pch=4) +
  geom_point(aes(x=daysampletoavg, y=sqm_flowm3s), color="blue", size=3, pch=1) +
  #geom_line(data=pad.sampling.flow, aes(x = timeplot, y = flow_m3s), color="black", alpha=.1) +
  #geom_line(data=chk.sampling.flow, aes(x = timeplot, y = flow_m3s), color="red", alpha=.1) +
  geom_line(data=sqm.sampling.flow, aes(x = timeplot, y = flow_m3s), color="blue", alpha=.5) +
  #geom_point(data=pad.discrete, aes(x=timeplot, y=flow_m3s), size=3, bg="black", pch=1) +
  theme_bw() +
  labs(y=bquote('Average Discharge '(m^3/s)), x="Date of Sampling", title="Squalicum Predictions")

pad_cf_plot <- ggplot(compare_chk_sqm) + 
  #geom_point(aes(x=daysampletoavg, y=flow_chk_by_padden), color="red", size=3, pch=4) +
  #geom_point(aes(x=daysampletoavg, y=chk_flowm3s), color="red", size=3, pch=1) + 
  #geom_point(aes(x=daysampletoavg, y=flow_sqm_by_padden), color="blue", size=3, pch=4) +
  #geom_point(aes(x=daysampletoavg, y=sqm_flowm3s), color="blue", size=3, pch=1) +
  geom_line(data=pad.sampling.flow, aes(x = timeplot, y = flow_m3s), color="black", alpha=.5) +
  #geom_line(data=chk.sampling.flow, aes(x = timeplot, y = flow_m3s), color="red", alpha=.1) +
  #geom_line(data=sqm.sampling.flow, aes(x = timeplot, y = flow_m3s), color="blue", alpha=.5) +
  geom_point(data=pad.discrete, aes(x=timeplot, y=flow_m3s), size=3, bg="black", pch=1) +
  theme_bw() +
  labs(y=bquote('Average Discharge '(m^3/s)), x="Date of Sampling", title="Padden (Used for Predictions)")


flowratestouse <- compare_chk_sqm %>% 
  dplyr::select(c(daysampletoavg, pad_flowm3s, flow_chk_by_padden,flow_sqm_by_padden)) %>% 
  mutate(flow_brn_by_padden = pad_flowm3s*1) %>%  #for now, just make something up.... 
  mutate(flow_prt_by_padden = pad_flowm3s*1) %>%  #for now, just make something up.... 
  distinct()

# write.csv(flowratestouse, here("Output","qpcr","20221121_flowrates_touse.csv"), row.names=FALSE)


####################################################################
# Plot and save  
####################################################################

#ggsave(pad.ts.sampling.discrete.plot, file=here("Output","Figures","padden_flow_sampling.png"))

g <- arrangeGrob(pad.ts.plot, chk.ts.plot, sqm.ts.plot, nrow=3) #generates g
ggsave(file=here("Output","SupplementalFigures","historical_flow.png"),g)

g <- arrangeGrob(pad.ts.facet.plot, chk.ts.facet.plot, sqm.ts.facet.plot, ncol=3) #generates g
ggsave(file=here("Output","SupplementalFigures","historical_flow_faceted.png"),g)

g <- arrangeGrob(pad.ts.year.plot, chk.ts.year.plot, sqm.ts.year.plot, nrow=3) #generates g
ggsave(file=here("Output","SupplementalFigures","historical_flow_year_avg.png"),g)

g <- arrangeGrob(pad.ts.year.discrete.plot, chk.ts.year.discrete.plot, sqm.ts.year.discrete.plot, nrow=3) #generates g
ggsave(file=here("Output","SupplementalFigures","historical_flow_year_avg_sampling.png"),g)

ggsave(flow_cf_plot, file=here("Output","SupplementalFigures","flow_chk_sqm_from_padden.png"))

g <- arrangeGrob(pad_cf_plot, chk_cf_plot, sqm_cf_plot, nrow=3) #generates g
ggsave(file=here("Output","SupplementalFigures","flow_by_cfs.png"),g)


