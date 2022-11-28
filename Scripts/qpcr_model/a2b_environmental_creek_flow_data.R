


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

pad.flow2 <- format_flow(pad.flow)
pad.year.avg <- yearly_avg_flow(pad.flow2)
pad.discrete <- flow_discrete(filtermetaplot, "Padden", pad.flow2)
pad.ts.plot <- plot_ts(pad.flow2, "Padden")
pad.ts.facet.plot <- plot_facet_year(pad.flow2, "Padden")
pad.ts.year.plot <- plot_year_avg(pad.year.avg, "Padden")
pad.ts.year.discrete.plot <- plot_year_avg_discrete(pad.year.avg, pad.discrete, "Padden")

sqm.flow2 <- format_flow(sqm.flow)
sqm.year.avg <- yearly_avg_flow(sqm.flow2)
sqm.discrete <- flow_discrete(filtermetaplot, "Squalicum", sqm.flow2)
sqm.ts.plot <- plot_ts(sqm.flow2, "Squalicum")
sqm.ts.facet.plot <- plot_facet_year(sqm.flow2, "Squalicum")
sqm.ts.year.plot <- plot_year_avg(sqm.year.avg, "Squalicum")
sqm.ts.year.discrete.plot <- plot_year_avg_discrete(sqm.year.avg, sqm.discrete, "Squalicum")

chk.flow2 <- format_flow(chk.flow)
chk.year.avg <- yearly_avg_flow(chk.flow2)
chk.discrete <- flow_discrete(filtermetaplot, "Chuckanut", chk.flow2)
chk.ts.plot <- plot_ts(chk.flow2, "Chuckanut")
chk.ts.facet.plot <- plot_facet_year(chk.flow2, "Chuckanut")
chk.ts.year.plot <- plot_year_avg(chk.year.avg, "Chuckanut")
chk.ts.year.discrete.plot <- plot_year_avg_discrete(chk.year.avg, chk.discrete, "Chuckanut")

# ggsave(pad.ts.year.discrete.plot, here("Output","padavgflow.png"))
# ggsave(sqm.ts.year.discrete.plot, here("Output","sqmavgflow.png"))
# ggsave(chk.ts.year.discrete.plot, here("Output","chkavgflow.png"))

#save
g <- arrangeGrob(pad.ts.year.discrete.plot, chk.ts.year.discrete.plot, sqm.ts.year.discrete.plot, nrow=3) #generates g
#ggsave(file=here("Output","SupplementalFigures","yearavg_flow_gauges.png"),g)



####################################################################
# Write function to convert to m3/s and manipulate dates 
####################################################################

format_flow <- function(flow_df){
  require(tidyverse)
  require(lubridate)
  
newflow <- flow_df %>% 
  mutate(flow_m3s = flow_cfs*0.028316847) %>% 
  mutate(timeplot = parse_date_time(TimeStamp, "mdy_HM", tz = "US/Pacific")) %>% 
  mutate(timeplot = with_tz(timeplot, "UTC")) %>% 
  mutate(daysampletoavg = floor_date(timeplot, unit="day")) %>% 
  mutate(yearplot = year(timeplot)) %>% 
  mutate(monthplot = month(timeplot)) %>% 
  mutate(dayplot = day(timeplot)) %>% 
  mutate(hourplot = hour(timeplot)) %>% 
  mutate(minplot = minute(timeplot)) %>% 
  mutate(datetime = lubridate::make_datetime(2022, monthplot, dayplot, hourplot, minplot)) %>% 
  drop_na() 

return(newflow)
}


####################################################################
# Write function to convert to m3/s and manipulate dates 
####################################################################

yearly_avg_flow <- function(newflow){
  require(tidyverse)
  require(lubridate)
  
  newflow2 <- newflow %>% 
  group_by(monthplot,dayplot,hourplot,minplot) %>% 
  summarise(val = mean(flow_m3s, na.rm = T)) %>%
  mutate(datetime = make_date(2022, monthplot,dayplot)) 
  
  return(newflow2)
}

####################################################################
# Write function to find closest discrete timepoint 
####################################################################

flow_discrete <- function(filtermetaplot, creekname, newflow){
  require(tidyverse)
  
  flow.discrete <- filtermetaplot %>% 
  filter(creek==creekname) %>% 
  left_join(newflow %>% dplyr::select(c(timeplot, flow_m3s)), by = "timeplot") %>% 
  mutate(yearplot = year(timeplot)) %>% 
  mutate(monthplot = month(timeplot)) %>% 
  mutate(dayplot = day(timeplot)) %>% 
  mutate(datetime = make_date(2022, monthplot,dayplot))
  
  return(flow.discrete)
}

####################################################################
# Plot by year, year averaged, year averaged with discrete points 
####################################################################

plot_ts <- function(newflow, creek){
  require(ggplot2)
  
  p1 <- ggplot(newflow, aes(x=timeplot, y=flow_m3s)) + 
  geom_line() +
  theme_bw() +
  labs(x="Date", y="Flow Rate (m3/s)", title = creek)
  
  return(p1)
}


plot_facet_year <- function(newflow, creek){
  require(ggplot2)
  
  p2 <- ggplot(newflow) +
  geom_line(aes(x = datetime, y = flow_m3s, color = factor(yearplot))) +
  labs(y="Flow Rate (m3/s)", x="Date", title = creek, colour = "Year") + 
  facet_wrap(~yearplot, nrow=7) +
  theme_bw()
  
  return(p2)
}

plot_year_avg <- function(newflow2, creek){
  require(ggplot2)
  
p3 <- ggplot(newflow2, aes(x = datetime, y = val)) + 
  geom_line() + 
  scale_x_date(date_labels = "%b") +
  #scale_x_datetime(breaks = pad.flow$datetime, labels = monthplot) +
  labs(y="Average Flow Rate (m3/s)", x="Month", title = paste0(creek, " Average Flow Rate (m3/s): 2015-2022"), colour = "Year") + 
  theme_bw()

return(p3)
}

plot_year_avg_discrete <- function(newflow2, flow.discrete, creek) {
  
p4 <- ggplot(newflow2, aes(x = datetime, y = val)) + 
  geom_line() + 
  scale_x_date(date_labels = "%b") +
  #scale_x_datetime(breaks = pad.flow$datetime, labels = monthplot) +
  labs(y="Average Flow Rate (m3/s)", x="Month", title = paste0(creek, " Average Flow Rate (m3/s): 2015-2022"), colour = "Year") + 
  theme_bw() +
  geom_point(data=flow.discrete, aes(x=datetime, y=flow_m3s), size=3, bg="blue", pch=21) 

return(p4)
}

