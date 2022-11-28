
### FLOW functions 


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
    mutate(datetime = lubridate::make_datetime(yearplot, monthplot, dayplot, hourplot, minplot)) %>% 
    filter(datetime < "2021-03-01") %>% 
    mutate(datetime = lubridate::make_datetime(2022, monthplot, dayplot, hourplot, minplot)) %>%
    drop_na() 
  
  return(newflow)
}

format_flow2 <- function(flow_df){
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
    mutate(datetime = lubridate::make_datetime(yearplot, monthplot, dayplot, hourplot, minplot)) %>% 
    drop_na() %>% 
    filter(datetime > "2021-03-01" & datetime <"2022-03-01") 
  
  return(newflow)
}


####################################################################
# Write function to find average daily flow over a few years
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
# Write function to find average monthly flow over a few years
####################################################################

monthly_avg_flow <- function(newflow){
  require(tidyverse)
  require(lubridate)
  
  newflow2 <- newflow %>% 
    group_by(monthplot) %>% 
    summarise(val = mean(flow_m3s, na.rm = T)) %>%
    mutate(datetime = make_date(2022, monthplot,15)) 
  
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
    labs(x="Date", y=bquote('Discharge '(m^3/s)), title = creek)
  
  return(p1)
}


plot_facet_year <- function(newflow, creek){
  require(ggplot2)
  
  p2 <- ggplot(newflow) +
    geom_line(aes(x = datetime, y = flow_m3s, color = factor(yearplot))) +
    labs(y=bquote('Discharge '(m^3/s)), x="Date", title = creek, colour = "Year") + 
    facet_wrap(~yearplot, nrow=7) +
    guides(color='none') +
    theme_bw()
  
  return(p2)
}

plot_year_avg <- function(newflow2, creek){
  require(ggplot2)
  
  p3 <- ggplot(newflow2, aes(x = datetime, y = val)) + 
    geom_line() + 
    scale_x_date(date_labels = "%b") +
    #scale_x_datetime(breaks = pad.flow$datetime, labels = monthplot) +
    labs(y=bquote('Average Discharge '(m^3/s)), x="Month", title = paste0(creek, " Average Discharge: 2015-2021"), colour = "Year") + 
    theme_bw()
  
  return(p3)
}

plot_year_avg_discrete <- function(newflow2, flow.discrete, creek) {
  
  p4 <- ggplot(newflow2, aes(x = datetime, y = val)) + 
    geom_line() + 
    scale_x_date(date_labels = "%b") +
    #scale_x_datetime(breaks = pad.flow$datetime, labels = monthplot) +
    labs(y=bquote('Average Discharge '(m^3/s)), x="Month", title = paste0(creek, " Average Discharge: 2015-2021"), colour = "Year") + 
    theme_bw() +
    geom_point(data=flow.discrete, aes(x=datetime, y=flow_m3s), size=3, bg="blue", pch=21) 
  
  return(p4)
}

plot_sampling_discrete <- function(newflow2, flow.discrete, creek) {
  
  p5 <- ggplot(newflow2, aes(x = timeplot, y = flow_m3s)) + 
    geom_line() + 
    #scale_x_date(date_labels = "%b") +
    #scale_x_datetime(breaks = pad.flow$datetime, labels = monthplot) +
    labs(y=bquote('Average Discharge '(m^3/s)), x="Month", title = paste0(creek, " Discharge: March 2021-March 2022"), colour = "Year") + 
    theme_bw() +
    geom_point(data=flow.discrete, aes(x=timeplot, y=flow_m3s), size=3, bg="blue", pch=21) 
  
  return(p5)
}

