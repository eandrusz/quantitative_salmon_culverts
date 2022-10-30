## NGN Map for Fig 1
# Author: Eily Allan 
# Person running: Eily
# Last modified: 10/19/22 by Eily
# Date of run: 10/19/22 by Eily 

# Overview 
# This script 

# Inputs: 
# 1) 

# Outputs: 
# 1) map

####################################################################
# Set up
####################################################################

# Load packages
library(here)
library(tidyverse)
library(ggmap)
library(scales)

####################################################################
# Read in lat/lon for sampling locations and stream gauges
####################################################################


# sampling locations
creeks <- rep(c("Portage","Barnes","Chuckanut","Padden","Squalicum"), each=2)
stations <- rep(c("Down","Up"), times=5)
lats <- c(48.18282,48.183108,48.665383,48.665415,48.689576,48.690745,48.71499,48.713933,48.800163,48.799909)
lons <- c(-122.130534,-122.128588,-122.373456,-122.368946,-122.409068,-122.409447,-122.478996,-122.478387,-122.405012,-122.404188)

sample_locs <- data.frame(creeks,stations,lats,lons)
sample_locs$type = "Sampling Site"

# gauge locations
creeks <- c("Chuckanut","Padden","Squalicum")
stations <- rep("Flow Gauge", times=3)
lats <- c(48.7024722,48.71494,48.76606)
lons <- c(-122.4824,-122.4996, -122.4996)

gauge_locs <- data.frame(creeks,stations,lats,lons)
gauge_locs$type = "Gauge"

all_locs <- rbind(sample_locs, gauge_locs)
down_locs <- sample_locs %>% filter(stations=="Down")
padden_locs <- all_locs %>% filter(creeks=="Padden")

####################################################################
# Get maps of Seattle / Bellingham (3 zooms)
####################################################################

#ggmap API
register_google(key = "AIzaSyC4lEbtuwdIKbBJNKyfTN73sm3DzEKCYt0") ## THIS IS ABBYS
#register_google(key = "AIzaSyD-_Nnel-29IzTHJ7shlDv3yBu3h6WgA-8")

#background maps
general_loc <-  c(lon = mean(all_locs$lons), lat = mean(all_locs$lats)) 
allbutportage_loc <-  c(lon = mean(all_locs$lons[3:10]), lat = mean(all_locs$lats[3:10])) 
justbham_loc <- c(lon = mean(all_locs$lons[3:8]), lat = mean(all_locs$lats[3:8]))

# furthest zoomed out version 
general_map1 <- get_map(location=general_loc, 
                       source='google',
                       zoom = 6,
                       maptype = 'satellite',
                       crop=FALSE)

map_zoom1 <- ggmap(general_map1) +
  labs(x='Longitude',
       y='Latitude') +
  theme_minimal()

# middle zoomed out version 
general_map2 <- get_map(location=general_loc, 
                        source='google',
                        zoom = 9,
                        maptype = 'satellite',
                        crop=FALSE)

map_zoom2 <- ggmap(general_map2) +
  labs(x='Longitude',
       y='Latitude') +
  theme_minimal()

# all but portage zoomed out version 
general_map3 <- get_map(location=allbutportage_loc, 
                        source='google',
                        zoom = 11,
                        maptype = 'satellite',
                        crop=FALSE)

map_zoom3 <- ggmap(general_map3) +
  labs(x='Longitude',
       y='Latitude') +
  theme_minimal()

# just bellingham sites version 
general_map4 <- get_map(location=justbham_loc, 
                        source='google',
                        zoom = 12,
                        maptype = 'satellite',
                        crop=FALSE)

map_zoom4 <- ggmap(general_map4) +
  labs(x='Longitude',
       y='Latitude') +
  theme_minimal()

# just padden version 
general_map5 <- get_map(location=c(all_locs$lons[12],all_locs$lats[12]), 
                        source='google',
                        zoom = 14,
                        maptype = 'satellite',
                        crop=FALSE)

map_zoom5 <- ggmap(general_map5) +
  labs(x='Longitude',
       y='Latitude') +
  theme_minimal()

####################################################################
# Add points to maps 
####################################################################

map_zoom2 +
  geom_point(data=down_locs, #sites
             aes(x=lons,
                 y=lats,
                 fill=creeks),
             pch=21,
             color='black',
             alpha=0.6,
             size=5.25) 

map_zoom5 +
  geom_point(data=padden_locs, #sites
             aes(x=lons,
                 y=lats,
                 fill=stations),
             pch=21,
             color='black',
             alpha=0.6,
             size=5.25) 
