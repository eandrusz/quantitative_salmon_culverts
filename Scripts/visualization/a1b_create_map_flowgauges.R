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
library(raster)
library(rosm)
library(prettymapr)
library(sf)
library(ggspatial)
library(ggrepel)
library(cowplot)
library(geosphere)

####################################################################
# Read in lat/lon for sampling locations and stream gauges
####################################################################

# sampling locations
sitename <- rep(c("Chuckanut","Padden","Squalicum"), each=3)
station <- rep(c("Down","Up","Gauge"), times=3)
lats <- c(48.689576,48.690745,48.70247, 48.71499,48.713933,48.71494, 48.800163,48.799909,48.76606)
longs <- c(-122.409068,-122.409447,-122.4824,-122.478996,-122.478387,-122.4996, -122.405012,-122.404188,-122.4996)
sample_locs <- data.frame(sitename,station,longs,lats)

# # gauge locations
# sitename2 <- c("Chuckanut","Padden","Squalicum")
# lats2 <- c(48.70247,48.71494,48.76606)
# longs2 <- c(-122.4824,-122.4996,-122.4996)
# gauge_locs <- data.frame(sitename2,lats2,longs2)

map_names <- sample_locs %>%
  st_as_sf(coords = c("longs", "lats"), crs = 4326, remove = F)  #make into spatial object w package `sf`, which you'll need for all kinds of spatial analysis

map_labels <- sample_locs %>%
  group_by(sitename) %>% 
  summarise(avglat = mean(lats, na.rm = T), avglong = mean(longs, na.rm = T)) %>% 
  replace(2, c(48.69426,48.72,48.78871)) %>% 
  st_as_sf(coords = c("avglong", "avglat"), crs = 4326, remove = F)  #make into spatial object w package `sf`, which you'll need for all kinds of spatial analysis


####################################################################
# Haversine distances between things
####################################################################
dns <- sample_locs %>% filter(station=="Down") %>% dplyr::select(longs,lats)
ups <- sample_locs %>% filter(station=="Gauge") %>% dplyr::select(longs,lats)

distmat <- apply(dns, 1, FUN=function(X) distHaversine(X, ups))
distdowngauge <- diag(distmat)
names(distdowngauge) <- c("Chuckanut","Padden","Squalicum")

mean_distdowngauge <- mean(distdowngauge)

####################################################################
# MAP 3: All but portage - up and down 
####################################################################

#make bounding box
map_bounding <-makebbox(w = min(map_names$longs),
                         e = max(map_names$longs),
                         n = max(map_names$lats),
                         s = min(map_names$lats))

#download data for map
map_raster <- osm.raster(map_bounding, 
                          type = "cartolight",
                          projection=4326)  

ggplot() +
  layer_spatial(map_raster) +
  geom_point(aes(y = map_names$lats,
                 x = map_names$longs,
                 color= factor(map_names$sitename),
                 pch=factor(map_names$station)),
             size=3) +
  xlab("") + ylab("") +
  guides(color='none', pch=guide_legend(title="Station")) + #, override.aes = aes(label = ""))) +
  geom_label_repel(data = map_labels,
                   stat = "sf_coordinates",
                   aes(geometry = geometry,
                       label = sitename,
                       fill = factor(sitename)),
                   alpha = 0.7,
                   label.r = 0,
                   min.segment.length = 0,
                   segment.size = 0) +
  guides(fill=guide_legend(title="Creek", override.aes = aes(label = ""))) +
  theme_map() + 
  theme(panel.border = element_rect(color = "black", fill=NA, size = 1)) +
  annotation_scale(unit_category = "metric",
                   width_hint = 0.15,
                   style = "ticks",
                   location = "bl") +
  annotation_north_arrow(which_north = "true",
                         height = unit(1, "cm"),
                         location = "tl",
                         style = north_arrow_nautical(text_col = "grey20"))

ggsave(here("Output","SupplementalFigures","map_gauges.png"))
