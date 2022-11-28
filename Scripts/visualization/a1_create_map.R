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
sitename <- rep(c("Portage","Barnes","Chuckanut","Padden","Squalicum"), each=2)
station <- rep(c("Down","Up"), times=5)
lats <- c(48.18282,48.183108,48.665383,48.665415,48.689576,48.690745,48.71499,48.713933,48.800163,48.799909)
longs <- c(-122.130534,-122.128588,-122.373456,-122.368946,-122.409068,-122.409447,-122.478996,-122.478387,-122.405012,-122.404188)

sample_locs <- data.frame(sitename,station,lats,longs)
#sample_locs$type = "Sampling Site"

map2_names <- sample_locs %>%
  filter(station=="Up") %>% 
  st_as_sf(coords = c("longs", "lats"), crs = 4326, remove = F)  #make into spatial object w package `sf`, which you'll need for all kinds of spatial analysis

map3_names <- sample_locs %>%
  filter(! sitename %in% c("Seattle","Portage")) %>% 
  st_as_sf(coords = c("longs", "lats"), crs = 4326, remove = F)  #make into spatial object w package `sf`, which you'll need for all kinds of spatial analysis

map3_names2a <- sample_locs %>%
  filter(station=="Up") %>%
  filter(sitename != "Portage") %>% 
  st_as_sf(coords = c("longs", "lats"), crs = 4326, remove = F)  #make into spatial object w package `sf`, which you'll need for all kinds of spatial analysis

map3_names2b <- sample_locs %>%
  filter(station=="Down") %>%
  filter(sitename != "Portage") %>% 
  st_as_sf(coords = c("longs", "lats"), crs = 4326, remove = F)  #make into spatial object w package `sf`, which you'll need for all kinds of spatial analysis

         
Seattle <- data.frame(Site = "Seattle",
                      lats = 47.60518519980011,
                      longs = -122.33924347206299)

Seattle <- st_as_sf(Seattle, coords = c("longs", "lats"), crs = 4326, remove = F)

####################################################################
# Distance between up and down stations
####################################################################
dns <- sample_locs %>% filter(station=="Down") %>% dplyr::select(longs,lats)
ups <- sample_locs %>% filter(station=="Up") %>% dplyr::select(longs,lats)

distmat <- apply(dns, 1, FUN=function(X) distHaversine(X, ups))
rownames(distmat) <- c("Portage","Barnes","Chuckanut","Padden","Squalicum")
distupdown <- diag(distmat)
names(distupdown) <- c("Portage","Barnes","Chuckanut","Padden","Squalicum")

mean_distupdown <- mean(distupdown)

####################################################################
# MAP 2: Middle zoom, all sampling sites
####################################################################

#make bounding box
# map2_bounding <-makebbox(w = min(all_locs$longs),
#                          e = max(all_locs$longs),
#                          n = max(all_locs$lats),
#                          s = min(all_locs$lats))

map2_bounding <-makebbox(w = -123,
                         e = -121,
                         n = 49,
                         s = 47.7)


#download data for map
map2_raster <- osm.raster(map2_bounding, 
                          type = "cartolight",
                          projection=4326)  

map2 <- ggplot() +
  layer_spatial(map2_raster) +
  geom_point(aes(y = map2_names$lats,
                 x = map2_names$longs)) +
  xlab("") + ylab("") +
  geom_label_repel(data = map2_names, 
                   stat = "sf_coordinates",
                   aes(geometry = geometry, 
                       label = sitename,
                       pch=1),
                   #nudge_x = -0.2,
                   fill = "grey90",
                   alpha = 0.7,
                   label.r = 0,
                   min.segment.length = 0,
                   #direction    = "y",
                   #hjust        = 1,
                   segment.size = 0.4) +
  geom_sf(data = Seattle, color = "gray40", size = 3) +
  geom_text_repel(data = Seattle,
                  stat = "sf_coordinates",
                  aes(geometry = geometry,
                      pch=8,
                      label = Site),
                  color = "gray40",
                  direction = "both",
                  nudge_y = -0.03,
                  nudge_x = -0.05,
                  hjust = 0,
                  #alpha = 0.8,
                  min.segment.length = 1) +
  geom_rect(aes(xmin = map3_bounding[1,1], xmax = map3_bounding[1,2], ymin = map3_bounding[2,1], ymax = map3_bounding[2,2]), alpha = .005, color="red") +
  theme_bw() +
  annotation_scale(unit_category = "metric",
                   #bar_cols = grey.colors(6)[c(3,6)], 
                   #line_width = unit(2, "cm"),
                   width_hint = 0.15,
                   style = "ticks",
                   location = "bl") +
  annotation_north_arrow(which_north = "true",
                         height = unit(1, "cm"),
                         location = "tr",
                         style = north_arrow_nautical(text_col = "grey20"))


####################################################################
# MAP 3: All but portage - up and down 
####################################################################

#make bounding box
map3_bounding <-makebbox(w = min(map3_names$longs),
                         e = max(map3_names$longs),
                         n = max(map3_names$lats),
                         s = min(map3_names$lats))

#download data for map
map3_raster <- osm.raster(map3_bounding, 
                          type = "cartolight",
                          projection=4326)  

map3 <- ggplot() +
  layer_spatial(map3_raster) +
  geom_point(aes(y = map3_names2a$lats,
                 x = map3_names2a$longs),
            pch = 1,
             size=3) +
  geom_point(aes(y = map3_names2b$lats,
                 x = map3_names2b$longs), 
                 pch = 2,
             size=3) +
  xlab("") + ylab("") +
  geom_label_repel(data = map3_names2a, 
                   stat = "sf_coordinates",
                   aes(geometry = geometry, 
                       #fill = sitename,
                       label = sitename),
                   fill = "grey90",
                   alpha = 0.7,
                   label.r = 0,
                   min.segment.length = 0,
                   segment.size = 0) +
  #labs(pch="Sampling Location") + 
  #guides(fill=guide_legend(title="Creek", override.aes = aes(label = ""))) +
  guides(fill="none", pch="none") +
  theme_bw() +
  theme_map() + 
  theme(panel.border = element_rect(color = "red", fill=NA, size = 2))


####################################################################
# PUT THEM TOGETHER 
####################################################################

plot.with.inset <-
  ggdraw() +
  draw_plot(map2) +
  draw_plot(map3,
            x = .6, y = .05, #location relative to background plot
            width = 0.4, #have inset take up a fifth of the width of the main figure
            height = 0.4*1.545) #maintain aspect ratio for inset plot
plot.with.inset

#ggsave(here("Output","Figures","SiteMap.png"))

######################################################################################################
######################################################################################################
######################################################################################################
# 
# # gauge locations
# creeks <- c("Chuckanut","Padden","Squalicum")
# stations <- rep("Flow Gauge", times=3)
# lats <- c(48.7024722,48.71494,48.76606)
# longs <- c(-122.4824,-122.4996, -122.4996)
# 
# gauge_locs <- data.frame(creeks,stations,lats,longs)
# gauge_locs$type = "Gauge"
# 
# all_locs <- rbind(sample_locs, gauge_locs)
# 
# pad_only <- all_locs %>% 
#   filter(creeks == "Padden")
# 
# ### RYan's map 
# 
# load("/Users/elizabethandruszkiewicz/Downloads/LocatorMap.Rdata")
# #add border to locator map
# locatorMap <- locatorMap +
#   theme(panel.border = element_rect(color = "black", size = 2))
# 
# 
# plot.with.inset <-
#   ggdraw() +
#   draw_plot(map2) +
#   draw_plot(locatorMap,
#             x = .6, y = .35, #location relative to background plot
#             width = 0.2, #have inset take up a fifth of the width of the main figure
#             height = 0.2*1.545) #maintain aspect ratio for inset plot
# plot.with.inset
# 
# 
# ####################################################################
# # MAP 1 Most zoomed out: Seattle / Bellingham 
# ####################################################################
# 
# #make bounding box
# map1_bounding <-makebbox(w = -125,
#                          e = -121,
#                          n = 49,
#                          s = 47)
# 
# #download data for map
# map1_raster <- osm.raster(map1_bounding, 
#                           type = "cartolight",
#                           projection=4326)  
# 
# map1_names <- data.frame( lats = c(47.639464, 48.74109), 
#                           longs = c(-122.332371, -122.45056),
#                           sitename = c("Seattle", "Bellingham")) %>%
#   st_as_sf(coords = c("longs", "lats"), crs = 4326, remove = F)  #make into spatial object w package `sf`, which you'll need for all kinds of spatial analysis
# 
# map1 <- ggplot() +
#   layer_spatial(map1_raster) +
#   geom_point(aes(y = map1_names$lats,
#                  x = map1_names$longs)) +
#   xlab("") + ylab("") +
#   geom_label_repel(data = map1_names, 
#                    stat = "sf_coordinates",
#                    aes(geometry = geometry, 
#                        label = sitename),
#                    #nudge_x = -0.2,
#                    fill = "grey90",
#                    alpha = 0.7,
#                    label.r = 0,
#                    min.segment.length = 0,
#                    #direction    = "y",
#                    #hjust        = 1,
#                    segment.size = 0.4) +
#   geom_rect(aes(xmin = map2_bounding[1,1], xmax = map2_bounding[1,2], ymin = map2_bounding[2,1], ymax = map2_bounding[2,2]), alpha = .005, color="red") +
#   theme_bw() +
#   theme_map() 
# 
# 
# ####################################################################
# # MAP 4: Most zoom, only Padden
# ####################################################################
# 
# #make bounding box
# map4_bounding <-makebbox(w = min(pad_only$longs),
#                          e = max(pad_only$longs),
#                          n = max(pad_only$lats),
#                          s = min(pad_only$lats))
# 
# map4_bounding <-makebbox(w = -122.55,
#                          e = -122.45,
#                          n = 48.72,
#                          s = 48.7)
# 
# #download data for map
# map4_raster <- osm.raster(map4_bounding, 
#                           type = "cartolight",
#                           projection=4326)  
# 
# map4_names <- pad_only %>%
#   st_as_sf(coords = c("longs", "lats"), crs = 4326, remove = F)  #make into spatial object w package `sf`, which you'll need for all kinds of spatial analysis
# 
# map4 <- ggplot() +
#   layer_spatial(map4_raster) +
#   geom_point(aes(y = map4_names$lats,
#                  x = map4_names$longs, 
#                  pch = map4_names$stations),
#              size=3) +
#   xlab("") + ylab("") +
#   geom_label_repel(data = map4_names, 
#                    stat = "sf_coordinates",
#                    aes(geometry = geometry, 
#                        fill = creeks,
#                        label = stations),
#                    alpha = 0.7,
#                    label.r = 0,
#                    min.segment.length = 0,
#                    segment.size = 0.4) +
#   labs(pch="Sampling Location") + 
#   guides(fill=guide_legend(title="Creek", override.aes = aes(label = ""))) +
#   theme_bw()
# 
