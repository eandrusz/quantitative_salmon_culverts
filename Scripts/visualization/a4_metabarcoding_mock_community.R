## NGN Visualizing metabarcoding results before / after Stan model
# Author: Eily Allan 
# Person running: Eily
# Last modified: 10/19/22 by Eily
# Date of run: 10/19/22 by Eily 

# Overview 
# This script 

# Inputs: 
# 1) 

# Outputs: 
# 1) 

####################################################################
# Set up
####################################################################

# Load packages
library(here)
library(tidyverse)
library(unikn)

metabeforeqm <- read.csv(here("Output","metabarcoding", "taxa_table.csv"))
#bayes_out <- readRDS("/Users/elizabethandruszkiewicz/Desktop/20221019-ngn-model-output/bayes_out_3spp.RDS")
bayes_out_prt <- readRDS(here("Output","metabarcoding","bayes_out_4salmonids_1Prt.RDS"))
bayes_out_brn <- readRDS(here("Output","metabarcoding","bayes_out_4salmonids_2Brn.RDS"))
bayes_out_chk <- readRDS(here("Output","metabarcoding","bayes_out_4salmonids_3Chk.RDS"))
bayes_out_pad <- readRDS(here("Output","metabarcoding","bayes_out_4salmonids_4Pad.RDS"))
bayes_out_sqm <- readRDS(here("Output","metabarcoding","bayes_out_4salmonids_5Sqm.RDS"))


####################################################################
# Set color palette
####################################################################

# Set color palette so don't change when plot
o_i_colors <- c(#rgb(230, 159,   0, maxColorValue = 255),  # orange
  rgb( 86, 180, 233, maxColorValue = 255),  # skyblue
  rgb(  0, 158, 115, maxColorValue = 255),  # green
  #rgb(240, 228,  66, maxColorValue = 255),  # yellow
  rgb(  0, 114, 178, maxColorValue = 255),  # blue
  rgb(204, 121, 167, maxColorValue = 255)   # purple
)
pal_okabe_ito <- newpal(col = o_i_colors,
                        names = c("Oncorhynchus clarkii","Oncorhynchus kisutch","Oncorhynchus mykiss","Oncorhynchus nerka"))

####################################################################
# Plot proportions of salmonids before QM correction
####################################################################

metabeforeqmplot <- metabeforeqm %>% 
  filter(! str_detect(Sample_name, "MC")) %>% 
  filter(! str_detect(Sample_name, "Kangaroo")) %>% 
  filter(! str_detect(Sample_name, "Up5")) %>% 
  filter(species %in% c("Oncorhynchus clarkii","Oncorhynchus kisutch","Oncorhynchus mykiss","Oncorhynchus nerka")) %>%  ### ONLY FOUR
  # group_by(Sample_name) %>% 
  # mutate(ReadDepth = sum(totalReads)) %>% 
  # mutate(propReads = totalReads/ReadDepth) %>% 
  separate(Sample_name, c("marker","time","creek","station","bio", "tech")) %>% 
  dplyr::select(-marker) %>% 
  mutate(station = case_when(station == "Up11" ~ "Up",
                             TRUE ~ station)) %>%
  mutate(station = case_when(station == "Dn" ~ "Down",
                             TRUE ~ station)) %>%
  mutate(creek = case_when(creek == "1Prt" ~ "Portage",
                           creek == "2Brn" ~ "Barnes",
                           creek == "3Chk" ~ "Chuckanut",
                           creek == "4Pad" ~ "Padden",
                           creek == "5Sqm" ~ "Squalicum",
                           TRUE ~ creek))  %>% 
  separate(time, into = c("month","year"), sep = 2, remove=FALSE) %>% 
  unite(newtime, c(year,month), sep="-", remove=FALSE) %>% 
  filter(tech==1) %>%
  unite(creekstntime, c(creek,station,time), remove=FALSE) %>% 
  group_by(creekstntime) %>% 
  mutate(ReadDepth = sum(totalReads)) %>% 
  mutate(propReads = totalReads/ReadDepth) %>%
  mutate(facetorder = factor(creek, levels=c('Padden','Portage','Chuckanut','Squalicum', 'Barnes')))

ggplot(metabeforeqmplot, aes(x=newtime, y=propReads, fill=species, color=species)) +
  geom_col() +
  facet_grid(~station~facetorder) +
  scale_fill_manual(values = pal_okabe_ito) +
  scale_color_manual(values = pal_okabe_ito) +
  guides(color= 'none') + 
  labs(x="Date (YY-MM)", y= "Proportion of DNA before QM correction") + 
  theme_bw() + 
  scale_x_discrete(guide = guide_axis(angle = -45))

# ggsave(here("Output","Figures","20221123_proportions_before_qm.png"), units="in", width=12, height=5)


####################################################################
# Plot proportions of salmonids after QM correction
####################################################################

prtpost <- bayes_out_prt$Bayes_estimates %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") 

brnpost <- bayes_out_brn$Bayes_estimates %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") 

chkpost <- bayes_out_chk$Bayes_estimates %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") 

padpost <- bayes_out_pad$Bayes_estimates %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") 

sqmpost <- bayes_out_sqm$Bayes_estimates %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") 

allpost <- rbind(prtpost,brnpost,chkpost,padpost,sqmpost)

postplot <- allpost %>%
  separate(col = sample, into = c("time", "creek", "station", "biol"), remove = FALSE) %>% 
  mutate(station = case_when(station == 1 ~ "Down",
                             station == 2 ~ "Up")) %>%
  mutate(creek = case_when(creek == "1Prt" ~ "Portage",
                           creek == "2Brn" ~ "Barnes",
                           creek == "3Chk" ~ "Chuckanut",
                           creek == "4Pad" ~ "Padden",
                           creek == "5Sqm" ~ "Squalicum",
                           TRUE ~ creek)) %>% 
  unite(creekstntime, c(creek,station,time), remove=FALSE) %>% 
  group_by(creekstntime) %>% 
  mutate(sumprop = sum(value)) %>% 
  mutate(avgprop = value/sumprop) %>%
  mutate(facetorder = factor(creek, levels=c('Padden','Portage','Chuckanut','Squalicum', 'Barnes'))) %>% 
  filter(value > 0.001) %>%
  separate(time, into = c("month","year"), sep = 2, remove=FALSE) %>% 
  unite(newtime, c(year,month), sep="-", remove=FALSE)

ggplot(postplot, aes(x = newtime, fill = species, color = species, y = avgprop)) +
  geom_col() +
  facet_grid(~station ~facetorder) +
  scale_fill_manual(values = pal_okabe_ito) + 
  scale_color_manual(values = pal_okabe_ito) +
  guides(color= 'none') +
  ylab("Proportion of DNA after QM correction") +
  labs(x="Date (YY-MM)", fill="Species") %>% 
  scale_x_discrete(guide = guide_axis(angle = -45)) +
  theme_bw()

ggsave(here("Output","Figures","20221123_proportions_after_qm.png"), units="in", width=12, height=5)




# bayes_out$Bayes_estimates %>%
#   rownames_to_column("sample") %>%
#   filter(! str_detect(sample, "Pad5")) %>% 
#   pivot_longer(-sample, names_to = "species") %>%
#   separate(col = sample, into = c("time", "creek", "station", "biol"), remove = FALSE) %>% 
#   mutate(station = case_when(station == 1 ~ "Down",
#                              station == 2 ~ "Up")) %>%
#   mutate(creek = case_when(creek == "1Prt" ~ "Portage",
#                            creek == "2Brn" ~ "Barnes",
#                            creek == "3Chk" ~ "Chuckanut",
#                            creek == "4Pad11" ~ "Padden",
#                            creek == "5Sqm" ~ "Squalicum",
#                            TRUE ~ creek)) %>% 
#   unite(creekstntime, c(creek,station,time), remove=FALSE) %>% 
#   group_by(creekstntime) %>% 
#   mutate(sumprop = sum(value)) %>% 
#   mutate(avgprop = value/sumprop) %>%
#   mutate(facetorder = factor(creek, levels=c('Padden','Portage','Chuckanut','Squalicum', 'Barnes'))) %>% 
#   filter(value > 0.001) %>%
#   separate(time, into = c("month","year"), sep = 2, remove=FALSE) %>% 
#   unite(newtime, c(year,month), sep="-", remove=FALSE) %>% 
#   ggplot(aes(x = newtime, fill = species, color = species, y = avgprop)) +
#   geom_col() +
#   #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
#   facet_grid(~facetorder ~station) +
#   scale_fill_manual(values = pal_okabe_ito) + 
#   scale_color_manual(values = pal_okabe_ito) +
#   labs(x="Date (YY-MM)", y="Proportion of DNA after QM correction", fill="Species", color="") %>% 
#   scale_x_discrete(guide = guide_axis(angle = -45)) +
#   theme_bw()
# 
# # ggsave(here("Output","Figures","20221120_proportions_after_qm.png"))
# 
