## NGN Run QM Model
# Author: Eily Allan - but really all from Ryan Kelly and Ole Shelton 
# Person running: Eily
# Last modified: 10/17/22 by Eily
# Date of run: 10/17/22 by Eily 

# Overview 
# 

# Inputs: 
# 1)  

# Outputs: 
# 1) 

####################################################################
# Set up
####################################################################

## Load packages
library(tidyverse)
library(rstan)
library(MCMCpack) #for rdirichelet function
library(here)
library(gridExtra)
library(unikn)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source(here("Scripts","functions", "calibrate_metabarcoding.R")) 

# read in taxa table 
taxa_table <- read.csv(here("Output","metabarcoding", "taxa_table.csv"))
mock <- readRDS(here("Output","metabarcoding","20221118_mockdatatocalibrate.RDS"))

enviro <- taxa_table %>% 
  filter(! str_detect(Sample_name, "MC")) %>% 
  separate(Sample_name, into=c("marker", "time", "creek", "station", "biol", "tech"), remove=TRUE) %>% 
  filter(station != "Up5") %>% 
  mutate(station = case_when(station == "Up11" ~ "Up", 
                              TRUE ~ station)) %>% 
  dplyr::rename(Nreads = totalReads) %>% 
  dplyr::select(-marker)

# only focus on four salmonids
mock <- mock %>%
  filter(species %in% c("Oncorhynchus clarkii","Oncorhynchus kisutch","Oncorhynchus mykiss","Oncorhynchus nerka"))

# Prepare for stan model 
qmdata <- format_metabarcoding_data(enviro, mock)
stan_metabarcoding_data <- makeDesign(qmdata, N_pcr_cycles = 43)

# Set color palette so don't change when plot
o_i_colors <- c(#rgb(230, 159,   0, maxColorValue = 255),  # orange
                rgb( 86, 180, 233, maxColorValue = 255),  # skyblue
                rgb(  0, 158, 115, maxColorValue = 255),  # green
                #rgb(240, 228,  66, maxColorValue = 255),  # yellow
                rgb(  0, 114, 178, maxColorValue = 255),  # blue
                rgb(204, 121, 167, maxColorValue = 255)   # purple
)
#o_i_colors = scales::hue_pal()(length(unique(mock$species)))
pal_okabe_ito <- newpal(col = o_i_colors,
                        names = unique(mock$species))

####################################################################
# Run ML QM model
####################################################################

ML_out <- QM_likelihood(here("Scripts", "functions", "quant_metabar_rosetta_noSampleEta.stan"), stan_metabarcoding_data)

#write_rds(ML_out, here("Output","metabarcoding","ML_out.RDS"))

ML_out$ML_estimates %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") %>%
  separate(col = sample, into = c("time", "creek", "station", "biol"), remove = FALSE) %>% 
  mutate(station = case_when(station == 1 ~ "Down",
                             station == 2 ~ "Up")) %>% 
  filter(value > 0.001) %>%
  ggplot(aes(x = biol, fill = species, y = value)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(~creek ~station ~time) +
  scale_fill_manual(values = pal_okabe_ito)

#ggsave(here("Output","Figures","ML_add_samples","allspecies","n200.png"))

####################################################################
# Run Bayesian QM model
####################################################################

bayes_out <- QM_bayes(here("Scripts", "functions", "quant_metabar_rosetta_noSampleEta.stan"), stan_metabarcoding_data)
#write_rds(bayes_out, "/Users/elizabethandruszkiewicz/Desktop/20221121_bayesout_4spp.RDS")

summaryout <- summary(bayes_out$Bayes_modelfit)$summary
#write.csv(summaryout, "/Users/elizabethandruszkiewicz/Desktop/20221121_bayesout_4spp_summary.csv")


bayes_out$Bayes_estimates %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") %>%
  separate(col = sample, into = c("time", "creek", "station", "biol"), remove = FALSE) %>% 
  mutate(station = case_when(station == 1 ~ "Down",
                             station == 2 ~ "Up")) %>% 
  separate(time, into = c("month","year"), sep = 2, remove=FALSE) %>% 
  unite(newtime, c(year,month), sep="-", remove=FALSE) %>% 
  filter(value > 0.001) %>%
  mutate(creek = case_when(creek == "1Prt" ~ "Portage",
                           creek == "2Brn" ~ "Barnes",
                           creek == "3Chk" ~ "Chuckanut",
                           creek == "4Pad" ~ "Padden",
                           creek == "5Sqm" ~ "Squalicum",
                           TRUE ~ creek))  %>% 
  mutate(facetorder = factor(creek, levels=c('Padden','Portage','Chuckanut','Squalicum', 'Barnes'))) %>% 
  unite(creekstntime, c(creek,station,time), remove=FALSE) %>% 
  group_by(creekstntime) %>% 
  mutate(sumprop = sum(value)) %>% 
  mutate(avgprop = value/sumprop) %>%
  ggplot(aes(x = newtime, fill = species, y = avgprop)) +
  geom_col() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
  facet_grid(~facetorder ~station) +
  scale_fill_manual(values = pal_okabe_ito) + 
  scale_color_manual(values = pal_okabe_ito) +
  ylab("Proportion of DNA after QM correction") +
  labs(x="Date (YY-MM)", y="Proportion of DNA after QM correction", fill="Species", color="") %>% 
  scale_x_discrete(guide = guide_axis(angle = -45)) +
  theme_bw()

#ggsave(here("Output","Figures","20221121_proportions_after_qm.png"))
