## NGN Merging QM results with QPCR results -- FLOW CORRECTED 
# Author: Eily Allan and Ryan Kelly
# Person running: Eily
# Last modified: 10/20/22 by Eily
# Date of run: 10/20/22 by Eily 

# Overview 
# This script 

# Inputs: 
# 1) 

# Outputs: 
# 1) 

####################################################################
# Set up
####################################################################

## Load packages
library(tidyverse)
library(here)
library(gtools)

source(here("Scripts", "functions","sample_posteriors.R"))

## Pull posteriors of Bayesian model for QM
bayes_out_prt <- readRDS(here("Output","metabarcoding","bayes_out_4salmonids_1Prt.RDS"))
bayes_out_brn <- readRDS(here("Output","metabarcoding","bayes_out_4salmonids_2Brn.RDS"))
bayes_out_chk <- readRDS(here("Output","metabarcoding","bayes_out_4salmonids_3Chk.RDS"))
bayes_out_pad <- readRDS(here("Output","metabarcoding","bayes_out_4salmonids_4Pad.RDS"))
bayes_out_sqm <- readRDS(here("Output","metabarcoding","bayes_out_4salmonids_5Sqm.RDS"))

## Pull posteriors of Bayesian model for qPCR
qMod_out <- readRDS("/Users/elizabethandruszkiewicz/Desktop/20221129_model_output/cut_qMod_out.RDS")
cut_names <- readRDS(here("Output","qpcr","cut_modeled_conc.RDS"))

## Flow data to correct
flow <- read.csv(here("Output","qpcr","20221121_flowrates_touse.csv"))

####################################################################
# Deal with QM model output
####################################################################

b_prt <- qm_pull_posteriors(bayes_out_prt)
b_brn <- qm_pull_posteriors(bayes_out_brn)
b_chk <- qm_pull_posteriors(bayes_out_chk)
b_pad <- qm_pull_posteriors(bayes_out_pad)
b_sqm <- qm_pull_posteriors(bayes_out_sqm)

b <- rbind(b_prt, b_brn, b_chk, b_pad, b_sqm)

####################################################################
# Deal with qPCR model output 
####################################################################

cut_names <- cut_names %>%
  mutate(creek = case_when(creek == "4Pad" & station == "Dn" ~ "4Pad11",
                           creek == "4Pad" & station == "Up11" ~ "4Pad11",
                           TRUE ~ creek)) %>% 
  mutate(creek = case_when(creek == "4Pad" & station == "Up5" ~ "4Pad5",
                           TRUE ~ creek)) %>% 
  mutate(station = case_when(creek == "4Pad11" & station == "Up11" ~ "Up",
                             TRUE ~ station)) %>% 
  mutate(station = case_when(creek == "4Pad5" & station == "Up5" ~ "Up",
                             TRUE ~ station)) %>% 
  mutate(station = case_when(station == "Dn" ~ 1,
                             station == "Up" ~ 2)) %>%
  unite(Sample, c("time","creek","station","biorep"), sep="_") %>% 
  dplyr::select(c(Sample,dilution,Adj_Vol))

qgrabthese <- grep("envir_concentration", names(qMod_out$qMod@sim$samples[[1]]))
qpostSamples <- qMod_out$qMod@sim$samples[[1]][qgrabthese]


####################################################################
# Deal with flow
####################################################################

flowtojoin <- flow %>% 
  pivot_longer(!daysampletoavg, names_to="creek", values_to="flow_m3s") %>% 
  mutate(creek = case_when(creek == "pad_flowm3s" ~ "4Pad",
                           creek == "flow_chk_by_padden" ~ "3Chk",
                           creek == "flow_sqm_by_padden" ~ "5Sqm",
                           creek == "flow_prt_by_padden" ~ "1Prt",
                           creek == "flow_brn_by_padden" ~ "2Brn",
                           TRUE ~ creek)) %>% 
  separate(daysampletoavg, c("year","month","day")) %>% 
  mutate(year= str_sub(year, 3, 4))  %>%
  unite(time, c(month, year), sep="") %>% 
  unite(timecreek, c(time, creek)) 

qa <- as.data.frame(qpostSamples) 
qb <- t(qa) %>% 
  as_tibble() %>% 
  add_column(cut_names) %>% 
  #add_column(coho_names) %>% 
  group_by(Sample, dilution, Adj_Vol) %>% 
  nest() %>% 
  rename(cutqpcr = data) %>% 
  rename(bottle = Sample) %>% 
  # for flow correction
  separate(bottle, into=c("time","creek","bio","tech"), remove=FALSE) %>% 
  filter(creek != "4Pad5") %>% 
  mutate(creek = case_when(creek == "4Pad11" ~ "4Pad",
                           TRUE ~ creek)) %>% 
  unite(timecreek, c(time,creek), remove=FALSE) %>%  
  unite(bottle, c("time","creek","bio","tech"), remove=FALSE) %>% 
  left_join(flowtojoin, by="timecreek") %>% 
  select(-c(timecreek,bio,tech, bottle)) 

####################################################################
# Actually merge
####################################################################

mergedb <- b %>% 
  left_join(qb, by="bottle") %>% 
  #filter(species == "Oncorhynchus clarkii") %>% 
  drop_na() %>% 
  #replace_na(list(cutqpcr = list(-0.30103), dilution=1, Adj_Vol = 2)) %>% ## back calculate -0.30103 to get final of 50 copies / L water for samples where there are reads but no qpcr data 
  mutate(cut_dnacopy_uL = map2(.x = cutqpcr, .y = dilution, .f = function(x, y) y*(10^x))) %>% 
  mutate(cut_dnacopy_L = map2(.x = cut_dnacopy_uL, .y = Adj_Vol, .f = function(x,y) x*(200/y))) %>% 
  ungroup()

mergedb <- mergedb %>%
  group_by(bottle, species) %>% 
  mutate(meanpropreads = mean(unlist(data))) %>%
  mutate(meancutqpcr = mean(unlist(cut_dnacopy_uL))) %>% 
  mutate(meancutdna = mean(unlist(cut_dnacopy_L))) %>% 
  separate(bottle, into=c("time","creek","station", "bio"), remove=FALSE)


#calculate mean total DNA, only using cutthroat
meantotdna <- mergedb %>% 
  filter(species == "Oncorhynchus clarkii") %>% 
  mutate(meantotdna = meancutdna/meanpropreads) %>% 
  dplyr::select(bottle, meantotdna)

#going old-school to get this to work
mergedc <- mergedb
mergedc$meantotdna <- meantotdna$meantotdna[match(mergedc$bottle, meantotdna$bottle)]

#check; are all values for total dna the same?
mergedc %>% 
  dplyr::select(bottle, species, meanpropreads, meancutdna, meantotdna) %>% 
  filter(bottle == "0421_4Pad_2_1")

mergedc <- mergedc %>% 
  ungroup() %>% 
  mutate(meandnaconc = meanpropreads*meantotdna)

#check; do concentration values make sense?
mergedc %>% 
  dplyr::select(bottle, species, meanpropreads, meancutdna, meantotdna, meandnaconc) %>% 
  filter(bottle == "0421_4Pad_1_1")


quants_to_plot <- mergedc %>% 
  mutate(station = case_when(station == 1 ~ "Down",
                             station == 2 ~ "Up")) %>%
  separate(time, into = c("month","year"), sep = 2) %>% 
  unite(newtime, c(year,month), sep="-") %>% 
  mutate(creek = case_when(creek == "1Prt" ~ "Portage",
                           creek == "2Brn" ~ "Barnes",
                           creek == "3Chk" ~ "Chuckanut",
                           creek == "4Pad" ~ "Padden",
                           creek == "5Sqm" ~ "Squalicum",
                           TRUE ~ creek)) %>% 
  # add on flow multiplier
  mutate(meandnaconcflow = meandnaconc*flow_m3s*1000) # convert 1000 L = m3

ggplot(quants_to_plot, aes(x=newtime, y=log10(meandnaconcflow), color=species)) +
  geom_point() +
  facet_grid(~creek ~station ~species) +
  # scale_color_manual(values = pal_okabe_ito) +
  guides(color= 'none') + 
  labs(x="Date (YY-MM)", y= "Log10 copies/s") + 
  theme_bw() + 
  scale_x_discrete(guide = guide_axis(angle = -45))

ggplot(quants_to_plot, aes(x=newtime, y=log10(meantotdna))) +
  geom_point() +
  facet_grid(~creek ~station) +
  labs(x="Date (YY-MM)", y= "Log10 copies total salmonid DNA/s") + 
  theme_bw() + 
  scale_x_discrete(guide = guide_axis(angle = -45))

#ggsave(here("Output","SupplementalFigures","totaldna_ts_after_qm.png"))


simple <- quants_to_plot %>% 
  select(c(bottle, newtime, creek, station, bio, species, meanpropreads, meancutdna, meantotdna, meandnaconc, meandnaconcflow))

#write_rds(simple, here("Output","20221129_abundance_flowcorrected.RDS"))



# 
# ## filter out <1% reads and BLOQ 
# 
# key <- mergedb %>% 
#   ungroup() %>% 
#   #dplyr::select(c(bottle, total_dnacopy_L))
#   dplyr::select(c(bottle, meantotdna))
# 
# all_salmonids <- key %>% 
#   left_join(b, by="bottle") %>%
#   mutate(absconc = map2(.x = meantotdna, .y = data, .f = function(x, y) x*y)) %>% 
#   select(c(bottle,species,absconc)) %>% 
#   ungroup() %>% 
#   mutate(meanconcspecies = mean(unlist(absconc)))
# 
# ## try to streamline a little 
# all_salmonids_by_species <- all_salmonids %>% 
#   group_by(species) %>% 
#   unnest() %>% 
#   separate(bottle, into=c("time","creek","station","bio")) 
# 
# 
# all_salmonids_by_species_wbios <- all_salmonids %>% 
#   separate(bottle, into=c("time","creek","station","bio")) %>% 
#   group_by(time, creek, station, bio, species) %>% 
#   mutate(meanabsconc = mean(unlist(absconc)))
# 
# ### EXPORT THIS TO MANIPULATE 
# write_rds(all_salmonids_by_species_wbios, here("Output","salmonids_abs_abundance_bio_posterior.RDS"))
# 
# check_cut <- all_salmonids_by_species_wbios %>%
#   filter(str_detect(species, "clarkii")) %>%
#   select(-species) %>%
#   unite(bottle, c(time,creek,station,bio))
# 
# check_cut_qpcr <- all_salmonids %>%
#   filter(str_detect(species, "clarkii")) %>%
#   group_by(bottle, species) %>%
#   mutate(meancutabs = mean(unlist(absconc))) %>%
#   left_join(check_cut, by = "bottle") %>%
#   mutate(checkdiff = (meancutabs - meanabsconc))
# 
# # now average biological replciates
# all_salmonids_by_species2 <- all_salmonids_by_species %>% 
#   group_by(time, creek, station, species) %>% 
#   summarise_at(.vars = colnames(all_salmonids_by_species)[6:1505],
#                .funs = c(mean="mean")) 
# 
# ### EXPORT THIS TO MANIPULATE BEFORE DOING UP / DOWN FOLD CHANGE
# write_rds(all_salmonids_by_species2, here("Output","salmonids_abs_abundance_posterior.RDS"))
# 
# 
# ## NOW GO ON TO DO FOLD CHANGES (per species, down / up)
# all_salmonids_by_species3 <- all_salmonids_by_species2 %>% 
#   group_by(time, creek, station, species) %>% 
#   nest() %>% 
#   pivot_wider(names_from=station, values_from = data) %>% 
#   rename(Dn = '1') %>% 
#   rename(Up = '2') %>% 
#   mutate(fc_d_u = map2(.x = Dn, .y = Up, .f = foldchange)) %>% 
#   dplyr::select(-c(Dn,Up)) %>% 
#   filter(is.logical(unlist(fc_d_u))==FALSE) 
# 
# alpha = 0.25
# degree.freedom = 1400
# t.score = qt(p=alpha/2, df=degree.freedom, lower.tail=F)
# 
# mean_magnitudes <- all_salmonids_by_species2 %>% 
#   group_by(time, creek, station, species) %>% 
#   nest() %>% 
#   pivot_wider(names_from=station, values_from = data) %>% 
#   rename(Dn = '1') %>% 
#   rename(Up = '2') %>% 
#   mutate(mean_down = mean(unlist(Dn))) %>% 
#   mutate(mean_up = mean(unlist(Up))) %>% 
#   dplyr::select(-c(Dn,Up)) %>% 
#   unite(timecreekspecies, c(time,creek,species), remove=TRUE)
#   
# all_salmonids_by_species4 <- all_salmonids_by_species3 %>%
#   mutate(mean_fc = mean(unlist(fc_d_u))) %>% 
#   mutate(sd_fc = sd(unlist(fc_d_u))) %>% 
#   dplyr::select(c(time,creek,species,mean_fc, sd_fc)) %>% 
#   mutate(Sign = sign(mean_fc)) %>% 
#   mutate(se_fc = sd_fc/sqrt(1500)) %>% 
#   mutate(me2575_fc = se_fc*t.score) %>% 
#   mutate(lb25_fc = mean_fc - me2575_fc) %>% 
#   mutate(ub75_fc = mean_fc + me2575_fc) %>% 
#   unite(timecreekspecies, c(time,creek,species), remove=FALSE) %>% 
#   left_join(mean_magnitudes, by="timecreekspecies") %>% 
#   dplyr::select(-timecreekspecies)
# 
# write_rds(all_salmonids_by_species4, here("Output", "salmonids_fc.RDS"))  
# 
# ggplot(all_salmonids_by_species4, aes(x=time, y=log10(abs(mean_fc)), color=as.character(Sign), size=log10(mean_down))) +
#   geom_point() +
#   geom_hline(yintercept=0) +
#   #geom_segment(aes(x = time, xend = time, y = log10(abs(lb25_fc))*Sign, yend = log10(abs(ub75_fc))*Sign)) +
#   facet_grid(~creek ~species) +
#   scale_color_discrete(name = "Legend", labels = c("Negative (Up > Down)", "Positive (Down > Up)")) +
#   labs(y="Log10(asb(mean fold change))", x="Month")
