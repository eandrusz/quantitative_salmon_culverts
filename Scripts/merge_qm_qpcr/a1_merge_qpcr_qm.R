## NGN Merging QM results with QPCR results
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

## Pull posteriors of Bayesian model for QM

# Ryan's latest

#bayes_out <- readRDS(here("Output/metabarcoding/bayes_out.RDS"))

# My latest
bayes_out <- readRDS("/Users/elizabethandruszkiewicz/Desktop/20221019-ngn-model-output/bayes_out_3spp.RDS")

grabthese <- grep("int_samp_small", names(bayes_out$Bayes_modelfit@sim$samples[[1]]))
postSamples <- bayes_out$Bayes_modelfit@sim$samples[[1]][grabthese]

a <- as.data.frame(postSamples) 
mynames <- names(a)
b <- t(a) %>% 
  as_tibble() %>% 
  rownames_to_column("samp") %>% 
  mutate(samp = mynames) %>% 
  mutate(samp = str_remove_all(samp,  "int_samp_small\\.")) %>% 
  separate(samp, into = c("bottle", "species", "x"), "\\.") %>% 
  dplyr::select(-x)

b <- b %>% 
  group_by(bottle, species) %>% 
  nest()

b$bottle <- rep(row.names(bayes_out$Bayes_estimates),
                times = ncol(bayes_out$Bayes_estimates))
b$species <- rep(colnames(bayes_out$Bayes_estimates),
                 each = nrow(bayes_out$Bayes_estimates))


## Pull posteriors for qPCR data for cutthroat 

qMod_out <- readRDS("/Users/elizabethandruszkiewicz/Desktop/20221019-ngn-model-output/cut_qMod_out.RDS")
# qMod_out <- readRDS(here("Output","qpcr","cut_qMod_out.RDS"))

cut_names <- readRDS(here("Output","qpcr","cut_modeled_conc.RDS")) %>%
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

qa <- as.data.frame(qpostSamples) 
qb <- t(qa) %>% 
  as_tibble() %>% 
  add_column(cut_names) %>% 
  group_by(Sample, dilution, Adj_Vol) %>% 
  nest() %>% 
  rename(cutqpcr = data) %>% 
  rename(bottle = Sample) 

mergedb <- b %>% 
  left_join(qb, by="bottle") %>% 
  filter(species == "Oncorhynchus clarkii") %>% 
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
  filter(bottle == "0321_1Prt_2_1")

mergedc <- mergedc %>% 
  ungroup() %>% 
  mutate(meandnaconc = meanpropreads*meantotdna)

#check; do concentration values make sense?
mergedc %>% 
  dplyr::select(bottle, species, meanpropreads, meancutdna, meantotdna, meandnaconc) %>% 
  filter(bottle == "0321_1Prt_2_1")


quants_to_plot <- mergedc %>% 
  filter(creek != "4Pad5") %>% 
  mutate(station = case_when(station == 1 ~ "Down",
                             station == 2 ~ "Up")) %>%
  separate(time, into = c("month","year"), sep = 2) %>% 
  unite(newtime, c(year,month), sep="-") %>% 
  mutate(creek = case_when(creek == "1Prt" ~ "Portage",
                           creek == "2Brn" ~ "Barnes",
                           creek == "3Chk" ~ "Chuckanut",
                           creek == "4Pad11" ~ "Padden",
                           creek == "5Sqm" ~ "Squalicum",
                           TRUE ~ creek))  

# ggplot(quants_to_plot, aes(x=newtime, y=log10(meandnaconc), color=species)) +
#   geom_point() +
#   facet_grid(~creek ~station ~species) +
#   # scale_color_manual(values = pal_okabe_ito) +
#   guides(color= 'none') + 
#   labs(x="Date (YY-MM)", y= "Log10 copies/L water") + 
#   theme_bw() + 
#   scale_x_discrete(guide = guide_axis(angle = -45))
# 
# # ggsave(here("Output","Figures","quant_ts_after_qm_EArun20221020.png"), units="in", width=12, height=8)
# 
# ggplot(quants_to_plot, aes(x=newtime, y=log10(meantotdna))) +
#   geom_point() +
#   facet_grid(~creek ~station) +
#   labs(x="Date (YY-MM)", y= "Log10 copies total salmonid DNA/L water") + 
#   theme_bw() + 
#   scale_x_discrete(guide = guide_axis(angle = -45))
# 
# # ggsave(here("Output","Figures","totaldna_ts_after_qm.png"))


simple <- quants_to_plot %>% 
  select(c(bottle, newtime, creek, station, bio, species, meanpropreads, meancutdna, meantotdna, meandnaconc))

# write_rds(simple, here("Output","salmonids_abs_abundance_notflowcorrected.RDS"))


## filter out <1% reads and BLOQ 
more1percent <- simple %>% 
  filter(meanpropreads > 0.01)

# write_rds(simple, here("Output","salmonids_abs_abundance_notflowcorrected_greater1percent.RDS"))




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
