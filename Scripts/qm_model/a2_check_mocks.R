## NGN Check mock communities
# Author: Eily Allan - from Ryan Kelly and Ole Shelton 
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
input.metadata <- read.csv(here("Input","metabarcoding", "mock_community","startMC.csv"))  


####################################################################
# Keep and format mock community data
####################################################################

# keep only mock community data 
mock <- taxa_table %>% 
  filter(str_detect(Sample_name, "MC")) %>% 
  filter(str_detect(Sample_name, "Tissue")) %>%
  filter(str_detect(species, "Oncorhynchus")) %>%    #### IF ONLY SALMONIDS
  group_by(Sample_name) %>% 
  mutate(ReadDepth = sum(totalReads)) %>% 
  mutate(propReads = totalReads/ReadDepth) %>% 
  separate(Sample_name, c("Community","Proportion","Type","GB","Tech_Rep"), remove=FALSE) %>%
  mutate(Tech_Rep = as.numeric(Tech_Rep)) %>% 
  mutate(Tech_Rep = case_when(GB == "GB" & Tech_Rep == 1 ~ 4, 
                              GB == "GB" & Tech_Rep == 2 ~ 5, 
                              GB == "GB" & Tech_Rep == 3 ~ 6,
                              TRUE ~ Tech_Rep))  %>% 
  unite(SimpleName,c(Community, Proportion, Tech_Rep), remove=FALSE) %>% 
  mutate(Tech_Rep = as.integer(Tech_Rep)) %>% 
  dplyr::rename(Species = species)

ggplot(mock, aes(x=Tech_Rep, y=propReads, fill=Species)) +
  geom_col() + 
  guides(fill = "none") +
  facet_grid(~Community ~ Proportion, scales="free") +
  labs(x="Technical Replicate", y="Proportion of Reads")

start.metadata <- expand_grid(input.metadata,Tech_Rep = 1:3)

start.metadata <- start.metadata %>% 
  filter(! str_detect(Community, "Amplicon")) %>% 
  filter(str_detect(Species, "Oncorhynchus")) %>%    #### IF ONLY SALMONIDS
  rename(Sample_name = Community) %>% 
  group_by(Sample_name, Tech_Rep) %>% 
  mutate(total_input = sum(Start)) %>% 
  mutate(start_prop = Start/total_input) %>% 
  separate(Sample_name, c("Community","Proportion","Type","GB"), remove=FALSE) %>% 
  unite(Sample_name, c("Community","Proportion","Type","GB","Tech_Rep"), remove=FALSE) %>% 
  mutate(Tech_Rep = as.numeric(Tech_Rep)) %>% 
  mutate(Tech_Rep = case_when(GB == "GB" & Tech_Rep == 1 ~ 4, 
                              GB == "GB" & Tech_Rep == 2 ~ 5, 
                              GB == "GB" & Tech_Rep == 3 ~ 6,
                              TRUE ~ Tech_Rep)) %>% 
  unite(SimpleName,c(Community, Proportion, Tech_Rep), remove=FALSE) %>% 
  mutate(., Species = case_when(Species == "Oncorhynchus tschawtscha" ~ "Oncorhynchus tshawytscha",
                                TRUE ~ Species))


ggplot(start.metadata, aes(x=SimpleName, y=start_prop, fill=Species)) +
  geom_col() + 
  guides(fill = "none") +
  facet_wrap(~Community ~ Type, scales="free", ncol=2) + 
  ggtitle("Starting Proportions")

meta.merge <- start.metadata %>% 
  ungroup() %>% 
  dplyr::select(c(SimpleName, Species, start_prop)) %>% 
  unite(samplespec, c("SimpleName","Species"))

miseq.merge <- mock %>% 
  ungroup() %>% 
  dplyr::select(-Sample_name) %>% 
  dplyr::select(c(SimpleName, Species, totalReads, ReadDepth, propReads)) %>% 
  unite(samplespec, c("SimpleName","Species")) 

combined.data <- miseq.merge %>%
  full_join(meta.merge, by= "samplespec") %>% 
  mutate(propReads = if_else(is.na(propReads), 0, propReads)) %>% 
  mutate(start_prop = if_else(is.na(start_prop), 0, start_prop)) %>% 
  mutate(N_pcr_mock = 43) %>%  
  separate(samplespec, c("Community","CommType","Tech_Rep","Genus","Species")) %>%
  dplyr::rename(nReads = totalReads) %>% 
  dplyr::rename(totReads = ReadDepth) %>% 
  unite("species", c(Genus, Species), sep = " ") %>% 
  filter(start_prop !=0)

# write output to file 
# write_rds(combined.data, file=here("Output","metabarcoding","20221018_mockdatatocalibrate.RDS"))

# write output to file - salmonids only 
# write_rds(combined.data, file=here("Output","metabarcoding","20221018_mockdatatocalibrate_salmonidonly.RDS"))

####################################################################
# Plot to see what they look like
####################################################################

plot(combined.data$start_prop, combined.data$propReads)

ggplot(combined.data, aes(x = Tech_Rep, y = propReads, fill = species)) +
  geom_col() +  
  guides(fill = "none") +
  facet_grid(~Community ~CommType, scales="free_x") +
  theme_bw() +
  labs(x="Technical Replicate", y="Proportion of Reads") 

plotmocks_obs <- combined.data %>% 
  dplyr::select(-c(nReads,totReads,N_pcr_mock, start_prop)) %>% 
  dplyr::rename(prop=propReads) %>% 
  mutate(typedata = "Observed") 

plotmocks_exp <- combined.data %>% 
  dplyr::select(-c(nReads,totReads,N_pcr_mock, Tech_Rep, propReads)) %>% 
  #filter(Tech_Rep==1) %>% 
  distinct() %>% 
  dplyr::rename(prop=start_prop) %>% 
  mutate(typedata = "Expected") %>% 
  mutate(Tech_Rep = 0)

plotmocks <- rbind(plotmocks_obs, plotmocks_exp)

####################################################################
# Prep to use even to calibrate skew, vice versa to check 
####################################################################

# ## If only salmonids 
# combined.data <- combined.data %>%
#   filter(str_detect(species, "Oncorhynchus"))

evenmocks <- combined.data %>% 
  filter(CommType == "Even") %>% 
  replace(is.na(.), 0) %>% 
  filter(!is.na(totReads))

skewmocks <- combined.data %>% 
  filter(CommType == "Skew") %>% 
  replace(is.na(.), 0)

evenmocktrue <- evenmocks 
evenmockunk <- evenmocks %>% 
  rename(time = Community) %>% 
  rename(station = CommType) %>% 
  rename(tech = Tech_Rep) %>% 
  mutate(creek=1, biol=1) %>% 
  filter(! is.na(species)) %>% 
  rename(Nreads = nReads) %>% 
  dplyr::select(c(time, creek, station, biol, tech, species, Nreads))

skewmocktrue <- skewmocks 
skewmockunk <- skewmocks %>%
  rename(time = Community) %>% 
  rename(station = CommType) %>% 
  rename(tech = Tech_Rep) %>% 
  mutate(creek=1, biol=1) %>% 
  filter(! is.na(species)) %>% 
  rename(Nreads = nReads) %>% 
  dplyr::select(c(time, creek, station, biol, tech, species, Nreads))

# Prepare for stan model 
eventrue <- format_metabarcoding_data(skewmockunk, evenmocktrue)
eventrue_stan <- makeDesign(eventrue, N_pcr_cycles = 43)

skewtrue <- format_metabarcoding_data(evenmockunk, skewmocktrue)
skewtrue_stan <- makeDesign(skewtrue, N_pcr_cycles = 43)

####################################################################
# Run ML to use even to calibrate skew, vice versa to check 
####################################################################

ML_out_eventrue <- QM_likelihood(here("Scripts", "functions", "quant_metabar_rosetta_noSampleEta.stan"), eventrue_stan)
ML_alpha_eventrue <- ML_out_eventrue$ML_alpha_est %>% 
  mutate(calibration="Even")

ML_out_skewtrue <- QM_likelihood(here("Scripts", "functions", "quant_metabar_rosetta_noSampleEta.stan"), skewtrue_stan)
ML_alpha_skewtrue <- ML_out_skewtrue$ML_alpha_est %>% 
  mutate(calibration="Skew")

#compare_alphas <- rbind(ML_alpha_eventrue, ML_alpha_skewtrue) 

mockmodeleven <- ML_out_skewtrue$ML_estimates %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") %>%
  separate(col = sample, into = c("time", "creek", "station", "biol"), remove = FALSE) %>% 
  filter(value > 0.001) %>% 
  dplyr::rename(Community = time) %>% 
  dplyr::select(-c(sample, creek, station, biol)) %>% 
  mutate(Tech_Rep = -1) %>% 
  mutate(typedata = "Modeled") %>% 
  mutate(CommType = "Even") %>% 
  dplyr::rename(prop = value) 

mockmodelskew <- ML_out_eventrue$ML_estimates %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") %>%
  separate(col = sample, into = c("time", "creek", "station", "biol"), remove = FALSE) %>% 
  filter(value > 0.001) %>% 
  dplyr::rename(Community = time) %>% 
  dplyr::select(-c(sample, creek, station, biol)) %>% 
  mutate(Tech_Rep = -1) %>% 
  mutate(typedata = "Modeled") %>% 
  mutate(CommType = "Skew") %>% 
  dplyr::rename(prop = value)

plotmocks <- rbind(plotmocks, mockmodeleven, mockmodelskew)

####################################################################
# Plot input, observed, modeled for each community  
####################################################################

# Set color palette so don't change when plot
o_i_colors <- c(rgb(230, 159,   0, maxColorValue = 255),  # orange
                rgb( 86, 180, 233, maxColorValue = 255),  # skyblue
                rgb(  0, 158, 115, maxColorValue = 255),  # green
                rgb(240, 228,  66, maxColorValue = 255),  # yellow
                rgb(  0, 114, 178, maxColorValue = 255),  # blue
                rgb(204, 121, 167, maxColorValue = 255)   # purple
)
#o_i_colors = scales::hue_pal()(length(unique(mock$species)))
pal_okabe_ito <- newpal(col = o_i_colors,
                        names = unique(plotmocks$species))

mc1_output_even <- plotmocks %>% 
  filter(Community == "MC1") %>% 
  filter(CommType == "Even") %>% 
  ggplot(aes(x = Tech_Rep, y = prop, fill = species)) +
  geom_col() +  
  guides(fill = "none") +
  theme_bw() +
  scale_fill_manual(values = pal_okabe_ito) +
  labs(x=" ", y="Proportion", title = "MC1 Even") +
  scale_x_discrete(breaks=c("-1","0","1","2","3","4","5","6"), labels=c("Modeled", "Expected", "Observed", "Observed", "Observed", "Observed", "Observed", "Observed"), guide = guide_axis(angle = -45))

mc1_output_skew <- plotmocks %>% 
  filter(Community == "MC1") %>% 
  filter(CommType == "Skew") %>% 
  ggplot(aes(x = Tech_Rep, y = prop, fill = species)) +
  geom_col() +  
  guides(fill = "none") +
  theme_bw() +
  scale_fill_manual(values = pal_okabe_ito) +
  labs(x=" ", y=" ", title = "MC1 Skew")  +
  scale_x_discrete(breaks=c("-1","0","1","2","3"), labels=c("Modeled", "Expected", "Observed", "Observed", "Observed"), guide = guide_axis(angle = -45))

mc2_output_even <- plotmocks %>% 
  filter(Community == "MC2") %>% 
  filter(CommType == "Even") %>% 
  ggplot(aes(x = Tech_Rep, y = prop, fill = species)) +
  geom_col() +  
  guides(fill = "none") +
  theme_bw() +
  scale_fill_manual(values = pal_okabe_ito) +
  labs(x=" ", y="Proportion", title = "MC2 Even") +
  scale_x_discrete(breaks=c("-1","0","1","2","3","4","5","6"), labels=c("Modeled", "Expected", "Observed", "Observed", "Observed", "Observed", "Observed", "Observed"), guide = guide_axis(angle = -45))

mc2_output_skew <- plotmocks %>% 
  filter(Community == "MC2") %>% 
  filter(CommType == "Skew") %>% 
  ggplot(aes(x = Tech_Rep, y = prop, fill = species)) +
  geom_col() +  
  guides(fill = "none") +
  theme_bw() +
  scale_fill_manual(values = pal_okabe_ito) +
  labs(x=" ", y=" ", title = "MC2 Skew")  +
  scale_x_discrete(breaks=c("-1","0","1","2","3"), labels=c("Modeled", "Expected", "Observed", "Observed", "Observed"), guide = guide_axis(angle = -45))

mc3_output_even <- plotmocks %>% 
  filter(Community == "MC3") %>% 
  filter(CommType == "Even") %>% 
  ggplot(aes(x = Tech_Rep, y = prop, fill = species)) +
  geom_col() +  
  guides(fill = "none") +
  theme_bw() +
  scale_fill_manual(values = pal_okabe_ito) +
  labs(x=" ", y="Proportion", title = "MC3 Even") +
  scale_x_discrete(breaks=c("-1","0","1","2","3","4","5","6"), labels=c("Modeled", "Expected", "Observed", "Observed", "Observed", "Observed", "Observed", "Observed"), guide = guide_axis(angle = -45))

mc3_output_skew <- plotmocks %>% 
  filter(Community == "MC3") %>% 
  filter(CommType == "Skew") %>% 
  ggplot(aes(x = Tech_Rep, y = prop, fill = species)) +
  geom_col() +  
  guides(fill = "none") +
  theme_bw() +
  scale_fill_manual(values = pal_okabe_ito) +
  labs(x=" ", y=" ", title = "MC3 Skew") +
  scale_x_discrete(breaks=c("-1","0","1","2","3"), labels=c("Modeled", "Expected", "Observed", "Observed", "Observed"), guide = guide_axis(angle = -45))

#save
grid.arrange(mc1_output_even, mc1_output_skew, mc2_output_even, mc2_output_skew, mc3_output_even, mc3_output_skew,  nrow=3, ncol=2)
g <- arrangeGrob(mc1_output_even, mc1_output_skew, mc2_output_even, mc2_output_skew, mc3_output_even, mc3_output_skew,  nrow=3, ncol=2)
ggsave(file=here("Output","Figures","mock_internal_calibration.png"),g)


####################################################################
# Prep to use different communities as true (both even and skew)
####################################################################

MC1 <- combined.data %>% 
  filter(Community == "MC1") %>% 
  replace(is.na(.), 0) %>% 
  filter(!is.na(totReads))

MC2 <- combined.data %>% 
  filter(Community == "MC2") %>% 
  replace(is.na(.), 0) %>% 
  filter(!is.na(totReads))

MC3 <- combined.data %>% 
  filter(Community == "MC3") %>% 
  replace(is.na(.), 0) %>% 
  filter(!is.na(totReads))

MC1true <- MC1
MC2true <- MC2
MC3true <- MC3

MC23unk <- rbind(MC2, MC3)
MC23unk <- MC23unk %>% 
  rename(time = Community) %>% 
  rename(station = CommType) %>% 
  rename(tech = Tech_Rep) %>% 
  mutate(biol=1) %>% 
  filter(! is.na(species)) %>% 
  rename(Nreads = nReads) %>% 
  mutate(creek = case_when(station == "Even" ~ 1,
                           station == "Skew" ~ 2)) %>% 
  dplyr::select(c(time, creek, station, biol, tech, species, Nreads))

MC13unk <- rbind(MC1, MC3)
MC13unk <- MC13unk %>% 
  rename(time = Community) %>% 
  rename(station = CommType) %>% 
  rename(tech = Tech_Rep) %>% 
  mutate(biol=1) %>% 
  filter(! is.na(species)) %>% 
  rename(Nreads = nReads) %>% 
  mutate(creek = case_when(station == "Even" ~ 1,
                           station == "Skew" ~ 2)) %>% 
  dplyr::select(c(time, creek, station, biol, tech, species, Nreads))


MC12unk <- rbind(MC1, MC2)
MC12unk <- MC12unk %>% 
  rename(time = Community) %>% 
  rename(station = CommType) %>% 
  rename(tech = Tech_Rep) %>% 
  mutate(biol=1) %>% 
  filter(! is.na(species)) %>% 
  rename(Nreads = nReads) %>% 
  mutate(creek = case_when(station == "Even" ~ 1,
                           station == "Skew" ~ 2)) %>% 
  dplyr::select(c(time, creek, station, biol, tech, species, Nreads))

# Prepare for stan model 
MC1true<- format_metabarcoding_data(MC23unk, MC1true)
MC1true_stan <- makeDesign(MC1true, N_pcr_cycles = 43)

MC2true<- format_metabarcoding_data(MC13unk, MC2true)
MC2true_stan <- makeDesign(MC2true, N_pcr_cycles = 43)

MC3true<- format_metabarcoding_data(MC12unk, MC3true)
MC3true_stan <- makeDesign(MC3true, N_pcr_cycles = 43)

####################################################################
# Run ML to use 1 to calibrate 2/3 and every combo
####################################################################

ML_out_MC1true <- QM_likelihood(here("Scripts", "functions", "quant_metabar_rosetta_noSampleEta.stan"), MC1true_stan)
ML_alpha_MC1true <- ML_out_MC1true$ML_alpha_est %>% 
  mutate(calibration="MC1")

ML_out_MC2true <- QM_likelihood(here("Scripts", "functions", "quant_metabar_rosetta_noSampleEta.stan"), MC2true_stan)
ML_alpha_MC2true <- ML_out_MC2true$ML_alpha_est %>% 
  mutate(calibration="MC2")

ML_out_MC3true <- QM_likelihood(here("Scripts", "functions", "quant_metabar_rosetta_noSampleEta.stan"), MC3true_stan)
ML_alpha_MC3true <- ML_out_MC3true$ML_alpha_est %>% 
  mutate(calibration="MC3")


####################################################################
# Plot comparison of alphas depending on what is "true"   
####################################################################

compare_alphas <- rbind(ML_alpha_eventrue, ML_alpha_skewtrue, ML_alpha_MC1true,  ML_alpha_MC2true,  ML_alpha_MC3true) 

ggplot(compare_alphas, aes(y=species, x=alpha_est, col=calibration)) + 
  geom_point() +
  labs(title = "Tissue-Derived Mock Communities") +
  theme_bw() + 
  labs(x="Estimated Alpha", y="Species", color="Community used as true")

ggsave(file=here("Output","Figures","mock_internal_calibration_compare_alphas.png"))


####################################################################
# Plot 1 to calibrate 2/3 etc   
####################################################################


mockmodel23 <- ML_out_MC1true$ML_estimates %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") %>%
  separate(col = sample, into = c("time", "creek", "station", "biol"), remove = FALSE) %>% 
  filter(value > 0.001) %>% 
  dplyr::rename(Community = time, CommType = creek) %>% 
  dplyr::select(-c(sample, station, biol)) %>% 
  #mutate(Tech_Rep = -1) %>% 
  #mutate(typedata = "Modeled") %>% 
  #mutate(CommType = "Even") %>% 
  dplyr::rename(prop = value) 

mockmodelskew <- ML_out_eventrue$ML_estimates %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "species") %>%
  separate(col = sample, into = c("time", "creek", "station", "biol"), remove = FALSE) %>% 
  filter(value > 0.001) %>% 
  dplyr::rename(Community = time) %>% 
  dplyr::select(-c(sample, creek, station, biol)) %>% 
  mutate(Tech_Rep = -1) %>% 
  mutate(typedata = "Modeled") %>% 
  mutate(CommType = "Skew") %>% 
  dplyr::rename(prop = value)

plotmocks <- rbind(plotmocks, mockmodeleven, mockmodelskew)

####################################################################
# Plot input, observed, modeled for each community  
####################################################################

# Set color palette so don't change when plot
o_i_colors <- c(rgb(230, 159,   0, maxColorValue = 255),  # orange
                rgb( 86, 180, 233, maxColorValue = 255),  # skyblue
                rgb(  0, 158, 115, maxColorValue = 255),  # green
                rgb(240, 228,  66, maxColorValue = 255),  # yellow
                rgb(  0, 114, 178, maxColorValue = 255),  # blue
                rgb(204, 121, 167, maxColorValue = 255)   # purple
)
#o_i_colors = scales::hue_pal()(length(unique(mock$species)))
pal_okabe_ito <- newpal(col = o_i_colors,
                        names = unique(plotmocks$species))

mc1_output_even <- plotmocks %>% 
  filter(Community == "MC1") %>% 
  filter(CommType == "Even") %>% 
  ggplot(aes(x = Tech_Rep, y = prop, fill = species)) +
  geom_col() +  
  guides(fill = "none") +
  theme_bw() +
  scale_fill_manual(values = pal_okabe_ito) +
  labs(x=" ", y="Proportion", title = "MC1 Even") +
  scale_x_discrete(breaks=c("-1","0","1","2","3","4","5","6"), labels=c("Modeled", "Expected", "Observed", "Observed", "Observed", "Observed", "Observed", "Observed"), guide = guide_axis(angle = -45))

mc1_output_skew <- plotmocks %>% 
  filter(Community == "MC1") %>% 
  filter(CommType == "Skew") %>% 
  ggplot(aes(x = Tech_Rep, y = prop, fill = species)) +
  geom_col() +  
  guides(fill = "none") +
  theme_bw() +
  scale_fill_manual(values = pal_okabe_ito) +
  labs(x=" ", y=" ", title = "MC1 Skew")  +
  scale_x_discrete(breaks=c("-1","0","1","2","3"), labels=c("Modeled", "Expected", "Observed", "Observed", "Observed"), guide = guide_axis(angle = -45))

mc2_output_even <- plotmocks %>% 
  filter(Community == "MC2") %>% 
  filter(CommType == "Even") %>% 
  ggplot(aes(x = Tech_Rep, y = prop, fill = species)) +
  geom_col() +  
  guides(fill = "none") +
  theme_bw() +
  scale_fill_manual(values = pal_okabe_ito) +
  labs(x=" ", y="Proportion", title = "MC2 Even") +
  scale_x_discrete(breaks=c("-1","0","1","2","3","4","5","6"), labels=c("Modeled", "Expected", "Observed", "Observed", "Observed", "Observed", "Observed", "Observed"), guide = guide_axis(angle = -45))

mc2_output_skew <- plotmocks %>% 
  filter(Community == "MC2") %>% 
  filter(CommType == "Skew") %>% 
  ggplot(aes(x = Tech_Rep, y = prop, fill = species)) +
  geom_col() +  
  guides(fill = "none") +
  theme_bw() +
  scale_fill_manual(values = pal_okabe_ito) +
  labs(x=" ", y=" ", title = "MC2 Skew")  +
  scale_x_discrete(breaks=c("-1","0","1","2","3"), labels=c("Modeled", "Expected", "Observed", "Observed", "Observed"), guide = guide_axis(angle = -45))

mc3_output_even <- plotmocks %>% 
  filter(Community == "MC3") %>% 
  filter(CommType == "Even") %>% 
  ggplot(aes(x = Tech_Rep, y = prop, fill = species)) +
  geom_col() +  
  guides(fill = "none") +
  theme_bw() +
  scale_fill_manual(values = pal_okabe_ito) +
  labs(x=" ", y="Proportion", title = "MC3 Even") +
  scale_x_discrete(breaks=c("-1","0","1","2","3","4","5","6"), labels=c("Modeled", "Expected", "Observed", "Observed", "Observed", "Observed", "Observed", "Observed"), guide = guide_axis(angle = -45))

mc3_output_skew <- plotmocks %>% 
  filter(Community == "MC3") %>% 
  filter(CommType == "Skew") %>% 
  ggplot(aes(x = Tech_Rep, y = prop, fill = species)) +
  geom_col() +  
  guides(fill = "none") +
  theme_bw() +
  scale_fill_manual(values = pal_okabe_ito) +
  labs(x=" ", y=" ", title = "MC3 Skew") +
  scale_x_discrete(breaks=c("-1","0","1","2","3"), labels=c("Modeled", "Expected", "Observed", "Observed", "Observed"), guide = guide_axis(angle = -45))

#save
grid.arrange(mc1_output_even, mc1_output_skew, mc2_output_even, mc2_output_skew, mc3_output_even, mc3_output_skew,  nrow=3, ncol=2)
g <- arrangeGrob(mc1_output_even, mc1_output_skew, mc2_output_even, mc2_output_skew, mc3_output_even, mc3_output_skew,  nrow=3, ncol=2)
ggsave(file=here("Output","Figures","mock_internal_calibration.png"),g)



