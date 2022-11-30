

####################################################################
# Set up
####################################################################

## Load packages
library(tidyverse)
library(here)

# read in hash key and asv table
hash_key <- read.csv(here("Input","metabarcoding","combined_MiFish_hash_key.csv"))
ASVs <- read.csv(here("Input","metabarcoding","combined_MiFish_ASV_table.csv"))
taxonomy <- read.csv(here("Output","metabarcoding","taxonomy_hash_key.csv"))
meta_metadata <- read.csv(here("Input","metabarcoding","all_mifish_metadata.csv"))

finalasvtotaxa <- ASVs %>% 
  left_join(taxonomy, by="Hash") %>% 
  group_by(Sample_name) %>% 
  mutate(ReadDepth = sum(nReads)) %>% 
  mutate(TotalASVs = length(unique(Hash)))

envsamples <- finalasvtotaxa %>% 
  filter(! str_detect(Sample_name, "MC")) %>% 
  filter(! str_detect(Sample_name, "Kangaroo")) %>% 
  filter(! str_detect(Sample_name, "Up5"))

envsamplesannotated <- envsamples %>% 
  filter(! is.na(species)) 

mcsamples <- finalasvtotaxa %>% 
  filter(str_detect(Sample_name, "MC")) 

mcsamplesannotated <- mcsamples %>% 
  filter(! is.na(species)) 

persample_annotation <- finalasvtotaxa %>% 
  filter(is.na(species)) %>% 
  group_by(Sample_name) %>% 
  mutate(NAreads = sum(nReads)) %>% 
  mutate(NAasvs = length(unique(Hash))) %>% 
  select(c(Sample_name, ReadDepth, TotalASVs, NAreads, NAasvs)) %>%
  distinct() %>% 
  mutate(Areads = ReadDepth-NAreads) %>% 
  mutate(Aasvs = TotalASVs-NAasvs) %>% 
  mutate(percAreads = Areads/ReadDepth*100) %>% 
  mutate(percAasvs = Aasvs/TotalASVs*100)
   
 
env_annotation <- persample_annotation %>% 
  filter(! str_detect(Sample_name, "MC")) %>% 
  filter(! str_detect(Sample_name, "Kangaroo")) %>% 
  filter(! str_detect(Sample_name, "Up5")) %>% 
  separate(Sample_name, into=c("marker","time","creek","station","bio","tech"), remove=FALSE) %>% 
  mutate(., station = case_when(station == "Up11" ~ "Up",
                                TRUE ~ station)) %>% 
  mutate(., creek = case_when(creek == "1Prt" ~ "Portage Creek",
                              creek == "2Brn" ~ "Barnes Creek",
                              creek == "3Chk" ~ "Chuckanut Creek",
                              creek == "4Pad" ~ "Padden Creek",
                              creek == "5Sqm" ~ "Squalicum Creek",
                              TRUE ~ creek)) 

env_annotation_stats <- env_annotation %>% 
  ungroup() %>% 
  select(c(ReadDepth, TotalASVs, percAreads, percAasvs)) %>% 
  summarize(meanpercAreads = mean(percAreads), meanpercAasvs = mean(percAasvs),
            medianpercAreads = median(percAreads), medianpercAasvs = median(percAasvs),
            minpercreads = min(percAreads), minpercasvs = min(percAasvs),
            maxpercreads = max(percAreads), maxpercasvs = max(percAasvs),
            meanTOTALreads = mean(ReadDepth), meanTOTALasvs = mean(TotalASVs),
            minTOTALreads = min(ReadDepth), minTOTALasvs = min(TotalASVs),
            maxTOTALreads = max(ReadDepth), maxTOTALasvs = max(TotalASVs))
  
env_annotation  %>% 
  ggplot(aes(x=ReadDepth, y=percAreads, color=time)) +
  geom_point() +
  scale_x_log10(limits = c(1E3, 1E6), guide = guide_axis(angle = -45)) +
  facet_grid(~station~creek) +
  labs(x="Total number of reads", y= "Percent of reads annotated") + 
  theme_bw() 

#ggsave(here("Output","SupplementalFigures","percentannotatedbytime.png"), units="in", width=10, height=4)

meta_to_join <- meta_metadata %>% 
  select(-Sample_name) %>% 
  dplyr::rename(Sample_name = Sample_ID) %>% 
  separate(Sample_name, into=c("marker","time","creek","station","bio")) %>% 
  mutate(tech=1) %>% 
  unite(Sample_name, c("marker","time","creek","station","bio", "tech"), sep=".")

  
env_annotation %>% 
  left_join(meta_to_join, by="Sample_name") %>% 
  ggplot(aes(x=ReadDepth, y=percAreads, color=factor(Sequencing.run))) +
  geom_point() +
  scale_x_log10(limits = c(1E3, 1E6), guide = guide_axis(angle = -45)) +
  facet_grid(~station~creek) +
  #scale_fill_manual(values = pal_okabe_ito) +
  #scale_color_manual(values = pal_okabe_ito) +
  #guides(color= 'none') + 
  labs(x="Total number of reads", y= "Percent of reads annotated", color="Sequencing Run") + 
  theme_bw() 

ggsave(here("Output","SupplementalFigures","percentannotatedbyrun.png"), units="in", width=10, height=4)


salmonidsused <- envsamples %>% 
  filter(species %in% c("Oncorhynchus kisutch", "Oncorhynchus nerka", "Oncorhynchus mykiss", "Oncorhynchus clarkii")) 


enviroreads <- sum(envsamples$nReads)
enviroasvs <- length(unique(envsamples$Hash))

enviroannotatedreads <- sum(envsamplesannotated$nReads)
enviroannotatedasvs <- length(unique(envsamplesannotated$Hash))

enviropercreadsannotated <- enviroannotatedreads/enviroreads*100
enviropercasvsannotated<- enviroannotatedasvs/enviroasvs*100

totalsalmonidreads <- sum(salmonidsused$nReads)
percentsalmonidreadstotal <- totalsalmonidreads/enviroreads*100
percentsalmonidreadsannotated <- totalsalmonidreads/enviroannotatedreads*100

mcreads <- sum(mcsamples$nReads)
mcannotatedreads <- sum(mcsamplesannotated$nReads)
mcpercannotated <- mcannotatedreads/mcreads*100


salmonidspersample <- salmonidsused %>% 
  select(c(Sample_name,species, nReads)) %>% 
  group_by(Sample_name,species) %>% 
  summarize(salmonreads = sum(nReads)) %>% 
  distinct()
  
salmonstats <- salmonidspersample %>% 
  group_by(species) %>% 
  summarize(numsamples = n()) %>% 
  mutate(percentsamples = numsamples/363*100)

oclarkiistats <- salmonidspersample %>% 
  group_by(Sample_name) %>% 
  mutate(foursalmonreads=sum(salmonreads)) %>% 
  mutate(percentofeach = salmonreads/foursalmonreads*100) %>% 
  filter(species == "Oncorhynchus clarkii")

overhalfclarkii <- oclarkiistats %>% filter(percentofeach >50)
percentoverhalfclarkii = nrow(overhalfclarkii)/nrow(env_annotation)*100
