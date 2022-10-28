## Prepare all MiSeq data for QM model 
## EAA 
## 8/19/22

## Load packages
library(here)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(DECIPHER)
select <- dplyr::select

## Point to files for MiSeq files 
all.enviro.asv.table <- read_csv(here("Input","metabarcoding", "combined.MiFish.ASV.table.csv"))
all.enviro.hash.key <- read_csv(here("Input","metabarcoding", "MiFish.all.previous.hashes.annotated.csv"))

## Select only salmonids from hash key 
enviro.salmon.hashes <- all.enviro.hash.key %>% 
  filter(family == "Salmonidae")

## Annotate environmental ASV table with only salmonid hashes 
enviro.salmon.asv.table <- all.enviro.asv.table %>% 
  filter(Hash %in% enviro.salmon.hashes$representative) %>% 
  dplyr::rename(representative = Hash) %>% 
  left_join(enviro.salmon.hashes, by = "representative")

## Read in fasta for salmonids from mock community 
salmonids <- readDNAStringSet(here("Input","mockcommunity","salmon_only_seqs.fasta"))
salmon.mc.hashes <- names(salmonids)
rename <- dplyr::rename

## Keep only the ASVs found in MC from environmental samples (unambiguous)
enviro.salmon.asv.MChashesonly <- enviro.salmon.asv.table %>% 
  filter(representative %in% salmon.mc.hashes) 

## Compare number of reads in all environmental salmon ASVs vs. just MC salmonid ASVs
enviro.all.salmon.reads <- sum(enviro.salmon.asv.table$nReads) # 17694962
enviro.salmon.MChashes.reads <- sum(enviro.salmon.asv.MChashesonly$nReads) # 16892531
percent = enviro.salmon.MChashes.reads/enviro.all.salmon.reads*100 # 95% - awesome!! 

## Format this ASV table to be easily used with the QM model
asv.table.for.qm <- enviro.salmon.asv.MChashesonly %>% 
  filter(! str_detect(Sample_name, "Kangaroo")) %>% 
  select(c(Sample_name, nReads, taxon)) %>% 
  group_by(Sample_name, taxon) %>% 
  mutate(totReads = sum(nReads)) %>% 
  select(-nReads) %>% 
  distinct() %>% 
  separate(Sample_name, into=c("marker", "time", "creek", "station", "biol", "tech"), remove=TRUE) %>% 
  select(-marker)

write_rds(asv.table.for.qm, file=paste0(here("Output","metabarcoding"),"/envirodata.salmonids.for.qm.RDS"))

