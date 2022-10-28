## NGN Classifying ASVs: MiFish taxonomic assignment 
# Author: Eily Allan - modified from Erin D'Agnese and Ramon Gallego 
# Person running: Eily
# Last modified: 10/17/22 by Eily
# Date of run: 10/17/22 by Eily 

# Overview 
# This script takes the combined ASV table from all Miseq runs and assigns taxonomy to each ASV in a multi-step process, as follows: 
# 1) use a tree-based algorithm using the R package "insect" to assign taxonomy (note: the tree is too large to be stored on Github so it is stored locally.) 
# 2) take any ASVs not assigned taxonomy or not recieving a species level assigment by insect and blast
# 3) read in the blast file and keep any species level assignments with 100% percent identity, or >95% identity but unambiguous
# 4) combine insect and blast assignments for species level assignments 
# 5) sum reads from all ASVs that assign to the exact same species 

# Inputs: 
# 1) combined ASV table for environmental samples and mock community samples 
# 2) combined hash key (hash ID, sequence) for environmental samples and mock community samples 
# 3) insect classifier (stored locally) 

# Outputs: 
# 1) taxonomy key with each ASV and species level assigments (no info on what is found in what samples)
# 2) table of samples with taxonomic assignment to ANY taxonomic level 
# 3) table of samples with taxonomic assignment to species level (without ASVs summed for the same species)
# 4) table of samples with summed reads for all species found (and number of reads not assigned to species level)

####################################################################
# Set up
####################################################################

## Load packages
library(tidyverse)
library(insect)
library(seqinr)
library(here)
library(taxonomizr)
library(DECIPHER)

# read in hash key and asv table
hash_key <- here("Input","metabarcoding","combined_MiFish_hash_key.csv")
ASVs <- here("Input","metabarcoding","combined_MiFish_ASV_table.csv")

## HARD CODE BECAUSE MUST BE STORED LOCALLY
classifier <- "/Users/elizabethandruszkiewicz/GoogleDrive/UW/GitHub/NextGenNEPA_LOCAL/Input/classifiers/classifier_12S_v1.rds"

## Read in files 
hash_key <- read_csv(hash_key)
Hash <- hash_key %>% 
  select(Hash, Sequence) %>% 
  distinct()
ALL.ASVs <- read_csv(ASVs)
tree <- read_rds(classifier)

####################################################################
# Classify
## Note that the classifier is maintained by insect. We are using version 1, which was last updated on 20181111. 
## Check here to see if there is a newer version: https://cran.r-project.org/web/packages/insect/vignettes/insect-vignette.html
####################################################################

## Prepare and run insect

## FOR 12S REMOVE ANYTHING TOO BIG (aka bacteria)
hash.length <- nchar(Hash$Sequence)
hash.keep <- hash.length < 200
Hash <- Hash[hash.keep,]
# going from 24577 hashes of any length to 1977 hashes of length <200 
# COME BACK AND CHECK NUMBER OF READS 

# format for insect
hashes.insect <- char2dna(Hash$Sequence)
names(hashes.insect) <- Hash$Hash

# just check it out 
hashes.insect

# do the classification! this takes a while
clasif.hashes <- classify (x = hashes.insect, tree = tree, cores = 4)

# rename columns to be useful
names(clasif.hashes) <- c('representative', 'taxID', 'taxon', 'rank', 'score', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')

####################################################################
# Save only species level assignments 
####################################################################

# see table of ranks 
ranks <- clasif.hashes %>% dplyr::count (rank) %>% arrange(desc(n))
# so of the 1977 hashes, only 90 are to species level, but this might be fine if its the majority of reads

# keep only ASVs assigned to species level 
speciesASVs <- clasif.hashes %>%
  filter(rank == "species")

# keep other ASVs to blast 
toblastASVs <- clasif.hashes %>%
  filter(rank != "species")

toblastseqs <- hash_key %>% 
  filter(Hash %in% toblastASVs$representative)

toblast <- DNAStringSet(toblastseqs$Sequence)
names(toblast) <- toblastseqs$Hash

writeXStringSet(toblast, file=here("Output","metabarcoding","ASVs_to_blast.fasta"))


####################################################################
# Read in BLAST results for non-species level IDs 
####################################################################

blast_results <- read_delim((here("Input","metabarcoding","ASVs_to_blast.txt")), 
                            col_names = c("qseqid", "sseqid", "sacc", "pident", "length", "mismatch", "gapopen", "qcovus", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "staxid", "qlen", "sscinames", "sseq"), 
                            delim = "\t" )

# only keep things that hit at 100% id to ONLY ONE species
blast_100only <- blast_results %>% 
  filter(pident == 100) %>% 
  select(c(qseqid, sscinames)) %>% 
  distinct() %>% 
  mutate(., sscinames = case_when(sscinames == "Dinornis robustus" ~ "Dinornis giganteus",
                                  sscinames == "Macropus cf. giganteus TL-2021" ~ "Macropus giganteus",
                                  sscinames == "Macropus fuliginosus" ~ "Macropus giganteus",
                                  sscinames == "Oncorhynchus clarkii lewisi" ~ "Oncorhynchus clarkii", 
                                  TRUE ~ sscinames)) %>% 
  distinct() %>% 
  group_by(qseqid) %>%
  filter(n() == 1)

# now look at things that were >95% ID but only hit to ONE thing and see if we can squeeze a few more annotations out
blast_unambiguous <- blast_results %>%
  filter(! qseqid %in% blast_100only$qseqid) %>% 
  select(c(qseqid, pident, sscinames)) %>% 
  distinct() %>% 
  group_by(qseqid) %>%
  filter(n() == 1) %>% 
  filter(pident > 95) %>% 
  select(c(qseqid, sscinames)) %>% 
  distinct() 

blast_both <- rbind(blast_100only, blast_unambiguous)

# didn't save hash ids before blasting so need to make a key to add them back in
toblastseqs <- toblastseqs %>% 
  mutate(id = row_number()) %>% 
  mutate(qseqid = paste0("Query_",id)) %>% 
  left_join(blast_both, by= "qseqid")

toaddblast <- toblastseqs %>% 
  filter(! is.na(sscinames)) %>% 
  select(c(Hash, sscinames)) %>% 
  dplyr::rename(species = sscinames)

####################################################################
# Combine insect and BLAST species level IDs and write taxonomy key
####################################################################

toaddinsect <- speciesASVs %>% 
  select(c(representative, species)) %>% 
  dplyr::rename(Hash = representative)

final_taxonomy_key <- rbind(toaddinsect, toaddblast)  

#
asvtotaxa <- ALL.ASVs %>% 
  left_join(final_taxonomy_key, by="Hash") %>% 
  group_by(Sample_name) %>% 
  mutate(ReadDepth = sum(nReads)) 

# add a few more by hand with manual blast that didn't get picked up somehow
checkmanually <- asvtotaxa %>% 
  ungroup() %>% 
  select(c(Hash,nReads,species)) %>% 
  group_by(Hash) %>% 
  mutate(totalReads = sum(nReads)) %>% 
  select(-nReads) %>% 
  filter(is.na(species)) %>% 
  distinct()

final_taxonomy_key <- final_taxonomy_key %>% 
  add_row(Hash="580a64b0134da56e49bf06482ed6ee164d8a6ec5", species="Dicamptodon tenebrosus") %>% 
  #add_row(Hash="49b916c007d3dbb75da7ceb655dd8cefb4e15ac4", species="Cottus") %>% 
  #add_row(Hash="23bf5f5d24f589cba120a8635ee13e9817d90a52", species="Oncorhynchus mykiss or Oncorhynchus clarkii") %>% 
  #add_row(Hash="be1bd85cbedd796b3b70618ed8248f5c66768c49", species="Oncorhynchus kisutch or Oncorhynchus keta but probably kisutch") %>% 
  add_row(Hash="7ed525639c40b9418b4743e59f833fa5b3bdf355", species="Pimephales promelas")


####################################################################
# Check how many annotated (ASVs and reads)
####################################################################

annotatedASVs <- asvtotaxa %>% 
  filter(!is.na(species)) %>% 
  group_by(Sample_name) %>% 
  mutate(AnnotatedReads = sum(nReads)) %>% 
  mutate(AnnotatedASVs = length(unique(Hash))) %>% 
  mutate(PercAnnoReads = AnnotatedReads/ReadDepth*100)%>% 
  mutate(PercAnnoASVs = AnnotatedASVs/TotalASVs*100)

annotation_stats <- annotatedASVs %>% 
  dplyr::select(c(Sample_name, ReadDepth, TotalASVs,AnnotatedReads,AnnotatedASVs,PercAnnoReads,PercAnnoASVs)) %>% 
  distinct()

meanpercentreadsannotated <- mean(annotation_stats$PercAnnoReads)

####################################################################
# Separate out mock community samples to see what is found
####################################################################

mocks <- asvtotaxa %>% 
  filter(str_detect(Sample_name, "MC"))

## MC1 
# 70e8e188ebdd7ea3d78d9c745f097cb48e934f95 - either Entosphenus tridentatus or Lampetra ayresii
# a2f2121e42c56aa5bd2b1fca2bd1b29869d2400a - either Entosphenus tridentatus or Lampetra ayresii
# 8b14c93f1397b9332dad055f0a0fc595f3646681 - either Cottus marginatus or Cottus asper
# 49b916c007d3dbb75da7ceb655dd8cefb4e15ac4 - either Cottus marginatus or Cottus asper
# 0a3f5fde1bbb570e4705a194d016c59a57f82481 - either Cottus marginatus or Cottus asper
# 1c6d757a385993942374fd4c94f562e4dda52ae5 - either Cottus marginatus or Cottus asper
# 13febd03296481e037b148ca82f8a9a59c7c61eb - either Cottus marginatus or Cottus asper
# 0caa98cce420d7464ced1558f1635a8ae531104a - either Cottus marginatus or Cottus asper
# 631402d8ed10c8b197e1d67c208b29b9ba254a64 -- HERE IS Salvelinus malma


## MC2 
# 49b916c007d3dbb75da7ceb655dd8cefb4e15ac4 -- HERE IS Cottus asper
# 70e8e188ebdd7ea3d78d9c745f097cb48e934f95 -- HERE IS Entosphenus tridentatus
# 2888591f23daf76ef3ae53a7d4bd6ffe13361bb8 -- Muskrat - Ondatra zibethicus 
# 2c79a691f61d8877b4c2f3953447d6e033680d33 -- says Salvelinus alpinus but CHANGE TO S. confluentus?


## MC 3
# 8b14c93f1397b9332dad055f0a0fc595f3646681 -- HERE IS Cottus marginatus 
# 0a3f5fde1bbb570e4705a194d016c59a57f82481 -- HERE IS Cottus marginatus
# a2f2121e42c56aa5bd2b1fca2bd1b29869d2400a -- HERE IS Lampetra ayresii  
# 631402d8ed10c8b197e1d67c208b29b9ba254a64 -- either Salvelinus malma or Salvelinus confluentus 
# 38e693e5c3a9002b3c7fe03ba0f5d3e6440d86ea -- Novumbra hubbsi
# bc392f2ab2d28f1c9aa2a98ebf78c0d707871993 -- Bufo exo CHANGE TO anaxyrus boreas?
# 1c6d757a385993942374fd4c94f562e4dda52ae5 -- HERE IS Cottus marginatus
# 2eb18eaf3053e4bb306b47b42bd9c34863aef20b -- Lithobates catesbeianus CHANGE TO Rana catesbeiana?
# 6f93809b6bc63b61b0d350733addd08a6ab86f2e -- Ardea cinerea CHANGE TO Ardea herodias? 
# 2c79a691f61d8877b4c2f3953447d6e033680d33 -- says Salvelinus alpinus but put in S. malma AND S. confluentus 
# cce45cc7e4ea0594ae1859ee98f0420718bcd348 -- Bos taurus 


## ROUND UP 
# 631402d8ed10c8b197e1d67c208b29b9ba254a64 -- HERE IS Salvelinus malma
# 49b916c007d3dbb75da7ceb655dd8cefb4e15ac4 -- HERE IS Cottus asper
# 70e8e188ebdd7ea3d78d9c745f097cb48e934f95 -- HERE IS Entosphenus tridentatus
# 2888591f23daf76ef3ae53a7d4bd6ffe13361bb8 -- Muskrat - Ondatra zibethicus 
# 2c79a691f61d8877b4c2f3953447d6e033680d33 -- says Salvelinus alpinus but CHANGE TO S. confluentus?
# 8b14c93f1397b9332dad055f0a0fc595f3646681 -- HERE IS Cottus marginatus 
# 0a3f5fde1bbb570e4705a194d016c59a57f82481 -- HERE IS Cottus marginatus
# a2f2121e42c56aa5bd2b1fca2bd1b29869d2400a -- HERE IS Lampetra ayresii 
# 38e693e5c3a9002b3c7fe03ba0f5d3e6440d86ea -- Novumbra hubbsi
# bc392f2ab2d28f1c9aa2a98ebf78c0d707871993 -- Bufo exo CHANGE TO anaxyrus boreas?
# 1c6d757a385993942374fd4c94f562e4dda52ae5 -- HERE IS Cottus marginatus
# 2eb18eaf3053e4bb306b47b42bd9c34863aef20b -- Lithobates catesbeianus CHANGE TO Rana catesbeiana?
# 6f93809b6bc63b61b0d350733addd08a6ab86f2e -- Ardea cinerea CHANGE TO Ardea herodias? 
# cce45cc7e4ea0594ae1859ee98f0420718bcd348 -- Bos taurus 



final_taxonomy_key <- final_taxonomy_key %>% 
  add_row(Hash="631402d8ed10c8b197e1d67c208b29b9ba254a64", species="Salvelinus malma") %>% 
  filter(Hash != "49b916c007d3dbb75da7ceb655dd8cefb4e15ac4") %>% 
  add_row(Hash="49b916c007d3dbb75da7ceb655dd8cefb4e15ac4", species="Cottus asper") %>% 
  filter(Hash != "70e8e188ebdd7ea3d78d9c745f097cb48e934f95") %>% 
  add_row(Hash="70e8e188ebdd7ea3d78d9c745f097cb48e934f95", species="Entosphenus tridentatus") %>% 
  add_row(Hash="2888591f23daf76ef3ae53a7d4bd6ffe13361bb8", species="Ondatra zibethicus") %>% 
  filter(Hash != "2c79a691f61d8877b4c2f3953447d6e033680d33") %>% 
  add_row(Hash="2c79a691f61d8877b4c2f3953447d6e033680d33", species="Salvelinus confluentus") %>% 
  add_row(Hash="8b14c93f1397b9332dad055f0a0fc595f3646681", species="Cottus marginatus") %>% 
  add_row(Hash="0a3f5fde1bbb570e4705a194d016c59a57f82481", species="Cottus marginatus") %>% 
  filter(Hash != "a2f2121e42c56aa5bd2b1fca2bd1b29869d2400a") %>% 
  add_row(Hash="a2f2121e42c56aa5bd2b1fca2bd1b29869d2400a", species="Lampetra ayresii") %>% 
  add_row(Hash="38e693e5c3a9002b3c7fe03ba0f5d3e6440d86ea", species="Novumbra hubbsi") %>% 
  filter(Hash != "bc392f2ab2d28f1c9aa2a98ebf78c0d707871993") %>% 
  add_row(Hash="bc392f2ab2d28f1c9aa2a98ebf78c0d707871993", species="Anaxyrus boreas") %>% 
  add_row(Hash="1c6d757a385993942374fd4c94f562e4dda52ae5", species="Cottus marginatus") %>% 
  filter(Hash != "2eb18eaf3053e4bb306b47b42bd9c34863aef20b") %>% 
  add_row(Hash="2eb18eaf3053e4bb306b47b42bd9c34863aef20b", species="Rana catesbeiana") %>% 
  filter(Hash != "6f93809b6bc63b61b0d350733addd08a6ab86f2e") %>% 
  add_row(Hash="6f93809b6bc63b61b0d350733addd08a6ab86f2e", species="Ardea herodias") %>% 
  add_row(Hash="cce45cc7e4ea0594ae1859ee98f0420718bcd348", species="Bos taurus") %>% 
  mutate(., species = case_when(species == "Canis lupus familiaris" ~ "Canis lupus",
                                species == "Canis lupus lycaon" ~ "Canis lupus",
                                species == "Gasterosteus aculeatus aculeatus" ~ "Gasterosteus aculeatus",
                                species == "Oncorhynchus clarkii lewisi" ~ "Oncorhynchus clarkii", 
                                  TRUE ~ species))

####################################################################
# Write final taxonomy file 
####################################################################

write.csv(final_taxonomy_key, here("Output","metabarcoding","taxonomy_hash_key.csv"), row.names=FALSE)


####################################################################
# Sum reads to same species for enviro samples and mocks 
####################################################################

finalasvtotaxa <- ALL.ASVs %>% 
  left_join(final_taxonomy_key, by="Hash") %>% 
  group_by(Sample_name) %>% 
  mutate(ReadDepth = sum(nReads)) %>% 
  mutate(TotalASVs = length(unique(Hash)))

annotatedASVs <- finalasvtotaxa %>% 
  filter(!is.na(species)) %>% 
  group_by(Sample_name) %>% 
  mutate(AnnotatedReads = sum(nReads)) %>% 
  mutate(AnnotatedASVs = length(unique(Hash))) %>% 
  mutate(PercAnnoReads = AnnotatedReads/ReadDepth*100)%>% 
  mutate(PercAnnoASVs = AnnotatedASVs/TotalASVs*100)

annotation_stats <- annotatedASVs %>% 
  dplyr::select(c(Sample_name, ReadDepth, TotalASVs,AnnotatedReads,AnnotatedASVs,PercAnnoReads,PercAnnoASVs)) %>% 
  distinct()

meanpercentreadsannotated <- mean(annotation_stats$PercAnnoReads)

taxa_table <- annotatedASVs %>% 
  select(c(Sample_name, nReads, species)) %>% 
  group_by(Sample_name, species) %>% 
  mutate(totalReads = sum(nReads)) %>% 
  select(-nReads) %>% 
  distinct()

write.csv(taxa_table, here("Output","metabarcoding", "taxa_table.csv"),row.names=FALSE)
