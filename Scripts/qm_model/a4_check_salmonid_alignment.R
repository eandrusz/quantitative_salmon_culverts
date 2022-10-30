## NGN plot alignment of salmonids with mifish primer
# Author: Eily Allan 
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

# Load packages
library(here)
library(tidyverse)
library(DECIPHER)

# read in fasta file of the seven salmonid species of interest and the F mifish primer (made in geneious)
fastafile <- here("Input","metabarcoding", "salmonids_to_align","MiFish-Salmonids.fasta")

# Oncorhynchus clarkii henshawi AY886762
# Oncorhynchus gorbuscha EF455489
# Oncorhynchus keta AP010773
# Oncorhynchus kisutch EF126369
# Oncorhynchus mykiss L29771
# Oncorhynchus nerka EF055889
# Oncorhynchus tshawytscha AF392054

####################################################################
# Use DECIPHER to align and plot the alignment 
####################################################################

seqs <- readDNAStringSet(fastafile)
dna <- DNAStringSet(seqs)
names(dna) <- c("Consensus","tshawytscha AF3920540","keta AP010773","clarkii AY886762","nerka EF055889","kisutch EF126369","gorbuscha EF455489","mykiss L29771","MiFish F")
aligned_DNA <- AlignSeqs(dna)


