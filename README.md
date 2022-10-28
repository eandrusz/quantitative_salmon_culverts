# quantitative_salmon_culverts
metabarcoding + qPCR for salmonids up/downstream of culverts

# Introduction
This repository includes all the scripts, fuctions, and files associated with the XXX manuscript. The work here takes

Metabarcoding data is generated from 13 Miseq runs (including mock community samples) and ~50 plates of quantitative PCR were run. The general order of operations is as follows:

1) Process qPCR data
  a) check for inhibition (cutthroat assay is multiplexed with an internal positive control)
  b) use Bayesian model to convert Ct values to copies/uL of extract for each sample
  c) account for dilution for inhibition and total volume of water filtered to get copies/L water filtered 
  d) repeat steps b and c for coho qPCR assay
  e) NOT YET DONE - CORRECT FOR FLOW
  
2) Process metabarcoding data 
  a) take ASVs from all environmental samples and assign taxonomy and format
  b) prepare mock community samples with known proportions before PCR and percent of reads after PCR 
  c) subset data for only species of interest
  d) run quantitative metabarcoding model (QM) on environmental samples using mock community to get corrected proportions of starting DNA concentrations from environmental samples
  
3) Merge qPCR and metabarcoding data 
  a) using the absolute concentrations of cutthroat DNA from qPCR, divide by the proportion of cutthroat DNA (of the given subset of species used in step 2c) to solve for the total concentration of DNA (of the given subset of species used in 2c) in each water sample
  b) using the total concentration of DNA (subset, see above) in each water sample, multiply the remaining species proportions by the total DNA to get quantitative (copies/L water) estimates of all species in each water sample
  
4) Run time series model for each species 
  a) for each species, run time series model with all data - extra term gives effect of culvert per species, per creek, per time
  b) to evaluate impact of construction, remove the time points in Padden Creek during construction and run the time series model with those as missing data, allowing for estimate of the contrafactual scenario that no construction occurred
  c) NOT DONE YET - compare the observed patterns (with construction) to the expected patterns without construction occuring in Padden Creek
  
5) Visualization 
  a) visualization of any/all raw, intermediate, final data generated in 1-4. 

# Repo Structure

## Input
This contains all data that are required for any of the analyses (qPCR, metabarcoding, etc.). The Input folder only contains raw data - any intermediate data goes to Output. For example, the only data in Input for qPCR are the files containing Ct values for each sample and the conversion to copies/uL are then in the Output folder. 

## Scripts
The scripts folder has several subfolders for the 

##

