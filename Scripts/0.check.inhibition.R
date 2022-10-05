## Check for inhibition using IPC kit (multiplexed with cutthroat assay)
## EAA 
## 8/18/22

## Load packages
library(tidyverse)
library(here)
library(readxl)

## Point to files for cutthroat and inhibition (IPC)
qPCRfiles <- list.files(path = here("Input","qpcr","Results","CUT"), pattern = "*data", recursive = T, full.names = T)
qPCRmeta <- read_excel(path=here("Input","qpcr","qPCR_samples_ALL.xlsx"), sheet = "info_samples")

## Write function to simplify output to be just IPC data 
simple_IPC_plate <- function(plate_num) {
  # first just get file open - they are listed by date so plate number is not in order
  filename = qPCRfiles[plate_num]
  filename2 = unlist(strsplit(filename, "_"))[4]
  filename3 = unlist(strsplit(filename2,"-"))[3]
  filename4 = str_sub(filename3, start=1, end =3)
  plate <- read_excel(qPCRfiles[plate_num], sheet = "Results", skip=6)
  colnames(plate)[7] <- "Ct"  
  colnames(plate)[2] <- "Sample"
  colnames(plate)[3] <- "Assay"
  plate$Plate <- as.numeric(filename4) # now right the actual plate number that it is (i.e., what Megan called it)
  
  # now only look at IPC data and Ct values 
  IPCdata <- plate %>% 
    filter(Assay == "IPC") %>% 
    select(c("Plate", "Well","Sample", "Ct")) %>% 
    filter(Ct != "Undetermined") %>% 
    mutate(Ct = as.numeric(Ct)) 

  return(IPCdata)
}

## Write function to read in plate, only look at IPC data, look at the IPC value in the NTC and the standards (which should not be inhibited), determine threshold for what is deemed inhibited 
determine_thresh <- function(IPCdata) {

# using ntcs, find mean and 2 sds, ct threshold = mean + 2 sds
IPCntc <- IPCdata %>% 
  filter(str_detect(Sample,"NTC")) %>% 
  filter(! is.na(Ct)) 
IPCntc_mean <- mean(IPCntc$Ct)
IPCntc_sd <- sd(IPCntc$Ct)
IPCntcthresh <- IPCntc_mean + 2*IPCntc_sd

# standards are also basically NTCs for this purpose because they should have NO inhibitors (they are just gBlocks)
IPCstds <- IPCdata %>% 
  filter(str_detect(Sample,"St")) %>% 
  filter(! is.na(Ct))
IPCstds_mean <- mean(IPCstds$Ct)
IPCstds_sd <- sd(IPCstds$Ct)
IPCstdsthresh <- IPCstds_mean + 2*IPCstds_sd

# take whichever is lower (ntc or stds)
thresh <- min(IPCntcthresh, IPCstdsthresh)
return(thresh)
}

## Write a function to select which samples are in fact inhibited
inhib_samp <- function(IPCdata,thresh) {
  inhibited <- IPCdata %>% 
    filter(!str_detect(Sample,"St")) %>% 
    filter(! is.na(Ct)) %>% 
    group_by(Sample) %>% 
    mutate(meanCt = mean(Ct)) %>% 
    mutate(inhib = meanCt > thresh) %>% 
    filter(inhib == TRUE) %>% 
    select(c(Plate,Sample,meanCt,inhib)) %>% 
    distinct()

  return(inhibited)
}


## OK TEST IT 
#plate_num = 1
#plate <- simple_IPC_plate(plate_num)
#thresh <- determine_thresh(plate)
#inhibited <- inhib_samp(plate, thresh) 
# works! 

# now write a function to take only samples in metadata that we are actually using 
still_inhib_samp <- function(plate_num, inhibited) {
  filename = qPCRfiles[plate_num]
  filename2 = unlist(strsplit(filename, "_"))[4]
  
  touse <- qPCRmeta %>% 
    filter(qPCR_no == filename2) %>% 
    filter(target == "cutt") %>%  # all the inhibition is on cutthroat plates
    filter(use == 1) %>% # sample by sample, plate by plate (only use good SC plates, and only use the samples we THOUGHT were not inhibited using the 3 ct threshold)
    filter(! str_detect(sample, "St")) %>% 
    filter(! str_detect(sample, "NTC")) %>% 
    filter(! str_detect(sample, "EB")) 
  
  stillinhibited <- inhibited %>% 
    filter(Sample %in% touse$sample)
  
  return(stillinhibited)
}

#stillinhibited <- still_inhib_samp(plate_num,inhibited)
# works again! 


#### RUN FOR ALL PLATES 
redoinhib <- list()
thresh <- vector()
for (i in 1:length(qPCRfiles)) {
  plate_num <- i
  plate <- simple_IPC_plate(plate_num)
  thresh[i] <- determine_thresh(plate)
  inhibited <- inhib_samp(plate, thresh[i]) 
  redoinhib[[i]] <- still_inhib_samp(plate_num,inhibited)
}

# plate_num = 18 --> plate 23 does not have IPC (only redoing standard deviations)
# plate_num = 24 --> plate 25 does not have IPC (only redoing standard deviations)

redoinhib <- do.call(rbind, redoinhib)  
write.csv(redoinhib, here("Output","qpcr","cut_fourth_inhibition.csv"), row.names=FALSE)



##### now check for any left overs with high standard deviations
bad_sd_samps <- function(plate_num, IPCdata, thresh) {
  filename = qPCRfiles[plate_num]
  filename2 = unlist(strsplit(filename, "_"))[4]
  
  tousestd <- qPCRmeta %>% 
    filter(qPCR_no == filename2) %>% 
    filter(target == "cutt") %>%  # all the inhibition is on cutthroat plates
    filter(use == 1) %>% # sample by sample, plate by plate (only use good SC plates, and only use the samples we THOUGHT were not inhibited using the 3 ct threshold)
    filter(str_detect(sample, "St")) 
  
  tousesamp <- qPCRmeta %>% 
    filter(qPCR_no == filename2) %>% 
    filter(target == "cutt") %>%  # all the inhibition is on cutthroat plates
    filter(use == 1) %>% # sample by sample, plate by plate (only use good SC plates, and only use the samples we THOUGHT were not inhibited using the 3 ct threshold)
    filter(! str_detect(sample, "St")) %>% 
    filter(! str_detect(sample, "NTC")) %>% 
    filter(! str_detect(sample, "EB")) 
  
  notinhib <- IPCdata %>% 
    filter(!str_detect(Sample,"St")) %>% 
    filter(! is.na(Ct)) %>% 
    group_by(Sample) %>% 
    mutate(meanCt = mean(Ct)) %>% 
    mutate(sdCt = sd(Ct)) %>% 
    mutate(inhib = meanCt > thresh) %>% 
    filter(inhib == FALSE) %>% 
    select(c(Plate,Sample,meanCt,sdCt,inhib)) %>% 
    distinct() %>% 
    filter(Sample %in% tousesamp$sample)
  
  std_sd_thresh <- IPCdata %>% 
    filter(str_detect(Sample,"St")) %>% 
    filter(! is.na(Ct)) %>% 
    group_by(Sample) %>% 
    mutate(meanCt = mean(Ct)) %>% 
    mutate(sdCt = sd(Ct)) %>% 
    select(c(Plate,Sample,meanCt,sdCt)) %>% 
    distinct() 
  
  std_sd_cutoff <- mean(2*std_sd_thresh$sdCt)
  
  badsd <- notinhib %>% 
    mutate(std_sd_cutoff=std_sd_cutoff) %>% 
    mutate(good_sd = sdCt < std_sd_cutoff) %>% 
    filter(good_sd == FALSE)
  
  return(badsd)
}

#### RUN FOR ALL PLATES 
badsd <- list()
for (i in 1:length(qPCRfiles)) {
  plate_num <- i
  plate <- simple_IPC_plate(plate_num)
  thresh[i] <- determine_thresh(plate)
  badsd[[i]] <- bad_sd_samps(plate_num,plate, thresh[i])
}
badsd <- do.call(rbind, badsd)  

## seems like a lot of samples and really low threshold -- use 1 instead
# 1121.2Brn.Up.1 is only rerun


#write.csv(redoinhib, here("Output","qpcr","cut_fourth_inhibition.csv"), row.names=FALSE)

