## Look at fold change in down / up over time by creek
## EAA
## 9/9/22 

library(tidyverse)
library(here)


df <- readRDS(here("Output", "salmonids_fc.RDS"))  

df <- df %>% 
  mutate(treat_control = case_when(creek == "4Pad11" ~ "treatment", 
                                   creek !="4Pad11" ~ "control"))  


ggplot(df, aes(x=time, y=log10(abs(mean_fc)), shape = creek, color=as.character(Sign), size=4)) +
  geom_point() +
  geom_hline(yintercept=0) +
  #geom_segment(aes(x = time, xend = time, y = log10(abs(lb25_fc))*Sign, yend = log10(abs(ub75_fc))*Sign)) +
  facet_grid(~species ~treat_control, scales="free") +
  scale_shape_manual(values=c(0,1,3,17,16,4)) + 
  scale_color_manual(values=c('#999999','#56B4E9')) +
  #scale_color_discrete(name = "Legend", labels = c("Negative (Up > Down)", "Positive (Down > Up)")) +
  labs(y="Log10(asb(mean fold change))", x="Month")
