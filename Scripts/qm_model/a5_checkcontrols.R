
library(tidyverse)
library(ggplot2)
library(here)
library(gridExtra)

taxatable <- read.csv(here("Output","metabarcoding", "taxa_table.csv"))

nonkangaroo <- taxatable %>% 
  filter(! str_detect(Sample_name, "MC")) %>% 
  filter(! str_detect(Sample_name, "Kangaroo")) %>% 
  group_by(Sample_name) %>% 
  mutate(ReadDepth = sum(totalReads)) %>% 
  mutate(propReads = totalReads/ReadDepth) 

nonkangaroo_withkangaroo <- nonkangaroo %>% 
  filter(str_detect(species, "Osphranter"))

nonkangaroo <- nonkangaroo %>% 
  filter(Sample_name %in% nonkangaroo_withkangaroo$Sample_name)

kangaroo <- taxatable %>% 
  filter(str_detect(Sample_name, "Kangaroo")) %>% 
  group_by(Sample_name) %>% 
  mutate(ReadDepth = sum(totalReads)) %>% 
  mutate(propReads = totalReads/ReadDepth) %>% 
  mutate(Sample_name = str_replace(Sample_name, "MiFish.Kangaroo.", ""))

kangaroo %>% 
  ggplot(aes(x=Sample_name, y=propReads, fill=species)) +
  geom_col() + 
  scale_fill_manual(values=c("#E69F00","#E69F00","#E69F00","#E69F00","#999999","#E69F00","#E69F00","#E69F00", "#999999","#999999")) + 
  labs(y="Proportion of Reads", x="Sample") +
  ggtitle("Positive Control Samples")

ggsave(file=here("Output","Figures","check_controls.png"))

nonkangaroo %>% 
  ggplot(aes(x=Sample_name, y=propReads, fill=species)) +
  geom_col() + 
  scale_fill_manual(values=c("#999999", "#999999","#999999","#999999","#999999","#999999","#999999","#999999","#999999","#999999","#999999","#E69F00", "#999999","#999999")) + 
  labs(y="Proportion of Reads", x="Sample") +
  ggtitle("Environmental Samples with Positive Control")

ggsave(file=here("Output","Figures","check_controls2.png"))
