## Impact of JUST construction (ignoring culverts for now)
## EAA 
## 9/9/22

library(tidyverse)
library(here)
library(lme4)
#library(epower)

df <- readRDS(here("Output","salmonids_abs_abundance_posterior.RDS"))

df <- df %>% 
  group_by(time, creek, station, species) %>% 
  nest() 

alpha = 0.25
degree.freedom = 1499
t.score = qt(p=alpha/2, df=degree.freedom, lower.tail=F)

down.data <- df %>% 
  filter(station == 1) %>% 
  mutate(mean_down = mean(unlist(data))) %>% 
  mutate(sd_down= sd(unlist(data))) %>% 
  mutate(se_down = sd_down/sqrt(1500)) %>% 
  mutate(me2575_down = se_down*t.score) %>% 
  mutate(lb25_down = mean_down - me2575_down) %>% 
  mutate(ub75_down = mean_down + me2575_down) %>% 
  mutate(treat_control = case_when(creek == "4Pad11" ~ "treatment", 
                                   creek !="4Pad11" ~ "control"))  %>% 
  group_by(creek, species) %>% 
  mutate(creek_species_max = max(mean_down)) %>% 
  separate(time, into=c("month","year"), sep=2, remove=FALSE) %>% 
  mutate(scaled_mean = mean_down / creek_species_max) %>% 
  filter(str_detect(species, "clarkii"))
  

ggplot(down.data, aes(x=time, y=log10(mean_down), shape = creek, size=4)) +
#ggplot(down.data, aes(x=time, y=log10(scaled_mean), shape = creek)) +
  geom_point(aes(color= treat_control)) +
  #scale_x_date(date_labels = "%m-%Y") + 
  #geom_segment(aes(x = time, xend = time, y = log10(lb25_down), yend = log10(ub75_down))) +
  facet_wrap(~species) +
  theme_bw() +
  scale_shape_manual(values=c(0,1,3,17,16,4)) + 
  scale_color_manual(values=c('#999999','#56B4E9')) +
  #scale_fill_manual(values=c('#999999','#999999','#999999','#56B4E9','#999999')) +
  labs(y="Log10 Downstream Concentration (copies / L water)", x="Month")


testdf <- down.data %>% 
  select(c(mean_down, time, creek, station, species, treat_control))

testlm <- lmer(mean_down ~ time + (1 | creek) + (1+ treat_control | creek), testdf)
testlm <- lmer(mean_down ~ (treat_control | creek) + (1 | creek) + (1 | time), testdf)



### example
 
lm <- lmer(Reaction ~ Days + (Days | Subject), sleep)
parsedFormula <- lFormula(formula = Reaction ~ Days + (Days | Subject),
                          data = sleepstudy)
devianceFunction <- do.call(mkLmerDevfun, parsedFormula)
optimizerOutput <- optimizeLmer(devianceFunction)
mkMerMod( rho = environment(devianceFunction),
          opt = optimizerOutput,
          reTrms = parsedFormula$reTrms,
          fr = parsedFormula$fr)
