## Metadata from backpack 

backpackpath <- here("Input","backpack")
files <- list.files(path = backpackpath, pattern = "*.csv", recursive = T, full.names = T)

file_vol <- matrix(NA, nrow=length(files), ncol=2)
for (i in 1:length(files)) {
  filename = files[i]
  
  lognum = unlist(strsplit(filename, "/"))[10]
  lognum = str_sub(lognum, start=1, end=5)
  lognum = as.numeric(lognum)
  
  data = read_csv(filename, n_max = 10, col_names = FALSE)
  vol = as.numeric(data[5,2])
  
  file_vol[i,1]=lognum
  file_vol[i,2]=vol
}
file_vol <- as_tibble(file_vol) %>% 
  rename(LogID = V1, Volume = V2)
  

backpackmeta <- read_csv(here("Input","vol_filtered.csv"))

backpackmeta <- backpackmeta %>% 
  left_join(file_vol, by="LogID") %>% 
  mutate(Volume = case_when(is.na(Volume.x) ~ Volume.y,
                           !is.na(Volume.x) ~ Volume.x)) %>% 
  dplyr::select(-c(Volume.x, Volume.y)) %>% 
  mutate(Adj_Vol = case_when(Trident == 0 ~ Volume,
                             Trident == 1 ~ Volume/3)) %>% 
  dplyr::select(c(Sample,Adj_Vol))

write_rds(backpackmeta, here("Input", "adj_vol_filtered.RDS"))
