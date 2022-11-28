
qm_pull_posteriors <- function(bayes_out){
  require(tidyverse)
  
  grabthese <- grep("int_samp_small", names(bayes_out$Bayes_modelfit@sim$samples[[1]]))
  postSamples <- bayes_out$Bayes_modelfit@sim$samples[[1]][grabthese]
  
  a <- as.data.frame(postSamples) 
  mynames <- names(a)
  b <- t(a) %>% 
    as_tibble() %>% 
    rownames_to_column("samp") %>% 
    mutate(samp = mynames) %>% 
    mutate(samp = str_remove_all(samp,  "int_samp_small\\.")) %>% 
    separate(samp, into = c("bottle", "species", "x"), "\\.") %>% 
    dplyr::select(-x)
  
  b <- b %>% 
    group_by(bottle, species) %>% 
    nest()
  
  b$bottle <- rep(row.names(bayes_out$Bayes_estimates),
                  times = ncol(bayes_out$Bayes_estimates))
  b$species <- rep(colnames(bayes_out$Bayes_estimates),
                   each = nrow(bayes_out$Bayes_estimates))
  
  return(b)
}
