## code to prepare `recaptures` dataset goes here
library(dplyr)
recaptures = read.csv('data-raw/siskin03_sim_001.csv') %>% mutate(id = factor(id))
usethis::use_data(recaptures, overwrite = TRUE)
