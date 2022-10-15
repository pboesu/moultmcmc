## code to prepare `recaptures` dataset goes here
library(dplyr)
recaptures = read.csv('data-raw/siskin03_sim_001.csv') %>% mutate(id = factor(id))
usethis::use_data(recaptures, overwrite = TRUE)
recaptures2 = read.csv('data-raw/siskin37_sim_001.csv') %>% mutate(id = factor(id))
usethis::use_data(recaptures2, overwrite = TRUE)


siskins = readRDS('data-raw/siskin_anon.rds')
usethis::use_data(siskins, overwrite = TRUE)
