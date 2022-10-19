## code to prepare `recaptures` dataset goes here
library(dplyr)
recaptures = read.csv('data-raw/siskin03_sim_001.csv') %>% mutate(id = factor(id))
usethis::use_data(recaptures, overwrite = TRUE)
recaptures2 = read.csv('data-raw/siskin37_sim_001.csv') %>% mutate(id = factor(id))
usethis::use_data(recaptures2, overwrite = TRUE)


siskins = readRDS('data-raw/siskin_sim.rds') %>%
  select(date_sampled, pfmg_sampled, id) %>%
  rename(yday = date_sampled, pfmg = pfmg_sampled)
usethis::use_data(siskins, overwrite = TRUE)

#reexport weavers data
library(moult)
data(weavers)

if (is.numeric(weavers$Moult)) {
  scores <- format(weavers$Moult, scientific = FALSE, trim = TRUE)
} else {
  scores <- weavers$Moult
}
mscores <- gsub('8', '0', substr(scores, 1, 9))#really this should distinguish position of the 8, if smaller numbers precede it should be a zero, if 5s or 4s follow it should really be 5

feather.mass <- c(10.4, 10.8, 11.5, 12.8, 14.4, 15.6, 16.3, 15.7, 15.7)
weavers$pfmg <- ms2pfmg(mscores, feather.mass)
weavers$day <- date2days(weavers$RDate, dateformat = "yyyy-mm-dd",
                         startmonth = 8)
weavers_processed <- weavers
usethis::use_data(weavers_processed, overwrite = TRUE)
