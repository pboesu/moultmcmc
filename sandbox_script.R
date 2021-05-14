#test type 3-5 models
library(moult)
library(moultmcmc)
library(dplyr)
library(ggplot2)
library(rstan)
data("sanderlings")
sanderlings$MCat <- case_when(sanderlings$MIndex == 0 ~ 1,
                              sanderlings$MIndex == 1 ~3,
                              TRUE ~ 2)

m1 = moult(MIndex ~ Day,data = sanderlings, type = 1)
uz1 = uz1_linpred(moult_cat_column = "MCat",
            date_column = "Day",
            data = sanderlings,
            log_lik = FALSE)
summary_table(m1)
summary_table(uz1)
compare_plot(m1, uz1)

m2 = moult(MIndex ~ Day,data = sanderlings, type = 2)
uz2 = uz2_linpred("MIndex",
                  date_column = "Day",
                  data = sanderlings,
                  log_lik = FALSE)
compare_plot(m2, uz2)

#compare type 1/2
compare_plot(m1,m2)
compare_plot(uz1,uz2)

m3 = moult(MIndex ~ Day,data = sanderlings, type = 3)
uz3 = uz3_linpred("MIndex",
                  date_column = "Day",
                  data = sanderlings,
                  log_lik = FALSE)
compare_plot(m3, uz3)
stan_trace(uz3)

m4 = moult(MIndex ~ Day,data = sanderlings, type = 4)
uz4 = uz4_linpred("MIndex",
                  date_column = "Day",
                  data = sanderlings,
                  log_lik = FALSE)
compare_plot(m4, uz4)

m5 = moult(MIndex ~ Day,data = sanderlings, type = 5)
uz5 = uz5_linpred("MIndex",
                  date_column = "Day",
                  data = sanderlings,
                  log_lik = FALSE)
compare_plot(m5,uz5)

data(weavers)
if (is.numeric(weavers$Moult)) {
  scores <- format(weavers$Moult, scientific = FALSE, trim = TRUE) } else {
    scores = weavers$Moult }
mscores = substr(scores, 1, 9)
feather.mass <- c(10.4, 10.8, 11.5, 12.8, 14.4, 15.6, 16.3, 15.7, 15.7)
weavers$pfmg <- ms2pfmg(mscores, feather.mass)
weavers$day <- date2days(weavers$RDate, dateformat = "yyyy-mm-dd",
                            startmonth = 8)
weavers <- na.omit(weavers)
plot(pfmg ~ day, weavers)
weavers_trunc <- filter(weavers, day > 50 & day <320)
points(pfmg ~ day, weavers_trunc, col = 'green')
 m88.3 <- moult(pfmg ~ day, data = weavers_trunc, type=3)
 uz88.3 <- uz3_linpred('pfmg',
                       'day',
                       log_lik = FALSE,
                       data = weavers_trunc)#type 2 fails on start values. possibly numerical underflow with bad start values?? - but ML converges, so a smarter way of estimating start values necessary? Type 3 works, but highly heterogeneous sampling times for each chain - use non-flat priors??
summary(m88.3)
compare_plot(m88.3,uz88.3)
weavers$MCat <- case_when(weavers$pfmg == 0 ~ 1,
                          weavers$pfmg == 1 ~3,
                              TRUE ~ 2)
table(weavers$MCat)
m88.5 <- moult(pfmg ~ day, data = weavers, type = 5)
m88.4 <- moult(pfmg ~ day, data = weavers, type = 4)


##combined model tests
set.seed(1234)
sanderlings_combined <- get(data("sanderlings"))
sanderlings_combined$MCat <- case_when(sanderlings$MIndex == 0 ~ 1,
                              sanderlings$MIndex == 1 ~3,
                              TRUE ~ 2)
remove_m_scores <- function(data, prop = 0.5){
  #randomly subsample 1/3 of rows and set moult score to NA
  m_birds = which(data$MCat == 2)
  cats <- sample(m_birds,size = floor(length(m_birds)*prop))
  data$MIndex[cats] <- NA
  return(data)}
table(sanderlings_combined$MCat)

uz12_25 <- uz12_linpred("MIndex",
             "MCat",
             date_column = "Day",
             data = remove_m_scores(sanderlings_combined, 0.25),
             log_lik = FALSE)

uz12_50 <- uz12_linpred("MIndex",
                        "MCat",
                        date_column = "Day",
                        data = remove_m_scores(sanderlings_combined, 0.5),
                        log_lik = FALSE)
uz12_75 <- uz12_linpred("MIndex",
                        "MCat",
                        date_column = "Day",
                        data = remove_m_scores(sanderlings_combined, 0.75),
                        log_lik = FALSE)
uz12_80 <- uz12_linpred("MIndex",
                        "MCat",
                        date_column = "Day",
                        data = remove_m_scores(sanderlings_combined, 0.8),
                        log_lik = FALSE)
uz12_90 <- uz12_linpred("MIndex",
                        "MCat",
                        date_column = "Day",
                        data = remove_m_scores(sanderlings_combined, 0.9),
                        log_lik = FALSE)

compare_plot(m1, uz12)
compare_plot(m2, uz12)
compare_plot(uz1,uz12)
compare_plot(uz2,uz12)

bind_rows(summary_table(uz1),
          #summary_table(uz12_90),
          #summary_table(uz12_80),
          summary_table(uz12_75),
          summary_table(uz12_50),
          summary_table(uz12_25),
          summary_table(uz2),
          .id='prop_cat') %>%
  left_join(tibble(prop_cat = as.character(1:5), prop_score = c(0,25,50,75,100))) %>%
  filter(!(parameter %in% c('lp__', 'log_sd_(Intercept)'))) %>%
  ggplot(aes(x = prop_score, y = estimate, ymin = lci, ymax = uci)) + geom_pointrange() + facet_wrap(~parameter, scales = 'free_y') + xlab("Proportion of scored birds in active moult category") + theme_classic()