#sandbox type 5 model with simulated data
library(moult)
library(moultmcmc)
library(dplyr)
library(ggplot2)

sim_data <- readRDS('../BTO94_moult_phenology/data/simulated_datasets/simulated_no_recaptures_start_duration_temp_no_no_linear_set99.rds') %>% ungroup()
set.seed(456)
sim_data_sub <- slice_sample(sim_data, n = 20)

m5 = moult(moult_score ~ yday,data = sim_data_sub, type = 5)
uz5 = uz5_linpred("moult_score",
                  date_column = "yday",
                  data = sim_data_sub,
                  log_lik = FALSE,
                  cores = 3,
                  chains = 3,
                  control = list(adapt_delta = 0.9),
                  iter = 3000)
compare_plot(m5,uz5, names = c('ML','MCMC')) + geom_hline(yintercept = c(unique(sim_data$pop_duration),unique(sim_data$pop_start_date),unique(sim_data$pop_start_date_sd)))
rstan::stan_trace(uz5$stanfit)
pairs(uz5$stanfit)
fixef(uz5)
moult_plot(m5)
moult_plot(m5, data = sim_data_sub)
moult_plot(uz5)
moult_plot(uz5, data = sim_data_sub)

post_samples <- fixef(uz5, summary = FALSE)

plot(0:365, dnorm(0:365, mean = post_samples[1,1], sd = post_samples[1,'sd_(Intercept)']),type = 'l')
for (i in 1:nrow(post_samples)) lines(0:365,dnorm(0:365, mean = post_samples[i,1], sd = post_samples[i,'sd_(Intercept)']))

plot(c(post_samples[1,1],post_samples[1,1]+post_samples[1,'duration_(Intercept)']), c(0,1),type = 'l', xlim = c(100,365))
for (i in 1:nrow(post_samples)) lines(c(post_samples[i,1],post_samples[i,1]+post_samples[i,'duration_(Intercept)']),c(0,1))
for (i in 1:nrow(post_samples)) lines(c(post_samples[i,1],post_samples[1,1]+post_samples[i,'duration_(Intercept)']+1.96*post_samples[i,'sd_(Intercept)']),c(0,1), col = 'red')



plot(moult_score ~ yday, data = sim_data_sub, col = moult_cat)
segments(coef(m5)[2], 0, coef(m5)[2]+coef(m5)[1], 1)
segments(fixef(uz5)[1,1], 0, fixef(uz5)[1,1]+fixef(uz5)[2,1], 1, lty = 2)
lines(1:365, dnorm(1:365, mean = coef(m5)[2], sd = coef(m5)[3]))
lines(1:365, dnorm(1:365, mean = fixef(uz5)[1,1], sd = fixef(uz5)[4,1]), lty = 2)

m4 = moult(moult_score ~ yday,data = sim_data_sub, type = 4)
uz4 = uz4_linpred("moult_score",
                  date_column = "yday",
                  data = sim_data_sub,
                  log_lik = FALSE,
                  cores = 4,
                  control = list(adapt_delta = 0.9))
compare_plot(m4,uz4, names = c('ML','MCMC')) + geom_hline(yintercept = c(unique(sim_data$pop_duration),unique(sim_data$pop_start_date),unique(sim_data$pop_start_date_sd)))
pairs(uz4$stanfit)
