#sandbox type 5 model with simulated data
sim_data <- readRDS('../BTO94_moult_phenology/data/simulated_datasets/simulated_no_recaptures_start_duration_temp_no_no_linear_set99.rds') %>% ungroup()
set.seed(456)
sim_data_sub <- slice_sample(sim_data, n = 20)

m5 = moult(moult_score ~ yday,data = sim_data_sub, type = 5)
uz5 = uz5_linpred("moult_score",
                  date_column = "yday",
                  data = sim_data_sub,
                  log_lik = FALSE,
                  cores = 3,
                  chains = 3)
compare_plot(m5,uz5, names = c('ML','MCMC')) + geom_hline(yintercept = c(unique(sim_data$pop_duration),unique(sim_data$pop_start_date),unique(sim_data$pop_start_date_sd)))
stan_trace(uz5$stanfit)

m4 = moult(moult_score ~ yday,data = sim_data_sub, type = 4)
uz4 = uz4_linpred("moult_score",
                  date_column = "yday",
                  data = sim_data_sub,
                  log_lik = FALSE,
                  cores = 4)
compare_plot(m4,uz4, names = c('ML','MCMC')) + geom_hline(yintercept = c(unique(sim_data$pop_duration),unique(sim_data$pop_start_date),unique(sim_data$pop_start_date_sd)))
