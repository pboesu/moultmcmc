#sandbox script for visualising recaptures data
library(moultmcmc)
library(dplyr)
library(ggplot2)

#data(recaptures)
recaptures = read.csv('data-raw/siskin03_sim_001.csv') %>% mutate(id = factor(id))
summary(recaptures)


recaptures %>% group_by(id) %>% summarize(n = n(), n_active_moult = sum(pfmg_sampled != 0 & pfmg_sampled != 1), start_date = unique(start_date), duration = unique(duration)) %>% arrange(desc(n_active_moult)) -> recapture_statistics
summary(recapture_statistics)


left_join(recaptures, recapture_statistics) -> data_joined
data_joined %>% ggplot(aes(x = date_sampled, y = pfmg_sampled, col = as.factor(n), groups = interaction(id,as.factor(n)))) + geom_point() + geom_line(lty = 2) + geom_line(data = filter(data_joined, n_active_moult > 1 & pfmg_sampled != 1 & pfmg_sampled != 0), lwd = 1)

fit_r <- uz2_linpred_recap('pfmg_sampled', 'date_sampled', id_column = 'id', data = recaptures, flat_prior = FALSE, chains = 2, log_lik = FALSE, iter = 4000, control=list(adapt_delta=0.99), cores=2)
rstan::traceplot(fit_r$stanfit)
summary_table(fit_r)

fit_ar <- uz2_linpred_recap('pfmg_sampled', 'date_sampled', id_column = 'id', data = recaptures, flat_prior = FALSE, chains = 2, cores=2, log_lik = FALSE, iter = 3000, active_moult_recaps_only = TRUE)

fit_3r_test <- uz3_linpred_recap('pfmg_sampled', 'date_sampled', id_column = 'id', data = recaptures, flat_prior = FALSE, chains = 2,cores=2, log_lik = FALSE, iter = 200)

fit_3r <- uz3_linpred_recap('pfmg_sampled', 'date_sampled', id_column = 'id', data = recaptures, flat_prior = FALSE, chains = 2,cores=2, log_lik = FALSE, iter = 2000)
rstan::traceplot(fit_3r$stanfit)
s3r = summary_table(fit_3r)
ind_intercepts = data.frame(start_date_est = colMeans(rstan::extract(fit_3r$stanfit, pars='mu_ind_out')$mu_ind_out),
                            start_date_lci = matrixStats::colQuantiles(rstan::extract(fit_3r$stanfit, pars='mu_ind_out')$mu_ind_out, probs = 0.025),
                            start_date_uci = matrixStats::colQuantiles(rstan::extract(fit_3r$stanfit, pars='mu_ind_out')$mu_ind_out, probs = 0.975),
                            id = unique(fit_3r$terms$id_original)) %>% #this only works because data happen to be ordered by individual!!
  left_join(recapture_statistics)


intercept_lookup_2r <- data.frame(index = unique(fit_r$individual_ids$index), id = unique(fit_r$individual_ids$id))

ind_intercepts_2r = data.frame(start_date_est = colMeans(rstan::extract(fit_r$stanfit, pars='mu_ind_out')$mu_ind_out),
                            start_date_lci = matrixStats::colQuantiles(rstan::extract(fit_r$stanfit, pars='mu_ind_out')$mu_ind_out, probs = 0.025),
                            start_date_uci = matrixStats::colQuantiles(rstan::extract(fit_r$stanfit, pars='mu_ind_out')$mu_ind_out, probs = 0.975)) %>% tibble::rownames_to_column(var='index') %>% mutate(index= as.numeric(index)) %>% left_join(intercept_lookup_2r) %>% left_join(recapture_statistics)

ggplot(ind_intercepts, aes(x = start_date, y = start_date_est, ymin = start_date_lci, ymax = start_date_uci, pch = n_active_moult >= 2, fill = duration - mean(recapture_statistics$duration))) + geom_pointrange()+geom_abline(slope = 1, interept = 0) + scale_fill_gradient2() + scale_shape_manual(values = c(21,24)) + theme_classic()

ggplot(filter(ind_intercepts_2r, between(start_date_est,0,365)), aes(x = start_date, y = start_date_est, pch = as.factor(n_active_moult),  fill = duration - mean(recapture_statistics$duration))) + geom_point()+geom_abline(intercept = 0, slope=1)+scale_shape_manual(values = c(21,22,24)) + scale_fill_gradient2()



fit_naive <- moultmcmc('pfmg_sampled', 'date_sampled', data = recaptures, flat_prior = FALSE, chains = 2, log_lik = FALSE)
fit_3naive <- moultmcmc('pfmg_sampled', 'date_sampled', data = recaptures, flat_prior = FALSE, chains = 2, log_lik = FALSE, type=3)
compare_plot(fit_naive, fit_3r, fit_3naive)

moult_plot(fit_naive)
moult_plot(fit_ar)
