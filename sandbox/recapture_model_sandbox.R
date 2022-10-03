#sandbox script for visualising recaptures data
library(moultmcmc)
library(dplyr)
library(ggplot2)

data(recaptures)
summary(recaptures)


recaptures %>% group_by(id) %>% summarize(n = n(), n_active_moult = sum(pfmg_sampled != 0 & pfmg_sampled != 1), start_date = unique(start_date), duration = unique(duration)) %>% arrange(desc(n_active_moult)) -> recapture_statistics
summary(recapture_statistics)


left_join(recaptures, recapture_statistics) -> data_joined
data_joined %>% ggplot(aes(x = date_sampled, y = pfmg_sampled, col = as.factor(n), groups = interaction(id,as.factor(n)))) + geom_point() + geom_line(lty = 2) + geom_line(data = filter(data_joined, n_active_moult > 1 & pfmg_sampled != 1 & pfmg_sampled != 0), lwd = 1)

fit_r <- uz2_linpred_recap('pfmg_sampled', 'date_sampled', id_column = 'id', data = recaptures, flat_prior = FALSE, chains = 2, log_lik = FALSE, iter = 4000, control=list(adapt_delta=0.99), cores=2)
rstan::traceplot(fit_r$stanfit)
summary_table(fit_r)

fit_ar <- uz2_linpred_recap('pfmg_sampled', 'date_sampled', id_column = 'id', data = recaptures, flat_prior = FALSE, chains = 2, cores=2, log_lik = FALSE, iter = 3000, active_moult_recaps_only = TRUE)

fit_3r <- uz3_linpred_recap('pfmg_sampled', 'date_sampled', id_column = 'id', data = recaptures, flat_prior = FALSE, chains = 2,cores=2, log_lik = FALSE, iter = 2000)
rstan::traceplot(fit_3r$stanfit)
s3r = summary_table(fit_3r)
ind_intercepts = data.frame(start_date_est = colMeans(rstan::extract(fit_3r$stanfit, pars='mu_ind_out')$mu_ind_out),
                            start_date_lci = matrixStats::colQuantiles(rstan::extract(fit_3r$stanfit, pars='mu_ind_out')$mu_ind_out, probs = 0.025),
                            start_date_uci = matrixStats::colQuantiles(rstan::extract(fit_3r$stanfit, pars='mu_ind_out')$mu_ind_out, probs = 0.975),
                            id = unique(fit_3r$terms$id_original)) %>% left_join(recapture_statistics)
ggplot(ind_intercepts, aes(x = start_date, y = start_date_est, ymin = start_date_lci, ymax = start_date_uci, pch = n_active_moult >= 2, fill = duration - mean(recapture_statistics$duration))) + geom_pointrange()+geom_abline(slope = 1, interept = 0) + scale_fill_gradient2() + scale_shape_manual(values = c(21,24)) + theme_classic()



fit_naive <- moultmcmc('pfmg_sampled', 'date_sampled', data = recaptures, flat_prior = FALSE, chains = 2, log_lik = FALSE)
fit_3naive <- moultmcmc('pfmg_sampled', 'date_sampled', data = recaptures, flat_prior = FALSE, chains = 2, log_lik = FALSE, type=3)
compare_plot(fit_naive, fit_3r, fit_3naive)

moult_plot(fit_naive)
moult_plot(fit_ar)
