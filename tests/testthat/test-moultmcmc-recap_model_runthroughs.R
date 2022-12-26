library(dplyr)
data(recaptures)

#recaptures %>% group_by(id) %>% summarize(start_date = unique(start_date))
#mcmc_iter = 200
#mcmc_refresh = max(mcmc_iter/4,1)
test_that("moultmcmc uz2_recap works", {
  m2r = moultmcmc(moult_column = "pfmg_sampled",
                  date_column = "date_sampled",
                  id_column = "id",
                  data = recaptures,
                  flat_prior = FALSE,
                  type = 2,
                  log_lik = FALSE,
                  chains = 1,
                  cores = 2,
                  control = list(adapt_delta = 0.99, max_treedepth = 11),
                  iter = 200,
                  refresh = 50,
                  open_progress = FALSE)
expect_s3_class(m2r, "moultmcmc")


as.data.frame(ranef(m2r)) %>% tibble::rownames_to_column('id') %>% left_join(recaptures %>% group_by(id) %>% summarize(start_date = unique(start_date), n_moult = sum(pfmg_sampled != 0 & pfmg_sampled != 1), n_total = n())) -> joined_df
#ggplot(joined_df, aes(x = start_date-196.83, y = mean, ymin = `2.5%`, ymax = `97.5%`, col = factor(n_moult), pch = factor(n_total), label = id)) + geom_pointrange() + geom_text(col = 'black', nudge_y=0.1) + geom_abline(slope = 1, intercept = 0)

})


test_that("moultmcmc uz2l_recap works", {
  uz2rl = moultmcmc(moult_column = "pfmg_sampled",
                    date_column = "date_sampled",
                    id_column = "id",
                    lump_non_moult = TRUE,
                    type = 2,
                    data = recaptures,
                    log_lik = FALSE,
                    chains = 1,
                    iter = 200)
  expect_s3_class(uz2rl, "moultmcmc")
})
test_that("moultmcmc uz2r_active_moult_only", {
  uz2ram = moultmcmc(moult_column = "pfmg_sampled",
                     date_column = "date_sampled",
                     id_column = "id",
                     lump_non_moult = FALSE,
                     active_moult_recaps_only = TRUE,
                     type = 2,
                     data = recaptures,
                     log_lik = FALSE,
                     chains = 2,
                     cores = 2,
                     iter = 400)
  expect_s3_class(uz2ram, "moultmcmc")
})
test_that("moultmcmc uz2r_active_moult_only with unmodelled heterogeneity in tau", {
  uz2ram2 = moultmcmc(moult_column = "pfmg_sampled",
                     date_column = "date_sampled",
                     id_column = "id",
                     lump_non_moult = FALSE,
                     active_moult_recaps_only = TRUE,
                     type = 2,
                     data = recaptures2,
                     log_lik = FALSE,
                     chains = 2,
                     cores = 2,
                     iter = 400)
  expect_s3_class(uz2ram2, "moultmcmc")
})
test_that("uz2l_phi_approx", {
  uz2rapprox = moultmcmc(moult_column = "pfmg_sampled",
                             date_column = "date_sampled",
                             id_column = "id",
                             lump_non_moult = FALSE,
                             use_phi_approx = TRUE,
                             data = recaptures,
                             log_lik = FALSE,
                             chains = 1,
                             iter = 200)
  expect_s3_class(uz2rapprox, "moultmcmc")
})
test_that("uz3r works", {
  uz3 = moultmcmc(moult_column = "pfmg_sampled",
                  date_column = "date_sampled",
                  id_column = "id",
                  data = subset(recaptures, pfmg_sampled != 0 & pfmg_sampled != 1),
                  type = 3,
                  log_lik = FALSE,
                  chains = 2,
                  cores = 2,#should work on gh-actions
                  iter = 300)
  expect_s3_class(uz3, "moultmcmc")
})
test_that("uz4r works", {
  uz4 = moultmcmc(moult_column = "pfmg_sampled",
                  date_column = "date_sampled",
                  id_column = "id",
                  data = subset(recaptures, pfmg_sampled != 0),
                  type = 4,
                  log_lik = FALSE,
                  chains = 1,
                  iter = 200)
  expect_s3_class(uz4, "moultmcmc")
})
test_that("uz5r works", {
  uz5 = moultmcmc(moult_column = "pfmg_sampled",
                          date_column = "date_sampled",
                          id_column = "id",
                          data = subset(recaptures, pfmg_sampled != 1),
                  type = 5,
                    log_lik = FALSE,
                    chains = 1,
                    iter = 200)
  expect_s3_class(uz5, "moultmcmc")
})

