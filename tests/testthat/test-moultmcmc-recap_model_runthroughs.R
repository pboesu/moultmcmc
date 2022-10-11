library(dplyr)
data(recaptures)

#recaptures %>% group_by(id) %>% summarize(start_date = unique(start_date))

test_that("moultmcmc uz2_recap works", {
  m2r = moultmcmc(moult_column = "pfmg_sampled",
                  date_column = "date_sampled",
                  id_column = "id",
                  data = recaptures,
                  type = 2,
                  log_lik = FALSE,
                  chains = 2,
                  cores = 2,
                  control = list(adapt_delta = 0.99, max_treedepth = 11),
                  iter = 2000)
expect_s3_class(m2r, "moultmcmc")
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
                            chains = 1,
                            iter = 200)
  expect_s3_class(uz2ram, "moultmcmc")
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

