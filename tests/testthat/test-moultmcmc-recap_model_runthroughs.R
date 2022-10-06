data(recaptures)

test_that("uz2_recap works", {
  uz2r = moultmcmc(moult_column = "pfmg_sampled",
                  date_column = "date_sampled",
                  id_column = "id",
                  data = recaptures,
                  log_lik = FALSE,
                  chains = 1,
                  iter = 200)
expect_s3_class(uz2r, "moultmcmc")
})
test_that("uz2l_recap works", {
  uz2rl = moultmcmc(moult_column = "pfmg_sampled",
                          date_column = "date_sampled",
                          id_column = "id",
                          lump_non_moult = TRUE,
                          data = recaptures,
                          log_lik = FALSE,
                          chains = 1,
                          iter = 200)
  expect_s3_class(uz2rl, "moultmcmc")
})
test_that("uz2l_active_moult_only", {
  uz2ram = moultmcmc(moult_column = "pfmg_sampled",
                            date_column = "date_sampled",
                            id_column = "id",
                            lump_non_moult = FALSE,
                            active_moult_recaps_only = TRUE,
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
                    log_lik = FALSE,
                    chains = 1,
                    iter = 200)
  expect_s3_class(uz3, "moultmcmc")
})
test_that("uz5r works", {
  uz5 = moultmcmc(moult_column = "pfmg_sampled",
                          date_column = "date_sampled",
                          id_column = "id",
                          data = subset(recaptures, pfmg_sampled != 1),
                    log_lik = FALSE,
                    chains = 1,
                    iter = 200)
  expect_s3_class(uz5, "moultmcmc")
})

