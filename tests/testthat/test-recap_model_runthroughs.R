data(recaptures)

test_that("uz2_recap works", {
  uz2r = uz2_linpred_recap(moult_index_column = "pfmg_sampled",
                  date_column = "date_sampled",
                  id_column = "id",
                  data = recaptures,
                  log_lik = FALSE,
                  chains = 1,
                  iter = 200)
expect_s3_class(uz2r, "moultmcmc")
})
test_that("uz2l_recap works", {
  uz2rl = uz2_linpred_recap(moult_index_column = "pfmg_sampled",
                          date_column = "date_sampled",
                          id_column = "id",
                          lump_non_moult = TRUE,
                          data = recaptures,
                          log_lik = FALSE,
                          chains = 1,
                          iter = 200)
  expect_s3_class(uz2rl, "moultmcmc")
})
test_that("uz3 works", {
  uz3 = uz3_linpred_recap(moult_index_column = "pfmg_sampled",
                          date_column = "date_sampled",
                    id_column = "id",
                    data = subset(recaptures, pfmg_sampled != 0 & pfmg_sampled != 1),
                    log_lik = FALSE,
                    chains = 1,
                    iter = 200)
  expect_s3_class(uz3, "moultmcmc")
})
test_that("uz5 works", {
  uz5 = uz5_linpred_recap(moult_index_column = "pfmg_sampled",
                          date_column = "date_sampled",
                          id_column = "id",
                          data = subset(recaptures, pfmg_sampled != 1),
                    log_lik = FALSE,
                    chains = 1,
                    iter = 200)
  expect_s3_class(uz5, "moultmcmc")
})

