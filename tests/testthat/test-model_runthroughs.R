sanderlings <- moult::sanderlings
sanderlings$MCat[sanderlings$MIndex == 0] <- 1
sanderlings$MCat[sanderlings$MIndex == 1] <- 3
sanderlings$MCat[is.na(sanderlings$MCat)] <- 2

test_that("uz1 works", {
  uz1 = uz1_linpred(moult_cat_column = "MCat",
                  date_column = "Day",
                  data = sanderlings,
                  log_lik = FALSE,
                  chains = 1)
expect_s3_class(uz1, "moultmcmc")
})
test_that("uz2 works", {
  uz2 = uz2_linpred(moult_index_column = "MIndex",
                    date_column = "Day",
                    data = sanderlings,
                    log_lik = FALSE,
                    chains = 1)
  expect_s3_class(uz2, "moultmcmc")
})
test_that("uz2l standalone works", {
  uz2ls = uz2_lumped(moult_index_column = "MIndex",
                    date_column = "Day",
                    data = sanderlings,
                    log_lik = FALSE,
                    chains = 1)
  expect_s3_class(uz2ls, "moultmcmc")
})
test_that("uz2l works", {
  uz2l = uz2_linpred(moult_index_column = "MIndex",
                   date_column = "Day",
                   lump_non_moult = TRUE,
                   data = sanderlings,
                   log_lik = FALSE,
                   chains = 1)
  expect_s3_class(uz2l, "moultmcmc")
})
test_that("uz3 works", {
  uz3 = uz3_linpred(moult_index_column = "MIndex",
                    date_column = "Day",
                    data = sanderlings,
                    log_lik = FALSE,
                    chains = 1)
  expect_s3_class(uz3, "moultmcmc")
})
test_that("uz4 works", {
  uz4 = uz4_linpred(moult_index_column = "MIndex",
                    date_column = "Day",
                    data = sanderlings,
                    log_lik = FALSE,
                    chains = 1)
  expect_s3_class(uz4, "moultmcmc")
})
test_that("uz5 works", {
  uz5 = uz5_linpred(moult_index_column = "MIndex",
                    date_column = "Day",
                    data = sanderlings,
                    log_lik = FALSE,
                    chains = 1)
  expect_s3_class(uz5, "moultmcmc")
})
test_that("uz12 works", {
  uz12 = uz12_linpred(moult_cat_column = "MCat",
                      moult_index_column = "MIndex",
                    date_column = "Day",
                    data = sanderlings,
                    log_lik = FALSE,
                    chains = 1)
  expect_s3_class(uz12, "moultmcmc")
})

