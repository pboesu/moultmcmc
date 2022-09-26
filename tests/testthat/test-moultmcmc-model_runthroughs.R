sanderlings <- moult::sanderlings
sanderlings$MCat[sanderlings$MIndex == 0] <- 1
sanderlings$MCat[sanderlings$MIndex == 1] <- 3
sanderlings$MCat[is.na(sanderlings$MCat)] <- 2

test_that("uz1 works", {
  uz1 = moultmcmc(moult_column = "MCat",
                  date_column = "Day",
                  data = sanderlings,
                  type = 1,
                  log_lik = FALSE,
                  chains = 1)
expect_s3_class(uz1, "moultmcmc")
})
test_that("uz1l works", {
  uz1l = moultmcmc(moult_column = "MCat",
                    date_column = "Day",
                    lump_non_moult = TRUE,
                    type = 1,
                    data = sanderlings,
                    log_lik = FALSE,
                    chains = 1)
  expect_s3_class(uz1l, "moultmcmc")
})
test_that("uz2 works", {
  uz2 = moultmcmc(moult_column = "MIndex",
                    date_column = "Day",
                    type = 2,
                    data = sanderlings,
                    log_lik = FALSE,
                    chains = 1)
  expect_s3_class(uz2, "moultmcmc")
})
test_that("uz2l works", {
  uz2l = moultmcmc(moult_column = "MIndex",
                   date_column = "Day",
                   type = 2,
                   lump_non_moult = TRUE,
                   data = sanderlings,
                   log_lik = FALSE,
                   chains = 1)
  expect_s3_class(uz2l, "moultmcmc")
})
test_that("uz3 works", {
  uz3 = moultmcmc(moult_column = "MIndex",
                  date_column = "Day",
                  type = 3,
                  data = sanderlings,
                  log_lik = FALSE,
                  chains = 1)
  expect_s3_class(uz3, "moultmcmc")
})
test_that("uz4 works", {
  uz4 = moultmcmc(moult_column = "MIndex",
                  date_column = "Day",
                  data = sanderlings,
                  type = 4,
                  log_lik = FALSE,
                  chains = 1)
  expect_s3_class(uz4, "moultmcmc")
})
test_that("uz5 works", {
  uz5 = moultmcmc(moult_column = "MIndex",
                  date_column = "Day",
                  data = sanderlings,
                  type = 5,
                  log_lik = FALSE,
                  chains = 1)
  expect_s3_class(uz5, "moultmcmc")
})
test_that("uz12 works", {
  uz12 = moultmcmc(moult_column = "MCat",
                   date_column = "Day",
                   data = sanderlings,
                   type = 12,
                   log_lik = FALSE,
                   chains = 1)
  expect_s3_class(uz12, "moultmcmc")
})

