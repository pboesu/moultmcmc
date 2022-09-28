sanderlings <- moult::sanderlings
sanderlings$MCat[sanderlings$MIndex == 0] <- 1
sanderlings$MCat[sanderlings$MIndex == 1] <- 3
sanderlings$MCat[is.na(sanderlings$MCat)] <- 2
#introduce a couple of NAs
sanderlings$MCat[12] <- NA
sanderlings$MIndex[15] <- NA
sanderlings$Day[105] <- NA

test_that("uz1 works with 1,2,3 inputs", {
  uz1 = moultmcmc(moult_column = "MCat",
                  date_column = "Day",
                  data = sanderlings,
                  type = 1,
                  log_lik = FALSE,
                  chains = 1)
expect_s3_class(uz1, "moultmcmc")
expect_error(moultmcmc(moult_column = "MCat",
                date_column = "Day",
                data = sanderlings,
                type = 14,
                log_lik = FALSE,
                chains = 1))

})
test_that("uz1 works with [0,1] inputs", {
  uz1sc = moultmcmc(moult_column = "MIndex",
                  date_column = "Day",
                  data = sanderlings,
                  type = 1,
                  log_lik = FALSE,
                  chains = 1)
  expect_s3_class(uz1sc, "moultmcmc")
})
test_that("uz1l works with 1,2,3 inputs", {
  uz1l = moultmcmc(moult_column = "MCat",
                    date_column = "Day",
                    lump_non_moult = TRUE,
                    type = 1,
                    data = sanderlings,
                    log_lik = FALSE,
                    chains = 1)
  expect_s3_class(uz1l, "moultmcmc")
})
test_that("uz1l works with [0,1] inputs", {
  uz1lsc = moultmcmc(moult_column = "MCat",
                   date_column = "Day",
                   lump_non_moult = TRUE,
                   type = 1,
                   data = sanderlings,
                   log_lik = FALSE,
                   chains = 1)
  expect_s3_class(uz1lsc, "moultmcmc")
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
test_that("uz12 data prep works", {
  sanderlings$MIndexCat <- consolidate_moult_records(sanderlings$MIndex, sanderlings$MCat)
  expect_true(all((sanderlings$MIndexCat %in% c(2,NA)) | (sanderlings$MIndexCat >= 0 & sanderlings$MIndexCat <= 1)))
})

test_that("uz12 works", {
  sanderlings$MIndex[100:110] <- NA
  sanderlings$MIndexCat <- consolidate_moult_records(sanderlings$MIndex, sanderlings$MCat)
  uz12 = moultmcmc(moult_column = "MIndexCat",
                   date_column = "Day",
                   data = sanderlings,
                   type = 12,
                   log_lik = FALSE,
                   chains = 1)
  expect_s3_class(uz12, "moultmcmc")
})

