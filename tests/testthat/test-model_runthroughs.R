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
