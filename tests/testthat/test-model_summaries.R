sanderlings <- moult::sanderlings
sanderlings$MCat[sanderlings$MIndex == 0] <- 1
sanderlings$MCat[sanderlings$MIndex == 1] <- 3
sanderlings$MCat[is.na(sanderlings$MCat)] <- 2
library(moult)
test_that("lci is smaller than uci - moult", {
m1 =  moult(MIndex ~ Day,data = sanderlings, type = 2)
testthat::expect_gte(summary_table(m1)$uci[1],summary_table(m1)$lci[1])
})
