library(moult)
sanderlings <- moult::sanderlings

#generate random covars to
sanderlings$noise <- rnorm(nrow(sanderlings))

m1 =  moult(MIndex ~ Day|noise|noise,data = sanderlings, type = 2)

test_that("lci is smaller than uci - moult", {
testthat::expect_gte(summary_table(m1)$uci[1],summary_table(m1)$lci[1])
})

test_that("moult_plot.moult runs through", {
  testthat::expect_error(moult_plot(m1), NA)
})


test_that("residual_plot.moult runs through", {
  testthat::expect_error(residual_plot(m1), NA)
})
