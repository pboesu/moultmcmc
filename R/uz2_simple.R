#' Bayesian inference for the Type 2 moult model
#'
#' @export
#' @param dates Numeric vector of input capture dates.
#' @param moult_indices Numeric vector of numerical index of moult progression in [0,1].
#' @param init Specification of initial values for all or some parameters. Can be the string "auto" for an automatic guess based on the data, or any of the permitted rstan options: the digit 0, the strings "0" or "random", or a function. See the detailed documentation for the init argument in ?rstan::stan.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
uz2_simple <- function(dates, moult_indices, init = "auto",...) {
  stopifnot(all(moult_indices >= 0 & moult_indices <= 1))
  #order data by moult index
  data <- data[order(data$moult_indices),]
  #prepare data structure for stan
  standata <- list(old_dates = dates[moult_indices==0],
                   N_old = length(dates[moult_indices==0]),
                   moult_dates  = dates[(moult_indices > 0 & moult_indices < 1)],
                   moult_indices = moult_indices[(moult_indices > 0 & moult_indices < 1)],
                   N_moult = length(dates[(moult_indices > 0 & moult_indices < 1)]),
                   new_dates = dates[moult_indices==1],
                   N_new = length(dates[moult_indices==1]))
  #guess initial values
  if(init == "auto"){
    mu_start = mean(max(standata$old_dates),min(standata$moult_dates))
    tau_start = mean(max(standata$moult_dates),min(standata$new_dates)) - mu_start
    sigma_start = min(10,sd(standata$moult_dates))
    initfunc <- function(chain_id = 1) {
      # cat("chain_id =", chain_id, "\n")
      list(mu = mu_start, tau = tau_start, sigma = sigma_start)
    }
    out <- rstan::sampling(stanmodels$uz2_simple, data = standata, init = initfunc, ...)
  } else {
    out <- rstan::sampling(stanmodels$uz2_simple, data = standata, init = init, ...)
  }

  return(out)
}
