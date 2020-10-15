#' Bayesian inference for the Type 2 moult model
#'
#' @export
#' @param start_formula model formula for start date
#' @param duration_formula model formula for duration
#' @param sigma_formula model formula for start date sigma
#' @param data Input data frame must contain a numeric column "date" and a column "moult_cat" which is a numeric vector of categorical moult codes (1 = old plumage,2 = moulting,3 = new plumage).
#' @param init Specification of initial values for all or some parameters. Can be the string "auto" for an automatic guess based on the data, or any of the permitted rstan options: the digit 0, the strings "0" or "random", or a function. See the detailed documentation for the init argument in ?rstan::stan.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
#TODO: implement an input data class which ensures column names and correct encoding for categorical variables
uz2_linpred <- function(start_formula = ~1, duration_formula = ~1, sigma_formula = ~1, data, init = "auto",...) {
  stopifnot(all(data$moult_indices >= 0 & data$moult_indices <= 1))

  #order data by moult category
  data <- data[order(data$moult_indices),]
  #setup model matrices
  X_mu <- model.matrix(start_formula, data)
  X_tau <- model.matrix(duration_formula, data)
  X_sigma <- model.matrix(sigma_formula, data)
  #prepare data structure for stan
  standata <- list(old_dates = dates[moult_indices==0],
                   N_old = length(dates[moult_indices==0]),
                   moult_dates  = dates[(moult_indices > 0 & moult_indices < 1)],
                   moult_indices = moult_indices[(moult_indices > 0 & moult_indices < 1)],
                   N_moult = length(dates[(moult_indices > 0 & moult_indices < 1)]),
                   new_dates = dates[moult_indices==1],
                   N_new = length(dates[moult_indices==1]),
                   X_mu = X_mu,
                   N_pred_mu = ncol(X_mu),
                   X_tau = X_tau,
                   N_pred_tau = ncol(X_tau),
                   X_sigma = X_sigma,
                   N_pred_sigma = ncol(X_sigma))
  #guess initial values
  if(init == "auto"){
    mu_start = mean(max(standata$old_dates),min(standata$moult_dates))
    tau_start = mean(max(standata$moult_dates),min(standata$new_dates)) - mu_start
    sigma_start = min(10,sd(standata$moult_dates))
    initfunc <- function(chain_id = 1) {
      # cat("chain_id =", chain_id, "\n")
      list(beta_mu = as.array(c(mu_start,rep(0, standata$N_pred_mu - 1))),
           beta_tau = as.array(c(tau_start, rep(0, standata$N_pred_tau - 1))),
           sigma = sigma_start)
    }
    out <- rstan::sampling(stanmodels$uz2_linpred, data = standata, init = initfunc, ...)
  } else {
    out <- rstan::sampling(stanmodels$uz2_linpred, data = standata, init = init, ...)
  }

  return(out)
}
