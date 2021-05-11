#' Bayesian inference for the Type 1 moult model
#'
#' @export
#' @param start_formula model formula for start date
#' @param duration_formula model formula for duration
#' @param sigma_formula model formula for start date sigma
#' @param data Input data frame
#' @param date_column the name the column in `data` containing sampling dates, encoded as days since an arbitrary reference date, i.e. a numeric vector
#' @param moult_cat_column the name the column in `data` containing moult categories, i.e. a numeric vector of categorical moult codes (1 = old plumage,2 = moulting,3 = new plumage)
#' @param init Specification of initial values for all or some parameters. Can be the string "auto" for an automatic guess based on the data, or any of the permitted rstan options: the digit 0, the strings "0" or "random", or a function. See the detailed documentation for the init argument in ?rstan::stan.
#' @param log_lik boolean retain pointwise log-likelihood in output? This enables model assessment and selection via the loo package. Defaults to true, can lead to very large output arrays if sample size is large.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
#TODO: implement an input data class which ensures column names and correct encoding for categorical variables
uz1_linpred <- function(moult_cat_column, date_column, start_formula = ~1, duration_formula = ~1, sigma_formula = ~1, data, init = "auto", log_lik = TRUE, ...) {
  stopifnot(all(data[[moult_cat_column]] %in% c(1,2,3)))
  stopifnot(is.numeric(data[[date_column]]))
  stopifnot(is.data.frame(data))
  #order data by moult category
  data <- data[order(data[[moult_cat_column]]),]
  #setup model matrices
  X_mu <- model.matrix(start_formula, data)
  X_tau <- model.matrix(duration_formula, data)
  X_sigma <- model.matrix(sigma_formula, data)
  #prepare data structure for stan
  standata <- list(old_dates = data[[date_column]][data[[moult_cat_column]]==1],
                   N_old = length(data[[date_column]][data[[moult_cat_column]]==1]),
                   moult_dates  = data[[date_column]][data[[moult_cat_column]]==2],
                   N_moult = length(data[[date_column]][data[[moult_cat_column]]==2]),
                   new_dates = data[[date_column]][data[[moult_cat_column]]==3],
                   N_new = length(data[[date_column]][data[[moult_cat_column]]==3]),
                   X_mu = X_mu,
                   N_pred_mu = ncol(X_mu),
                   X_tau = X_tau,
                   N_pred_tau = ncol(X_tau),
                   X_sigma = X_sigma,
                   N_pred_sigma = ncol(X_sigma))
  #include pointwise log_lik matrix  in output?
  if(log_lik){
    outpars <- c('beta_mu','beta_tau','beta_sigma', 'sigma_intercept', 'log_lik')
  } else {
    outpars <- c('beta_mu','beta_tau','beta_sigma', 'sigma_intercept')
  }
  #guess initial values
  if(init == "auto"){
    mu_start = mean(c(max(standata$old_dates),min(standata$moult_dates)))
    tau_start = max(10,mean(c(max(standata$moult_dates),min(standata$new_dates))) - mu_start)
    sigma_start = min(10,sd(standata$moult_dates))
    initfunc <- function(chain_id = 1) {
      # cat("chain_id =", chain_id, "\n")
      list(beta_mu = as.array(c(mu_start,rep(0, standata$N_pred_mu - 1))), #initialize intercept term from data, set inits for all other effects to 0
           beta_tau = as.array(c(tau_start, rep(0, standata$N_pred_tau - 1))),
           beta_sigma = as.array(c(log(sigma_start)), rep(0, standata$N_pred_tau - 1)))#NB this is on log link scale
    }

    out <- rstan::sampling(stanmodels$uz1_linpred, data = standata, init = initfunc, pars = outpars, ...)
  } else {
    out <- rstan::sampling(stanmodels$uz1_linpred, data = standata, init = init, pars = outpars, ...)
  }
  #rename regression coefficients for output
  names(out)[grep('beta_mu', names(out))] <- paste('mean',colnames(X_mu), sep = '_')
  names(out)[grep('beta_tau', names(out))] <- paste('duration',colnames(X_tau), sep = '_')
  names(out)[grep('beta_sigma', names(out))] <- paste('log_sd',colnames(X_sigma), sep = '_')
  names(out)[grep('sigma_intercept', names(out))] <- 'sd_(Intercept)'
  return(out)
}
