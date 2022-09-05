#' Bayesian inference for the combined Type 1 and Type 2 moult model
#'
#' @export
#' @param date_column the name the column in `data` containing sampling dates, encoded as days since an arbitrary reference date, i.e. a numeric vector
#' @param moult_index_column the name the column in `data` containing moult indices, i.e. a numeric vector of (linearized) moult scores (0 = old plumage,1 = new plumage).
#' @param moult_cat_column the name the column in `data` containing moult categories, i.e. a numeric vector of categorical moult codes (1 = old plumage,2 = moulting,3 = new plumage)
#' @param start_formula model formula for start date
#' @param duration_formula model formula for duration
#' @param sigma_formula model formula for start date sigma
#' @param data Input data frame must contain a numeric column "date" and a column "moult_cat" which is a numeric vector of categorical moult codes (1 = old plumage,2 = moulting,3 = new plumage).
#' @param init Specification of initial values for all or some parameters. Can be the string "auto" for an automatic guess based on the data, or any of the permitted rstan options: the digit 0, the strings "0" or "random", or a function. See the detailed documentation for the init argument in ?rstan::stan.
#' @param log_lik boolean retain pointwise log-likelihood in output? This enables model assessment and selection via the loo package. Defaults to true, can lead to very large output arrays if sample size is large.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
#' @importFrom stats sd model.matrix
#'
#TODO: implement an input data class which ensures column names and correct encoding for categorical variables
uz12_linpred <- function(moult_index_column, moult_cat_column, date_column, start_formula = ~1, duration_formula = ~1, sigma_formula = ~1, data, init = "auto", log_lik = TRUE,...) {
  stopifnot(all((data[[moult_index_column]] >= 0 & data[[moult_index_column]] <= 1) | is.na(data[[moult_index_column]])))#need to handle NAs
  stopifnot(all(data[[moult_cat_column]] %in% c(1,2,3,NA)))
  stopifnot(is.numeric(data[[date_column]]))
  stopifnot(is.data.frame(data))
  #fill in NA scores for O and N birds
  data[[moult_index_column]][is.na(data[[moult_index_column]]) & data[[moult_cat_column]] == 1] <- 0
  data[[moult_index_column]][is.na(data[[moult_index_column]]) & data[[moult_cat_column]] == 3] <- 1
  #create subset of categorical observations of active moult
  moult_cat_data <- subset(data, is.na(get(moult_index_column)) & get(moult_cat_column) == 2)
  #order data by moult category
  data <- data[order(data[[moult_index_column]]),]
  #split up scoren data by state, then reassemble with the categorical data to get the model matrix into the right shape for stan
  #TODO: test that this works if any of the components is of zero length!
  #TODO: this is a really memory expensive way of doing this once data frames get large!
  old_data <- subset(data, get(moult_index_column) == 0)
  moult_score_data <- subset(data, get(moult_index_column) != 0 & get(moult_index_column) != 1)
  new_data <- subset(data, get(moult_index_column) == 1)
  data <- bind_rows(old_data, moult_score_data, moult_cat_data, new_data)
  #setup model matrices
  X_mu <- model.matrix(start_formula, data)
  X_tau <- model.matrix(duration_formula, data)
  X_sigma <- model.matrix(sigma_formula, data)
  #prepare data structure for stan
  standata <- list(old_dates = old_data[[date_column]][old_data[[moult_index_column]]==0],
                   N_old = nrow(old_data),
                   moult_dates  = moult_score_data[[date_column]][(moult_score_data[[moult_index_column]] > 0 & moult_score_data[[moult_index_column]] < 1)],
                   moult_indices = moult_score_data[[moult_index_column]][(moult_score_data[[moult_index_column]] > 0 & moult_score_data[[moult_index_column]] < 1)],
                   N_moult = nrow(moult_score_data),
                   N_moult_cat = nrow(moult_cat_data),
                   moult_cat_dates = moult_cat_data[[date_column]],
                   new_dates = new_data[[date_column]][new_data[[moult_index_column]]==1],
                   N_new = nrow(new_data),
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
    mu_start = mean(c(max(standata$old_dates),min(c(standata$moult_dates,standata$moult_cat_dates))))
    tau_start = max(10, mean(c(max(c(standata$moult_dates,standata$moult_cat_dates)),min(standata$new_dates))) - mu_start)
    sigma_start = min(10,sd(c(standata$moult_dates,standata$moult_cat_dates)))
    initfunc <- function(chain_id = 1) {
      # cat("chain_id =", chain_id, "\n")
      list(beta_mu = as.array(c(mu_start,rep(0, standata$N_pred_mu - 1))), #initialize intercept term from data, set inits for all other effects to 0
           beta_tau = as.array(c(tau_start, rep(0, standata$N_pred_tau - 1))),
           beta_sigma = as.array(c(log(sigma_start), rep(0, standata$N_pred_sigma - 1))))#NB this is on log link scale
    }
    out <- rstan::sampling(stanmodels$uz12_linpred, data = standata, init = initfunc, pars = outpars, ...)
  } else {
    out <- rstan::sampling(stanmodels$uz12_linpred, data = standata, init = init, pars = outpars, ...)
  }
  #TODO: implement exception handling so function errors gracefully when sampling fails
  #rename regression coefficients for output
  names(out)[grep('beta_mu', names(out))] <- paste('mean',colnames(X_mu), sep = '_')
  names(out)[grep('beta_tau', names(out))] <- paste('duration',colnames(X_tau), sep = '_')
  names(out)[grep('beta_sigma', names(out))] <- paste('log_sd',colnames(X_sigma), sep = '_')
  names(out)[grep('sigma_intercept', names(out))] <- 'sd_(Intercept)'
  out_struc <- list()
  out_struc$stanfit <- out
  out_struc$terms$start_formula <- start_formula
  out_struc$terms$duration_formula <- duration_formula
  out_struc$terms$sigma_formula <- sigma_formula
  out_struc$data <- data
  class(out_struc) <- 'moultmcmc'
  return(out_struc)
}
