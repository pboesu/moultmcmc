#' Bayesian inference for the Type 5 moult model with repeat measures
#'
#' @export
#' @param date_column the name the column in `data` containing sampling dates, encoded as days since an arbitrary reference date, i.e. a numeric vector
#' @param moult_index_column the name the column in `data` containing moult indices, i.e. a numeric vector of (linearized) moult scores (0 = old plumage,1 = new plumage).
#' @param id_column factor identifier. Usually a season-individual combination to encode within-season recaptures
#' @param start_formula model formula for start date
#' @param duration_formula model formula for duration
#' @param sigma_formula model formula for start date sigma
#' @param data Input data frame
#' @param init Specification of initial values for all or some parameters. Can be the string "auto" for an automatic guess based on the data, or any of the permitted rstan options: the digit 0, the strings "0" or "random", or a function. See the detailed documentation for the init argument in ?rstan::stan.
#' @param log_lik boolean retain pointwise log-likelihood in output? This enables model assessment and selection via the loo package. Defaults to true, can lead to very large output arrays if sample size is large.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
#TODO: implement an input data class which ensures column names and correct encoding for categorical variables
uz5_linpred_recap <- function(moult_index_column, date_column, id_column, start_formula = ~1, duration_formula = ~1, sigma_formula = ~1, data, init = "auto", log_lik = TRUE,...) {
  stopifnot(all(data[[moult_index_column]] >= 0 & data[[moult_index_column]] <= 1))
  #TODO: Assess the relative amount of old vs moult data, and especially fail/warn when there is no data of either category
  stopifnot(is.numeric(data[[date_column]]))
  stopifnot(is.factor(data[[id_column]]))
  stopifnot(is.data.frame(data))
  #remove new data by moult category
  data <- droplevels(subset(data, get(moult_index_column) != 1))#TODO:record / notice of how many obs removed missing
  #TODO: remove consecutive old captures within a season ?
  #order data by moult category
  data <- data[order(data[[moult_index_column]]),]
  #setup model matrices
  #TODO: check for consistency of covariates among recaptures?!
  #TODO: do these need to be reduced to N_ind?? or does there need to be checking that they are identical for recaptures?
  X_mu <- model.matrix(start_formula, data)#TODO: na.action missing
  X_tau <- model.matrix(duration_formula, data)
  X_sigma <- model.matrix(sigma_formula, data)
  #create vector of first occurrence of each individual in the model matrix to index unique linear predictor values for each individual
  id_first <- match(unique(data[[id_column]]), data[[id_column]])
  #create vectors of indices of replicated/non_replicated observations
  N_moult = length(data[[date_column]][(data[[moult_index_column]] > 0 & data[[moult_index_column]] < 1)])
  N_old = length(data[[date_column]][data[[moult_index_column]]==0])#TODO: this is not robust to NA's being thrown out by model.matrix/model.frame
  replicated <- which(data[[id_column]] %in% names(table(data[[id_column]])[table(data[[id_column]])>1]))
  not_replicated <- which(data[[id_column]] %in% names(table(data[[id_column]])[table(data[[id_column]])==1]))
  is_replicated <- ifelse(table(data[[id_column]])>1, 1,0)
  #prepare data structure for stan
  #TODO: all of these need to be as.array if they are not generally a scalar, otherwise stan indexing falls over for length-one vectors
  standata <- list(moult_dates  = as.array(data[[date_column]][(data[[moult_index_column]] > 0 & data[[moult_index_column]] < 1)]),
                   moult_indices = as.array(data[[moult_index_column]][(data[[moult_index_column]] > 0 & data[[moult_index_column]] < 1)]),
                   N_moult = N_moult,
                   old_dates = as.array(data[[date_column]][data[[moult_index_column]]==0]),
                   N_old = N_old,#TODO: this is not robust to NA's being thrown out by model.matrix/model.frame
                   N_ind = length(unique(data[[id_column]])),

                   individual = as.numeric(data[[id_column]]),#TODO: the resultant individual intercepts are hard to map onto original factor levels - this should be handled in the postprocessing of the model output
                   individual_first_index = as.array(id_first),
                   replicated = as.array(replicated),
                   not_replicated = as.array(not_replicated),
                   not_replicated_old = as.array(not_replicated[not_replicated <= N_old]),
                   not_replicated_moult = as.array(not_replicated[not_replicated > N_old] - N_old),
                   is_replicated = as.array(is_replicated),
                   Nobs_replicated = length(replicated),
                   Nobs_not_replicated_old = length(not_replicated[not_replicated <= N_old]),
                   Nobs_not_replicated_moult = length(not_replicated[not_replicated > N_old]),
                   X_mu = X_mu,
                   N_pred_mu = ncol(X_mu),
                   X_tau = X_tau,
                   N_pred_tau = ncol(X_tau),
                   X_sigma = X_sigma,
                   N_pred_sigma = ncol(X_sigma))
  #include pointwise log_lik matrix  in output?
  if(log_lik){
    outpars <- c('beta_mu','beta_tau','beta_sigma', 'sigma_intercept', 'sigma_mu_ind','beta_star','finite_sd', 'mu_ind_star', 'mu_ind', 'log_lik')
  } else {
    outpars <- c('beta_mu','beta_tau','beta_sigma', 'sigma_intercept', 'sigma_mu_ind','beta_star','finite_sd', 'mu_ind_star', 'mu_ind')
  }
  #guess initial values
  if(init == "auto"){
    mu_start = mean(c(min(standata$moult_dates),max(standata$old_dates)))
    tau_start = max(10, max(standata$moult_dates)-max(standata$old_dates))
    sigma_start = min(10,sd(standata$moult_dates))
    initfunc <- function(chain_id = 1) {
      # cat("chain_id =", chain_id, "\n")
      list(beta_mu = as.array(c(mu_start,rep(0, standata$N_pred_mu - 1))), #initialize intercept term from data, set inits for all other effects to 0
           beta_tau = as.array(c(tau_start, rep(0, standata$N_pred_tau - 1))),
           beta_sigma = as.array(c(log(sigma_start), rep(0, standata$N_pred_sigma - 1))),#NB this is on log link scale
           mu_ind = as.array(rep(0, standata$N_ind)))
    }
    out <- rstan::sampling(stanmodels$uz5_recap, data = standata, init = initfunc, pars = outpars, ...)
  } else {
    out <- rstan::sampling(stanmodels$uz5_recap, data = standata, init = init, pars = outpars, ...)
  }
  #rename regression coefficients for output
  names(out)[grep('beta_mu', names(out))] <- paste('mean',colnames(X_mu), sep = '_')
  names(out)[grep('beta_tau', names(out))] <- paste('duration',colnames(X_tau), sep = '_')
  names(out)[grep('beta_sigma', names(out))] <- paste('log_sd',colnames(X_sigma), sep = '_')
  names(out)[grep('sigma_intercept', names(out))] <- 'sd_(Intercept)'
  out_struc <- list()
  out_struc$stanfit <- out
  out_struc$terms$date_column <- date_column
  out_struc$terms$moult_index_column <- moult_index_column
  out_struc$terms$moult_cat_column <- NA
  out_struc$terms$id_column <- id_column
  class(out_struc) <- 'moultmcmc'
  return(out_struc)
}
