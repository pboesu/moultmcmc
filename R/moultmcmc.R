#' Bayesian inference for Underhill-Zucchini moult models and expansions
#'
#' @export
#' @param date_column the name the column in `data` containing sampling dates, encoded as days since an arbitrary reference date, i.e. a numeric vector
#' @param moult_column the name the column in `data` containing moult indices, i.e. a numeric vector of (linearized) moult scores (0 = old plumage,1 = new plumage) or numerical moult codes (1 = old plumage,2 = moulting, 3 = new plumage), depending on the model type.
#' @param id_column (optional) factor identifier. Usually a season-individual combination to encode within-season recaptures, defaults to NULL. When provided moultmcmc will attempt to fit the relevant recaptures model.
#' @param start_formula model formula for start date
#' @param duration_formula model formula for duration
#' @param sigma_formula model formula for start date sigma
#' @param type integer (one of 1,2,3,4,5,12) referring to type of moult data and consequently model to be fitted (see details)
#' @param lump_non_moult logical; should pre- and post-moult observations be treated as indistinguishable? if TRUE and type %in% c(1,2,12), the relevant lumped model will be fitted (see details).
#' @param data Input data frame must contain a numeric column "date" and a column "moult_cat" which is a numeric vector of categorical moult codes (1 = old plumage,2 = moulting,3 = new plumage).
#' @param init Specification of initial values for all or some parameters. Can be the string "auto" for an automatic guess based on the data, or any of the permitted rstan options: the digit 0, the strings "0" or "random", or a function. See the detailed documentation for the init argument in ?rstan::stan.
#' @param flat_prior use uniform prior on start date and duration (TRUE) or vaguely informative truncated normal prior (FALSE). Defaults to TRUE.
#' @param beta_sd use zero-centred normal priors for regression coefficients other than intercepts? If <= 0 the stan default of improper flat priors is used.
#' @param log_lik boolean retain pointwise log-likelihood in output? This enables model assessment and selection via the loo package. Defaults to true, can lead to very large output arrays if sample size is large.
#' @param use_phi_approx logical flag whether to use stan's Phi_approx function to calculate the "old" likelihoods
#' @param active_moult_recaps_only logical flag whether to ignore repeated observations outside the active moult phase
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#'
#' @details
#' type refers to the type of moult data available (see Underhill and Zucchini (1998) and Underhill, Zucchini and Summers (1990)).
#'
#' type = 1
#' sample is representative of entire population (not yet moulting, in moult, and birds which have completed moult). For type 1 data, any value between 0 and 1 (> 0 and < 1) can be used as the moult index for birds in active moult. The value used does not matter, only the fact that they are in moult.
#' type = 2
#' (default) sample is representative of entire population (not yet moulting, in moult, and birds which have completed moult). Moult scores are required.
#'
#' type = 3
#' sample is representative only of birds in moult. Individuals with moult scores 0 or 1 are ignored.
#'
#' type = 4
#' sample is representative only of birds in moult and those that have completed moult. Individuals with moult scores 0 are ignored.
#'
#' type = 5
#' sample is representative only of birds that have not started moult or that are in moult. Individuals with moult scores 1 are ignored.
#'
#' @return An object of class `stanfit` returned by `rstan::sampling`
#' @importFrom nlme asOneFormula

#TODO: implement an input data class which ensures column names and correct encoding for categorical variables
moultmcmc <- function(moult_column,
                        date_column,
                        id_column = NULL,
                        start_formula = ~1,
                        duration_formula = ~1,
                        sigma_formula = ~1,
                        type = 2,
                        lump_non_moult = FALSE,
                        data,
                        init = "auto",
                        flat_prior = TRUE,
                        beta_sd = 0,
                        log_lik = TRUE,
                        use_phi_approx = FALSE,
                        active_moult_recaps_only = FALSE,
                        ...) {
  #check input data are as expected
  if(type %in% c(2:5)) {
    stopifnot(all(data[[moult_column]] >= 0 & data[[moult_column]] <= 1))
  } else {
      if (type == 1) {
        stopifnot(all(data[[moult_column]] %in% c(1,2,3)))
      } else {
        if (type == 12) {
          stopifnot(all(data[[moult_column]] >= 0 & data[[moult_column]] <= 1 | data[[moult_column]]==2))
        } else {
          stop('Unsupported model type. Argument type must be one of 1,2,3,4,5,12')
        }
      }
    }
  stopifnot(is.numeric(data[[date_column]]))
  stopifnot(is.data.frame(data))
  #prepare model.frame and handle missing values
  data <- model.frame(nlme::asOneFormula(start_formula, duration_formula, sigma_formula, as.formula(paste(moult_column, '~ ', date_column))), data = data)
  #order data by moult category
  data <- data[order(data[[moult_column]]),]
  #setup model matrices
  X_mu <- model.matrix(start_formula, data)
  X_tau <- model.matrix(duration_formula, data)
  X_sigma <- model.matrix(sigma_formula, data)
  #prepare data structure for stan
  standata <- list(old_dates = data[[date_column]][data[[moult_index_column]]==0],
                   N_old = length(data[[date_column]][data[[moult_index_column]]==0]),
                   moult_dates  = data[[date_column]][(data[[moult_index_column]] > 0 & data[[moult_index_column]] < 1)],
                   moult_indices = data[[moult_index_column]][(data[[moult_index_column]] > 0 & data[[moult_index_column]] < 1)],
                   N_moult = length(data[[date_column]][(data[[moult_index_column]] > 0 & data[[moult_index_column]] < 1)]),
                   new_dates = data[[date_column]][data[[moult_index_column]]==1],
                   N_new = length(data[[date_column]][data[[moult_index_column]]==1]),
                   X_mu = X_mu,
                   N_pred_mu = ncol(X_mu),
                   X_tau = X_tau,
                   N_pred_tau = ncol(X_tau),
                   X_sigma = X_sigma,
                   N_pred_sigma = ncol(X_sigma),
                   lumped = as.numeric(lump_non_moult),
                   llik = as.numeric(log_lik))
  #include pointwise log_lik matrix  in output?
  if(log_lik){
    outpars <- c('beta_mu','beta_tau','beta_sigma', 'sigma_intercept', 'log_lik')
  } else {
    outpars <- c('beta_mu','beta_tau','beta_sigma', 'sigma_intercept')
  }
  #guess initial values
  if(init == "auto"){
    mu_start = mean(c(max(standata$old_dates),min(standata$moult_dates)))
    tau_start = max(10, mean(c(max(standata$moult_dates),min(standata$new_dates))) - mu_start)
    sigma_start = min(40,sd(standata$moult_dates))
    initfunc <- function(chain_id = 1) {
      # cat("chain_id =", chain_id, "\n")
      list(beta_mu = as.array(c(mu_start,rep(0, standata$N_pred_mu - 1))), #initialize intercept term from data, set inits for all other effects to 0
           beta_tau = as.array(c(tau_start, rep(0, standata$N_pred_tau - 1))),
           beta_sigma = as.array(c(log(sigma_start), rep(0, standata$N_pred_sigma - 1))))#NB this is on log link scale
    }
    out <- rstan::sampling(stanmodels$uz2_linpred, data = standata, init = initfunc, pars = outpars, ...)
  } else {
    out <- rstan::sampling(stanmodels$uz2_linpred, data = standata, init = init, pars = outpars, ...)
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
  out_struc$terms$start_formula <- start_formula
  out_struc$terms$duration_formula <- duration_formula
  out_struc$terms$sigma_formula <- sigma_formula
  out_struc$data <- data
  out_struc$na.action <- attr(data, "na.action")
  class(out_struc) <- 'moultmcmc'
  return(out_struc)
}
