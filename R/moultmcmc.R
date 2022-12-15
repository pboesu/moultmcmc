#' Bayesian inference for Underhill-Zucchini moult models and expansions
#'
#' @export
#' @param moult_column the name the column in `data` containing moult indices, i.e. a numeric vector of (linearized) moult scores in \[0,1\] (0 = old plumage, 1 = new plumage; for model types 1-5),   numerical moult codes (1 = old plumage, 2 = moulting, 3 = new plumage; for model type 1), or a mixed column created by \code{\link{consolidate_moult_records}} for model type 12.
#' @param date_column the name the column in `data` containing sampling dates, encoded as days since an arbitrary reference date, i.e. a numeric vector
#' @param id_column (optional) factor identifier. Usually a season-individual combination to encode within-season recaptures, defaults to NULL. When provided moultmcmc will attempt to fit the relevant recaptures model.
#' @param start_formula model formula for start date
#' @param duration_formula model formula for duration
#' @param sigma_formula model formula for start date sigma
#' @param type integer (one of 1,2,3,4,5,12) referring to type of moult data and consequently model to be fitted (see details)
#' @param lump_non_moult logical; should pre- and post-moult observations be treated as indistinguishable? if TRUE and type %in% c(1,2,12), the relevant lumped model will be fitted (see details).
#' @param data Input data frame must contain a numeric column "date" and a column "moult_cat" which is a numeric vector of categorical moult codes (1 = old plumage,2 = moulting,3 = new plumage).
#' @param init Specification of initial values for all or some parameters. Can be the string "auto" for an automatic guess based on the data, or any of the permitted \code{rstan} options: the digit 0, the strings "0" or "random", or a function. See the detailed documentation for the init argument in \code{\link[rstan:stan]{rstan::stan}}.
#' @param flat_prior use uniform prior on start date and duration (TRUE) or vaguely informative truncated normal prior (FALSE). Defaults to TRUE.
#' @param beta_sd use zero-centred normal priors for regression coefficients other than intercepts? If <= 0 the stan default of improper flat priors is used.
#' @param log_lik boolean retain pointwise log-likelihood in output? This enables model assessment and selection via the loo package. Defaults to FALSE, can lead to very large output arrays when sample size is large.
#' @param use_phi_approx logical flag whether to use stan's Phi_approx function to calculate the "old" likelihoods
#' @param active_moult_recaps_only logical flag whether to ignore repeated observations outside the active moult phase
#' @param same_sigma logical flag, currently unused
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
#' @importFrom stats lm model.frame as.formula

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
                        log_lik = FALSE,
                        use_phi_approx = FALSE,
                        active_moult_recaps_only = FALSE,
                        same_sigma = FALSE,
                        ...) {
  #check input data are as expected
  stopifnot(is.data.frame(data))
  #prepare model.frame and handle missing values
  #subset data to relevant columns
  if (is.null(id_column)){
    implicit_vars_formula <- as.formula(paste(moult_column, '~ ', date_column))
  } else {
    implicit_vars_formula <- as.formula(paste(moult_column, '~ ', date_column, '+', id_column))
  }
  data <- model.frame(nlme::asOneFormula(start_formula, duration_formula, sigma_formula, implicit_vars_formula), data = data)

  #check data encoding is as expected
  stopifnot(is.numeric(data[[date_column]]))
  if(type %in% c(2:5)) {
    stopifnot(all(data[[moult_column]] >= 0 & data[[moult_column]] <= 1))
  } else {
      if (type == 1) {
        stopifnot(all(data[[moult_column]] %in% c(1,2,3))|all(data[[moult_column]] >= 0 & data[[moult_column]] <= 1))
      } else {
        if (type == 12) {
          stopifnot(all(data[[moult_column]] >= 0 & data[[moult_column]] <= 1 | data[[moult_column]]==2))
        } else {
          stop('Unsupported model type. Argument `type` must be one of 1,2,3,4,5,12')
        }
      }
  }
  #delete irrelevant observations depending on model type
  data <- switch(as.character(type),
                 "3"=droplevels(data[data[[moult_column]] > 0 & data[[moult_column]] < 1,]),
                 "4"=droplevels(data[data[[moult_column]] > 0,]),
                 "5"=droplevels(data[data[[moult_column]] < 1,]),
                 droplevels(data))
  #recode type 1 data in 1,2,3 format
  #TODO: this overwrites the input data, so decision to be made whether model slot data returns data as fitted or data as supplied - latter makes more sense in terms of values, but how to deal with deleted values in type 3,4,5 models? could pass them through and change indexing in stan model?!
  if (type == 1 & all(data[[moult_column]] %in% c(1,2,3))){
    data[[moult_column]][data[[moult_column]]==1] <- 0
    data[[moult_column]][data[[moult_column]]==2] <- 0.5
    data[[moult_column]][data[[moult_column]]==3] <- 1
  }


  #order data by moult category
  data <- data[order(data[[moult_column]]),]
  #TODO: ensure moult data is correctly coded for type 1 and type 12 model


  #setup model matrices
  X_mu <- model.matrix(start_formula, data)
  X_tau <- model.matrix(duration_formula, data)
  X_sigma <- model.matrix(sigma_formula, data)

    #prepare data structure for stan
  standata <- list(old_dates = data[[date_column]][data[[moult_column]]==0],
                   N_old = length(data[[date_column]][data[[moult_column]]==0]),
                   moult_dates  = data[[date_column]][(data[[moult_column]] > 0 & data[[moult_column]] < 1)],
                   moult_indices = data[[moult_column]][(data[[moult_column]] > 0 & data[[moult_column]] < 1)],
                   N_moult = length(data[[date_column]][(data[[moult_column]] > 0 & data[[moult_column]] < 1)]),
                   new_dates = data[[date_column]][data[[moult_column]]==1],
                   N_new = length(data[[date_column]][data[[moult_column]]==1]),
                   X_mu = X_mu,
                   N_pred_mu = ncol(X_mu),
                   X_tau = X_tau,
                   N_pred_tau = ncol(X_tau),
                   X_sigma = X_sigma,
                   N_pred_sigma = ncol(X_sigma),
                   lumped = as.numeric(lump_non_moult),
                   llik = as.numeric(log_lik),
                   use_phi_approx = as.numeric(use_phi_approx),
                   active_moult_recaps_only = as.numeric(active_moult_recaps_only),
                   same_sigma = as.numeric(same_sigma))
  # setup replication information
  if (!is.null(id_column)){
    stopifnot(is.factor(data[[id_column]]))
    #create vector of first occurrence of each individual in the model matrix to index unique linear predictor values for each individual
    id_first <- match(unique(data[[id_column]]), data[[id_column]])
    #create indicator  variable to identify active_moult records
    active_moult <- data[[moult_column]] != 1 & data[[moult_column]] != 0

    if(active_moult_recaps_only==TRUE){
      replicated <- which(data[[id_column]] %in% rownames(table(data[[id_column]],active_moult))[table(data[[id_column]],active_moult)[,'TRUE']>0 & rowSums(table(data[[id_column]],active_moult))>1])
      is_replicated <- ifelse(table(data[[id_column]],active_moult)[,'TRUE']>0 & rowSums(table(data[[id_column]],active_moult))>1, 1,0)#vector of length N_ind
    } else {
      replicated <- which(data[[id_column]] %in% names(table(data[[id_column]])[table(data[[id_column]])>1]))#indices of vector of length N_obs
      is_replicated <- ifelse(table(data[[id_column]])>1, 1,0)#vector of length N_ind
    }
    not_replicated <- setdiff(seq_along(data[[id_column]]),replicated)

    standata$N_ind = length(unique(data[[id_column]]))
    standata$N_ind_rep = length(unique(as.numeric(data[[id_column]])[replicated]))
    standata$individual = as.numeric(data[[id_column]])#TODO: the resultant individual intercepts are hard to map onto original factor levels - this should be handled in the postprocessing of the model output
    standata$individual_first_index = as.array(id_first)
    standata$replicated = as.array(replicated)
    standata$not_replicated = as.array(not_replicated)
    if (type %in% c(2,5)){
      standata$not_replicated_old = as.array(not_replicated[not_replicated <= standata$N_old])
      standata$not_replicated_moult = as.array(not_replicated[not_replicated > standata$N_old] - standata$N_old)
    }
    if (type == 4){
      standata$not_replicated_moult = as.array(not_replicated[not_replicated <= standata$N_moult])
      standata$not_replicated_new = as.array(not_replicated[not_replicated > standata$N_moult] - standata$N_moult)
    }
    standata$is_replicated = as.array(is_replicated)
    standata$replicated_individuals = as.array(unique(as.numeric(data[[id_column]])[replicated]))
    standata$Nobs_replicated = length(replicated)
    if (type %in% c(2,5)){
      standata$Nobs_not_replicated_old = length(standata$not_replicated_old)
      standata$Nobs_not_replicated_moult = length(standata$not_replicated_moult)
    }
    if (type == 4){
      standata$Nobs_not_replicated_new = length(standata$not_replicated_new)
      standata$Nobs_not_replicated_moult = length(standata$not_replicated_moult)
    }
    if(standata$N_ind_rep < 10) warning(paste0('There are only ', standata$N_ind_rep, ' individuals with recapture information.\n This dataset may not be suitable for the recaptures moult model.'))
  }
  #additional input data for type 12 model
  if (type == 12){
    standata$N_moult_cat = length(data[[date_column]][data[[moult_column]]==2])
    standata$moult_cat_dates =data[[date_column]][data[[moult_column]]==2]
  }
  #include pointwise log_lik matrix  in output?
  if(log_lik){
    outpars <- c('beta_mu','beta_tau','beta_sigma', 'sigma_intercept', 'log_lik')
  } else {
    outpars <- c('beta_mu','beta_tau','beta_sigma', 'sigma_intercept')
  }
  #get name of relevant model object
  if (is.null(id_column)){
    stan_model_name <- paste0('uz',type,'_linpred')
  } else {
    stan_model_name <- paste0('uz',type,'_recap')
    outpars <- c(outpars, 'mu_ind_star', 'mu_ind_out','sigma_mu_ind')
    outpars <- gsub('beta_mu','beta_mu_out', outpars)
  }

  #guess initial values and sample
  if(init == "auto"){
    date_on_score_lm <- lm(standata$moult_dates ~ standata$moult_indices)
    mu_start = coef(date_on_score_lm)[1]
    tau_start = max(coef(date_on_score_lm)[2],2*sd(standata$moult_dates), na.rm=T)
    sigma_start = sd(standata$moult_dates)
    sigma_mu_ind_start = jitter(20)
    initfunc <- function(chain_id = 1) {
      # cat("chain_id =", chain_id, "\n")
      list(beta_mu = as.array(c(mu_start,rep(0, standata$N_pred_mu - 1))), #initialize intercept term from data, set inits for all other effects to 0
           beta_tau = as.array(c(tau_start, rep(0, standata$N_pred_tau - 1))),
           beta_sigma = as.array(c(log(sigma_start), rep(0, standata$N_pred_sigma - 1))),#NB on log link scale
           sigma_mu_ind = sigma_mu_ind_start)
    }
    out <- rstan::sampling(stanmodels[[stan_model_name]], data = standata, init = initfunc, pars = outpars, ...)
  } else {
    out <- rstan::sampling(stanmodels[[stan_model_name]], data = standata, init = init, pars = outpars, ...)
  }
  #rename regression coefficients for output
  names(out)[grep('beta_mu', names(out))] <- paste('mean',colnames(X_mu), sep = '_')
  names(out)[grep('beta_tau', names(out))] <- paste('duration',colnames(X_tau), sep = '_')
  names(out)[grep('beta_sigma', names(out))] <- paste('log_sd',colnames(X_sigma), sep = '_')
  names(out)[grep('sigma_intercept', names(out))] <- 'sd_(Intercept)'
  out_struc <- list()
  out_struc$stanfit <- out
  out_struc$terms$date_column <- date_column
  out_struc$terms$moult_index_column <- moult_column
  out_struc$terms$moult_cat_column <- NA
  out_struc$terms$id_column <- id_column
  out_struc$terms$start_formula <- start_formula
  out_struc$terms$duration_formula <- duration_formula
  out_struc$terms$sigma_formula <- sigma_formula
  out_struc$data <- data
  out_struc$na.action <- attr(data, "na.action")
  out_struc$type = paste0(type,ifelse(lump_non_moult,'L', ''), ifelse(is.null(id_column), '','R'))
  out_struc$individual_ids <- if (is.null(id_column)) { NA } else { data.frame(index = as.numeric(unique(data[[id_column]])), id = unique(data[[id_column]])) }
  out_struc$replicated_ids <- if (is.null(id_column)) { NA } else { data.frame(index = 1:standata$N_ind_rep, id = levels(data[[id_column]])[standata$replicated_individuals]) }
  class(out_struc) <- 'moultmcmc'
  return(out_struc)
}
