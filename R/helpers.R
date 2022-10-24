#helper functions to plot moult and moultmcmc models

#' Consolidate a mixture of continuous and categorical moult records
#'
#' This is a helper function to format input data for the Type 1+2 moult model.
#'
#' @param moult_score a numeric vector of (linearized) moult scores in \[0,1\] (0 = old plumage,1 = new plumage).
#' @param moult_cat a numeric vector of categorical moult codes (1 = old plumage,2 = moulting, 3 = new plumage)
#'
#' @return a numeric vector of scores and categorical records with values \[0,1\] for old, new and continuous active moult scores, and value 2 for categorical active moult records.
#' @export
#'
#'
consolidate_moult_records <- function(moult_score, moult_cat){
  stopifnot(all((moult_score >= 0 & moult_score <= 1) | is.na(moult_score)))
  stopifnot(all(moult_cat %in% c(1,2,3,NA)))
  stopifnot(length(moult_score)==length(moult_cat))
  moult_score[is.na(moult_score) & moult_cat == 1] <- 0
  moult_score[is.na(moult_score) & moult_cat == 2] <- 2
  moult_score[is.na(moult_score) & moult_cat == 3] <- 1
  return(moult_score)
}


#' Extract Population-Level Estimates
#'
#' Extract the population-level ('fixed') effects
#' from a \code{moultmcmc} object.
#'
#' @aliases fixef
#'
#' @param object a moultmcmc model
#' @param summary logical, should posterior samples be summarised for each parameter
#' @param probs numeric, desired quantiles for summary statistics
#' @param pars Optional names of coefficients to extract.
#'   By default, all coefficients are extracted.
#' @param ... Currently ignored.
#'
#' @return If \code{summary} is \code{TRUE}, a matrix for the population-level effects.
#'   If \code{summary} is \code{FALSE}, a matrix with one row per
#'   posterior draw and one column per population-level effect.
#'
#'
#' @method fixef moultmcmc
#' @export
#' @export fixef
#' @importFrom nlme fixef
fixef.moultmcmc <-  function(object, summary = TRUE,
                           probs = c(0.025, 0.975), pars = NULL, ...) {
  fpars <- names(object$stanfit)[!grepl('lp__|log_lik|mu_ind',names(object$stanfit))] #TODO: need to filter out hierarchical pars here and things like lp__ and log_lik
  if (!is.null(pars)) {
    fpars <- as.character(pars)
  }
  if (!length(fpars)) {
    return(NULL)
  }
  out <- as.matrix(object$stanfit, pars = fpars)
  if (summary) {
    out <- summary(object$stanfit, probs = probs, pars = fpars)$summary
  }
  out
}

#' Extract Individual-Level Estimates
#'
#' Extract the individual-level ('random') effects
#' from a \code{moultmcmc} object.
#'
#' @aliases ranef
#'
#' @param object a moultmcmc model
#' @param summary logical, should posterior samples be summarised for each parameter
#' @param probs numeric, desired quantiles for summary statistics
#' @param pars Optional names of coefficients to extract.
#'   By default, all coefficients are extracted.
#' @param ... Currently ignored.
#'
#' @return If \code{summary} is \code{TRUE}, a matrix for the individual-level effects.
#'   If \code{summary} is \code{FALSE}, a matrix with one row per
#'   posterior draw and one column per individual-level effect.
#'
#'
#' @method ranef moultmcmc
#' @export
#' @export ranef
#' @importFrom nlme ranef
ranef.moultmcmc <-  function(object, summary = TRUE,
                             probs = c(0.025, 0.975), pars = NULL, ...) {
  fpars <- names(object$stanfit)[grepl('mu_ind_star',names(object$stanfit))]
  if (!is.null(pars)) {
    fpars <- as.character(pars)
  }
  if (!length(fpars)) {
    return(NULL)
  }
  out <- as.matrix(object$stanfit, pars = fpars)
  if (summary) {
    out <- summary(object$stanfit, probs = probs, pars = fpars)$summary
    rownames(out)<-object$replicated_ids$id
  }
  out
}


#need to define generic function
#' Summary table generic
#'
#' @param ... ...
#'
#' @return ...
#' @export
#'
summary_table <- function(...){UseMethod("summary_table")}

#' Summary Table  for moult model
#'
#' @param x moult model object created with moult::moult
#' @param prob nominal coverage probability of confidence interval
#' @param tidy_names adjust default parameter names from moult() to follow consistent nomenclature of model.matrix()
#' @param ... not currently used
#'
#' @return a tibble
#'
#' @importFrom stats coef qnorm
#' @importFrom rlang .data
#' @export
#'
#' @examples \dontrun{m1 <- moult::moult(Mindex ~ Day, data = sanderlings)}
#'\dontrun{summary_table(m1)}
#'
summary_table.moult <-function(x, prob = 0.95, tidy_names = TRUE, ...){
  probs = c((1-prob)/2, 1 -(1-prob)/2)
  plotdata <- tibble::tibble(parameter = names(coef(x)), estimate = coef(x), stderr = x$standard.errors, lci = coef(x) + qnorm(probs[1])*x$standard.errors, uci = coef(x) + qnorm(probs[2])*x$standard.errors, prob = prob, method = 'ML',type = as.character(x$type), Rhat = NA)
  if (tidy_names)
  plotdata <- mutate(plotdata, parameter = stringr::str_replace_all(.data$parameter, "intercept.1", "(Intercept)"))
  plotdata <- mutate(plotdata, parameter = stringr::str_replace_all(.data$parameter, "duration_duration", "duration_(Intercept)"))
  plotdata <- mutate(plotdata, parameter = stringr::str_replace_all(.data$parameter, "mean-start-day", "(Intercept)"))
  plotdata <- mutate(plotdata, parameter = stringr::str_replace_all(.data$parameter, "SD-start-day", "(Intercept)"))
  plotdata <- mutate(plotdata, parameter = stringr::str_replace_all(.data$parameter, "intercept", "(Intercept)"))
  return(plotdata)
}


#' Summary Table
#'
#' @param x moultmcmc fit object
#' @param pars A character vector of parameter names. The default is all parameters for which samples are saved. If include = FALSE, then the specified parameters are excluded from the printed summary.
#' @param prob nominal coverage probability of credible interval
#' @param include Logical scalar (defaulting to TRUE) indicating whether to include or exclude the parameters named by the pars argument.
#' @param ... Additional arguments passed to the summary method for stanfit objects.
#' @importFrom rstan summary
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#' @return a tibble
#' @export
#'
summary_table.moultmcmc <- function (x, pars = x$stanfit@sim$pars_oi, prob = 0.95, include = TRUE, ...){
  probs = c((1-prob)/2, 1 -(1-prob)/2)

  if (x$stanfit@mode == 1L) {
    cat("Stan model '", x$stanfit@model_name, "' is of mode 'test_grad';\n",
        "sampling is not conducted.\n", sep = "")
    return(invisible(NULL))
  }
  else if (x$stanfit@mode == 2L) {
    cat("Stan model '", x$stanfit@model_name, "' does not contain samples.\n",
        sep = "")
    return(invisible(NULL))
  }
  if (!include)
    pars <- setdiff(x$stanfit@sim$pars_oi, pars)
  s <- tibble::as_tibble(summary(x$stanfit, pars, probs, ...)$summary, rownames = "parameter") %>%
    dplyr::rename(estimate = mean) %>% dplyr::rename(lci = .data$`2.5%`, uci = .data$`97.5%`) %>% #TODO: this falls over when prob != 0.95. need to select by position or by assembling column names from prob or the stanfit summary
    mutate(prob = 0.95, method = 'MCMC', type = as.character(x$type))
  return(s)
}

#' Summary method for moultmcmc models
#'
#' @method summary moultmcmc
#' @param object a moultmcmc object
#' @param ... passed to summary_table
#'
#' @return a summary table
#' @export
#'
summary.moultmcmc <- function(object, ...){
  summary_table(object, ...)
}

#' internal helper function to create a named list
#'
#' original version by Ben Bolker: https://stackoverflow.com/a/16951524
#'
#' @param ... args
#' @importFrom stats setNames
#' @return a named list
#'
namedList <- function(...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)),deparse)[-1]
  if (is.null(nm <- names(L))) nm <- snm
  if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
  setNames(L,nm)
}



#' Efficient approximate leave-one-out cross-validation (LOO)
#'
#' @param x A moultmcmc model object
#' @param cores The number of cores to use for parallelization. This defaults to the option mc.cores which can be set for an entire R session by options(mc.cores = NUMBER). The old option loo.cores is now deprecated but will be given precedence over mc.cores until loo.cores is removed in a future release. As of version 2.0.0 the default is now 1 core if mc.cores is not set, but we recommend using as many (or close to as many) cores as possible. Note for Windows 10 users: it is strongly recommended to avoid using the .Rprofile file to set mc.cores (using the cores argument or setting mc.cores interactively or in a script is fine).
#'
#' @importFrom loo extract_log_lik relative_eff loo
#' @return The loo() methods return a named list with class c("psis_loo", "loo"). See ?loo::loo
#' @export
#'
loo.moultmcmc <- function(x, cores = getOption("mc.cores", 1)){
  # Extract pointwise log-likelihood
  # using merge_chains=FALSE returns an array, which is easier to
  # use with relative_eff()
 log_lik_1 <- loo::extract_log_lik(x, merge_chains = FALSE)

  # as of loo v2.0.0 we can optionally provide relative effective sample sizes
  # when calling loo, which allows for better estimates of the PSIS effective
  # sample sizes and Monte Carlo error
  r_eff <- loo::relative_eff(exp(log_lik_1), cores = cores)

  # preferably use more than 2 cores (as many cores as possible)
  # will use value of 'mc.cores' option if cores is not specified
  loo_1 <- loo::loo(log_lik_1, r_eff = r_eff, cores = cores)
  return(loo_1)
}




