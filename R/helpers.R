#helper functions to plot moult and moultmcmc models

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
  fpars <- names(object$stanfit) #TODO: need to filter out hierarchical pars here and things like lp__ and log_lik
  if (!is.null(pars)) {
    fpars <- as.character(pars)
  }
  if (!length(fpars)) {
    return(NULL)
  }
  out <- as.matrix(object$stanfit, pars = fpars)
  if (summary) {
    out <- summary(object$stanfit, probs = probs)$summary
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
  plotdata <- tibble::tibble(parameter = names(coef(x)), estimate = coef(x), stderr = x$standard.errors, lci = coef(x) + qnorm(probs[1])*x$standard.errors, uci = coef(x) + qnorm(probs[2])*x$standard.errors, prob = prob, method = 'ML', Rhat = NA)
  if (tidy_names)
  plotdata <- mutate(plotdata, parameter = stringr::str_replace_all(.data$parameter, "intercept.1", "(Intercept)"))
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
    mutate(prob = 0.95, method = 'MCMC')
  return(s)
}

#' Visual comparison of moult models
#'
#' @param ... two or a moult or moultmcmc model
#' @param names optional character vector of model names
#'
#'
#' @return a plot comparing parameter estimates and their uncertainties
#'
#' @importFrom dplyr bind_rows mutate filter
#' @importFrom ggplot2 ggplot geom_pointrange aes position_dodge
#' @importFrom rlang .data
#' @export
#'
compare_plot <- function(...,names = NULL){
  parlist <- list(...)
  stopifnot(all(sapply(parlist, class) %in% c('moult','moultmcmc')))
  #TODO:type checking etc
  #TODO:import necessary dplyr and ggplot components
  if(is.null(names)) names = as.character(seq(1, length(parlist), by = 1))
  names(parlist) <- names

  plotdata <- dplyr::bind_rows(lapply(parlist, function(x){summary_table(x)}), .id = 'model')
  plotdata$not_converged <- ifelse(plotdata$Rhat > 1.05 & !is.na(plotdata$Rhat), TRUE, FALSE)
  dplyr::filter(plotdata, !grepl("lp__|log_sd_\\(Intercept\\)|\\blp\\b|log_lik[[0-9]+]|mu_ind[[0-9]+]|mu_ind_star", .data$parameter)) %>%
    ggplot(aes(x = .data$model, y = .data$estimate, col = .data$model, ymin = .data$lci, ymax = .data$uci, shape = .data$not_converged)) + geom_pointrange(position = position_dodge(0.1)) + facet_wrap(~ .data$parameter, scales = 'free')
}

#' Plot method for moult models
#'
#' @param x maximum likelihood fit of a moult model from the moult package (class moult)
#' @param prob coverage probability of the confidence interval
#' @param scales argument to facet_wrap. default is "free_x", "free" may enable a better comparison if parameter values are disparate
#' @param ... not currently used
#'
#' @return a plot
#' @importFrom ggplot2 facet_wrap
#' @importFrom rlang .data
#'
#' @export
#'
plot.moult <- function(x, prob = 0.95,scales = "free_x", ...){
  plotdata <- tibble::tibble(param = names(coef(x)), estimate = coef(x), stderr = x$standard.errors)
  #plot(as.factor(names(coef(x))), coef(x))
  #dotchart(coef(x))
  ggplot(plotdata, aes(x = .data$param, y = .data$estimate, ymin = .data$estimate - qnorm(prob)*stderr, ymax = .data$estimate + qnorm(prob)*stderr)) + geom_pointrange() + facet_wrap(~.data$param, scales = scales)
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


#' moult_plot generic
#'
#' @param ... ...
#'
#' @return ...
#' @export
#'
moult_plot <- function(...){UseMethod("moult_plot")}

#' Moult plot for moult model
#'
#' @param x moult model object created with moult::moult
#' @param prob coverage probability of active moult in the sampled population
#' @param data optional data.frame of observations, column names need to match model date and moult index variables
#' @param plot logical, if TRUE (default) return a plot, else return a dataframe with start and end dates of moult based on model intercepts
#' @param ... not currently used
#'
#' @return a plot or data.frame
#'
#' @importFrom stats coef qnorm
#' @importFrom ggplot2 geom_segment aes scale_linetype_manual ggplot theme_classic xlab ylab geom_point
#' @export
#'
#'
moult_plot.moult <-function(x, prob = 0.95, data = NULL, plot = TRUE, ...){
  probs = c((1-prob)/2, 1 -(1-prob)/2)
  quantile_name <- paste(prob*100, '% quantile')
  data_x <- strsplit(as.character(x$terms$full)[3], split = ' ')[[1]][1]
  data_y <- as.character(x$terms$full)[2]
  plotdata <- tibble::tibble(start_date = c(coef(x)['mean_intercept'],
                                            qnorm(probs[1])*coef(x)['sd_(Intercept)']+coef(x)['mean_intercept'],
                                            qnorm(probs[2])*coef(x)['sd_(Intercept)']+coef(x)['mean_intercept']),
                             end_date = c(coef(x)['mean_intercept']+coef(x)['duration_intercept.1'],
                                            qnorm(probs[1])*coef(x)['sd_(Intercept)']+coef(x)['mean_intercept']+coef(x)['duration_intercept.1'],
                                            qnorm(probs[2])*coef(x)['sd_(Intercept)']+coef(x)['mean_intercept']+coef(x)['duration_intercept.1']),
                             line_type = c('Population mean', quantile_name, quantile_name))
  if (plot) {
      mplot <- ggplot(plotdata, aes(x = start_date, xend = end_date, y = 0, yend = 1, lty = line_type)) + scale_linetype_manual(values = c(3,1), name = '') + geom_segment() + theme_classic() + xlab('Date') + ylab('Moult Index')
      if(!is.null(data)){
      #check data matches model
      if(data_x %in% names(data) & data_y %in% names(data)){
        mplot <- ggplot(data = data, aes(x = get(data_x), y = get(data_y))) + scale_linetype_manual(values = c(3,1), name = '') + geom_point(col = 'darkgrey') + geom_segment(data = plotdata, aes(x = start_date, xend = end_date, y = 0, yend = 1, lty = line_type)) + theme_classic() + xlab('Date') + ylab('Moult Index')}
      else {warning(paste('data does not contain model date and moult variables:',data_x, data_y, ' - plotting model only'))}
    }
    return(mplot) } else {
      return(plotdata)
    }
}

#' Moult plot for moultmcmc model
#'
#' @param x moultmcmc model object
#' @param prob coverage probability of active moult in the sampled population
#' @param prob_ci coverage probability of credible interval
#' @param plot logical, if TRUE (default) return a plot, else return a dataframe with start and end dates of moult based on model intercepts
#' @param ... not currently used
#'
#' @return a plot or data.frame
#'
#' @importFrom stats coef qnorm
#' @importFrom rlang .data
#' @importFrom ggplot2 geom_segment aes scale_linetype_manual ggplot theme_classic xlab ylab geom_point
#' @export
#'
#'
moult_plot.moultmcmc <-function(x, prob = 0.95, prob_ci = NULL, data = NULL, plot= TRUE, ...){
  probs = c((1-prob)/2, 1 -(1-prob)/2)
  if(!is.null(prob_ci)) probs_ci = c((1-prob_ci)/2, 1 -(1-prob_ci)/2)
  quantile_name <- paste(prob*100, '% quantile')
  data_x <- x$terms$date_column
  data_y <- ifelse(!is.na(x$terms$moult_index_column), x$terms$moult_index_column, x$terms$moult_cat_column)
  plotdata <- tibble::tibble(start_date = c(fixef(x)['mean_(Intercept)',1],
                                            qnorm(probs[1])*fixef(x)['sd_(Intercept)',1]+fixef(x)['mean_(Intercept)',1],
                                            qnorm(probs[2])*fixef(x)['sd_(Intercept)',1]+fixef(x)['mean_(Intercept)',1]),
                             end_date = c(fixef(x)['mean_(Intercept)',1]+fixef(x)['duration_(Intercept)',1],
                                          qnorm(probs[1])*fixef(x)['sd_(Intercept)',1]+fixef(x)['mean_(Intercept)',1]+fixef(x)['duration_(Intercept)',1],
                                          qnorm(probs[2])*fixef(x)['sd_(Intercept)',1]+fixef(x)['mean_(Intercept)',1]+fixef(x)['duration_(Intercept)',1]),
                             line_type = c('Population mean', quantile_name, quantile_name))
  if(!is.null(prob_ci)){
  }
if (plot) {
    mplot <- ggplot(plotdata, aes(x = start_date, xend = end_date, y = 0, yend = 1, lty = line_type)) + scale_linetype_manual(values = c(3,1), name = '') + geom_segment() + theme_classic() + xlab('Date') + ylab('Moult Index')
    if(!is.null(data)){
      #check data matches model
      if(data_x %in% names(data) & data_y %in% names(data)){
        mplot <- ggplot(data = data, aes(x = get(data_x), y = get(data_y))) + scale_linetype_manual(values = c(3,1), name = '') + geom_point(col = 'darkgrey') + geom_segment(data = plotdata, aes(x = start_date, xend = end_date, y = 0, yend = 1, lty = line_type)) + theme_classic() + xlab('Date') + ylab('Moult Index')}
      else {warning(paste('data does not contain model date and moult variables:',data_x, data_y, ' - plotting model only'))}
    }
    return(mplot) } else {
      return(plotdata)
    }
}
