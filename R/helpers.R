#helper functions to plot moult and moultmcmc models

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
  plotdata <- tibble::tibble(parameter = names(coef(x)), estimate = coef(x), stderr = x$standard.errors, lci = coef(x) - qnorm(probs[1])*x$standard.errors, uci = coef(x) + qnorm(probs[1])*x$standard.errors, prob = prob, method = 'ML')
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

#' Visual comparison of ML and MCMC fits
#'
#' @param m1 a moult or moultmcmc model
#' @param m2 a moult or moultmcmc model
#' @param names optional character vector of model names
#'
#' @return a plot
#'
#' @importFrom dplyr bind_rows mutate filter
#' @importFrom ggplot2 ggplot geom_pointrange aes position_dodge
#' @importFrom rlang .data
#' @export
#'
compare_plot <- function(m1,m2,names = NULL){
  #TODO:type checking etc
  #TODO:import necessary dplyr and ggplot components
  if(is.null(names)) names = as.character(1:2)

  plotdata <- dplyr::bind_rows(summary_table(m1) %>% dplyr::mutate(model = names[1]),
                               summary_table(m2) %>% dplyr::mutate(model = names[2]))
  dplyr::filter(plotdata, !(.data$parameter %in% c('lp__', 'log_sd_(Intercept)'))) %>%
    ggplot(aes(x = .data$parameter, y = .data$estimate, col = .data$model, ymin = .data$lci, ymax = .data$uci)) + geom_pointrange(position = position_dodge(0.1))
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


