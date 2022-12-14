
### Markov chain traceplots - too confused by S3/S4 method dispatch

# #' Draw the traceplot corresponding to one or more Markov chains, providing a visual way to inspect sampling behavior and assess
# #' mixing across chains and convergence.
# #' A simple wrapper for `rstan::traceplot`
# #' @importMethodsFrom rstan traceplot
# #' @export
# traceplot.moultmcmc <- function(object, ...){
#   rstan::traceplot(object$stanfit, ...)
# }

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
#' @param plot_data logical, if TRUE plot observations
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
moult_plot.moult <-function(x, prob = 0.95, plot_data = TRUE, plot = TRUE, ...){
  probs = c((1-prob)/2, 1 -(1-prob)/2)
  quantile_name <- paste(prob*100, '% quantile')
  data_x <- strsplit(as.character(x$terms$full)[3], split = ' ')[[1]][1]
  data_y <- as.character(x$terms$full)[2]
  plotdata <- tibble::tibble(start_date = c(x$coefficients$mean[1],
                                            qnorm(probs[1])*x$coefficients$sd[1]+x$coefficients$mean[1],
                                            qnorm(probs[2])*x$coefficients$sd[1]+x$coefficients$mean[1]),
                             end_date = c(x$coefficients$mean[1]+x$coefficients$duration[1],
                                          qnorm(probs[1])*x$coefficients$sd[1]+x$coefficients$mean[1]+x$coefficients$duration[1],
                                          qnorm(probs[2])*x$coefficients$sd[1]+x$coefficients$mean[1]+x$coefficients$duration[1]),
                             line_type = c('Population mean', quantile_name, quantile_name))
  if (plot) {
    mplot <- ggplot(plotdata, aes(x = .data$start_date, xend = .data$end_date, y = 0, yend = 1, lty = .data$line_type)) + scale_linetype_manual(values = c(3,1), name = '') + geom_segment() + theme_classic() + xlab('Date') + ylab('Moult Index')
    if(plot_data){
      #check data matches model
      if(data_x %in% names(x$X) & data_y %in% names(x$X)){
        mplot <- ggplot(data = x$X, aes(x = get(data_x), y = get(data_y))) + scale_linetype_manual(values = c(3,1), name = '') + geom_point(col = 'darkgrey') + geom_segment(data = plotdata, aes(x = .data$start_date, xend = .data$end_date, y = 0, yend = 1, lty = .data$line_type)) + theme_classic() + xlab('Date') + ylab('Moult Index')}
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
#' @param plot_data logical, if TRUE plot observations
#' @param newdata optional data.frame of conditions to plot (passed to predict.moultmcmc)
#' @param col.data optional, name of column in model data to colour data points
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
moult_plot.moultmcmc <-function(x, prob = 0.95, prob_ci = NULL, plot_data = TRUE, plot= TRUE, newdata = NULL, col.data = NULL, ...){
  probs = c((1-prob)/2, 1 -(1-prob)/2)
  if(!is.null(prob_ci)) probs_ci = c((1-prob_ci)/2, 1 -(1-prob_ci)/2)
  quantile_name <- paste(prob*100, '% quantile')
  data_x <- x$terms$date_column
  data_y <- ifelse(!is.na(x$terms$moult_index_column), x$terms$moult_index_column, x$terms$moult_cat_column)
  if (is.null(newdata)){
    plotdata <- tibble::tibble(start_date = c(fixef(x)['mean_(Intercept)',1],
                                              qnorm(probs[1])*fixef(x)['sd_(Intercept)',1]+fixef(x)['mean_(Intercept)',1],
                                              qnorm(probs[2])*fixef(x)['sd_(Intercept)',1]+fixef(x)['mean_(Intercept)',1]),
                               end_date = c(fixef(x)['mean_(Intercept)',1]+fixef(x)['duration_(Intercept)',1],
                                            qnorm(probs[1])*fixef(x)['sd_(Intercept)',1]+fixef(x)['mean_(Intercept)',1]+fixef(x)['duration_(Intercept)',1],
                                            qnorm(probs[2])*fixef(x)['sd_(Intercept)',1]+fixef(x)['mean_(Intercept)',1]+fixef(x)['duration_(Intercept)',1]),
                               line_type = c('Population mean', quantile_name, quantile_name),
                               scenario = "1")

  } else {
    plotdata <- bind_rows(predict(x, predict.type = 'parameters', newdata = newdata) %>%
                            mutate(line_type = "Population mean", scenario = rownames(newdata)),
                          predict(x, predict.type = 'parameters', newdata = newdata) %>%
                            mutate(line_type = quantile_name,
                                   scenario = rownames(newdata),
                                   start_date = .data$start_date + qnorm(probs[1])*.data$start_sd,
                                   end_date = .data$end_date + qnorm(probs[1])*.data$start_sd),
                          predict(x, predict.type = 'parameters', newdata = newdata) %>%
                            mutate(line_type = quantile_name,
                                   scenario = rownames(newdata),
                                   start_date = .data$start_date - qnorm(probs[1])*.data$start_sd,
                                   end_date = .data$end_date - qnorm(probs[1])*.data$start_sd))
    newdata$scenario = rownames(newdata)
    plotdata = left_join(plotdata,newdata, by = 'scenario')
  }
  if(!is.null(prob_ci)){
  }
  if (plot) {
    mplot <- ggplot(plotdata, aes(x = .data$start_date, xend = .data$end_date, y = 0, yend = 1, lty = .data$line_type, col = .data$scenario)) + scale_linetype_manual(values = c(3,1), name = '') + geom_segment() + theme_classic() + xlab('Date') + ylab('Moult Index')
    if(plot_data){
      #check data matches model
      if(data_x %in% names(x$data) & data_y %in% names(x$data)){
        if (!is.null(col.data)){
          mplot <- ggplot(data = x$data, aes(x = get(data_x), y = get(data_y))) + scale_linetype_manual(values = c(3,1), name = '') + geom_point(aes(col = get(col.data)), alpha = 0.5) + geom_segment(data = plotdata, aes(x = .data$start_date, xend = .data$end_date, y = 0, yend = 1, lty = .data$line_type, col = .data$scenario), lwd = 1) + theme_classic() + xlab('Date') + ylab('Moult Index')
        } else {
          mplot <- ggplot(data = x$data, aes(x = get(data_x), y = get(data_y))) + scale_linetype_manual(values = c(3,1), name = '') + geom_point(col='darkgrey') + geom_segment(data = plotdata, aes(x = .data$start_date, xend = .data$end_date, y = 0, yend = 1, lty = .data$line_type, col = .data$scenario), lwd = 1) + theme_classic() + xlab('Date') + ylab('Moult Index')
        }
        }
      else {warning(paste('data does not contain model date and moult variables:',data_x, data_y, ' - plotting model only'))}
    }
    return(mplot) } else {
      return(plotdata)
    }
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
  parlist <- moultmcmc:::namedList(...)
  stopifnot(all(sapply(parlist, class) %in% c('moult','moultmcmc')))
  #TODO:type checking etc
  #TODO:import necessary dplyr and ggplot components
  if(is.null(names)) names = names(parlist)
  names(parlist) <- names

  plotdata <- dplyr::bind_rows(lapply(parlist, function(x){summary_table(x)}), .id = 'model')
  plotdata$not_converged <- ifelse(plotdata$Rhat > 1.05 & !is.na(plotdata$Rhat), TRUE, FALSE)
  dplyr::filter(plotdata, !grepl("lp__|log_sd_\\(Intercept\\)|\\blp\\b|log_lik[[0-9]+]|mu_ind[[0-9]+]|mu_ind_star|mu_ind_out|tau_ind_out", .data$parameter)) %>%
    ggplot(aes(x = .data$model, y = .data$estimate, col = .data$model, ymin = .data$lci, ymax = .data$uci, shape = .data$not_converged)) + geom_pointrange(position = position_dodge(0.1)) + facet_wrap(~ .data$parameter, scales = 'free')
}

#' residual_plot generic
#'
#' @param ... ...
#'
#' @return ...
#' @export
#'
residual_plot <- function(...){UseMethod("residual_plot")}

#' Residual plot for moult model
#'
#' This function displays a plot that shows where observed dates and moult scores fall relative to the predictions of the fitted model. Active moult observations should generally fall within +/- 3 start-date standard deviations of the regression line.
#'
#' @param x moult model object created with moult::moult
#' @param plot logical, if TRUE (default) return a plot, else return a data.frame with calculated quantities
#' @param ... not currently used
#'
#' @return a plot
#'
#' @importFrom stats coef qnorm
#' @importFrom ggplot2 geom_hline aes scale_shape_manual ggplot theme theme_classic xlab ylab ggtitle geom_point
#' @export
#'
#'
residual_plot.moult <-function(x, plot = TRUE, ...){
  model_data <- x$X
  model_data$duration_pred <- predict(x, newdata = x$X, predict.type = 'duration')$duration
  model_data$start_pred <- predict(x, newdata = x$X, predict.type = 'start')$mean.start
  #calculate sd predictions by hand
  p.sd <- unlist(x$coefficients$sd)
  mm <- model.matrix(x$terms$sd, x$X)
  model_data$sd_pred <- mm %*% p.sd

  model_data$delta_day <-  x$Day - (x$y[[1]]*model_data$duration_pred+model_data$start_pred)
  model_data$delta_day_zscore <- model_data$delta_day / model_data$sd_pred
  model_data$active_moult <- x$y[[1]] > 0 & x$y[[1]] < 1

  if(plot == TRUE){
    ggplot(model_data, aes(x = get(names(model_data)[1]), y = .data$delta_day_zscore, pch = .data$active_moult)) + geom_point() + geom_hline(yintercept = -3:3, lty = 2, col = 'grey') + geom_hline(yintercept = 0) + scale_shape_manual(values = c(1,16)) + xlab('PFMG_observed') + ylab('Population SD of start date') + theme_classic() + theme(legend.position = 'bottom') + ggtitle('Experimental "residual" plot')
  } else {
    return(model_data)
  }

}
