#' Predict method for moultmcmc models
#'
#' @param object a fitted moultmcmc model
#' @param newdata data.frame with explanatory variables for which to make predictions
#' @param predict.type specifies form of predictions, see details.
#' @param summary logical, if TRUE (default) return a data.frame of posterior means, otherwise return a list of arrays of the full posterior sample of the predicted quantity (with one list element per predicted quantity and array dimensions nrow(newdata) *number or posterior draws)
#' @param intervals not currently used
#' @param ... further arguments
#'
#' @return a data.frame or list, depending on input arguments
#' @importFrom stats predict
#' @export
#'
predict.moultmcmc <- function(object, newdata = NULL, predict.type = "parameters", summary = TRUE, intervals = 0.1, ...){
  #if (!is.element(predict.type, c("start","duration","end"))) stop("Type of prediction not available")

  #set up data
  if(is.null(newdata)){
    df_pred <- object$data
  } else {
    df_pred <- newdata
  }

  #extract model coefficients
  #names(out)[grep('beta_mu', names(out))] <- paste('mean',colnames(X_mu), sep = '_')
  beta_mu <- do.call(cbind, rstan::extract(object$stanfit, pars = names(object$stanfit)[grepl('mean_', names(object$stanfit))]))
  beta_tau <- do.call(cbind,rstan::extract(object$stanfit, pars = names(object$stanfit)[grepl('duration_', names(object$stanfit))]))
  beta_sigma <- do.call(cbind,rstan::extract(object$stanfit, pars = names(object$stanfit)[grepl('log_sd_', names(object$stanfit))]))

  #must coerce the lists returned by rstan::extract into a matrix - use do.call(cbind, list)?

  switch(predict.type,

         parameters = {
           X_mu <- model.matrix(object$terms$start_formula, df_pred)#TODO: na.action missing
           X_tau <- model.matrix(object$terms$duration_formula, df_pred)
           X_sigma <- model.matrix(object$terms$sigma_formula, df_pred)

           start_date <- X_mu %*% t(beta_mu)
           duration <- X_tau %*% t(beta_tau)
           start_sd <- exp(X_sigma %*% t(beta_sigma))

           if (summary){
             return(data.frame(start_date = rowMeans(start_date), duration = rowMeans(duration), start_sd = rowMeans(start_sd), end_date= rowMeans(start_date + duration)%%365))
           } else {
             return(list(start_date = start_date, duration = duration, start_sd = start_sd, end_date= (start_date + duration)%%365))
           }



         },

         stop('predict.type not valid'))
}
