#' Predict method for moultmcmc models
#'
#' @param object a fitted moultmcmc model
#' @param newdata dataframe with explanatory variables for which to make predictions
#' @param predict.type specifies form of predictions, see details.
#' @param intervals not currently used
#' @param ... further arguments
#'
#' @return a data.frame?!
#' @importFrom stats predict
#' @export
#'
predict.moultmcmc <- function(object, newdata = NULL, predict.type = "start", intervals = 0.1, ...){
  #if (!is.element(predict.type, c("start","duration","end"))) stop("Type of prediction not available")

  switch(predict.type,

         start = {

         },

         stop('predict.type not valid'))
}
