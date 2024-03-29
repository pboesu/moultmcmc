% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict.moultmcmc}
\alias{predict.moultmcmc}
\title{Predict method for moultmcmc models}
\usage{
\method{predict}{moultmcmc}(
  object,
  newdata = NULL,
  predict.type = "parameters",
  summary = TRUE,
  intervals = 0.1,
  prob = 0.95,
  ...
)
}
\arguments{
\item{object}{a fitted moultmcmc model}

\item{newdata}{data.frame with explanatory variables for which to make predictions}

\item{predict.type}{specifies form of predictions, see details.}

\item{summary}{logical, if TRUE (default) return a data.frame of posterior means, otherwise return a list of arrays of the full posterior sample of the predicted quantity (with one list element per predicted quantity and array dimensions nrow(newdata) *number or posterior draws)}

\item{intervals}{not currently used}

\item{prob}{nominal coverage probability of credible interval}

\item{...}{further arguments}
}
\value{
a data.frame or list, depending on input arguments
}
\description{
Predict method for moultmcmc models
}
