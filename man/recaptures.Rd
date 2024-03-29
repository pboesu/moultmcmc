% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{recaptures}
\alias{recaptures}
\title{Simulated moult scores for 235 individuals with no heterogeneity in duration}
\format{
A data frame with 299 rows and 13 variables:
\describe{
\item{id}{individual identifier}
\item{sample_size}{number of observations for the said individual}
\item{start_date}{individual moult initiation date }
\item{duration}{individual duration}
\item{set}{simulation set}
\item{date_sampled}{observation date (days since January 01)}
\item{pfmg_sampled}{true PFMG at sampling occasion}
\item{pfmg_with_error}{PFMG incorporating measurement error, if any}
\item{sampling_bias}{descriptor of sampling bias type, if any}
\item{prop_recaptured_individuals}{proportion of individuals with more than 1 capture}
\item{rel_tau_sd}{relative population standard deviation of the moult duration}
\item{obs_error}{measurement error}
\item{scenario}{simulation scenario}
}
}
\source{
P. Boersch-Supan
}
\usage{
recaptures
}
\description{
A dataset containing the moult scores and capture dates for a simulated population of passerines including recaptures
}
\keyword{datasets}
