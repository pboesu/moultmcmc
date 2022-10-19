#' Simulated moult scores for 235 individuals with no heterogeneity in duration
#'
#' A dataset containing the moult scores and capture dates for a simulated population of passerines including recaptures
#'
#' @format A data frame with 299 rows and 13 variables:
#' \describe{
#'   \item{id}{individual identifier}
#'   \item{sample_size}{number of observations for the said individual}
#'   \item{start_date}{individual moult initiation date }
#'   \item{duration}{individual duration}
#'   \item{set}{simulation set}
#'   \item{date_sampled}{observation date (days since January 01)}
#'   \item{pfmg_sampled}{true PFMG at sampling occasion}
#'   \item{pfmg_with_error}{PFMG incorporating measurement error, if any}
#'   \item{sampling_bias}{descriptor of sampling bias type, if any}
#'   \item{prop_recaptured_individuals}{proportion of individuals with more than 1 capture}
#'   \item{rel_tau_sd}{relative population standard deviation of the moult duration}
#'   \item{obs_error}{measurement error}
#'   \item{scenario}{simulation scenario}
#' }
#' @source P. Boersch-Supan
"recaptures"

#' Simulated moult scores for 235 individuals with heterogeneity in durations
#'
#' A dataset containing the moult scores and capture dates for a simulated population of passerines including recaptures. Simulation is based on an underlying population moult start date of 196.83, moult duration of 77.80 days and a population standard deviation of the start date of 8.1 days.
#'
#' @format A data frame with 299 rows and 13 variables:
#' \describe{
#'   \item{id}{individual identifier}
#'   \item{sample_size}{number of observations for the said individual}
#'   \item{start_date}{individual moult initiation date }
#'   \item{duration}{individual duration}
#'   \item{set}{simulation set}
#'   \item{date_sampled}{observation date (days since January 01)}
#'   \item{pfmg_sampled}{true PFMG at sampling occasion}
#'   \item{pfmg_with_error}{PFMG incorporating measurement error, if any}
#'   \item{sampling_bias}{descriptor of sampling bias type, if any}
#'   \item{prop_recaptured_individuals}{proportion of individuals with more than 1 capture}
#'   \item{rel_tau_sd}{relative population standard deviation of the moult duration}
#'   \item{obs_error}{measurement error}
#'   \item{scenario}{simulation scenario}
#' }
#' @source P. Boersch-Supan
"recaptures2"


#' Moult records of Eurasian Siskins
#'
#' A dataset containing simulated moult scores and capture dates of Eurasian Siskins (Spinus spinus), based on a sub-urban Scottish population.
#'
#' @format A data frame with 299 rows and 3 variables:
#' \describe{
#'   \item{yday}{sampling date (days since Jan 01)}
#'   \item{pfmg}{PFMG at sampling occasion}
#'   \item{id}{individual identifier}
#' }
#' @source Philipp Boersch-Supan
#' @references Insley et al. (in prep). Breeding and moult phenology of siskins in the Scottish Highlands.
"siskins"

#' Moult records of Southern Masked Weavers
#'
#' A dataset containing the moult scores and capture dates for Southern Masked Weavers (Ploceus velatus) from the Western Cape, South Africa.
#'
#' @format A data frame with 5282 rows and 6 variables:
#' \describe{
#'   \item{RDate}{sampling date (calendar date)}
#'   \item{Sex}{SAFRING sex code}
#'   \item{Year}{Sampling Year}
#'   \item{Moult}{feather scores for primaries}
#'   \item{pfmg}{Proportion of feather mass grown}
#'   \item{day}{sampling date (days since August 01 of each year)}
#' }
#' @source Erni et al. 2013
#' @references Oschadleus 2005, Patterns of primary moult in the weavers, Ploceidae. University of Cape Town.
"weavers_processed"
