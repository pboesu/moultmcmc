#' The 'moultmcmc' package.
#'
#' @description Precompiled Stan models to conduct full Bayesian inference in regression models for the phenology of avian primary moult based on the approach of Underhill & Zucchini (1988) <DOI:10.1111/j.1474-919X.1988.tb00993.x> and related hierarchical models to accommodate repeated-measures data (within-season recaptures of individuals).
#'
#' @docType package
#' @name moultmcmc-package
#' @aliases moultmcmc
#' @useDynLib moultmcmc, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' Boersch-SUpan et al. (2022) Bayesian inference for models of moult duration and timing in birds arXiv <DOI:10.48550/arXiv.2205.12120>
#'
#' Underhill & Zucchini (1988) A model for avian primary moult. Ibis 130:358 <DOI:10.1111/j.1474-919X.1988.tb00993.x>
#'
#' Stan Development Team (2022). RStan: the R interface to Stan. R package version 2.21.5. https://mc-stan.org
#'
NULL

## usethis namespace: start
#' @importFrom magrittr %>%
## usethis namespace: end
NULL
