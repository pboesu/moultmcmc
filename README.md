# moultmcmc: Bayesian inference for moult phenology models

<!-- badges: start -->
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![R-CMD-check](https://github.com/pboesu/moultmcmc/workflows/R-CMD-check/badge.svg)](https://github.com/pboesu/moultmcmc/actions)
[![codecov](https://codecov.io/gh/pboesu/moultmcmc/branch/master/graph/badge.svg?token=Y7PB0302FH)](https://codecov.io/gh/pboesu/moultmcmc)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/moultmcmc)]()
<!-- badges: end -->

In most free-living bird populations moult progression and duration in individuals can not be observed fully.
Instead snapshot measurements of (re)captured individuals are typically used to infer these parameters on a population level.
As an additional complication, recording of moult in the field may take various forms both in terms of the subset of the population that is sampled and whether moult is recorded as a categorical state, or a (semi-)continuous progression.

[Underhill & Zucchini (1989; Ibis 130:358)](https://doi.org/10.1111/j.1474-919X.1988.tb00993.x) proposed a general modelling framework to accommodate many of these features, implemented in the [R package `moult`](https://cran.r-project.org/package=moult) [(Erni et al. 2013; J Stat Soft 52:8)](http://dx.doi.org/10.18637/jss.v052.i08).
A related approach based on the probit GLM was suggested by [Rothery & Newton (2002; Ibis 144:526)](http://dx.doi.org/10.1046/j.1474-919X.2002.00072.x) for use with categorical moult data, allowing for separate variances on start and end dates.
Both models are special cases of more general categorical regression models.    

`moultmcmc` implements a Bayesian inference framework for this class of models with the aim of allowing the inclusion of hierarchical model structures to accommodate 
1) the integration of moult data sets using different modes of recording, 
2) individual heterogeneity in moult timing and progression, and 
3) hierarchical spatial/temporal effects for multi-site/multi-season data sets.

`moultmcmc` implements fast inference for these models using Hamiltonian Monte Carlo samplers from Stan. 

**Please note that `moultmcmc` is under active development and API changes may still occur.**
