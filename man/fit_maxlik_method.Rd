% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcs_fitting.R
\name{fit_maxlik_method}
\alias{fit_maxlik_method}
\title{Perform max-likelihood optimization fitting}
\usage{
fit_maxlik_method(
  data,
  model_formula,
  model_family = c("automatic", "poisson", "quasipoisson", "nb2"),
  fit_link,
  aberr_module = c("dicentrics", "translocations", "micronuclei")
)
}
\arguments{
\item{data}{Count data.}

\item{model_formula}{Model formula.}

\item{model_family}{Model family.}

\item{fit_link}{Family link.}

\item{aberr_module}{Aberration module.}
}
\value{
List object containing maxLik fit results.
}
\description{
Method based on the paper by Oliveira, M. et al. (2016). Zero-inflated
regression models for radiation-induced chromosome aberration data:
A comparative study. Biometrical Journal, 58(2), 259-279.
<doi:10.1002/bimj.201400233>.
}
