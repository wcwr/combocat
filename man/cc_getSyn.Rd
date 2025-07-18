% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cc_getSyn.R
\name{cc_getSyn}
\alias{cc_getSyn}
\title{Calculate Synergy Scores}
\usage{
cc_getSyn(norm_data, dr_data, conc_rounding_factor = 6)
}
\arguments{
\item{norm_data}{List. Output from `cc_norm`, containing normalized data and metadata.}

\item{dr_data}{List. Output from `cc_getDR`, containing dose-response and IC50 data.}

\item{conc_rounding_factor}{Numeric. Number of decimal places to round concentrations to. Default is `6`. 
Rounding is important for merging drugs with slightly different concentrations due to floating-point errors
This is particularly important due to floating point errors that may arise from `cc_makeMeta`}
}
\description{
This function calculates Bliss and Loewe synergy scores for both dense and sparse mode data, based on normalized input data.
The results are added to the reference data frame.
}
