% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p-curve.R
\name{p_value}
\alias{p_value}
\title{Calculate p-value}
\usage{
p_value(h, test_output, dataf, verbose = TRUE)
}
\arguments{
\item{h}{Bare (unquoted) expression for H0 (eg, \code{delta == 0.2})}

\item{test_output}{Bare (unquoted) expression for the test output (eg, \code{aic_comp == 'quadratic'})}

\item{dataf}{Data frame of simulation results, as returned by \code{many_metas()}}

\item{verbose}{Message \code{h} and \code{test_output}}
}
\value{
Dataframe with columns
\item{h}{Expression for H0, as a string}
\item{test_output}{Expression for the test output, as a string}
\item{n_false}{Number of rows post-filtering where \code{test_output} is true}
\item{n_true}{Number of rows post-filtering where \code{test_output} is false}
\item{p}{Share of rows \code{n_true/(n_false + n_true)}}
}
\description{
Calculate the p-value, as the share of simulation runs in which \code{test_output} is true, given the condition \code{h}
}
