% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p-curve.R
\name{young_curve}
\alias{young_curve}
\title{Young's p-value plot}
\usage{
young_curve(studies, alpha = 0.05, draw_alpha = TRUE, ...)
}
\arguments{
\item{studies}{A dataframe of studies to be plotted, as returned by \code{draw_studies()}. Minimally, a dataframe with columns \code{rank} and \code{p.value}.}

\item{alpha}{Threshold for statistical significance; default (.05)}

\item{draw_alpha}{Draw a horizontal line at \code{alpha}?}

\item{...}{Other aesthetics, passed at the plot level}
}
\value{
A \code{ggplot2} plot
}
\description{
S. Stanley Young's "p-value plot" or "p-curve." Equivalent to a QQ-plot against a uniform distribution with a rescaled x-axis.
}
\details{
Young's p-curve plots the p-value of each study against its rank, by ascending p-value (1 = smallest p-value). It is equivalent to a QQ-plot against a uniform distribution, with the x-axis running from 1 to the number of studies rather than the quantiles \code{(1:N)/N}.
}
