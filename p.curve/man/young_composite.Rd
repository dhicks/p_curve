% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p-curve.R
\name{young_composite}
\alias{young_composite}
\title{A "composite" version of Young's p-value plot}
\usage{
young_composite(sims, alpha = 0.5, ...)
}
\arguments{
\item{sims}{A dataframe of runs of the simulation, as returned by \code{many_metas}.}

\item{alpha}{Transparency for lines; defaults to .5}

\item{...}{Other aesthetics, passed to \code{aes()}}
}
\value{
A \code{ggplot2} plot
}
\description{
This version of Young's p-value plot allows us to visualize the plot across all runs of a simulation. Each simulation run is represented by a line in the plot rather than a sequence of points.
}
