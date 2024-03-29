% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p-curve.R
\name{p_gap}
\alias{p_gap}
\alias{p_gap_inference}
\title{Calculate the largest "gap" in a Young's p-value plot}
\usage{
p_gap(p_values)

p_gap_inference(studies, threshold = 0.125)
}
\arguments{
\item{p_values}{A vector of p-values}

\item{studies}{A dataframe of studies, as generated by \code{draw_studies}}

\item{threshold}{Threshold for concluding that a collection of p-values is "gappy"}
}
\value{
\itemize{
\item \code{p_gap()} returns a length-one numeric, the largest gap
}

\itemize{
\item \code{p_gap_inference()} returns a dataframe with columns:
\itemize{
\item \code{gap} Numeric; the largest gap
\item \code{gappy} Whether the plot is considered gappy; either \code{"gappy"} or \code{"not gappy"}.
}
}
}
\description{
The "gap" in Young's p-value plot is the greatest difference between sequential p-values.  A collection of p-values is considered "gappy" if the largest gap is larger than \code{threshold}.
}
