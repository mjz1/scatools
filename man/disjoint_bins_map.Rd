% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rebin_sce.R
\name{disjoint_bins_map}
\alias{disjoint_bins_map}
\title{Map between bins}
\usage{
disjoint_bins_map(x, y)
}
\arguments{
\item{x}{x Granges}

\item{y}{y Granges}
}
\value{
data frame with bin mapping
}
\description{
Will return a mapping between \code{x} and \code{y} bins using the \code{disjoin} function. The results can be used to map bins together taking weighted mean of the \code{disjoint_width}
}
