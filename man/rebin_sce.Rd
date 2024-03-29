% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rebin_sce.R
\name{rebin_sce}
\alias{rebin_sce}
\title{Rebin \code{SingleCellExperiment} object}
\usage{
rebin_sce(sce, assays, new_bins, cell_ids = NULL, ncores = 1)
}
\arguments{
\item{sce}{\code{SingleCellExperiment} object}

\item{assays}{Names of assays to include in the rebinning process}

\item{new_bins}{GenomicRanges of new bins}

\item{cell_ids}{Cells to include in rebinning. If \code{NULL} defaults to all cells}

\item{ncores}{Number of cores to speed up rebinning}
}
\value{
Rebinned \code{SingleCellExperiment} object
}
\description{
Will rebin data into larger bins. Have not tested with smaller bins.
}
