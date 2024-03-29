% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{scale_sub}
\alias{scale_sub}
\title{Scale and subset SCE assays}
\usage{
scale_sub(
  sce,
  assay_name = "counts",
  log2 = FALSE,
  scale = "none",
  verbose = FALSE,
  new_assay = NULL,
  center = FALSE
)
}
\arguments{
\item{sce}{SingleCellExperiment}

\item{assay_name}{Assay to transform}

\item{log2}{Logical: log transform}

\item{scale}{One of :c("none", "cells", "bins", "both"). Specifies how scaling should be dome}

\item{verbose}{Message verbosity}

\item{new_assay}{New assay name}

\item{center}{Center the matrix}
}
\value{
A single cell experiment object
}
\description{
Scale and subset SCE assays
}
