% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quality_control.R
\name{gc_modal_qc_filter}
\alias{gc_modal_qc_filter}
\title{GC modal QC filter}
\usage{
gc_modal_qc_filter(
  sce,
  assay = "counts_gc_modal",
  filter_prop = 0.05,
  verbose = TRUE
)
}
\arguments{
\item{sce}{SingleCellExperiment object}

\item{assay}{Assay containing GC modal corrected counts}

\item{filter_prop}{Filter proportion to keep cells with less then X\% NA}

\item{verbose}{Logical. Message verbosity}
}
\value{
Filtered SingleCellExperiment
}
\description{
Removes cells that have high numbers of NA bins after GC modal correction
}
