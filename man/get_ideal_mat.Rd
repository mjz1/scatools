% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bin_utils.R
\name{get_ideal_mat}
\alias{get_ideal_mat}
\alias{add_ideal_mat}
\title{Get ideal bin matrix}
\usage{
get_ideal_mat(
  mat,
  gc,
  n_freq,
  map,
  min_reads = 1,
  max_N_freq = 0.05,
  reads_outlier = 0.01,
  gc_outlier = 0.001,
  min_map = 0.9,
  ncores = 1,
  verbose = FALSE
)

add_ideal_mat(
  sce,
  assay_name = "counts",
  gc = rowData(sce)$gc,
  n_freq = rowData(sce)$n_freq,
  map = rowData(sce)$map,
  min_reads = 1,
  max_N_freq = 0.05,
  reads_outlier = 0.01,
  gc_outlier = 0.001,
  min_map = 0.9,
  ncores = 1,
  verbose = FALSE
)
}
\arguments{
\item{mat, sce}{A count matrix or SCE object depending on the function}

\item{gc}{Vector of gc content}

\item{n_freq}{Vector of bin N frequency (proportion of N bases in a bin)}

\item{map}{Vector of bin mappability}

\item{min_reads}{Minimum number of reads to consider a bin}

\item{max_N_freq}{Maximum allowable frequency of N bases to consider a bin. Range (0, 1)}

\item{reads_outlier}{Flag bins with reads in the top quantile given by this value. Range (0, 1)}

\item{gc_outlier}{Flag bins with GC content in the top and bottom quantule given by this value. Range (0, 1)}

\item{min_map}{Minimum allowable mappability score for a bin. Range (0, 1)}

\item{ncores}{number of cores for parallel evaluation (requires \code{pbmcapply} package)}

\item{verbose}{message verbosity}

\item{assay_name}{Name of assay}
}
\value{
Boolean matrices of ideal and valid bins

SCE object with ideal and valid boolean matrices
}
\description{
Given a matrix of bin counts, bin gc and N frequency, and filtering parameters, return a boolean matrix flagging ideal bins
}
