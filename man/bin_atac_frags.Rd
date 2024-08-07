% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bin_utils.R
\name{bin_atac_frags}
\alias{bin_atac_frags}
\title{Bin atac fragments}
\usage{
bin_atac_frags(
  sample_id,
  fragment_file,
  cells = NULL,
  bins,
  bin_name = NULL,
  blacklist = NULL,
  outdir,
  ncores = 1,
  bpparam = BiocParallel::SerialParam(),
  overwrite = FALSE,
  return_mat = FALSE
)
}
\arguments{
\item{sample_id}{Sample id}

\item{fragment_file}{Fragments file}

\item{cells}{A vector of cell barcodes. If not provided will use all barcodes}

\item{bins}{Genomic ranges object with bins to use}

\item{bin_name}{bin_name Name of the bins (e.g. \code{'10Mb'}, \code{'500Kb'}, \code{'chr_arm'}). If not provided is automatically detected based on binwidth.}

\item{blacklist}{Genome blacklist regions to filter against}

\item{outdir}{Output directory}

\item{ncores}{Number of cores}

\item{bpparam}{BiocParallel Parameters}

\item{overwrite}{Whether to overwrite previous binned counts}

\item{return_mat}{Logical. Return the binned depth matrix (default = FALSE)}
}
\value{
If \code{return_mat=TRUE}, returns a sparse binned depth matrix. Otherwise returns \code{NULL}, but will save the results to a 10X style market matrix \code{outdir}
}
\description{
Convenience wrapper to bin fragments of a given size and save them as \code{.mtx} files
}
