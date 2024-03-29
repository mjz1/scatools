% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bin_utils.R
\name{bin_frags_chr}
\alias{bin_frags_chr}
\title{Bin scATAC fragments}
\usage{
bin_frags_chr(chr, bins, blacklist = NULL, ArrowFile)
}
\arguments{
\item{chr}{A single chromsome to compute depth information}

\item{bins}{A list of bins (can include all chromosomes)}

\item{ArrowFile}{Path to an ArrowFile generated by \code{ArchR}}
}
\value{
Sparse matrix of binned fragment counts
}
\description{
\code{bin_frags_chr} computes the fragments across bins in a single chromosome from an ArchR ArrowFile
}
\examples{
\dontrun{
dp_mat <- bin_frags_chr(chr = "chr1", bins = get_chr_arm_bins("hg38"), ArrowFile = ArrowFile)
}
}
