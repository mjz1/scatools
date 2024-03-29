% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_gene_copy.R
\name{get_gene_copy}
\alias{get_gene_copy}
\title{Get Gene Copies}
\usage{
get_gene_copy(
  sce,
  assay_name,
  group_var = "all",
  gain_cutoff = 0.75,
  loss_cutoff = -gain_cutoff,
  clonal_thr = 0.75,
  gene_names = NULL,
  cancer_genes = TRUE
)
}
\arguments{
\item{sce}{Single cell experiment object}

\item{assay_name}{Name of the assay from which to pull copy number}

\item{group_var}{Grouping variable (defaults to all)}

\item{gain_cutoff}{Relative copy increase to be considered a copy gain}

\item{loss_cutoff}{Relative copy decrease to be considered a copy loss}

\item{clonal_thr}{Proportion of clones harboring a gain or loss to be considered clonal}

\item{gene_names}{List of genes to get copy number}

\item{cancer_genes}{Logical. Whether to include results for cancer genes. Default=TRUE.}
}
\value{
A list with three elements:
\describe{
\item{copy}{Wide form data frame with one row per  gene, and columns reflecting the average per clone copy, the relative copy increase or decrease versus clone average copy, clone average copy}
\item{copy_long}{A melted version of \code{copy}}
\item{clonal}{A filtered version of \code{copy_long} that classifies each mutation as private, shared, or clonal, and whether the type of gain or loss matches the gene's oncogene or tumor suppressor type.}
}
}
\description{
Calculate per group copy number for a list of genes.
}
