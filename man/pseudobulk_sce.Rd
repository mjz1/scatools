% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{pseudobulk_sce}
\alias{pseudobulk_sce}
\title{Pseudobulk cell CNA profiles}
\usage{
pseudobulk_sce(sce, assay_name, group_var = "all", statistics = "mean")
}
\arguments{
\item{sce}{SingleCellExperiment Object}

\item{assay_name}{Name of the assay}

\item{group_var}{Grouping variable to pseudobulk across. Default: "all"}

\item{statistics}{Character vector specifying the type of statistics to be computed, see \link[scuttle:summarizeAssayByGroup]{scuttle::summarizeAssayByGroup}.}
}
\value{
SingleCellExperiment object
}
\description{
Pseudobulk cell CNA profiles
}
