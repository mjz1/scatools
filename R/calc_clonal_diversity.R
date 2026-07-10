#' Calculate clonal diversity
#'
#' @param sce sce
#' @param sample_column sample
#' @param clone_name clone
#' @param index diversity index to use. One of `"shannon"`, `"simpson"`, or
#'   `"invsimpson"`.
#'
#' @return named list
#' @export
#'
calc_clonal_diversity <- function(sce, sample_column, clone_name, index = "shannon") {
  clone_counts <- table(colData(sce)[[sample_column]], colData(sce)[[clone_name]])

  res <- apply(clone_counts, MARGIN = 1, FUN = diversity_index, index = index)
}

#' Ecological diversity index of a vector of counts
#'
#' Base-R equivalent of `vegan::diversity()` for a single sample.
#'
#' @param x Numeric vector of per-category counts (or proportions)
#' @param index One of `"shannon"`, `"simpson"`, or `"invsimpson"`
#' @return A single numeric diversity value
#' @noRd
diversity_index <- function(x, index = c("shannon", "simpson", "invsimpson")) {
  index <- match.arg(index)
  x <- x[x > 0]
  p <- x / sum(x)
  switch(index,
    shannon = -sum(p * log(p)),
    simpson = 1 - sum(p^2),
    invsimpson = 1 / sum(p^2)
  )
}
