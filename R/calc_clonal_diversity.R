#' Calculate clonal diversity
#'
#' @param sce sce
#' @param sample_column sample
#' @param clone_name clone
#' @param index diversity index to use. See [vegan::diversity]
#'
#' @return named list
#' @export
#'
calc_clonal_diversity <- function(sce, sample_column, clone_name, index = "shannon") {
  clone_counts <- table(colData(sce)[[sample_column]], colData(sce)[[clone_name]])

  res <- vegan::diversity(clone_counts, index = index, MARGIN = 1)
}
