#' GC modal QC filter
#'
#' Removes cells that have high numbers of NA bins after GC modal correction
#'
#' @param sce SingleCellExperiment object
#' @param assay Assay containing GC modal corrected counts
#' @param filter_prop Filter proportion
#' @param verbose Logical. Message verbosity
#'
#' @return Filtered SingleCellExperiment
#' @export
#'
gc_modal_qc_filter <- function(sce, assay = "counts_gc_modal", filter_prop = 0.05, verbose = TRUE) {
  prop_na <- colSums(apply(assay(sce, assay), 2, is.na)) / nrow(sce)

  sce[["prop_na"]] <- prop_na

  keep_cells <- which(prop_na <= filter_prop)

  sce@metadata$cell_filter_info <- list(cells_in = ncol(sce), cells_kept = length(keep_cells), filter_threshold = filter_prop)

  logger::log_info("Keeping {length(keep_cells)} of {ncol(sce)} cells: {signif(length(keep_cells) / ncol(sce) * 100, 3)}% (filter threshold: {filter_prop})")

  return(sce[, keep_cells])
}
