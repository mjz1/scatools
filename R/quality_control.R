#' Bin and/or Cell-wise quality filtering
#'
#' Note -- this filters GC bins, then bins by counts, and then cells. So the denominator changes at each step.
#'
#' @param sce SingleCellExperiment object
#' @param assay_name Name of assay to pull counts from. Ideally raw counts
#' @param filt Filter cells or bins
#' @param min_counts Minimum counts per cell/bin
#' @param min_prop Minimum proportion of cells/bins with `min_counts`
#' @param flag_only Only flag cells/bins, do not remove
#' @param gc_range Min and max GC values per bin to keep

#'
#' @return SingleCellExperiment object
#' @export
#'
filter_sce <- function(sce, assay_name = "counts", which = c("bins", "cells"), min_counts = 1, min_prop = 0.95, flag_only = FALSE, gc_range = c(-Inf, Inf)) {

  which = match.arg(arg = which, choices = c("bins", "cells"), several.ok = TRUE)


  gc_keeps <- rowRanges(sce)$gc> gc_range[1] & rowRanges(sce)$gc < gc_range[2]
  logger::log_info("Removing bins outside of {paste(gc_range, collapse = ' - ')} gc proportion: {sum(!gc_keeps)} bins removed")

  sce <- sce[gc_keeps,]

  for (f in which) {
    if (f == "bins") {
      m = 1
      bin_bool <- (colSums(apply(X = assay(sce, assay_name), MARGIN = m, function(X) X >= min_counts)) / ncol(sce)) >= min_prop

      # Flag
      rowData(sce)[[paste0("keep_", f)]] <- bin_bool

      logger::log_info("Keeping {sum(bin_bool)} of {length(bin_bool)} bins")

      if (!flag_only) {
        sce <- sce[bin_bool,]
      }
    }

    if (f == "cells") {
      m = 2
      cell_bool <- (colSums(apply(X = assay(sce, assay_name), MARGIN = m, function(X) X >= min_counts)) / nrow(sce)) >= min_prop

      # Flag
      colData(sce)[[paste0("keep_", f)]] <- cell_bool

      logger::log_info("Keeping {sum(cell_bool)} of {length(cell_bool)} cells")

      if (!flag_only) {
        sce <- sce[,which(cell_bool)]
      }

    }


  }

  return(sce)

}


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

