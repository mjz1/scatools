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
  which <- match.arg(arg = which, choices = c("bins", "cells"), several.ok = TRUE)


  gc_keeps <- rowRanges(sce)$gc > gc_range[1] & rowRanges(sce)$gc < gc_range[2]
  logger::log_info("Removing bins outside of {paste(gc_range, collapse = ' - ')} gc proportion: {sum(!gc_keeps)} bins removed")

  sce <- sce[gc_keeps, ]

  for (f in which) {
    if (f == "bins") {
      m <- 1
      bin_bool <- (colSums(apply(X = assay(sce, assay_name), MARGIN = m, function(X) X >= min_counts)) / ncol(sce)) >= min_prop

      # Flag
      rowData(sce)[[paste0("keep_", f)]] <- bin_bool

      logger::log_info("Keeping {sum(bin_bool)} of {length(bin_bool)} bins")

      if (!flag_only) {
        sce <- sce[bin_bool, ]
      }
    }

    if (f == "cells") {
      m <- 2
      cell_bool <- (colSums(apply(X = assay(sce, assay_name), MARGIN = m, function(X) X >= min_counts)) / nrow(sce)) >= min_prop

      # Flag
      colData(sce)[[paste0("keep_", f)]] <- cell_bool

      logger::log_info("Keeping {sum(cell_bool)} of {length(cell_bool)} cells")

      if (!flag_only) {
        sce <- sce[, which(cell_bool)]
      }
    }
  }

  return(sce)
}


#' Flag doublets
#'
#' This function is modified from `ArchR::filterDoublets` to enable flagging of doublets prior to removal.
#'
#' @param sce A `SingleCellExperiment` object.
#' @param cutEnrich The minimum numeric cutoff for `DoubletEnrichment`. This number is equivalent to the number of simulated
#' doublets identified as a nearest neighbor to the cell divided by the expected number given a random uniform distribution.
#' @param cutScore The minimum numeric cutoff for `DoubletScore` which represents the `-log10(binomial adjusted p-value)` for the `DoubletEnrichment`.
#' @param filterRatio The maximum ratio of predicted doublets to filter based on the number of pass-filter cells.
#' For example, if there are 5000 cells, the maximum would be `filterRatio * 5000^2 / (100000)` (which simplifies to `filterRatio * 5000 * 0.05`).
#' This `filterRatio` allows you to apply a consistent filter across multiple different samples that may have different
#' percentages of doublets because they were run with different cell loading concentrations.
#' The higher the `filterRatio`, the greater the number of cells potentially removed as doublets.
#' @param remove Logical. Whether to remove doublets from the object.
#'
#'
#' @export
flagDoublets <- function(sce, cutEnrich = 1, cutScore = -Inf, filterRatio = 1, remove = FALSE) {
  sce$doublet <- NA

  df <- colData(sce)[, c("Sample", "DoubletEnrichment", "DoubletScore")]
  splitDF <- split(seq_len(nrow(df)), as.character(df$Sample))

  cellsFilter <- lapply(splitDF, function(y) {
    x <- df[y, , drop = FALSE]

    n <- nrow(x)

    x <- x[order(x$DoubletEnrichment, decreasing = TRUE), ]

    if (!is.null(cutEnrich)) {
      x <- x[which(x$DoubletEnrichment >= cutEnrich), ]
    }

    if (!is.null(cutScore)) {
      x <- x[which(x$DoubletScore >= cutScore), ]
    }

    if (nrow(x) > 0) {
      head(rownames(x), filterRatio * n * (n / 100000))
    } else {
      NULL
    }
  }) %>% unlist(use.names = FALSE)

  logger::log_info("Identified {length(cellsFilter)} doublet cells across all samples")
  tabRemove <- table(df[cellsFilter, ]$Sample)
  tabAll <- table(df$Sample)
  samples <- unique(df$Sample)
  for (i in seq_along(samples)) {
    if (!is.na(tabRemove[samples[i]])) {
      message("\t", samples[i], " : ", tabRemove[samples[i]], " of ", tabAll[samples[i]], " (", round(100 * tabRemove[samples[i]] / tabAll[samples[i]], 1), "%)")
    } else {
      message("\t", samples[i], " : ", 0, " of ", tabAll[samples[i]], " (0%)")
    }
  }

  if (length(cellsFilter) > 0) {
    sce[, colnames(sce) %in% cellsFilter]$doublet <- TRUE
    sce[, !colnames(sce) %in% cellsFilter]$doublet <- FALSE
  } else {
    sce$doublet <- FALSE
  }

  if (remove) {
    logger::log_info("Removing doublets!")
    sce <- sce[, !sce$doublet]
  } else {
    logger::log_warn("Doublets flagged but not removed!")
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
