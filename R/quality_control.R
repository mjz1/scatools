#' Bin and/or Cell-wise quality filtering
#'
#' Note -- this filters GC bins, then bins by counts, and then cells. So the denominator changes at each step.
#'
#' @param sce SingleCellExperiment object
#' @param assay_name Name of assay to pull counts from. Ideally raw counts
#' @param which Filter cells or bins
#' @param flag_only Only flag cells/bins, do not remove
#'
#' @inheritParams run_scatools
#' @return SingleCellExperiment object
#' @export
#'
filter_sce <- function(sce,
                       assay_name = "counts",
                       which = c("bins", "cells"),
                       min_bin_counts = 1,
                       min_bin_prop = 0.95,
                       min_cell_counts = min_bin_counts,
                       min_cell_prop = min_bin_prop,
                       flag_only = FALSE,
                       gc_range = c(-Inf, Inf)) {
  which <- match.arg(arg = which, choices = c("bins", "cells"), several.ok = TRUE)


  gc_keeps <- rowRanges(sce)$gc > gc_range[1] & rowRanges(sce)$gc < gc_range[2]
  cli::cli_alert_info("Removing bins outside of {paste(gc_range, collapse = ' - ')} gc proportion: {sum(!gc_keeps)} bins removed")

  sce <- sce[gc_keeps, ]

  for (f in which) {
    if (f == "bins") {
      # m <- 1
      # bin_bool <- (colSums(apply(X = assay(sce, assay_name), MARGIN = m, function(X) X >= min_counts)) / ncol(sce)) >= min_prop
      bin_bool <- (Matrix::rowSums(assay(sce, assay_name) >= min_bin_counts) / ncol(sce)) >= min_bin_prop

      # Flag
      rowData(sce)[[paste0("keep_", f)]] <- bin_bool

      cli::cli_alert_info("Keeping {sum(bin_bool)} of {length(bin_bool)} bins with at least {min_bin_counts} counts in {min_bin_prop*100}% of cells")

      if (!flag_only) {
        sce <- sce[bin_bool, ]
      }
    }

    if (f == "cells") {
      m <- 2
      # cell_bool <- (rowSums(apply(X = assay(sce, assay_name), MARGIN = m, function(X) X >= min_counts)) / nrow(sce)) >= min_prop
      cell_bool <- (Matrix::colSums(assay(sce, assay_name) >= min_cell_counts) / nrow(sce)) >= min_cell_prop

      # Flag
      colData(sce)[[paste0("keep_", f)]] <- cell_bool

      cli::cli_alert_info("Keeping {sum(cell_bool)} of {length(cell_bool)} cells with at least {min_cell_counts} counts in {min_cell_prop*100}% of bins")

      if (!flag_only) {
        sce <- sce[, which(cell_bool)]
      }
    }
  }

  # Record filter provenance so per-sample bin/cell drops are traceable
  # (which bins were dropped is also flagged in rowData(sce)$keep_bins).
  sce@metadata$filter_info <- list(
    gc_range = gc_range, gc_removed = sum(!gc_keeps),
    min_bin_counts = min_bin_counts, min_bin_prop = min_bin_prop,
    min_cell_counts = min_cell_counts, min_cell_prop = min_cell_prop,
    n_bins = nrow(sce), n_cells = ncol(sce)
  )

  return(sce)
}


#' Bin and cell data-quality metrics
#'
#' Computes per-cell and per-bin QC metrics (counts, median counts, mean
#' coverage, fraction of zero bins) from a count assay, stores a summary in
#' `metadata(sce)$qc_summary`, and optionally plots diagnostic histograms.
#'
#' @param sce SingleCellExperiment object
#' @param assay_name Name of the count assay to summarise
#' @param plot Logical: draw diagnostic histograms
#'
#' @return The `sce` with QC metrics added to `colData`/`rowRanges` and a
#'   summary in `metadata(sce)$qc_summary`.
#' @export
#'
do_qc <- function(sce, assay_name = "counts", plot = TRUE) {
  m <- assay(sce, assay_name)

  sce$counts_per_cell <- Matrix::colSums(m)
  sce$median_percell_counts <- MatrixGenerics::colMedians(m)
  rowRanges(sce)$counts_per_bin <- Matrix::rowSums(m)
  rowRanges(sce)$median_perbin_counts <- MatrixGenerics::rowMedians(m)
  rowRanges(sce)$mean_perbin_counts <- Matrix::rowMeans(m)
  # fraction of cells with zero counts in each bin (efficient on sparse input)
  rowRanges(sce)$frac_zero_perbin <- 1 - (Matrix::rowSums(m != 0) / ncol(sce))

  sce@metadata$qc_summary <- list(
    assay = assay_name,
    n_cells = ncol(sce),
    n_bins = nrow(sce),
    median_counts_per_cell = stats::median(sce$counts_per_cell),
    median_counts_per_bin = stats::median(rowRanges(sce)$counts_per_bin),
    mean_coverage_per_bin = mean(rowRanges(sce)$mean_perbin_counts),
    median_frac_zero_per_bin = stats::median(rowRanges(sce)$frac_zero_perbin)
  )

  if (plot) {
    p1 <- colData(sce) %>%
      as.data.frame() %>%
      ggplot(aes(x = counts_per_cell)) +
      geom_histogram(bins = 50)

    p2 <- colData(sce) %>%
      as.data.frame() %>%
      ggplot(aes(x = median_percell_counts)) +
      geom_histogram(bins = 50)

    p3 <- rowData(sce) %>%
      as.data.frame() %>%
      ggplot(aes(x = counts_per_bin)) +
      geom_histogram(bins = 50)

    p4 <- rowData(sce) %>%
      as.data.frame() %>%
      ggplot(aes(x = median_perbin_counts)) +
      geom_histogram(bins = 50)

    p5 <- rowData(sce) %>%
      as.data.frame() %>%
      ggplot(aes(x = mean_perbin_counts)) +
      geom_histogram(bins = 50)

    p6 <- rowData(sce) %>%
      as.data.frame() %>%
      ggplot(aes(x = frac_zero_perbin)) +
      geom_histogram(bins = 50)

    title <- if (!is.null(sce$Sample)) paste(unique(sce$Sample), collapse = ", ") else NULL
    bw <- rowRanges(sce)$binwidth
    subtitle <- if (!is.null(bw)) glue::glue("{prettyMb(getmode(bw))} bins") else glue::glue("{nrow(sce)} bins")

    pcomb <- patchwork::wrap_plots(list(p1, p2, p3, p4, p5, p6), ncol = 2) +
      patchwork::plot_annotation(title = title, subtitle = subtitle)
    print(pcomb)
  }

  return(sce)
}


#' GC modal QC filter
#'
#' Removes cells that have high numbers of NA bins after GC modal correction
#'
#' @param sce SingleCellExperiment object
#' @param assay Assay containing GC modal corrected counts
#' @param filter_prop Filter proportion to keep cells with less then X% NA
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

  cli::cli_alert_info("Keeping {length(keep_cells)} of {ncol(sce)} cells: {signif(length(keep_cells) / ncol(sce) * 100, 3)}% (filter threshold: {filter_prop})")

  return(sce[, keep_cells])
}
