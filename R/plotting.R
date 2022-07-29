
#' @export
plot_cell_cna <- function(sce, cell_id = NULL, assay_name = "counts") {
  if (is.null(cell_id)) {
    if (length(cell_id) > 20) {
      logger::log_warn("No cell ids provided and plotting many cells. Are you sure you want to do this!?")
    }
    cell_id = colnames(sce)
  }

  # TODO: Try to intelligently find the start positions from the provided sce
  bindat <- as.data.frame(SummarizedExperiment::rowRanges(sce))
  if (nrow(bindat) == 0) {
    logger::log_warn("Bin data not found in rowRanges(sce). Checking rowData(sce)")
    bindat <- as.data.frame(SummarizedExperiment::rowData(sce))
    if (nrow(bindat) == 0) {
      logger::log_error("Bin data not found in rowData(sce).")
    }
  } else {
    # This means we do have rowranges.
    # make sure seqnames are labelled as chr
    bindat$chr <- bindat$seqnames
  }
  if (is.null(rownames(sce))) {
    logger::log_error("rownames(sce) cannot be NULL.")
  }
  bindat$bin_id <- rownames(sce)

  # TODO: subset for main chromosomes and reorder)

  plot_dat <- makePerCellDF(sce[, cell_id], features = rownames(sce), assay.type = assay_name, use.coldata = F) %>%
    rownames_to_column(var = "barcode") %>%
    pivot_longer(cols = -dplyr::starts_with("barcode"), names_to = "bin_id", values_to = "counts") %>%
    left_join(as.data.frame(bindat))

  # Keep cell ordering as provided
  plot_dat$barcode <- factor(plot_dat$barcode, levels = cell_id)

  ggplot(plot_dat) +
    geom_point(aes(x = start, y = counts), size = 0.5) +
    facet_grid(barcode ~ chr, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = assay_name) +
    theme(
      panel.spacing = unit(0, "lines"),
      panel.border = element_rect(fill = NA),
      axis.text.x = element_blank()
    )
}
