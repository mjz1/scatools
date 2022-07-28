
#' @export
plot_cell_cna <- function(sce, cell_id = NULL, assay_name = "counts") {
  if (is.null(cell_id)) {
    logger::log_warn("No cell ids provided! Plotting all cells. Are you sure you want to do this!?")
    cell_id = colnames(sce)
  }

  # TODO: Try to intelligently find the start positions from the provided sce
  bindat <- as.data.frame(SummarizedExperiment::rowRanges(sce))
  bindat$bin_id <- rownames(sce)

  makePerCellDF(sce[, cell_id], features = rownames(sce), assay.type = assay_name, use.coldata = F) %>%
    rownames_to_column(var = "barcode") %>%
    pivot_longer(cols = dplyr::starts_with("chr"), names_to = "bin_id", values_to = "counts") %>%
    left_join(as.data.frame(bindat)) %>%
    ggplot() +
    geom_point(aes(x = start, y = counts), size = 0.5) +
    facet_grid(barcode ~ seqnames, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = assay_name) +
    theme(
      panel.spacing = unit(0, "lines"),
      panel.border = element_rect(fill = NA),
      axis.text.x = element_blank()
    )
}
