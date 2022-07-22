
#' @export
plot_cell_cna <- function(sce, cell_id, assay_name = "counts") {
  makePerCellDF(sce[, cell_id], features = rownames(sce), assay.type = assay_name, use.coldata = F) %>%
    rownames_to_column(var = "barcode") %>%
    pivot_longer(cols = starts_with("chr"), names_to = "bin_id", values_to = "counts") %>%
    left_join(as.data.frame(bins)) %>%
    ggplot() +
    geom_point(aes(x = start, y = counts), size = 0.5) +
    facet_grid(barcode ~ seqnames, scales = "free_x", space = "free_x") +
    theme(panel.spacing = unit(0, "lines"), 
          panel.border = element_rect(fill = NA), 
          axis.text.x = element_blank())
}
