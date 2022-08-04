#' Plot Cell Copy Number
#'
#' @param sce sce object
#' @param cell_id Vector of cell ids
#' @param assay_name Name of assay to plot
#' @param col_fun Color mapping function from [circlize::colorRamp2()]
#'
#' @return A ggplot object
#' @export
#'
plot_cell_cna <- function(sce, cell_id = NULL, assay_name = "counts", col_fun = NULL) {
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

  plot_dat <- scuttle::makePerCellDF(sce[, cell_id], features = rownames(sce), assay.type = assay_name, use.coldata = F) %>%
    tibble::rownames_to_column(var = "barcode") %>%
    tidyr::pivot_longer(cols = -dplyr::starts_with("barcode"), names_to = "bin_id", values_to = "counts") %>%
    dplyr::left_join(as.data.frame(bindat), by = "bin_id") %>%
    filter(!is.na(counts))

  plot_dat$chr_no <- gsub("chr", "", plot_dat$chr)
  plot_dat$chr_no <- factor(plot_dat$chr_no, levels = unique(chr_reorder(plot_dat$chr_no)))

  # Keep cell ordering as provided
  plot_dat$barcode <- factor(plot_dat$barcode, levels = cell_id)

  # check for assay type and add colors
  if (grepl("state", assay_name) ){
    cn_colors = state_cn_colors()
    # Clip counts to 11
    plot_dat[plot_dat$counts > 11,'counts'] <- 11
    plot_dat$counts <- as.factor(plot_dat$counts)
    plot_dat$counts <- dplyr::recode_factor(plot_dat$counts, "11" = "11+")

    p <- ggplot(plot_dat) +
      geom_point(aes(x = start, y = counts, color = counts), size = 0.5) +
      facet_grid(barcode ~ chr_no, scales = "free_x", space = "free_x") +
      labs(x = NULL, y = assay_name, color = assay_name) +
      theme(
        panel.spacing = unit(0, "lines"),
        panel.border = element_rect(fill = NA),
        axis.text.x = element_blank()
      ) +
      scale_color_manual(values = cn_colors) +
      guides(colour = guide_legend(override.aes = list(size=3)))

  } else if (!is.null(col_fun)) {
    cn_colors = col_fun
    # This is hacky -- maybe there is a more direct way to pass color pallete from col_fun

    p <- ggplot(plot_dat) +
      geom_point(aes(x = start, y = counts, color = counts), size = 0.5) +
      facet_grid(barcode ~ chr_no, scales = "free_x", space = "free_x") +
      labs(x = NULL, y = assay_name, color = assay_name) +
      theme(
        panel.spacing = unit(0, "lines"),
        panel.border = element_rect(fill = NA),
        axis.text.x = element_blank()
      ) +
      scale_color_gradient2(low = attr(cn_colors, "colors")[1],
                            mid = attr(cn_colors, "colors")[2],
                            high = attr(cn_colors, "colors")[3],
                            limits = c(attr(cn_colors, "breaks")[1], attr(cn_colors, "breaks")[3]))

  } else {
    # Attempt some intelligent color mapping
    cn_colors = NULL

    p <- ggplot(plot_dat) +
      geom_point(aes(x = start, y = counts), size = 0.5) +
      facet_grid(barcode ~ chr_no, scales = "free_x", space = "free_x") +
      labs(x = NULL, y = assay_name, color = assay_name) +
      theme(
        panel.spacing = unit(0, "lines"),
        panel.border = element_rect(fill = NA),
        axis.text.x = element_blank()
      )
  }

  return(p)
}

}
