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

#' Plot multiple cell assays together
#'
#' @param sce SCE object
#' @param cell_id Cell ids
#' @param assays Assays to plot
#'
#' @return List of ggplot objects
#' @export
#'
plot_cell_multi <- function(sce, cell_id, assays) {

  plots <- vector(mode = "list")

  for (cell in cell_id) {
    cell_ps <- vector(mode = "list")
    for (assay in assays) {
      cell_p <- plot_cell_cna(sce = sce, cell_id = cell, assay_name = assay)
      cell_ps[[assay]] <- cell_p
    }
    plots[[cell]] <- patchwork::wrap_plots(cell_ps, ncol = 1) + patchwork::plot_annotation(title = cell)
  }

  # Arrange output
  return(plots)

}


#' Plot copy number heatmap
#'
#' @param sce A SingleCellExperiment object
#' @param assay_name Name of the assay to plot
#' @param cell_order Optional: Order of the cells
#' @param log2 Logical: Log2 transform the matrix prior to plotting
#' @param scale One of `'cells', 'bins', 'both' or 'none'`. Determines what kind of scaling is done.
#' @param clustering_results Clustering results to provide to inform cell ordering and cluster labelling. From [perform_umap_clustering]
#' @param col_fun Color mapping function from [circlize::colorRamp2()]
#' @param legend_name Name of the legend
#' @param ... Additional parameters that can be passed to [ComplexHeatmap::Heatmap()]
#'
#' @return A heatmap
#' @export
#'
cnaHeatmap <- function(sce, assay_name = "state", cell_order = NULL, log2 = FALSE, scale = c("none", "cells", "bins", "both"), clustering_results = NULL, col_fun = NULL, legend_name = assay_name, ...) {

  # TODO: seperate out the clustering to be containing within the sce object and allow the user to pass specified ordering or clusters
  cn_mat <- as.matrix(assay(sce, assay_name))

  cn_mat <- scale_mat(cn_mat, log2 = log2, scale = scale)

  sce <- sce[rownames(cn_mat),colnames(cn_mat)]


  if (is.null(cell_order) & is.null(clustering_results)) {

    clustering_results <- perform_umap_clustering(cn_matrix = cn_mat)

    ordered_cell_ids <- clustering_results$clustering[order(clustering_results$clustering$clone_id), "cell_id"]
    cnv_clusters <- clustering_results$clustering[order(clustering_results$clustering$clone_id), "clone_id"]

  } else if (!is.null(clustering_results)) {

    ordered_cell_ids <- clustering_results$clustering[order(clustering_results$clustering$clone_id), "cell_id"]
    cnv_clusters <- clustering_results$clustering[order(clustering_results$clustering$clone_id), "clone_id"]

  } else {
    if (!all(cell_order %in% rownames(cn_mat))) {
      logger::log_error("Provided cell ordering contains cells not contained in the input")
    }
    ordered_cell_ids <- cell_order
  }

  # Reorder cells
  cn_mat <- cn_mat[,ordered_cell_ids]



  if (grepl("state", assay_name)) {
    cn_mat[cn_mat >= 11] <- "11+"
    cn_colors <- state_cn_colors()
  } else if (!is.null(col_fun)) {
    cn_colors <- col_fun
  } else {
    # Attempt some intelligent color mapping
    cn_colors <- NULL
  }

  # Split columns by chromosome
  chrs <- as.vector(gsub("chr", "", GenomeInfoDb::seqnames(SummarizedExperiment::rowRanges(sce))))
  col_split <- factor(chrs, levels = unique(gtools::mixedsort(chrs)))

  suppressMessages(ht_plot <- ComplexHeatmap::Heatmap(
    name = legend_name,
    matrix = t(cn_mat),
    show_row_names = FALSE,
    cluster_rows = FALSE,
    col = cn_colors,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    na_col = "white",
    use_raster = TRUE,
    raster_quality = 10,
    column_split = col_split,
    ...
  ))

  return(ht_plot)
}
