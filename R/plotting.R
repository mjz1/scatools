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
    cell_id <- colnames(sce)
  }

  if (is.numeric(cell_id)) {
    cell_id <- colnames(sce[,cell_id])
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
    rownames(sce) <- with(bindat, paste(chr, start, end, sep = "_"))
    # logger::log_error("rownames(sce) cannot be NULL.")
  }
  bindat$bin_id <- rownames(sce)

  # TODO: subset for main chromosomes and reorder)

  plot_dat <- scuttle::makePerCellDF(sce[, cell_id], features = rownames(sce), assay.type = assay_name, use.coldata = FALSE, use.dimred = FALSE, use.altexps = FALSE) %>%
    tibble::rownames_to_column(var = "barcode") %>%
    tidyr::pivot_longer(cols = -dplyr::starts_with("barcode"), names_to = "bin_id", values_to = "counts") %>%
    dplyr::left_join(as.data.frame(bindat), by = "bin_id") %>%
    dplyr::filter(!is.na(counts))

  plot_dat$chr_no <- gsub("chr", "", plot_dat$chr)
  plot_dat$chr_no <- factor(plot_dat$chr_no, levels = unique(chr_reorder(plot_dat$chr_no)))

  # Keep cell ordering as provided
  plot_dat$barcode <- factor(plot_dat$barcode, levels = cell_id)

  # Base plot
  base_p <- ggplot(plot_dat) +
    geom_segment(aes(x = start, xend = end, y = counts, yend = counts), size = 1) +
    facet_grid(barcode ~ chr_no, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = assay_name, color = assay_name) +
    theme_bw() +
    theme(
      panel.spacing = unit(0, "lines"),
      panel.border = element_rect(fill = NA),
      axis.text.x = element_blank()
    )

  # check for assay type and add colors
  if (grepl("state", assay_name)) {
    cn_colors <- state_cn_colors()
    # Clip counts to 11
    plot_dat[plot_dat$counts > 11, "counts"] <- 11
    plot_dat$counts <- as.factor(plot_dat$counts)
    plot_dat$counts <- dplyr::recode_factor(plot_dat$counts, "11" = "11+")

    p <- base_p +
      scale_color_manual(values = cn_colors) +
      guides(colour = guide_legend(override.aes = list(size = 3)))
  } else if (!is.null(col_fun)) {
    cn_colors <- col_fun
    # This is hacky -- maybe there is a more direct way to pass color pallete from col_fun

    p <- base_p  +
      scale_color_gradient2(
        low = attr(cn_colors, "colors")[1],
        mid = attr(cn_colors, "colors")[2],
        high = attr(cn_colors, "colors")[3],
        limits = c(attr(cn_colors, "breaks")[1], attr(cn_colors, "breaks")[3])
      )
  } else {
    # Attempt some intelligent color mapping
    cn_colors <- NULL

    p <- base_p
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

  if (is.numeric(cell_id)) {
    cell_id <- colnames(sce[,cell_id])
  }

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
#' @param clone_name Name of clone_id column in sce object
#' @param cell_order Optional: Order of the cells
#' @param log2 Logical: Log2 transform the matrix prior to plotting
#' @param scale One of `'cells', 'bins', 'both' or 'none'`. Determines what kind of scaling is done.
#' @param clustering_results Clustering results to provide to inform cell ordering and cluster labelling. From [perform_umap_clustering]
#' @param col_fun Color mapping function from [circlize::colorRamp2()]
#' @param col_clones Optional: A named vector of clone colors.
#' @param legend_name Name of the legend
#' @param clust_annot Annotate cluster and sample labels
#' @param verbose Logical: Message verbosity
#' @param ... Additional parameters that can be passed to [ComplexHeatmap::Heatmap()]
#'
#' @return A heatmap
#' @export
#'
cnaHeatmap <- function(sce,
                       assay_name = "state",
                       clone_name = NULL,
                       cell_order = NULL,
                       cluster_rows = NULL,
                       log2 = FALSE,
                       scale = c("none", "cells", "bins", "both"),
                       clustering_results = NULL,
                       col_fun = NULL,
                       col_clones = NULL,
                       legend_name = assay_name,
                       clust_annot = TRUE,
                       verbose = TRUE,
                       ...) {

  if (is.null(rownames(sce))) {
    rownames(sce) <- 1:nrow(sce)
  }

  cn_mat <- as.matrix(assay(sce, assay_name))

  cn_mat <- scale_mat(cn_mat, log2 = log2, scale = scale)

  sce <- sce[rownames(cn_mat), colnames(cn_mat)]

  if (!is.null(clone_name)) {
    if (clone_name %in% colnames(colData(sce))) {
      # If clone data is present use it
      logger::log_info("Using {clone_name} as clones...")
    } else {
      logger::log_warn("{clone_name} not found in sce object. Reperforming clustering...")
      clone_name = "clone"
    }
  } else if (is.null(cell_order) & is.null(clustering_results)) {
    clone_name = "clone"

    clustering_results <- perform_umap_clustering(cn_matrix = cn_mat, verbose = verbose)

    # ordered_cell_ids <- clustering_results$clustering[order(clustering_results$clustering$clone_size, decreasing = TRUE), "cell_id"]
    # cnv_clusters <- clustering_results$clustering[order(clustering_results$clustering$clone_size, decreasing = TRUE), "clone_id"]

    sce[[clone_name]] <- clustering_results$clustering$clone_id[match(colnames(sce), clustering_results$clustering$cell_id)]

  } else if (!is.null(clustering_results)) {
    clone_name = "clone"

    sce[[clone_name]] <- clustering_results$clustering$clone_id[match(colnames(sce), clustering_results$clustering$cell_id)]

    # ordered_cell_ids <- clustering_results$clustering[order(clustering_results$clustering$clone_size, decreasing = TRUE), "cell_id"]
    # cnv_clusters <- clustering_results$clustering[order(clustering_results$clustering$clone_size, decreasing = TRUE), "clone_id"]
  } else {
    if (!all(cell_order %in% colnames(cn_mat))) {
      logger::log_error("Provided cell ordering contains cells not contained in the input")
    }
    ordered_cell_ids <- cell_order
  }

  # TODO: Separate the below into a separate function
  # First get the average per clone signal to order the clones properly
  avg_exp <- scuttle::summarizeAssayByGroup(sce,
                                            assay.type = assay_name,
                                            ids = sce[[clone_name]],
                                            statistics = "mean")
  rowRanges(avg_exp) <- rowRanges(sce)

  # Compute the distance based on correlation
  d <- as.dist(1-cor(assay(avg_exp, 'mean')))
  hc <- hclust(d, method = "complete")

  # Pull out the order from the hc object
  clone_order <- hc$labels[hc$order]

  # clone_order <- names(sort(table(sce[[clone_name]]),
  #                           decreasing = TRUE))

  cell_order <- order(match(sce[[clone_name]], clone_order))

  ordered_cell_ids <- as.character(colData(sce)[cell_order, 'cell_id'])

  cnv_clusters <- colData(sce)[cell_order, clone_name]

  row_split <- factor(cnv_clusters, levels = as.character(hc$labels[hc$order]))
  row_title <- hc$labels[hc$order]

  # Reorder cells
  cn_mat <- cn_mat[, ordered_cell_ids]

  if (class(clust_annot) == "HeatmapAnnotation") {
    left_annot <- clust_annot
  } else if (clust_annot) {
    if (is.null(col_clones)) {
      col_clones <- scales::hue_pal()(length(unique(cnv_clusters)))
      names(col_clones) <- levels(factor(cnv_clusters))
    }

    left_annot <- ComplexHeatmap::HeatmapAnnotation(
      Clone = cnv_clusters,
      # Sample = sce[, ordered_cell_ids]$Sample,
      col = list(Clone = col_clones),
      which = "row",
      show_legend = c(TRUE, FALSE)
    )
  } else {
    left_annot <- NULL
  }



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

  if (is.null(cluster_rows)) {
    cluster_rows <- FALSE
  } else if (class(cluster_rows) == "dendrogram") {
    row_split = length(levels(cnv_clusters))
    row_title = hc$labels[hc$order]
  }

  suppressMessages(ht_plot <- ComplexHeatmap::Heatmap(
    name = legend_name,
    matrix = t(cn_mat),
    show_row_names = FALSE,
    col = cn_colors,
    cluster_columns = FALSE,
    cluster_rows = cluster_rows,
    show_column_names = FALSE,
    na_col = "white",
    use_raster = TRUE,
    raster_quality = 10,
    column_split = col_split,
    left_annotation = left_annot,
    row_split = row_split,
    row_title = row_title,
    ...
  ))

  return(ht_plot)
}

#' @export
cloneCnaHeatmap <- function(sce, assay_name = "counts", clone_name = NULL, scale = c("none", "cells", "bins", "both"), log2 = FALSE, clustering_results = NULL, clust_lab = TRUE, ...) {

  if (is.null(rownames(sce))) {
    rownames(sce) <- 1:nrow(sce)
  }

  orig_sce <- sce

  clust_mat <- scale_mat(assay(sce, assay_name), scale = scale, log2 = log2)

  # Subset the SCE to the same dimensions
  sce <- sce[rownames(clust_mat), colnames(clust_mat)]

  new_assay <- paste(assay_name, "avg", sep = "_")

  assay(sce, new_assay) <- clust_mat

  # if (!is.null(clone_name)) {
  #
  # }

  # if (is.null(clustering_results)) {
  #   clustering_results <- perform_umap_clustering(assay(sce, new_assay))
  # }

  # sce$clone_id <- clustering_results$clustering[match(sce$Barcode, clustering_results$clustering$cell_id), "clone_id"]

  avg_exp <- scuttle::summarizeAssayByGroup(sce, assay.type = new_assay, ids = sce[[clone_name]], statistics = "mean")
  rowRanges(avg_exp) <- rowRanges(sce)

  # Order by clone size
  avg_exp <- avg_exp[, order(avg_exp$ncells, decreasing = TRUE)]

  if (grepl("state", new_assay)) {
    # Round to integers
    assay(avg_exp, new_assay) <- round(assay(avg_exp, "mean"))
  } else {
    assay(avg_exp, new_assay) <- assay(avg_exp, "mean")
  }

  # if (clust_lab) {
  #   clust_lab <- ComplexHeatmap::row_anno_text(avg_exp$ids, rot = 90, just = "center")
  # } else {
  #   clust_lab <- NULL
  # }

  # left_annot <- ComplexHeatmap::HeatmapAnnotation(
  #   # ClusterLab = clust_lab,
  #   Clone = avg_exp$ids,
  #   which = "row",
  #   show_legend = c(FALSE, FALSE)
  # )


  metadata(orig_sce)[[new_assay]] <- avg_exp

  avg_exp$cell_id <- avg_exp$ids

  # TODO: Figure out how to handle returning the avg exp data and plot
  # Ideally want to replace the SCE in the parent environment while also returning the plot

  # row_split <- factor(sort(avg_exp$ids))

  # Get dendrogram
  # Compute the distance based on correlation
  d <- as.dist(1-cor(assay(avg_exp, 'mean')))
  hc <- hclust(d, method = "complete")

  tree <- as.dendrogram(hc)

  # Pull out the order from the hc object
  # clone_order <- hc$labels[hc$order]

  # clone_order <- names(sort(table(sce[[clone_name]]),
  #                           decreasing = TRUE))

  # cell_order <- order(match(sce[[clone_name]], clone_order))
  #
  # ordered_cell_ids <- as.character(colData(sce)[cell_order, 'cell_id'])
  #
  # cnv_clusters <- colData(sce)[cell_order, clone_name]

  # row_split <- factor(cnv_clusters, levels = hc$labels[hc$order])

  # left_annot <- ComplexHeatmap::HeatmapAnnotation(
  #   Clone = hc$labels[hc$order],
  #   # Sample = sce[, ordered_cell_ids]$Sample,
  #   col = list(Clone = col_clones),
  #   # pch = hc$labels[hc$order],
  #   which = "row",
  #   show_legend = c(TRUE, FALSE)
  # )

  # ht_plot <- cnaHeatmap(sce = avg_exp, assay_name = new_assay, scale = "none", clone_name = "ids", border = TRUE, cluster_rows = tree, ...)
  ht_plot <- cnaHeatmap(sce = avg_exp, assay_name = new_assay, scale = "none", clone_name = "ids", border = TRUE, ...)

  # print(ht_plot)
  return(list(plot = ht_plot, sce = orig_sce))
}
