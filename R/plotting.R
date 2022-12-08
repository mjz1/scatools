# TODO: Separate out the single cell plot generation from data preprocessing to allow for more flexibility


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
    cell_id <- colnames(sce[, cell_id])
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
    dplyr::left_join(as.data.frame(bindat), by = "bin_id") #%>%
    # dplyr::filter(!is.na(counts))

  plot_dat$chr_no <- gsub("chr", "", plot_dat$chr)
  plot_dat$chr_no <- factor(plot_dat$chr_no, levels = unique(chr_reorder(plot_dat$chr_no)))

  # Keep cell ordering as provided
  plot_dat$barcode <- factor(plot_dat$barcode, levels = cell_id)

  # Base plot
  base_p <- ggplot(plot_dat) +
    facet_grid(barcode ~ chr_no, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = assay_name, color = assay_name) +
    theme_bw() +
    theme(
      panel.spacing.x = unit(0, "lines"),
      panel.border = element_rect(fill = NA),
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  # check for assay type and add colors
  if ((!is.null(col_fun) && grepl("state", assay_name))) {
    cn_colors <- state_cn_colors()
    # Clip counts to 11
    plot_dat[plot_dat$counts > 11, "counts"] <- 11
    plot_dat$counts <- as.factor(plot_dat$counts)
    plot_dat$counts <- dplyr::recode_factor(plot_dat$counts, "11" = "11+")

    p <- base_p +
      geom_segment(aes(x = start, xend = end, y = counts, yend = counts, color = counts), size = 1) +
      scale_color_manual(values = cn_colors) +
      guides(colour = guide_legend(override.aes = list(size = 3)))
  } else if (!is.null(col_fun)) {
    cn_colors <- col_fun
    # This is hacky -- maybe there is a more direct way to pass color pallete from col_fun
    p <- base_p +
      geom_segment(aes(x = start, xend = end, y = counts, yend = counts, color = counts), linewidth = 1) +
      scale_color_gradient2(
        low = attr(cn_colors, "colors")[1],
        mid = attr(cn_colors, "colors")[2],
        high = attr(cn_colors, "colors")[3],
        midpoint = attr(cn_colors, "breaks")[2],
        limits = c(attr(cn_colors, "breaks")[1], attr(cn_colors, "breaks")[3])
      )
  } else {
    # Attempt some intelligent color mapping? For now, no.
    cn_colors <- NULL

    # base_p$mapping$colour <- NULL # remove color mapping

    p <- base_p + geom_segment(aes(x = start, xend = end, y = counts, yend = counts), size = 1)
  }

  return(p)
}


#' Plot Psuedobulk cell CNA profiles
#'
#' @inheritParams plot_cell_cna
#' @param col_fun Color mapping function from [circlize::colorRamp2()]
#' @param aggr_fun Function to use to psuedobulk data.
#'
#' @return ggplot
#' @export
#'
plot_cell_psuedobulk_cna <- function(sce, assay_name, group_var = "all", aggr_fun = mean, col_fun = NULL) {
  if (group_var == "all") {
    sce[[group_var]] <- "all"
  }

  # avg_exp <- pseudobulk_sce(sce = sce, assay_name = assay_name, group_var = group_var)

  avg_exp <- pseudo_groups(sce, assay_name = assay_name, ids = sce[[group_var]], FUN = aggr_fun, na.rm = TRUE)

  plot_cell_cna(sce = avg_exp, assay_name = "pseudo", col_fun = col_fun) + labs(title = group_var, y = assay_name)
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
  # TODO: Enable passing of assay specific plot options
  plots <- vector(mode = "list")

  if (is.numeric(cell_id)) {
    cell_id <- colnames(sce[, cell_id])
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
#' @param cluster_cells If the value is a logical, it controls whether to make cluster on rows. The value can also be a `stats::hclust` or a `stats::dendrogram` which already contains clustering. Note this will override ordering of specified clones. Check https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#clustering.
#' @param log2 Logical: Log2 transform the matrix prior to plotting
#' @param center logical: center the matrix prior to plotting
#' @param scale One of `'cells', 'bins', 'both' or 'none'`. Determines what kind of scaling is done.
#' @param label_genes Optional: Vector of gene names to label in the heatmap. Note: [overlap_genes] must be run prior to labelling genes.
#' @param col_fun Color mapping function from [circlize::colorRamp2()]
#' @param col_clones Optional: A named vector (or unnamed) of clone colors.
#' @param cluster_clones Logical: Whether or not to order the clones by clustering
#' @param legend_name Name of the legend
#' @param clust_annot Annotate cluster and sample labels
#' @param bulk_cn_col Name of column in `rowRanges(sce)` that contains bulk copy number data to plot on top of heatmap
#' @param verbose Logical: Message verbosity
#' @param ... Additional parameters that can be passed to [ComplexHeatmap::Heatmap()]
#'
#' @return A heatmap
#' @export
#'
cnaHeatmap <- function(sce,
                       assay_name = "state",
                       clone_name = NULL,
                       cluster_cells = FALSE,
                       log2 = FALSE,
                       center = FALSE,
                       scale = c("none", "cells", "bins", "both"),
                       label_genes = NULL,
                       col_fun = NULL,
                       col_clones = NULL,
                       cluster_clones = FALSE,
                       legend_name = assay_name,
                       clust_annot = TRUE,
                       verbose = TRUE,
                       top_annotation = NULL,
                       bulk_cn_col = NULL,
                       ...) {
  # if (is.null(rownames(sce))) {
  #   rownames(sce) <- 1:nrow(sce)
  # }

  sce <- scale_sub(
    sce = sce, assay_name = assay_name, log2 = log2,
    scale = scale, new_assay = "heatmat", center = center
  )

  orig_order <- colnames(sce)

  if (!is.null(clone_name)) {
    if (!clone_name %in% colnames(colData(sce))) {
      logger::log_warn("{clone_name} not found in sce object")
    }
  }

  if (!is.null(clone_name)) {
    if (cluster_clones) {
      # First get the average per clone signal to order the clones properly
      avg_exp <- scuttle::summarizeAssayByGroup(sce,
        assay.type = assay_name,
        ids = sce[[clone_name]],
        statistics = "mean"
      )
      rowRanges(avg_exp) <- rowRanges(sce)

      # Compute the distance based on correlation
      d <- as.dist(1 - cor(assay(avg_exp, "mean")))
      hc <- hclust(d, method = "complete")

      # Pull out the order from the hc object
      clone_order <- hc$labels[hc$order]
    } else {
      clone_order <- as.character(sort(unique(sce[[clone_name]])))
    }

    # Do we want to reorder the sce now?
    cell_order <- order(match(sce[[clone_name]], clone_order))

    sce <- sce[, cell_order]

    row_split <- factor(sce[[clone_name]], levels = clone_order)
    row_title <- clone_order
  } else {
    row_split <- NULL
    row_title <- NULL
  }

  if (class(clust_annot) == "HeatmapAnnotation") {
    left_annot <- clust_annot
  } else if (clust_annot & !is.null(clone_name)) {
    if (is.null(col_clones)) {
      col_clones <- scales::hue_pal()(length(unique(sce[[clone_name]])))
      names(col_clones) <- levels(factor(sce[[clone_name]]))
    }

    # Code to allow passing of any length vector of unnamed colors to make it
    # easier across plots/cluster lengths
    if (is.null(names(col_clones))) {
      col_clones <- col_clones[1:length(unique(sce[[clone_name]]))]
      names(col_clones) <- levels(factor(sce[[clone_name]]))
    }

    left_annot <- ComplexHeatmap::HeatmapAnnotation(
      Clone = sce[[clone_name]],
      col = list(Clone = col_clones),
      which = "row",
      show_legend = c(TRUE, FALSE)
    )
  } else {
    left_annot <- NULL
  }

  # Label genes
  if (!is.null(label_genes)) {
    if (is.null(sce@metadata$gene_overlap)) {
      # Attempt to perform the overlaps on the fly
      if (requireNamespace("EnsDb.Hsapiens.v86")) {
        logger::log_warn("No gene overlaps detected in SCE input. Performing overlaps now.")
        sce <- overlap_genes(sce = sce, ensDb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, gene_biotype = "protein_coding")
      } else {
        logger::log_error("No gene overlaps detected in SCE input. Please run 'overlap_genes' prior to labelling genes.")
      }
    }
    # Now create a heatmap where genes are highlighted in the bottom panel
    match_idx <- match(label_genes, sce@metadata$gene_overlap$symbol)

    # Remove and warn about missing genes
    na_idx <- which(is.na(match_idx))
    missing_g <- label_genes[na_idx]
    if (length(na_idx) > 0) {
      logger::log_warn("Genes not found: {paste(missing_g, collapse = '; ')}")
    }

    label_genes <- label_genes[!is.na(match_idx)]
    match_idx <- match_idx[!is.na(match_idx)]

    gene_bins <- mcols(sce@metadata$gene_overlap[match_idx])[["bin_id"]]
    gene_bins_idx <- match(gene_bins, get_bin_ids(rowRanges(sce)))

    # Label genes
    bottom_ha_genes <- HeatmapAnnotation(genes = anno_mark(at = c(gene_bins_idx), labels = label_genes, which = "column", side = "bottom", link_width = unit(3, "mm")))
  } else {
    # Make null so we can pass to heatmap function regardless
    bottom_ha_genes <- NULL
  }

  if (grepl("state", assay_name) || length(col_fun) == length(state_cn_colors())) {
    assay(sce, "heatmat")[assay(sce, "heatmat") >= 11] <- "11+"
    cn_colors <- state_cn_colors()
  } else if (!is.null(col_fun)) {
    cn_colors <- col_fun
  } else {
    # Attempt some intelligent color mapping?
    cn_colors <- NULL
  }

  # Split columns by chromosome
  chrs <- as.vector(gsub("chr", "", GenomeInfoDb::seqnames(SummarizedExperiment::rowRanges(sce))))
  col_split <- factor(chrs, levels = unique(gtools::mixedsort(chrs)))

  if (class(cluster_cells) %in% c("dendrogram", "hclust") || cluster_cells == TRUE) {
    # Revert to original ordering
    # Maybe there is a better way to handle this
    sce <- sce[, orig_order]

    if (!is.null(clone_name)) {
      row_split <- length(unique(sce[[clone_name]]))
      row_title <- NULL

      # Catch single row split
      if (row_split == 1) {
        row_split <- NULL
      }

      left_annot <- ComplexHeatmap::HeatmapAnnotation(
        Clone = sce[[clone_name]],
        col = list(Clone = col_clones),
        which = "row",
        show_legend = c(TRUE, FALSE)
      )
    }
  }

  # Add bulk CN annotation
  if (!is.null(bulk_cn_col)) {
    if (!bulk_cn_col %in% colnames(mcols(rowRanges(sce)))) {
      logger::log_warn("{bulk_cn_col} not found in provided object. Not plotting...")
      top_annotation <- NULL
    } else {
      cn_dat <- mcols(rowRanges(sce))[, bulk_cn_col]
      top_annotation <- ComplexHeatmap::HeatmapAnnotation(bulk_cn_col = ComplexHeatmap::anno_points(cn_dat, border = T, pch = 15, size = unit(1, "mm")))
      # Need to figure out how to rename the anno.
      # top_annotation@anno_list$bulk_cn_col@name <- bulk_cn_col
    }
  }

  suppressMessages(ht_plot <- ComplexHeatmap::Heatmap(
    name = legend_name,
    matrix = t(assay(sce, "heatmat")),
    show_row_names = FALSE,
    col = cn_colors,
    cluster_columns = FALSE,
    cluster_rows = cluster_cells,
    show_column_names = FALSE,
    na_col = "white",
    use_raster = TRUE,
    raster_quality = 10,
    column_split = col_split,
    left_annotation = left_annot,
    row_split = row_split,
    row_title = row_title,
    bottom_annotation = bottom_ha_genes,
    top_annotation = top_annotation,
    ...
  ))

  return(ht_plot)
}

#' @export
cloneCnaHeatmap <- function(sce, assay_name = "counts", clone_name = NULL, scale = c("none", "cells", "bins", "both"), log2 = FALSE, center = FALSE, clustering_results = NULL, clust_lab = TRUE, aggr_fun = mean, round = FALSE, ...) {
  if (is.null(rownames(sce))) {
    rownames(sce) <- 1:nrow(sce)
  }

  orig_sce <- sce

  clust_mat <- scale_mat(assay(sce, assay_name), scale = scale, log2 = log2, center = center)

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

  # avg_exp <- scuttle::summarizeAssayByGroup(sce, assay.type = new_assay, ids = sce[[clone_name]], statistics = "mean")
  # rowRanges(avg_exp) <- rowRanges(sce)

  avg_exp <- pseudo_groups(sce, assay_name = new_assay, ids = sce[[clone_name]], FUN = aggr_fun, na.rm = TRUE)

  # Order by clone size
  avg_exp <- avg_exp[, order(avg_exp$ncells, decreasing = TRUE)]

  if (round) {
    # Round to integers
    assay(avg_exp, new_assay) <- round(assay(avg_exp, "pseudo"))
  } else {
    assay(avg_exp, new_assay) <- assay(avg_exp, "pseudo")
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

  # If gene overlaps are present propagate to the pseudobulked
  if (!is.null(sce@metadata$gene_overlap)) {
    avg_exp@metadata$gene_overlap <- sce@metadata$gene_overlap
  }

  # TODO: Figure out how to handle returning the avg exp data and plot
  # Ideally want to replace the SCE in the parent environment while also returning the plot

  # row_split <- factor(sort(avg_exp$ids))

  # Get dendrogram
  # Compute the distance based on correlation
  # d <- as.dist(1 - cor(assay(avg_exp, "mean")))
  # hc <- hclust(d, method = "complete")
  #
  # tree <- as.dendrogram(hc)

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
  return(ht_plot)
}


#' Clone Copy Number Comparison Plot
#'
#' @param sce sce object
#' @param clone_column column name with clone information
#' @param subset_clones optional vector of clones to subset for comparison
#' @param subset_chr optional vector of chromosomes to subset for comparison
#' @param assay_name name of assay to use
#' @param center_point Center point of the assay for drawing dashed lines
#' @param pseudobulk_fun mean or median
#'
#' @return ggapply plot
#' @export
#'
plot_clone_comp <- function(sce,
                            clone_column = "clusters",
                            subset_clones = NULL,
                            subset_chr = NULL,
                            assay_name = "segment_merged_logratios",
                            center_point = 0,
                            pseudobulk_fun = mean) {
  # TODO: Allow for inversion of the plot to plot chromosomes on the facets by clones?


  avg_exp <- pseudo_groups(sce,
    assay_name = assay_name,
    ids = sce[[clone_column]],
    FUN = pseudobulk_fun,
    na.rm = TRUE
  ) %>%
    assay()

  # Get the range prior to subsetting
  xy_range <- range(avg_exp, na.rm = TRUE)

  # Subset if specified
  if (!is.null(subset_clones)) {
    avg_exp <- avg_exp[, as.character(subset_clones)]
  }

  if (!is.null(subset_chr)) {
    avg_exp <- avg_exp[grepl(paste(paste0(subset_chr, "_"), collapse = "|"), rownames(avg_exp)), ]
  }

  columns <- colnames(avg_exp)

  # Convert to data frame and join with range information
  avg_exp <- avg_exp %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID") %>%
    left_join(as.data.frame(rowRanges(sce)))

  # Plot
  pcomb <- GGally::ggpairs(
    data = avg_exp,
    columns = columns, legend = length(columns) + 1,
    lower = list(
      continuous = GGally::wrap(my_cont, limits = xy_range, size = 0.75),
      mapping = aes(color = seqnames)
    ),
    diag = list(continuous = GGally::wrap(my_dens, limits = xy_range, center_point = center_point))
  ) +
    theme_bw() +
    theme(legend.position = "right") + labs(color = "Chromosome")
  return(pcomb)
}


my_cont <- function(data, mapping, size = 1, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(size = size) +
    geom_abline(alpha = 0.75, linetype = "dashed") +
    scale_x_continuous(...) +
    scale_y_continuous(...)
}

my_dens <- function(data, mapping, center_point = 0, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_density() +
    scale_x_continuous(...) +
    geom_vline(xintercept = center_point, linetype = "dashed", alpha = 0.75)
}
