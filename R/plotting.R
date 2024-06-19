#' Plot Cell Copy Number
#'
#' @param sce sce object
#' @param cell_id Vector of cell ids
#' @param assay_name Name of assay to plot
#' @param col_fun Color mapping function from [circlize::colorRamp2()]
#' @param linewidth Width of the line segments in plot
#'
#' @return A ggplot object
#' @export
#'
plot_cell_cna <- function(sce,
                          cell_id = NULL,
                          assay_name = "counts",
                          col_fun = NULL,
                          linewidth = 1) {
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
  bindat$bin_id <- get_bin_ids(rowRanges(sce))

  # TODO: subset for main chromosomes and reorder

  plot_dat <- get_assay_dat(sce = sce, assay_names = assay_name, cell_id = cell_id)

  p <- plot_sc_track(plot_dat = plot_dat, assay_name = assay_name, col_fun = col_fun, linewidth = linewidth)

  return(p)
}


plot_sc_track <- function(plot_dat, assay_name, col_fun = NULL, linewidth = 1) {
  plot_dat$chr_no <- gsub("chr", "", plot_dat$seqnames)
  plot_dat$chr_no <- factor(plot_dat$chr_no, levels = unique(chr_reorder(plot_dat$chr_no)))

  # Grab y var
  y_var <- sym(assay_name)

  # Base plot
  base_p <- ggplot(plot_dat) +
    facet_grid(id ~ chr_no, scales = "free_x", space = "free_x") +
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
    plot_dat[plot_dat$assay_name > 11, assay_name] <- 11
    plot_dat$assay_name <- as.factor(plot_dat$assay_name)
    plot_dat$assay_name <- dplyr::recode_factor(plot_dat$assay_name, "11" = "11+")

    p <- base_p +
      geom_segment(aes(x = start, xend = end, y = !!y_var, yend = !!y_var, color = !!y_var), linewidth = linewidth) +
      scale_color_manual(values = cn_colors) +
      guides(colour = guide_legend(override.aes = list(size = 3)))
  } else if (!is.null(col_fun)) {
    cn_colors <- col_fun
    # This is hacky -- maybe there is a more direct way to pass color pallete from col_fun
    p <- base_p +
      geom_segment(aes(x = start, xend = end, y = !!y_var, yend = !!y_var, color = !!y_var), linewidth = linewidth) +
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

    p <- base_p + geom_segment(aes(x = start, xend = end, y = !!y_var, yend = !!y_var), linewidth = linewidth)
  }
  return(p)
}



#' Plot Psuedobulk cell CNA profiles
#'
#' @inheritParams plot_cell_cna
#' @param col_fun Color mapping function from [circlize::colorRamp2()]
#' @param aggr_fun Function to use to psuedobulk data
#' @param group_var Grouping variable. If `all` will pseudobulk all cells
#'
#' @return ggplot
#' @export
#'
plot_cell_psuedobulk_cna <- function(sce,
                                     assay_name,
                                     group_var = "all",
                                     aggr_fun = mean,
                                     col_fun = NULL) {
  if (group_var == "all") {
    sce[[group_var]] <- "all"
  }

  # avg_exp <- pseudobulk_sce(sce = sce, assay_name = assay_name, group_var = group_var)

  avg_exp <- pseudo_groups(sce, assay_name = assay_name, group_var = group_var, FUN = aggr_fun, na.rm = TRUE)

  plot_cell_cna(sce = avg_exp, assay_name = assay_name, col_fun = col_fun) + labs(title = group_var, y = assay_name)
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


#' Plot data with segments
#'
#' Will plot copy number tracks with raw input and overlaid segmented data
#'
#' @param sce `SingleCellExperiment` object with data
#' @param seg_assay segmented assay
#' @param input_assay input data assay (pre-segmented)
#' @param cell_id Cell to plot
#'
#' @return ggplot
#' @export
#'
plot_segs <- function(sce, seg_assay, input_assay, cell_id) {
  # Get consistent y axis scaling from the global assay object
  ylims <- range(assay(sce, input_assay))

  plot_dat <- get_assay_dat(sce, assay_names = c(seg_assay, input_assay), cell_id = cell_id)

  p <- plot_sc_track(plot_dat = plot_dat, assay_name = input_assay) +
    geom_segment(aes(x = start, xend = end, y = !!sym(seg_assay), yend = !!sym(seg_assay)), linewidth = 1, color = "red") +
    scale_y_continuous(limits = ylims)

  return(p)
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
#' @param top_annotation Top annotation as per [ComplexHeatmap::columnAnnotation()]
#' @param row_split Row split for [ComplexHeatmap::Heatmap()]
#' @param raster_quality Quality of raster (default: 10)
#' @param verbose Logical: Message verbosity
#' @param ... Additional parameters that can be passed to [ComplexHeatmap::Heatmap()]
#'
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
                       raster_quality = 10,
                       row_split = NULL,
                       ...) {
  # TODO: Enable multiple annotations
  # TODO: Enable multiple plots layered on top
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
      logger::log_error("{clone_name} not found in sce object")
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
      SummarizedExperiment::rowRanges(avg_exp) <- SummarizedExperiment::rowRanges(sce)

      # Compute the distance based on correlation
      d <- as.dist(1 - cor(assay(avg_exp, "mean")))
      hc <- hclust(d, method = "complete")

      # Pull out the order from the hc object
      clone_order <- hc$labels[hc$order]
    } else {
      # clone_order <- as.character(sort(unique(sce[[clone_name]])))
      clone_order <- as.character(unique(sce[[clone_name]]))
    }

    # Do we want to reorder the sce now?
    cell_order <- order(match(sce[[clone_name]], clone_order))

    sce <- sce[, cell_order]

    row_split <- factor(sce[[clone_name]], levels = clone_order)
    row_title <- clone_order
  } else {
    row_split <- row_split
    row_title <- NULL
  }

  if (is(clust_annot, "HeatmapAnnotation")) {
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
      if (!requireNamespace("EnsDb.Hsapiens.v86")) {
        logger::log_error("EnsDb.Hsapiens.v86 not installed. No gene overlaps detected in SCE input. Please run 'overlap_genes' prior to labelling genes.")
      } else {
        logger::log_warn("No gene overlaps detected in SCE input. Performing overlaps now.")
        sce <- overlap_genes(sce = sce, ensDb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, gene_biotype = "protein_coding")
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
    gene_bins_idx <- match(gene_bins, get_bin_ids(SummarizedExperiment::rowRanges(sce)))

    # Label genes
    bottom_ha_genes <- ComplexHeatmap::HeatmapAnnotation(genes = ComplexHeatmap::anno_mark(at = c(gene_bins_idx), labels = label_genes, which = "column", side = "bottom", link_width = unit(3, "mm")))
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
    if (!bulk_cn_col %in% colnames(mcols(SummarizedExperiment::rowRanges(sce)))) {
      logger::log_warn("{bulk_cn_col} not found in provided object. Not plotting...")
      top_annotation <- NULL
    } else {
      cn_dat <- data.frame(mcols(SummarizedExperiment::rowRanges(sce))[, bulk_cn_col])

      top_annotation <- ComplexHeatmap::HeatmapAnnotation(bulk_cn_col = ComplexHeatmap::anno_points(cn_dat, border = T, pch = 15, size = unit(1, "mm")))
      top_annotation@anno_list$bulk_cn_col@label <- bulk_cn_col
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
    raster_quality = raster_quality,
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


#' @rdname cnaHeatmap
#' @param aggr_fun Aggregation function
#' @param clust_lab Cluster labels
#' @param round Round pseudobulked values to integers
#' @export
#'
cloneCnaHeatmap <- function(sce, assay_name = "counts", clone_name = NULL, scale = c("none", "cells", "bins", "both"), log2 = FALSE, center = FALSE, clust_lab = TRUE, aggr_fun = mean, round = FALSE, ...) {
  if (is.null(rownames(sce))) {
    rownames(sce) <- 1:nrow(sce)
  }

  orig_sce <- sce

  # clust_mat <- scale_mat(assay(sce, assay_name), scale = scale, log2 = log2, center = center)

  # Subset the SCE to the same dimensions
  # sce <- sce[rownames(clust_mat), colnames(clust_mat)]

  # new_assay <- paste(assay_name, "avg", sep = "_")

  # assay(sce, new_assay) <- clust_mat

  # if (!is.null(clone_name)) {
  #
  # }

  # if (is.null(clustering_results)) {
  #   clustering_results <- perform_umap_clustering(assay(sce, new_assay))
  # }

  # sce$clone_id <- clustering_results$clustering[match(sce$Barcode, clustering_results$clustering$cell_id), "clone_id"]

  # avg_exp <- scuttle::summarizeAssayByGroup(sce, assay.type = new_assay, ids = sce[[clone_name]], statistics = "mean")
  # SummarizedExperiment::rowRanges(avg_exp) <- SummarizedExperiment::rowRanges(sce)

  avg_exp <- pseudo_groups(sce, assay_name = assay_name, group_var = clone_name, FUN = aggr_fun, na.rm = TRUE)

  # Order by clone size
  # avg_exp <- avg_exp[, order(avg_exp$ncells, decreasing = TRUE)]

  if (round) {
    # Round to integers
    assay(avg_exp, new_assay) <- round(assay(avg_exp, clone_name))
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


  # metadata(orig_sce)[[new_assay]] <- avg_exp

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
  ht_plot <- cnaHeatmap(sce = avg_exp, assay_name = assay_name, scale = scale, log2 = log2, center = center, clone_name = clone_name, border = TRUE, ...)

  # print(ht_plot)
  return(ht_plot)
}


#' Clone Copy Number Comparison Plot
#'
#' @param sce sce object
#' @param group_var column name with clone information
#' @param subset_clones optional vector of clones to subset for comparison
#' @param subset_chr optional vector of chromosomes to subset for comparison
#' @param assay_name name of assay to use
#' @param center_point Center point of the assay for drawing dashed lines
#' @param pseudobulk_fun mean or median
#' @param pt_size Size of the points in the plot
#' @param lg_pt_size Size of the points in the legend
#' @param xline Numeric for X intercept line
#' @param yline Numeric for Y intecept line
#'
#' @return ggapply plot
#' @export
#'
plot_clone_comp <- function(sce,
                            group_var = "clusters",
                            subset_clones = NULL,
                            subset_chr = NULL,
                            assay_name = "segment_merged_logratios",
                            center_point = NULL,
                            pt_size = 0.75,
                            xline = NULL,
                            yline = NULL,
                            lg_pt_size = 3,
                            pseudobulk_fun = mean) {
  # TODO: Allow for inversion of the plot to plot chromosomes on the facets by clones?

  if (!requireNamespace("GGally")) {
    logger::log_warn("The GGally package is required for the clone_comp_plot function. Please install")
    stop()
  }

  if (is.null(group_var)) {
    group_var <- "all"
    sce$all <- "all"
  }

  rownames(sce) <- rowRanges(sce)$bin_id <- get_bin_ids(rowRanges(sce))

  avg_exp <- pseudo_groups(sce,
    assay_name = assay_name,
    group_var = group_var,
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
    rownames_to_column(var = "bin_id") %>%
    left_join(as.data.frame(SummarizedExperiment::rowRanges(sce)))

  # Plot
  pcomb <- GGally::ggpairs(
    data = avg_exp,
    columns = columns, legend = length(columns) + 1,
    lower = list(
      continuous = GGally::wrap(my_cont, limits = xy_range, size = pt_size, lg_pt_size = lg_pt_size, xline = xline, yline = yline),
      mapping = aes(color = seqnames)
    ),
    diag = list(continuous = GGally::wrap(my_dens, limits = xy_range, center_point = center_point)), progress = FALSE
  ) +
    theme_bw() +
    theme(legend.position = "right") + labs(color = "Chromosome")
  return(pcomb)
}


my_cont <- function(data, mapping, size = 1, lg_pt_size = 3, xline = NULL, yline = NULL, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(size = size) +
    geom_abline(alpha = 0.75, linetype = "dashed") +
    scale_x_continuous(...) +
    scale_y_continuous(...) +
    scale_color_viridis_d() +
    guides(color = guide_legend(override.aes = list(size = lg_pt_size))) +
    geom_hline(yintercept = yline, linetype = "dashed", alpha = 0.75) +
    geom_vline(xintercept = xline, linetype = "dashed", alpha = 0.75)
}

my_dens <- function(data, mapping, center_point = 0, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_density() +
    scale_x_continuous(...) +
    geom_vline(xintercept = center_point, linetype = "dashed", alpha = 0.75)
}

get_label_centers <- function(obj, group_var = "clusters", reduced_dim = "UMAP") {
  # To any scatter of a umap can add + geom_label_repel(data = get_label_centers(sce), aes(x = x, y = y, label = clusters))

  if ("SingleCellExperiment" %in% class(obj)) {
    x_means <- lapply(split(reducedDim(obj, reduced_dim)[, 1], obj[[group_var]]), mean) %>% unlist()
    y_means <- lapply(split(reducedDim(obj, reduced_dim)[, 2], obj[[group_var]]), mean) %>% unlist()
    centers <- data.frame(x = x_means, y = y_means)
    centers[[group_var]] <- rownames(centers)
  }

  if ("Seurat" %in% class(obj)) {
    x_means <- lapply(split(Seurat::Embeddings(obj, reduced_dim)[, 1], obj[[group_var]]), mean) %>% unlist()
    y_means <- lapply(split(Seurat::Embeddings(obj, reduced_dim)[, 2], obj[[group_var]]), mean) %>% unlist()
    centers <- data.frame(x = x_means, y = y_means)
    centers[[group_var]] <- rownames(centers)
  }

  return(centers)
}


#' Plot Gene CNA
#'
#' @param sce sce
#' @param gene gene
#' @param assay_type assay
#' @param group_by grouping variable
#' @param color_by coloring variable
#' @param return_data If `TRUE`, will return plot data without plotting
#'
#' @return plot
#' @export
#'
plot_gene_cna <- function(sce,
                          gene,
                          assay_type = "counts",
                          group_by = NULL,
                          color_by = group_by,
                          return_data = FALSE) {
  # Check for the overlaps and create if possible
  if (is.null(sce@metadata$gene_overlap)) {
    # Attempt to perform the overlaps on the fly
    if (requireNamespace("EnsDb.Hsapiens.v86")) {
      logger::log_warn("No gene overlaps detected in SCE input. Performing overlaps now.")
      sce <- overlap_genes(sce = sce, ensDb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, gene_biotype = "protein_coding")
    } else {
      logger::log_error("No gene overlaps detected in SCE input. Please run 'overlap_genes' prior to labelling genes.")
    }
  }

  if (is.null(group_by)) {
    group_by <- "all"
    sce[[group_by]] <- group_by
  }

  sce[["group"]] <- sce[[group_by]]


  bin_id <- as.data.frame(sce@metadata$gene_overlap[which(sce@metadata$gene_overlap$gene_name == gene), ])[, "bin_id"]

  # Create plots
  plot_df <- scater::makePerCellDF(sce, assay.type = assay_type, use.coldata = TRUE, features = bin_id) %>%
    add_count(group) %>%
    mutate(sample_lab = glue("{group} n=({n})"))
  plot_df[[gene]] <- plot_df[[bin_id]]

  if (return_data) {
    return(plot_df)
  } else {
    plot_df %>%
      ggplot(aes(x = fct_reorder(sample_lab, .data[[bin_id]], median), y = .data[[bin_id]], color = .data[[color_by]])) +
      geom_pointrange(
        stat = "summary",
        fun.min = function(z) {
          quantile(z, 0.25)
        },
        fun.max = function(z) {
          quantile(z, 0.75)
        },
        fun = median
      ) +
      scale_y_continuous(trans = "log2", breaks = seq(0.5, 1.5, by = c(0.05))) +
      # facet_grid(. ~ Time.Point, scales = "free_x", space = "free_x") +
      guides(x = guide_axis(angle = 90)) +
      labs(x = group_by, y = glue("{gene} Bin Copy Number\n({bin_id})")) +
      theme(strip.background = element_blank(), axis.line = element_blank(), panel.border = element_rect(fill = NA, linewidth = 0.5))
  }
}
