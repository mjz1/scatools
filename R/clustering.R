
#' Cluster using Seurat
#'
#' @param sce SingleCellExperiment object
#' @param assay_name Assay name
#' @param raw_counts Assay containing raw counts
#' @param do.scale scale
#' @param do.center center
#' @param algorithm clustering algorithm
#' @param resolution clustering resolution
#' @param n.neighbors neighbors for umap
#' @param metric metric for umap
#' @param suffix Suffix name to add to the PCA, UMAP, and clusters
#' @param PCA_name Name to store PCA dimred
#' @param UMAP_name Name to store UMAP dimred
#' @param cluster_name Name to store seurat clusters
#' @param verbose Message verbosity
#'
#' @return SingleCellExperiment obj
#' @export
#'
cluster_seurat <- function(sce,
                           assay_name,
                           raw_counts = "counts",
                           do.scale = FALSE,
                           do.center = FALSE,
                           algorithm = 1,
                           resolution = 0.8,
                           n.neighbors = 10,
                           suffix = "",
                           PCA_name = paste0("PCA", suffix),
                           UMAP_name = paste0("UMAP", suffix),
                           cluster_name = paste0("clusters", suffix),
                           metric = "correlation",
                           verbose = TRUE) {

  # TODO: Improve documentation of this function.
  # TODO: Return neighbors object if possible to allow future umap projection
  # TODO:

  if (!requireNamespace("Seurat")) {
    logger::log_error("Seurat not installed. Please install Seurat to use this function.")
  }

  sce_orig <- sce
  # For safety, clear any reduced dims present before we process
  for (dim_name in reducedDimNames(sce)) {
    reducedDim(sce, dim_name) <- NULL
  }

  srt <- suppressWarnings(Seurat::as.Seurat(sce, counts = raw_counts, data = assay_name))

  srt <- Seurat::ScaleData(srt, do.scale = do.scale, do.center = do.center, verbose = verbose)

  srt <- Seurat::RunPCA(srt, features = rownames(srt), verbose = FALSE)

  srt <- Seurat::FindNeighbors(srt, dims = 1:50, verbose = verbose)
  srt <- Seurat::FindClusters(srt, resolution = resolution, algorithm = algorithm, verbose = verbose)

  srt <- HGC::FindClusteringTree(srt, graph.type = "SNN")

  srt <- Seurat::RunUMAP(srt, dims = 1:50, n.neighbors = n.neighbors, metric = metric, verbose = verbose)

  # Put the PCA, UMAP, and clustering results into the original SCE
  reducedDim(sce_orig, PCA_name) <- srt@reductions$pca@cell.embeddings
  reducedDim(sce_orig, UMAP_name) <- srt@reductions$umap@cell.embeddings
  sce_orig[[cluster_name]] <- srt[[]]$seurat_clusters

  # Keep the graphs stored in the metadata
  sce_orig@metadata[[paste("graphs", suffix, sep = "_")]] <- srt@graphs

  return(sce_orig)
}

#' Add UMAP clusters to SCE
#'
#' @param sce SCE object
#' @param assay_name Name of assay to use for clustering
#' @param name 	String specifying the name to be used to store the result in the [reducedDims] of the output.
#' @param clone_colname Name of clone column in resulting sce object.
#' @param force Logical. Overwrite identical column names in sce object
#' @inheritParams perform_umap_clustering
#'
#' @return An SCE object with the umap results and clones
#' @export
#'
add_umap_clusters <- function(sce,
                              assay_name,
                              n_neighbors = 10,
                              min_dist = 0.1,
                              minPts = 30,
                              scale = c("none", "cells", "bins", "both"),
                              name = "UMAP",
                              clone_colname = "clone_id",
                              log2 = FALSE,
                              seed = 3,
                              metric = "correlation",
                              verbose = TRUE,
                              force = TRUE) {
  clust_results <- perform_umap_clustering(
    cn_matrix = assay(sce, assay_name),
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    minPts = minPts,
    scale = scale,
    log2 = log2,
    seed = seed,
    metric = metric,
    verbose = verbose
  )

  # Add the UMAP results to the sce object
  reducedDim(sce, name, withDimnames = TRUE) <- clust_results$umapresults$embedding

  # Add the clone_ids
  # Non overwriting naming
  clone_colname_orig <- clone_colname
  n_idx <- 0
  while (clone_colname %in% colnames(colData(sce)) & !force) {
    if (verbose) {
      logger::log_warn("{clone_colname} already present as column in sce object. Picking alternate name")
    }
    clone_colname <- paste0(clone_colname_orig, "_", n_idx)
  }

  if (verbose) {
    logger::log_info("Saving clones in column: '{clone_colname}'")
  }
  sce[[clone_colname]] <- clust_results$clustering$clone_id[match(colnames(sce), clust_results$clustering$cell_id)]
  sce[[paste0(clone_colname, "_clonesize")]] <- clust_results$clustering$clone_size[match(colnames(sce), clust_results$clustering$cell_id)]

  return(sce)
}



#' Perform UMAP and clustering
#'
#' Clusters a copy number matrix using `hdbscan`.
#'
#' @param cn_matrix Copy number matrix with cells in columns and bins in rows
#'
#' @return A list with the clustering results
#' \describe{
#'   \item{clustering}{Data frame of cell and cluster identities}
#'   \item{hdbscanresults}{Results from [dbscan::hdbscan]}
#'   \item{umapresults}{Results from [uwot::umap]}
#'   \item{tree}{Results from [ape::as.phylo] on `hdbscanresults$hc`}
#' }
#' @export
#'
#' @inheritParams cnaHeatmap
#' @inheritParams uwot::umap
perform_umap_clustering <- function(cn_matrix,
                                    n_neighbors = 10,
                                    min_dist = 0.1,
                                    minPts = 30,
                                    scale = c("none", "cells", "bins", "both"),
                                    log2 = FALSE,
                                    seed = 3,
                                    metric = "correlation",
                                    verbose = TRUE) {

  # TODO: Integrate with the SCE object better -- in the reduced dim slot
  # TODO: Add clone id colors to the clustering object so that they are fixed across plots

  cn_matrix <- as.matrix(cn_matrix)

  # Scale and clean matrix
  # This function expects bins in cols and cells in rows
  cn_matrix <- t(scale_mat(cn_matrix, log2 = log2, scale = scale))


  if (ncol(cn_matrix) < n_neighbors) {
    n_neighbors <- ncol(cn_matrix) - 1
  }

  minPts <- max(minPts, 2)

  set.seed(seed)

  logger::log_info("Calculating UMAP in {nrow(cn_matrix)} cells and {ncol(cn_matrix)} bins")

  if (nrow(cn_matrix) > 500 & is.null(seed)) {
    pca <- min(50, ncol(cn_matrix))
    fast_sgd <- TRUE
  } else {
    pca <- NULL
    fast_sgd <- FALSE
  }
  umapresults <- uwot::umap(cn_matrix,
    metric = metric,
    n_neighbors = n_neighbors,
    n_components = 2,
    min_dist = min_dist,
    ret_model = TRUE,
    ret_nn = TRUE,
    pca = pca,
    fast_sgd = fast_sgd,
    verbose = verbose
  )

  dfumap <- data.frame(
    umap1 = umapresults$embedding[, 1],
    umap2 = umapresults$embedding[, 2],
    cell_id = row.names(cn_matrix)
  )

  dfumap$umap1 <- unlist(lapply(dfumap$umap1, function(x) ifelse(!is.finite(x), 0.0, x))) # remove non finite values
  dfumap$umap2 <- unlist(lapply(dfumap$umap2, function(x) ifelse(!is.finite(x), 0.0, x))) # remove non finite values
  rownames(dfumap) <- row.names(cn_matrix)

  logger::log_info("Clustering cells using hdbscan...")
  gentree <- FALSE
  while (gentree == FALSE) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    hdbscanresults <- try(dbscan::hdbscan(dfumap[, 1:2],
      minPts = minPts,
      gen_hdbscan_tree = FALSE,
      gen_simplified_tree = FALSE
    ))
    if (class(hdbscanresults) == "try-error") {
      logger::log_warn("Only 1 cluster found, reducing minPts size by 10...")
      minPts <- round(minPts - 10)
      if (minPts <= 0) {
        message("Only 1 cluster can be found")
        gentree <- TRUE
      }
      logger::log_info(paste0("Cluster size = ", minPts))
    } else {
      gentree <- TRUE
    }
  }

  # TODO: Clone names alphabetical in order of size
  clusterids <- hdbscanresults$cluster
  clusterids[clusterids == 0] <- 702
  LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
  dfumap$clone_id <- LETTERS702[clusterids]
  dfumap <- dfumap %>%
    dplyr::mutate(clone_id = ifelse(clone_id == "ZZ", "0", clone_id)) %>%
    dplyr::add_count(clone_id, name = "clone_size") %>%
    dplyr::arrange(dplyr::desc(clone_size))
  rownames(dfumap) <- dfumap$cell_id

  logger::log_info(paste0("Identified ", length(unique(dfumap$clone_id)), " clusters"))
  logger::log_info("Distribution of clusters:")
  f <- sort(table(dfumap$clone_id), decreasing = TRUE)
  for (cl in unique(dfumap$clone_id)) {
    message(paste0("  Cluster ", cl, ": ", f[[cl]]))
  }

  if (length(unique(dfumap$clone_id)) == 1) {
    dfumap$clone_id <- "A"
  }

  tree <- ape::as.phylo(hdbscanresults$hc, use.labels = TRUE)
  tree$tip.label <- row.names(cn_matrix)[as.numeric(tree$tip.label)]

  return(list(
    clustering = dfumap,
    hdbscanresults = hdbscanresults,
    umapresults = umapresults,
    tree = tree
  ))
}
