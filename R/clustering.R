#' Cluster using Seurat
#'
#' @param sce SingleCellExperiment object
#' @param assay_name Assay name. Can provide two assay names to perform joint clustering across both
#' @param do.scale scale
#' @param do.center center
#' @param algorithm clustering algorithm
#' @param resolution clustering resolution
#' @param n.neighbors neighbors for umap
#' @param npcs.pca Total Number of PCs to compute and store (50 by default)
#' @param features.pca One of 'all', 'variable', or a vector of features to include in dimensionality reduction. Defaults to 'all'.
#' @param nvar.features Number of variable features if `features.pca='variable'`
#' @param dims Number of reduced dimensions to use for FindNeighbors and UMAP
#' @param k.param Defines k for the k-nearest neighbor algorithm
#' @param annoy.metric Metric for [Seurat::FindNeighbors]
#' @param umap.metric Metric for [Seurat::RunUMAP]
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
                           do.scale = FALSE,
                           do.center = FALSE,
                           algorithm = 1,
                           resolution = 0.8,
                           n.neighbors = 10,
                           npcs.pca = 50,
                           features.pca = "all",
                           nvar.features = NULL,
                           dims = 1:npcs.pca,
                           k.param = 20,
                           suffix = "",
                           PCA_name = paste0("PCA", suffix),
                           UMAP_name = paste0("UMAP", suffix),
                           cluster_name = paste0("clusters", suffix),
                           umap.metric = "correlation",
                           annoy.metric = "cosine",
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

  # For joint clustering
  if (length(assay_name) == 2) {
    logger::log_info("Jointly clustering assays '{assay_name[1]}' and '{assay_name[2]}'")
    a1 <- scale(assay(sce, assay_name[1]))
    a2 <- scale(assay(sce, assay_name[2]))

    joint_data <- rbind(a1, a2)

    idx_1 <- 1:nrow(assay(sce, assay_name[1]))
    idx_2 <- (nrow(assay(sce, assay_name[1])) + 1):nrow(joint_data)

    rownames(joint_data)[idx_1] <- paste0(assay_name[1], "_", rownames(joint_data)[idx_1])
    rownames(joint_data)[idx_2] <- paste0(assay_name[2], "_", rownames(joint_data)[idx_2])

    srt <- Seurat::CreateSeuratObject(counts = joint_data)
    srt[["RNA"]]$data <- srt[["RNA"]]$counts
  } else {
    srt <- Seurat::CreateSeuratObject(counts = assay(sce, assay_name))
    srt[["RNA"]]$data <- srt[["RNA"]]$counts
  }

  if (features.pca == "all") {
    features.pca <- rownames(srt)
  } else if (features.pca == "variable") {
    if (is.null(nvar.features)) {
      logger::log_error("Variable features must provide 'nvar.features")
      stop()
    }
    srt <- Seurat::FindVariableFeatures(srt)
    # Will use Seurats find variable features
    features.pca <- NULL
  }

  srt <- Seurat::ScaleData(srt,
    do.scale = do.scale,
    do.center = do.center,
    verbose = verbose
  )

  if (length(features.pca) < npcs.pca) {
    logger::log_error("{length(features.pca)} features provided for PCA but requesting {npcs.pca} PCA dimensions. Please adjust.")
  }

  if (ncol(srt) < npcs.pca) {
    logger::log_warn("Not enough cells: {ncol(srt)} for requesting pcs: {npcs.pca}. Reducing to {ncol(srt)-1}")
    npcs.pca <- ncol(srt) - 1
  }

  srt <- Seurat::RunPCA(srt,
    features = features.pca,
    npcs = npcs.pca,
    verbose = FALSE
  )

  srt <- Seurat::FindNeighbors(srt,
    dims = dims,
    verbose = verbose,
    annoy.metric = annoy.metric,
    k.param = k.param
  )

  if (algorithm %in% c(4, "leiden")) {
    logger::log_info("Finding clusters using leiden algorithm")
    srt$seurat_clusters <- factor(leiden_wrapper(adj_mat = srt@graphs$RNA_snn, resolution = resolution))
    logger::log_success("Found ", length(unique(srt$seurat_clusters)), " communities")
  } else {
    srt <- Seurat::FindClusters(srt, resolution = resolution, algorithm = algorithm, verbose = verbose)
  }

  if (requireNamespace("HGC", quietly = TRUE)) {
    srt <- HGC::FindClusteringTree(srt, graph.type = "SNN")
  }

  srt <- Seurat::RunUMAP(srt, dims = dims, n.neighbors = n.neighbors, metric = umap.metric, verbose = verbose, seed.use = 3)

  # Put the PCA, UMAP, and clustering results into the original SCE
  reducedDim(sce_orig, PCA_name) <- srt@reductions$pca@cell.embeddings
  reducedDim(sce_orig, UMAP_name) <- srt@reductions$umap@cell.embeddings
  sce_orig[[cluster_name]] <- srt$seurat_clusters

  # Keep the graphs stored in the metadata
  sce_orig@metadata[[paste("graphs", suffix, sep = "_")]] <- srt@graphs

  return(sce_orig)
}


#' Wrapper for the Leiden Algorithm
#'
#' @param adj_mat Adjacency matrix
#' @param group_singletons Group singletons into nearest cluster. If FALSE, assign all singletons to a "singleton" group
#' @param resolution Resolution paramter
#'
#' @return cluster memberships
#' @export
#'
leiden_wrapper <- function(adj_mat, group_singletons = TRUE, resolution = 1) {
  if (!requireNamespace("igraph")) {
    logger::log_error("Package 'igraph' required for leiden clustering. Please install.")
    stop()
  }

  if (!requireNamespace("leidenbase")) {
    logger::log_error("Package 'leidenbase' required for leiden clustering. Please install.")
    stop()
  }

  # https://github.com/satijalab/seurat/discussions/6754?sort=top

  # Requires igraph, leidenbase
  graph_obj <- igraph::graph_from_adjacency_matrix(adj_mat, weighted = TRUE)

  res <- leidenbase::leiden_find_partition(
    graph_obj,
    partition_type = "RBConfigurationVertexPartition",
    resolution_parameter = resolution,
    num_iter = 10,
    seed = 3
  )

  ids <- res$membership
  names(ids) <- colnames(adj_mat)

  if (group_singletons) {
    ids <- Seurat:::GroupSingletons(ids, SNN = adj_mat, group.singletons = group_singletons, verbose = TRUE)
  }

  return(ids)
}
