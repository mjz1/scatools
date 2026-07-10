#' Cluster cells in a SingleCellExperiment
#'
#' SCE-native clustering workflow: PCA ([irlba::prcomp_irlba()]) ->
#' shared-nearest-neighbour graph ([BiocNeighbors::findKNN()] + Jaccard
#' weighting) -> Leiden community detection ([igraph::cluster_leiden()]) ->
#' UMAP embedding ([uwot::umap()]). This replaces the former Seurat-based
#' implementation and carries no Seurat dependency.
#'
#' @param sce SingleCellExperiment object
#' @param assay_name Assay name. Provide two names to jointly cluster across both.
#' @param do.scale,do.center Scale/center features before PCA
#' @param resolution Leiden resolution parameter
#' @param n.neighbors Neighbours for UMAP
#' @param npcs.pca Number of principal components to compute
#' @param features.pca One of `'all'`, `'variable'`, or a vector of features
#' @param nvar.features Number of variable features when `features.pca = 'variable'`
#' @param dims PCs to use for the graph and UMAP (default: all computed)
#' @param k.param k for the nearest-neighbour graph
#' @param suffix Suffix appended to the output PCA/UMAP/cluster names
#' @param PCA_name,UMAP_name,cluster_name Output names
#' @param umap.metric UMAP distance metric
#' @param run_umap Logical; compute a UMAP embedding (requires `uwot`)
#' @param verbose Message verbosity
#' @param ... (For the deprecated `cluster_seurat()` alias) further arguments
#'   passed to `cluster_sce()`.
#'
#' @return SingleCellExperiment with PCA/UMAP in `reducedDims` and cluster
#'   assignments in `colData`.
#' @export
cluster_sce <- function(sce,
                        assay_name,
                        do.scale = FALSE,
                        do.center = FALSE,
                        resolution = 0.8,
                        n.neighbors = 10,
                        npcs.pca = 50,
                        features.pca = "all",
                        nvar.features = NULL,
                        dims = NULL,
                        k.param = 20,
                        suffix = "",
                        PCA_name = paste0("PCA", suffix),
                        UMAP_name = paste0("UMAP", suffix),
                        cluster_name = paste0("clusters", suffix),
                        umap.metric = "cosine",
                        run_umap = TRUE,
                        verbose = TRUE) {
  # Clear any existing reduced dims
  for (dn in reducedDimNames(sce)) reducedDim(sce, dn) <- NULL

  # Build a feature x cell data matrix (optionally joint across two assays)
  if (length(assay_name) == 2) {
    if (verbose) cli::cli_alert_info("Jointly clustering assays '{assay_name[1]}' and '{assay_name[2]}'")
    a1 <- scale(as.matrix(assay(sce, assay_name[1])))
    a2 <- scale(as.matrix(assay(sce, assay_name[2])))
    mat <- rbind(a1, a2)
    rownames(mat) <- c(
      paste0(assay_name[1], "_", rownames(a1)),
      paste0(assay_name[2], "_", rownames(a2))
    )
  } else {
    mat <- as.matrix(assay(sce, assay_name))
  }

  # Feature selection
  if (identical(features.pca, "all")) {
    feats <- rownames(mat)
  } else if (identical(features.pca, "variable")) {
    if (is.null(nvar.features)) {
      cli::cli_abort("Provide 'nvar.features' when features.pca = 'variable'")
    }
    v <- apply(mat, 1, stats::var, na.rm = TRUE)
    feats <- rownames(mat)[order(v, decreasing = TRUE)][seq_len(min(nvar.features, nrow(mat)))]
  } else {
    feats <- features.pca
  }
  mat <- mat[feats, , drop = FALSE]
  mat[is.na(mat)] <- 0

  if (do.scale || do.center) {
    mat <- t(scale(t(mat), center = do.center, scale = do.scale))
    mat[is.na(mat)] <- 0
  }

  # PCA (cells as observations)
  npcs <- min(npcs.pca, nrow(mat) - 1, ncol(mat) - 1)
  if (npcs < 2) cli::cli_abort("Too few features/cells for PCA ({nrow(mat)} features, {ncol(mat)} cells)")
  emb <- withr::with_seed(3, irlba::prcomp_irlba(t(mat), n = npcs, center = TRUE, scale. = FALSE)$x)
  rownames(emb) <- colnames(sce)
  colnames(emb) <- paste0("PC_", seq_len(ncol(emb)))

  if (is.null(dims)) dims <- seq_len(ncol(emb))
  dims <- dims[dims <= ncol(emb)]
  emb_use <- emb[, dims, drop = FALSE]

  # SNN graph + Leiden clustering
  g <- build_snn_graph(emb_use, k = min(k.param, nrow(emb_use) - 1))
  memb <- cluster_leiden_graph(g, resolution = resolution)
  if (verbose) cli::cli_alert_success("Found {length(unique(memb))} clusters")

  # UMAP (optional)
  if (run_umap) {
    if (requireNamespace("uwot", quietly = TRUE)) {
      um <- withr::with_seed(3, uwot::umap(emb_use,
        n_neighbors = min(n.neighbors, nrow(emb_use) - 1),
        metric = umap.metric
      ))
      rownames(um) <- colnames(sce)
      colnames(um) <- c("UMAP_1", "UMAP_2")
      reducedDim(sce, UMAP_name) <- um
    } else {
      cli::cli_alert_warning("Package 'uwot' not installed; skipping UMAP")
    }
  }

  reducedDim(sce, PCA_name) <- emb
  sce[[cluster_name]] <- factor(memb)
  sce@metadata[[paste("snn_graph", suffix, sep = "_")]] <- g

  return(sce)
}

#' @rdname cluster_sce
#' @description `cluster_seurat()` is a deprecated alias for `cluster_sce()`.
#' @export
cluster_seurat <- function(sce, assay_name, ...) {
  .Deprecated("cluster_sce")
  cluster_sce(sce = sce, assay_name = assay_name, ...)
}

#' Build a shared-nearest-neighbour graph from a low-dimensional embedding
#'
#' @param emb Cell x dimension embedding matrix
#' @param k Number of nearest neighbours
#' @param prune Minimum Jaccard weight to keep an edge
#' @return An undirected weighted `igraph` graph
#' @noRd
build_snn_graph <- function(emb, k = 20, prune = 1 / 15) {
  knn <- BiocNeighbors::findKNN(emb, k = k, warn.ties = FALSE)$index
  n <- nrow(emb)
  nbr <- cbind(seq_len(n), knn) # include self
  kk <- ncol(nbr)
  A <- Matrix::sparseMatrix(
    i = rep(seq_len(n), times = kk),
    j = as.vector(nbr), x = 1, dims = c(n, n)
  )
  shared <- Matrix::tcrossprod(A) # shared-neighbour counts
  jac <- shared / (2 * kk - shared) # Jaccard weighting
  Matrix::diag(jac) <- 0
  jac <- methods::as(jac, "CsparseMatrix")
  jac@x[jac@x < prune] <- 0
  jac <- Matrix::drop0(jac)
  g <- igraph::graph_from_adjacency_matrix(jac, mode = "undirected", weighted = TRUE, diag = FALSE)
  igraph::V(g)$name <- rownames(emb)
  g
}

#' Leiden community detection on a weighted igraph
#' @noRd
cluster_leiden_graph <- function(g, resolution = 0.8, n_iterations = 10) {
  cl <- withr::with_seed(3, igraph::cluster_leiden(g,
    objective_function = "modularity",
    weights = igraph::E(g)$weight,
    resolution_parameter = resolution,
    n_iterations = n_iterations
  ))
  memb <- igraph::membership(cl)
  names(memb) <- igraph::V(g)$name
  memb
}

#' Wrapper for the Leiden Algorithm on an adjacency matrix
#'
#' @param adj_mat Adjacency (SNN) matrix
#' @param group_singletons Reassign singleton clusters to their most connected
#'   neighbouring cluster. If `FALSE`, singletons are left as-is.
#' @param resolution Resolution parameter
#'
#' @return Named vector of cluster memberships
#' @export
#'
leiden_wrapper <- function(adj_mat, group_singletons = TRUE, resolution = 1) {
  g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
  memb <- cluster_leiden_graph(g, resolution = resolution)
  if (group_singletons) {
    memb <- group_singleton_clusters(memb, adj_mat)
  }
  memb
}

#' Reassign singleton clusters to their most strongly connected neighbour cluster
#' @noRd
group_singleton_clusters <- function(ids, snn) {
  tab <- table(ids)
  singletons <- names(tab)[tab == 1]
  if (length(singletons) == 0) {
    return(ids)
  }
  snn <- as.matrix(snn)
  for (s in singletons) {
    cell <- which(ids == s)
    conn <- tapply(snn[cell, ], ids, sum)
    conn[as.character(s)] <- 0
    best <- names(which.max(conn))
    if (length(best) == 1 && conn[[best]] > 0) ids[cell] <- best
  }
  ids
}
