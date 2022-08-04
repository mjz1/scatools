#' @export
perform_umap_clustering <- function(cn_matrix,
                                    n_neighbors = 10,
                                    min_dist = 0.1,
                                    minPts = 30,
                                    scale = c("none", "cells", "bins", "both"),
                                    log2 = FALSE,
                                    seed = 3,
                                    umapmetric = "correlation") {

  # TODO: Integrate with the SCE object better -- in the reduced dim slot

  cn_matrix <- as.matrix(cn_matrix)

  # Scale and clean matrix
  # This function expects bins in cols and cells in rows
  cn_matrix <- t(scale_mat(cn_matrix, log2 = log2, scale = scale))


  if (ncol(cn_matrix) < n_neighbors) {
    n_neighbors <- ncol(cn_matrix) - 1
  }

  minPts <- max(minPts, 2)

  set.seed(seed)

  logger::log_info("Calculating UMAP dimensionality reduction in {nrow(cn_matrix)} cells")
  if (nrow(cn_matrix) > 500 & is.null(seed)) {
    pca <- min(50, ncol(cn_matrix))
    fast_sgd <- TRUE
  } else {
    pca <- NULL
    fast_sgd <- FALSE
  }
  umapresults <- uwot::umap(cn_matrix,
    metric = umapmetric,
    n_neighbors = n_neighbors,
    n_components = 2,
    min_dist = min_dist,
    ret_model = TRUE,
    ret_nn = TRUE,
    pca = pca,
    fast_sgd = fast_sgd
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
  clusterids <- hdbscanresults$cluster
  clusterids[clusterids == 0] <- 702
  LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
  dfumap$clone_id <- LETTERS702[clusterids]
  dfumap <- dfumap %>%
    dplyr::mutate(clone_id = ifelse(clone_id == "ZZ", "0", clone_id))

  logger::log_info(paste0("Identified ", length(unique(dfumap$clone_id)), " clusters"))
  logger::log_info("Distribution of clusters:")
  f <- table(dfumap$clone_id)
  for (cl in sort(unique(dfumap$clone_id))) {
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

