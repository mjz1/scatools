#' Smooth outlier copy number counts
#'
#' @param sce SingleCellExperiment object
#' @param assay_name Name of assay to smooth
#' @param ncores Number of cores to use
#' @param smooth_name Name of returned assay with smoothed counts
#' @param ...
#'
#' @return A SingleCellExperiument Object
#' @export
#'
smooth_counts <- function(sce, assay_name, ncores = 1, smooth_name = paste(assay_name, "smoothed", sep = "_"), ...) {
  chrs <- as.vector(seqnames(rowRanges(sce)))
  starts <- start(rowRanges(sce))
  sample_ids <- colnames(sce)
  logger::log_info("Smoothing {assay_name}")
  smoothed_counts <- pbmcapply::pbmclapply(1:ncol(sce), mc.cores = ncores, FUN = function(i) {
    x <- as.vector(assay(x = sce, assay_name)[,i])
    obj <- DNAcopy::CNA(genomdat = x, chrom = chrs, maploc = starts, data.type = "logratio", sampleid = sample_ids[i], presorted = T)
    res <- round(withr::with_seed(3, smoothed_CNA_counts <- DNAcopy::smooth.CNA(obj,
                                                                                smooth.region = 4,
                                                                                outlier.SD.scale = 4,
                                                                                smooth.SD.scale = 2,
                                                                                trim = 0.025))[,3], 2)
  })
  names(smoothed_counts) <- sample_ids
  smoothed_counts <- as.matrix(dplyr::bind_rows(smoothed_counts))
  rownames(smoothed_counts) <- rownames(sce)
  smoothed_counts[smoothed_counts <= 0 ] <- 1e-4 # negative or zero values post smoothed.

  assay(sce, smooth_name) <- smoothed_counts
  logger::log_success("Smoothing complete! Smoothed counts in assay '{smooth_name}'")
  return(sce)
}

#' @export
segment_cnv <- function(sce, assay_name, new_assay = paste(assay_name, "segment", sep = "_"), alpha = 0.2, nperm = 10, min.width = 2, undo.splits = "none", verbose = 0, ncores = 1, ...) {
  chrs <- as.vector(seqnames(rowRanges(sce)))
  starts <- start(rowRanges(sce))
  sample_ids <- colnames(sce)
  # Perform segmentation
  # TODO: Make this function chromosome arm aware (ie segment within arms rather than chrs)
  segmented_counts <- pbmcapply::pbmclapply(1:ncol(sce), mc.cores = ncores, FUN = function(i) {
    x <- as.vector(assay(sce, assay_name)[,i])
    obj <- DNAcopy::CNA(genomdat = x, chrom = chrs, maploc = starts, data.type = "logratio", sampleid = sample_ids[i], presorted = T)
    res <- withr::with_seed(3, smoothed_CNA_counts <- DNAcopy::segment(obj,
                                                                       alpha = alpha,
                                                                       nperm = nperm,
                                                                       min.width = min.width,
                                                                       undo.splits = undo.splits,
                                                                       verbose = verbose
    ))

    # test0 <- res$segRows[[2]] + 1 - res$segRows[[1]]

    df <- data.frame(idx = 1:length(x), seg.mean = NA)

    for (j in 1:nrow(res$segRows)) {
      # Fails if no counts on final segments so we put try
      try(df[res$segRows[j,1]:res$segRows[j,2], 'seg.mean'] <- res$output[j, 'seg.mean'])
    }
    return(df$seg.mean)
  })

  names(segmented_counts) <- sample_ids
  segmented_counts <- as.matrix(dplyr::bind_rows(segmented_counts))
  rownames(segmented_counts) <- rownames(sce)

  assay(sce, new_assay) <- segmented_counts
  return(sce)
}

#' @export
merge_segments <- function(sce, smooth_assay, segment_assay, new_assay = "segment_merged", ncores = 1) {
  smooth_counts <- log2(assay(sce, smooth_assay))

  segment_df <- as.data.frame(assay(sce, segment_assay))
  segment_df[segment_df == 0] <- 1e-4
  segment_df <- log2(segment_df)

  logger::log_info("Merging segments using {ncores} cores for {ncol(segment_df)} cells")

  seg_ml_list <- pbmcapply::pbmclapply(seq_along(segment_df), mc.cores = ncores, function(i) {
    cell_name <- names(segment_df)[i]
    seg_means_ml <- copykit::mergeLevels(vecObs = smooth_counts[, i],
                                         vecPred = segment_df[, i],
                                         verbose = 0,
                                         pv.thres = 1e-4
    )$vecMerged
  })

  names(seg_ml_list) <- names(segment_df)
  seg_ml_df <- dplyr::bind_cols(seg_ml_list)
  seg_ml_df <- round(2^seg_ml_df, 2)
  rownames(seg_ml_df) <- rownames(sce)
  assay(sce, new_assay) <- seg_ml_df

  # saving as segment ratios
  seg_ratios <- sweep(seg_ml_df, 2, apply(seg_ml_df, 2, mean, na.rm = T), "/")
  rownames(seg_ratios) <- rownames(sce)
  assay(sce, paste(new_assay, "ratios", sep = "_")) <- as.matrix(round(seg_ratios, 2))

  sce <- copykit::logNorm(sce, assay = paste(new_assay, "ratios", sep = "_"), name = paste(new_assay, "logratios", sep = "_"))

  logger::log_info("Merged segments in: {new_assay}")
  logger::log_info("Merged segments ratios in: {paste(new_assay, 'ratios', sep = '_')}")
  logger::log_info("Merged segments logratios in: {paste(new_assay, 'logratios', sep = '_')}")

  return(sce)
}

#' @export
identify_normal <- function(sce, assay_name, group_by = "clusters", n_normal_clusts = 1) {

  if (length(unique(sce[[group_by]])) == 1) {
    logger::log_error("Only one group detected")
    sce$tumor_cell <- NA
    return(sce)
  }

  sce$seg_sd <- colSdDiffs(assay(sce, assay_name), na.rm = TRUE)

  # Use pairwise t tests between groups?
  # comps <- colData(sce) %>%
  #   as.data.frame() %>%
  #   rstatix::pairwise_t_test(as.formula(paste0("seg_sd ~ ", group_by)), pool.sd = TRUE, p.adjust.method = "bonferroni")

  # Simply identify grou with lowest median sequential segmental difference
  normal_clust <- colData(sce) %>%
    as.data.frame() %>%
    dplyr::group_by(!!as.symbol(group_by)) %>%
    dplyr::summarise(med = median(seg_sd)) %>%
    dplyr::slice_min(med, n = n_normal_clusts) %>%
    pull(var = 1)

  sce$tumor_cell <- FALSE

  sce$tumor_cell[which(!sce[[group_by]] %in% normal_clust)] <- TRUE

  logger::log_info("{table(sce$tumor_cell)[[1]]} normal cells identified in clusters: {paste(normal_clust, collapse =', ')}")

  return(sce)
}

#' @export
calc_ratios <- function (sce, assay_name, fun = c("mean", "median"), new_assay = paste(assay_name, "ratios", sep = "_")) {
  fun <- match.arg(fun)

  ratios <- sweep(assay(sce, assay_name), 2, apply(assay(sce, assay_name), 2, fun, na.rm = T), "/")
  assay(sce, new_assay) <- round(ratios, 2)
  return(sce)
}
