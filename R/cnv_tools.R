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
segment_cnv <- function(sce, assay_name, new_assay = paste(assay_name, "segment", sep = "_"), alpha = 0.2, nperm = 10, min.width = 2, undo.splits = "none", verbose = 0, ...) {
  chrs <- as.vector(seqnames(rowRanges(sce)))
  starts <- start(rowRanges(sce))
  sample_ids <- colnames(sce)
  # Perform segmentation
  segmented_counts <- pbmcapply::pbmclapply(1:ncol(sce), mc.cores = 8, FUN = function(i) {
    x <- as.vector(assay(sce, assay_name)[,i])
    obj <- DNAcopy::CNA(genomdat = x, chrom = chrs, maploc = starts, data.type = "logratio", sampleid = sample_ids[i], presorted = T)
    res <- withr::with_seed(3, smoothed_CNA_counts <- DNAcopy::segment(obj,
                                                                       alpha = alpha,
                                                                       nperm = nperm,
                                                                       min.width = min.width,
                                                                       undo.splits = undo.splits,
                                                                       verbose = verbose
    ))

    test0 <- res$segRows[[2]] + 1 - res$segRows[[1]]

    df <- data.frame(idx = 1:length(x), seg.mean = NA)

    for (j in 1:nrow(res$segRows)) {
      df[res$segRows[j,1]:res$segRows[j,2], 'seg.mean'] <- res$output[j, 'seg.mean']
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
calc_ratios <- function (sce, assay_name, fun = c("mean", "median"), new_assay = paste(assay_name, "ratios", sep = "_")) {
  fun <- match.arg(fun)

  ratios <- sweep(assay(sce, assay_name), 2, apply(assay(sce, assay_name), 2, fun, na.rm = T), "/")
  assay(sce, new_assay) <- round(ratios, 2)
  return(sce)
}
