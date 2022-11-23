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

  assay(sce, "smoothed") <- smoothed_counts
  logger::log_success("Smoothing complete! Smoothed counts in assay '{smooth_name}'")
  return(sce)
}
