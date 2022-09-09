#' Summarise copu number signal over chromosome arms
#'
#' @param sce SingleCellExperiment objet
#' @param assay_name Name of assay
#' @param genome_ver Genome version ("hg38" or "hg19")
#'
#' @return A chromosome arm level binned SingleCellExperiment
#' @export
#'
# TODO: Make this more general to aggregate signal over defined bins
summarise_chr_arm <- function(sce, assay_name, genome_ver = "hg38") {
  chr_arm_granges <- get_chr_arm_bins(genome = genome_ver)

  y <- rowRanges(sce)
  mcols(y)[colnames(sce)] <- assay(sce, assay_name)

  chr_arm_granges <- integrate_segments(x = chr_arm_granges, y, granges_signal_colname = colnames(sce))

  arm_sce <- SingleCellExperiment(assays = list("counts" = as.matrix(mcols(chr_arm_granges)[colnames(sce)])), rowRanges = chr_arm_granges, colData = colData(sce))

  # Silly workaround
  assay(arm_sce, assay_name) <- assay(arm_sce, "counts")
  assay(arm_sce, "counts") <- NULL

  rownames(arm_sce) <- get_bin_ids(rowRanges(arm_sce))

  return(arm_sce)
}
