#' Summarise copy number signal over chromosome arms
#'
#' @param sce SingleCellExperiment objet
#' @param assay_name Name of assay
#' @param genome_ver Genome version ("hg38" or "hg19")
#' @param cn_granges Optional: Copy number granges to integrate into the arm level bins
#' @param granges_signal_colname Optional: Name of the column containing the CNA information in `cn_granges`
#'
#' @return A chromosome arm level binned SingleCellExperiment
#' @export
#'
summarise_chr_arm <- function(sce, assay_name, cn_granges = NULL, granges_signal_colname = NULL, genome_ver = "hg38") {
  # TODO: Make this more general to aggregate signal over defined bins

  chr_arm_granges <- get_chr_arm_bins(genome = genome_ver)

  y <- rowRanges(sce)
  mcols(y)[colnames(sce)] <- assay(sce, assay_name)

  chr_arm_granges <- integrate_segments(x = chr_arm_granges, y, granges_signal_colname = colnames(sce))

  # Split out the count matrix from the ranges data
  counts <- as.matrix(mcols(chr_arm_granges)[colnames(sce)])
  mcols(chr_arm_granges)[colnames(sce)] <- NULL # Remove all the matrix columns

  arm_sce <- SingleCellExperiment(assays = list("counts" = counts), rowRanges = chr_arm_granges, colData = colData(sce))

  if (!is.null(cn_granges)) {
    rowRanges(arm_sce) <- integrate_segments(x = rowRanges(arm_sce), y = cn_granges, granges_signal_colname = granges_signal_colname, drop_na = FALSE)
  }

  # Silly workaround
  assay(arm_sce, assay_name) <- assay(arm_sce, "counts")
  assay(arm_sce, "counts") <- NULL

  rownames(arm_sce) <- get_bin_ids(rowRanges(arm_sce))

  return(arm_sce)
}
