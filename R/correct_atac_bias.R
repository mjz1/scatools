#' Correct ATAC Bias
#'
#' Correct the ATAC signal bias in a cancer sample using a bulk copy number profile. First, the copy number profile of interest (in granges format) is integrated with the bin sizes provided in the input `SingleCellExperiment` object. The copy number profile is then divided out of the input counts and averaged across all cells to approximate the amount of atac specific bias. The input counts are then divided by this bias measure to yield corrected counts.
#'
#'
#'
#' @param sce Single cell experiment object
#' @param assay_name Name of assay to correct
#' @param corrected_name Name of corrected assay
#' @param cn_granges GRanges object containing bulk sample copy number information
#' @param granges_signal_colname Column in `cn_granges` that contains the copy number information.
#' @param drop_missing_bins Logical: Drop bins missing information in `cn_granges`.
#'
#' @return A `SingleCellExperiment` object containing slots `bias` and `corrected_counts`.
#' @export
#'
correct_atac_bias <- function(sce, assay_name, corrected_name = "corrected_counts", cn_granges, granges_signal_colname, drop_missing_bins = FALSE) {

  # Integrate the true CN data to our data
  rowRanges(sce) <- integrate_segments(x = rowRanges(sce), y = cn_granges, granges_signal_colname = granges_signal_colname, drop_na = drop_missing_bins)

  assay(sce, "bias") <- assay(sce, assay_name) / mcols(rowRanges(sce))[[granges_signal_colname]]

  mcols(rowRanges(sce))["mean_bias"] <- apply(assay(sce, "bias"), 1, mean)

  assay(sce, corrected_name) <- assay(sce, assay_name) / mean_bias

  if (drop_missing_bins) {
    sce <- sce[!is.na(mcols(rowRanges(sce))[[granges_signal_colname]]), ]
  }

  return(sce)
}
