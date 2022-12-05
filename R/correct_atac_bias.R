#' Correct ATAC Bias
#'
#' Correct the ATAC signal bias in a cancer sample using a bulk copy number profile. First, the copy number profile of interest (in granges format) is integrated with the bin sizes provided in the input `SingleCellExperiment` object. The copy number profile is then divided out of the input counts and averaged across all cells to approximate the amount of atac specific bias. The input counts are then divided by this bias measure to yield corrected counts.
#'
#'
#'
#' @param sce Single cell experiment object
#' @param assay_name Name of assay to correct
#' @param corrected_name Name of corrected assay
#' @param cn_granges GRanges object containing bulk sample copy number information. Only required if data was not already integrated.
#' @param granges_signal_colname Column in `cn_granges` that contains the copy number information.
#' @param drop_missing_bins Logical: Drop bins missing information in `cn_granges`.
#'
#' @return A `SingleCellExperiment` object containing slots `bias` and `counts_corrected`.
#' @export
#'
correct_atac_bias <- function(sce, assay_name, corrected_name = "counts_corrected", cn_granges = NULL, granges_signal_colname, drop_missing_bins = FALSE) {
  # Get bin_ids to keep track in case of dropped bins
  bin_ids <- get_bin_ids(rowRanges(sce))

  # Integrate the true CN data to our data
  if (!is.null(cn_granges)) {
    integ_ranges <- integrate_segments(x = rowRanges(sce), y = cn_granges, granges_signal_colname = granges_signal_colname, drop_na = drop_missing_bins)

    new_bin_ids <- get_bin_ids(integ_ranges)

    # Safely subset based on new bin ids
    sce <- sce[new_bin_ids, ]

    rowRanges(sce) <- integ_ranges
    rownames(sce) <- new_bin_ids
  }

  assay(sce, "bias") <- assay(sce, assay_name) / mcols(rowRanges(sce))[[granges_signal_colname]]

  mcols(rowRanges(sce))["mean_bias"] <- apply(assay(sce, "bias"), 1, mean, na.rm = TRUE)

  assay(sce, corrected_name) <- sweep(x = assay(sce, assay_name), MARGIN = 1, STATS = mcols(rowRanges(sce))[["mean_bias"]], FUN = "/")

  # Keep rownames
  # rownames(sce) <- get_bin_ids(rowRanges(sce))

  return(sce)
}
