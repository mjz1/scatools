#' Title
#'
#' @param snp_granges SNP granges object with GT column
#' @param binsize Size of bins
#' @param select_chrs Chromosomes to include
#' @param bins Optional override of bins
#'
#' @return SCE object with phased binned snps
#' @export
#'
bin_snp_data <- function(snp_sce, binsize = 500000, select_chrs = NULL, bins = NULL) {
  # THIS FUNCTION NEEDS WORK
  # TODO
  if (is.null(select_chrs)) {
    select_chrs <- paste("chr", c(1:22, "X"), sep = "")
  }

  # Create bins
  if (is.null(bins)) {
    bins <- get_tiled_bins(bs_genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, tilewidth = binsize, select_chrs = select_chrs)
  } else {
    bins <- bins[seqnames(bins) %in% select_chrs]
  }

  # Pull out snp granges
  snp_granges <- rowRanges(snp_sce)

  # Then overlap with findoverlaps
  hits <- GenomicRanges::findOverlaps(bins, snp_granges)

  # Add the bin index for aggregation
  snp_granges$bin_idx <- S4Vectors::queryHits(hits)

  # For each cell we want to aggregate the SNP depths per cell correcting for phasing
  # First get indices of the 0|1 vs 1|0 gts
  test <- snp_sce[, 1]

  colSums(assay(test, "ref"))

  which(snp_granges$gt == "0|1")
  which(snp_granges$gt == "1|0")

  # get the bin_ids
  bin_ids <- as.data.frame(bins) %>%
    select(seqnames, start, end) %>%
    unite("bin_id") %>%
    pull()

  # Have to aggregate per cell and end up with a matrix same shape as the bins
  snp_bins <- snp_granges %>%
    as_tibble() %>%
    filter(GT != "1|1") %>%
    # Adjust the phasing
    mutate(AD_phased = ifelse(GT == "1|0", DP - AD, AD)) %>%
    group_by(cell, bin_idx, GT) %>%
    summarise(
      DP = sum(DP),
      AD = sum(AD),
      AD_phased = sum(AD_phased),
      n_snps = n()
    ) %>%
    mutate(BAF = AD / DP)


  snp_bins$bin_id <- bin_ids[snp_bins$bin_idx]

  snp_bins <- snp_bins %>%
    separate(col = bin_id, into = c("chr", "start", "end"), sep = "_", remove = F)

  return(snp_bins)
}
