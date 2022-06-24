#' Get chromosome arm bins
#'
#' @param genome Genome version ('hg38', 'hg19')
#'
#' @return A GRanges object of chromosome arm bins
#' @export
#'
#' @examples
#' bins <- get_chr_arm_bins('hg38')
get_chr_arm_bins <- function(genome = 'hg38') {
  bins <- get_cytobands() %>%
    dplyr::group_by(CHROM, arm, genome) %>%
    dplyr::summarise(
      start = min(start),
      end = max(end)
    ) %>%
    dplyr::mutate(bin_id = paste(CHROM, start, end, sep = "_")) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  return(bins)
}

#' Get genome cytobands
#'
#' @param genome Genome version (hg38 or hg19)
#'
#' @return Dataframe of genome cytobands
#' @export
#'
#' @examples
#' hg38_cyto <- get_cytobands('hg38')
get_cytobands <- function(genome = 'hg38') {
  cyto_url <- paste0("http://hgdownload.cse.ucsc.edu/goldenpath/", genome, "/database/cytoBand.txt.gz")
  cyto <- vroom::vroom("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz",
              col_names = c("CHROM", "start", "end", "cytoband", "unsure")
  ) %>%
    dplyr::filter(!is.na(cytoband)) %>%
    dplyr::mutate(across(where(is.character), as.factor),
                  start = start + 1,
                  arm = factor(substr(cytoband, 0, 1)),
                  genome = genome)
  return(cyto)
}


#' Add GC frequency
#'
#' @param bs_genome BSGenome object
#' @param bins GRanges bins object
#'
#' @return GC frequency per bin
#' @export
#'
add_gc_freq <- function(bs_genome, bins) {
  freqs <- Biostrings::alphabetFrequency(Biostrings::getSeq(genome, bins))
  bins$gc <- (freqs[,'C'] + freqs[,'G'])/rowSums(freqs)
  return(bins)
}
