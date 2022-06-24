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

#' Get tiled bins
#'
#' @param bs_genome A BSgenome object
#' @param tilewidth Bin size
#' @param select_chrs Vector of chromosomes to include
#'
#' @return A GRanges object of bins
#' @export
#'
#' @examples
#' \dontrun{
#' bins <- get_tiled_bins(BSgenome.Hsapiens.UCSC.hg38, tilewidth = 500000)
#' }
get_tiled_bins <- function(bs_genome, tilewidth = 500000, select_chrs = NULL) {
  stopifnot(class(bs_genome) %in% "BSgenome")

  if (is.null(select_chrs)) {
    select_chrs <- paste("chr", c(1:22, "X"), sep = "")
  }

  bins <- GenomicRanges::tileGenome(BSgenome::seqinfo(bs_genome), tilewidth = tilewidth, cut.last.tile.in.chrom = TRUE)
  bins <- add_gc_freq(bs_genome, bins)
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
  stopifnot(class(bs_genome) %in% "BSgenome")
  freqs <- BSgenome::alphabetFrequency(BSgenome::getSeq(bs_genome, bins))
  bins$gc <- (freqs[,'C'] + freqs[,'G'])/rowSums(freqs)
  return(bins)
}
