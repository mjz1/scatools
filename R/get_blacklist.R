#' Get blacklist regions
#'
#' @param genome Genome version: hg38, hg19, mm10, ce10, or dm3
#'
#' @return Genomic ranges of blacklist regions
#' @export
#'
get_blacklist <- function(genome = NULL) {
  encodeBL <- c(
    "hg19" = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz",
    "hg38" = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz",
    "mm10" = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz",
    "mm9" = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/Blacklist_v1/mm9-blacklist.bed.gz",
    "ce10" = "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/ce10-C.elegans/ce10-blacklist.bed.gz",
    "dm3" = "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/dm3-D.melanogaster/dm3-blacklist.bed.gz"
  )

  if (tolower(genome) %in% names(encodeBL)) {
    bl <- tryCatch(
      {
        blacklist <- rtracklayer::import.bed(encodeBL[tolower(genome)])
      },
      error = function(x) {
        message("Blacklist not downloaded! Continuing without, be careful for downstream biases..")
        GRanges()
      }
    )
  } else {
    message("Blacklist not downloaded! Continuing without, be careful for downstream biases..")
    bl <- GRanges()
  }

  bl
}
