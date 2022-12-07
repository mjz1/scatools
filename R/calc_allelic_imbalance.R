
#' Calculate SNP allelic imbalance
#'
#'
#' @param snp SNP SingleCellExperiment object
#' @param ncores Number of cores to speed up computation
#'
#' @return a snp sce
#' @export
#'
calc_snp_ai <- function(snp, ncores = 1) {
  assay(snp, "total") <- assay(snp, "ref") + assay(snp, "alt")
  keep_snps <- names(which(Matrix::rowSums(assay(snp, "total")) > 0))
  snp <- snp[keep_snps, ]
  # To ensure memory efficiency, we need to encode a sparse binary matrix that
  # will mask any SNPs that are 0,0 in ref and alt in a cell. This will allow
  # us to convert NAs to zeroes in sparse form, and ensure calculations are
  # only performed on true data
  assay(snp, "cov") <- (assay(snp, "ref") != 0 & assay(snp, "alt") != 0)
  x <- assay(snp, "ref")
  y <- assay(snp, "alt")
  res <- pbmcapply::pbmclapply(X = 1:ncol(snp), mc.cores = ncores, FUN = function(i) {
    z <- calc_ai(ref_counts = x[, i], alt_counts = y[, i])
    z[is.na(z)] <- 0 # mask NAs as zero
    z <- as(z, "sparseMatrix")
    return(z)
  })
  res <- do.call(cbind, res)

  colnames(res) <- colnames(snp)
  rownames(res) <- rownames(snp)

  assay(snp, "ai") <- res
  return(snp)
}

#' Calculate allelic imbalance
#'
#' @param ref_counts vector of ref counts
#' @param alt_counts vector of alt counts
#'
#' @return measure of allelic imbalance
#' @export
#'
#' @examples
#' calc_ai(c(4, 2, 1), c(2, 3, 1))
calc_ai <- function(ref_counts, alt_counts) {

  # coerce to vectors
  ref_counts <- as.vector(ref_counts)
  alt_counts <- as.vector(alt_counts)

  x <- pmax(ref_counts, alt_counts)
  y <- pmin(ref_counts, alt_counts)
  ai <- (x - y) / x

  return(ai)
}


