
#' Calculate allelic information over cells and bins
#'
#'
#' @param snp SNP SingleCellExperiment object with slots 'ref' and 'alt' per cell counts
#' @param ncores Number of cores to speed up computation
#' @param bins Bins over which to aggregate allele information. If not provided will calculate for each SNP
#' @param group_var Column containing cell grouping information
#' @param FUN A function which can take two values (or vectors) of reference and alternate counts to produce another value (or vector)
#'
#' @return a snp sce
#' @export
#'
calc_allelic <- function(snp, ncores = 1, bins = NULL, group_var = NULL,
                         FUN = calc_ai) {

  if (is.null(bins)) {
    bins <- rowRanges(snp)
  }

  if (is.null(group_var)) {
    group_var = "all"
    snp$group_var <- "all"
  }

  # Get matching indices between snps and bins
  snp <- get_snp_bidx(snp, bins = bins)

  # Remove SNPs with no matching bin?
  snp <- snp[!is.na(rowRanges(snp)$bin_id),]

  # Get indices for each group split
  gr_split <- split(x = seq_along(snp[[group_var]]), f = snp[[group_var]])

  bin_split <- split(x = seq_along(rowRanges(snp)), f = rowRanges(snp)$bin_id)

  # Apply over cell groups
  res <- lapply(gr_split, FUN = function(gr) {
    ref = Matrix::rowSums(assay(snp[,gr], "ref"))
    alt = Matrix::rowSums(assay(snp[,gr], "alt"))

    # apply over bins multithreaded
    bin_res <- pbmcapply::pbmclapply(X = bin_split, mc.cores = ncores, FUN = function(b) {
      # Take the total count weighted mean across the bin
      z <- FUN(ref_counts = ref[b], alt_counts = alt[b])
      per_snp_cov <- ref[b] + alt[b]

      z_w <- weighted.mean(x = z, w = per_snp_cov, na.rm = TRUE)

      per_bin_cov <- sum(per_snp_cov)

      # bin_alt = sum(ref[b])
      # bin_ref = sum(alt[b])
      # z <- FUN(bin_ref, bin_alt)
      # z[is.na(z)] <- 0 # mask NAs as zero
      z_w <- as(z_w, "sparseMatrix")
      return(z_w)
    })
    bnames <- names(bin_res)
    bin_res <- do.call(rbind, bin_res)
    rownames(bin_res) <- bnames
    return(bin_res)
  })

  res <- do.call(cbind, res)
  colnames(res) <- names(gr_split)

  # Reorder results
  res <- res[gtools::mixedsort(rownames(res)),]

  sce <- SingleCellExperiment(list(res = res))

  rowRanges(sce) <- bins[match(rownames(sce), names(bins)),]
  return(sce)
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


