#' Calculate allelic information over cells and bins
#'
#'
#' @param snp SNP SingleCellExperiment object with slots 'ref' and 'alt' per cell counts
#' @param ncores Number of cores to speed up computation
#' @param bins Bins over which to aggregate allele information. If not provided will calculate for each SNP
#' @param group_var Column containing cell grouping information
#' @param FUN A function which can take two values (or vectors) of reference and alternate counts to produce another value (or vector)
#' @param min_cov Minumum coverage per SNP across pseudobulked cells to be included in the calculation
#'
#' @return a snp sce with slots:
#' * `snp` - Containing the main result
#' * `tot_cov` - Total SNP coverage in the bin for the set of cells after filtering for `min_cov`
#' * `n_snp` - Number of SNPs remaining in the bin after filtering for `min_cov`
#'
#' @export
#'
calc_allelic <- function(snp, ncores = 1,
                         bins = NULL,
                         group_var = NULL,
                         FUN = calc_ai,
                         min_cov = 1) {
  if (is.null(bins)) {
    bins <- rowRanges(snp)
  }

  if (is.null(group_var)) {
    group_var <- "all"
    snp$all <- "all"
  }

  # Get matching indices between snps and bins
  snp <- get_snp_bidx(snp, bins = bins)

  # Remove SNPs with no matching bin?
  snp <- snp[!is.na(rowRanges(snp)$bin_id), ]

  # Get indices for each group split
  gr_split <- split(x = seq_along(snp[[group_var]]), f = snp[[group_var]])

  bin_split <- split(x = seq_along(rowRanges(snp)), f = rowRanges(snp)$bin_id)

  # Apply over cell groups
  res <- lapply(seq_along(gr_split), FUN = function(i) {
    gr <- gr_split[[i]]
    gr_name <- names(gr_split[i])
    ref <- Matrix::rowSums(assay(snp[, gr], "ref"))
    alt <- Matrix::rowSums(assay(snp[, gr], "alt"))

    # apply over bins multithreaded
    bin_res <- pbmcapply::pbmclapply(X = bin_split, mc.cores = ncores, FUN = function(b) {
      # Compute the per SNP coverage within the bin
      per_snp_cov <- ref[b] + alt[b]

      s_idx <- per_snp_cov >= min_cov

      # Take the total count weighted mean across the bin
      z <- FUN(ref_counts = ref[b][s_idx], alt_counts = alt[b][s_idx])

      z_w <- weighted.mean(x = z, w = per_snp_cov[s_idx], na.rm = TRUE)

      per_bin_cov <- sum(per_snp_cov[s_idx])

      # z_w <- as(z_w, "sparseMatrix")
      return(list(mean = z_w, n_snp = sum(s_idx), tot_cov = per_bin_cov))
    })
    bin_res <- do.call(rbind.data.frame, bin_res)
    bin_res$bin_id <- rownames(bin_res)
    bin_res[[group_var]] <- gr_name
    rownames(bin_res) <- NULL
    return(bin_res)
  })

  res <- do.call(rbind.data.frame, res)

  # Format results
  mean_mat <- pivot_wider(res, id_cols = "bin_id", names_from = group_var, values_from = "mean") %>%
    column_to_rownames("bin_id")
  tot_cov_mat <- pivot_wider(res, id_cols = "bin_id", names_from = group_var, values_from = "tot_cov") %>%
    column_to_rownames("bin_id")
  n_snp_mat <- pivot_wider(res, id_cols = "bin_id", names_from = group_var, values_from = "n_snp") %>%
    column_to_rownames("bin_id")

  sce <- SingleCellExperiment(list("snp" = mean_mat, "tot_cov" = tot_cov_mat, "n_snp" = n_snp_mat))

  sce <- sce[gtools::mixedsort(rownames(sce)), ]

  metadata(sce)$calc_allelic <- list(FUN = FUN, min_cov = min_cov)
  rowRanges(sce) <- bins[match(rownames(sce), get_bin_ids(bins)), ]

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
