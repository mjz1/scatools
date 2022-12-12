#' Join SNP and scATAC depth
#'
#' @inheritParams calc_allelic
#' @inheritParams pseudo_groups
#' @param sce sce obj
#' @param snp snp obj
#' @param sce_assay sce_assay to aggr
#' @param sce_aggr_fun fx for aggr
#'
#' @return joint obj
#' @export
#'
pseudo_join <- function(sce, snp, sce_assay = "counts", ncores = 1, group_var = NULL, min_cov = 5, sce_aggr_fun = mean) {
  snp_ai <- calc_allelic(snp, ncores = 8, bins = rowRanges(sce), group_var = group_var, min_cov = 5)
  sce_cov <- pseudo_groups(sce = sce, assay_name = sce_assay, group_var = group_var, FUN = mean, na.rm = T)

  joint <- snp_ai
  colnames(joint) <- colnames(sce_cov)
  assay(joint, sce_assay) <- assay(sce_cov)

  joint$id <- colnames(joint)

  return(joint)
}
