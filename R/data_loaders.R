#' Read Vatrix Data
#'
#' Reads in `Vartix` scSNP pileups. Assumes Vartrix is run in `-s coverage` mode, outputting barcodes and variants table.
#'
#' If no `dir_path` is provided, then user must provide individual file paths for each file.
#'
#' @param dir_path Path to `Vartrix` outputs directory
#' @param mtx_ref Path to reference allele `.mtx` file
#' @param mtx_alt Path to alternate allele `.mtx` file
#' @param barcodes Path to barcodes file
#' @param variants Path to variants file
#'
#' @return A `SingleCellExperiment` object
#' @export
read_vartrix <- function(dir_path = NULL, mtx_ref, mtx_alt, barcodes, variants) {
  if (!is.null(dir_path)) {
    mtx_ref <- file.path(dir_path, "ref_matrix.mtx")
    mtx_alt <- file.path(dir_path, "alt_matrix.mtx")
    barcodes <- file.path(dir_path, "barcodes.tsv")
    variants <- file.path(dir_path, "variants.tsv")
  }

  ref <- as(Matrix::readMM(mtx_ref), "dgCMatrix")
  alt <- as(Matrix::readMM(mtx_alt), "dgCMatrix")
  cells <- read.table(barcodes, header = FALSE, col.names = "barcodes")
  snps <- read.table(variants, header = FALSE, col.names = "pos") %>%
    tidyr::separate(pos, into = c("chr", "pos"), sep = "_") %>%
    # Big memory savings plus we rebase the positions to match original VCF
    dplyr::mutate("chr" = as.factor(chr), "pos" = as.integer(as.integer(pos)+1))

  colnames(ref) <- cells$barcodes
  colnames(alt) <- cells$barcodes

  sce <- SingleCellExperiment(list(ref=ref, alt=alt), rowData = snps, colData = cells)
  return(sce)
}
