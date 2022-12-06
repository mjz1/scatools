#' Read Vatrix Data
#'
#' Reads in `Vartix` scSNP pileups. Assumes Vartrix is run in `-s coverage` mode, outputting barcodes and variants table. If no `dir_path` is provided, then user must provide individual file paths for each file. If phased VCF and input VCF is provided this function will also remove any non-heterozygous loci, and indels.
#'
#' @param dir_path Path to `Vartrix` outputs directory
#' @param mtx_ref Path to reference allele `.mtx` file
#' @param mtx_alt Path to alternate allele `.mtx` file
#' @param barcodes Path to barcodes file
#' @param variants Path to variants file
#' @param input_vcf Germline VCF before phasing
#' @param phased_vcf Phased VCF
#' @param verbose Verbosity
#'
#' @return A `SingleCellExperiment` object
#' @export
read_vartrix <- function(dir_path = NULL, mtx_ref = NULL, mtx_alt = NULL, barcodes = NULL, variants = NULL, input_vcf = NULL, phased_vcf = NULL, verbose = FALSE) {
  if (!is.null(dir_path)) {
    mtx_ref <- dir(dir_path, pattern = "ref", full.names = TRUE)
    mtx_alt <- dir(dir_path, pattern = "alt", full.names = TRUE)
    barcodes <- dir(dir_path, pattern = "_barcodes", full.names = TRUE)
    variants <- dir(dir_path, pattern = "_variants", full.names = TRUE)
  }

  logger::log_info("Loading SNP matrices")
  ref <- as(Matrix::readMM(mtx_ref), "CsparseMatrix")
  alt <- as(Matrix::readMM(mtx_alt), "CsparseMatrix")
  cells <- read.table(barcodes, header = FALSE, col.names = "barcodes")

  snps <- read.table(variants, header = FALSE)[[1]] %>%
    as.vector() %>%
    stringr::str_split_fixed(string = ., pattern = "_", n = 2) %>%
    as.data.frame()
  colnames(snps) <- c("chr", "pos")
  # Big memory savings plus we rebase the positions to match original VCF
  snps <- snps %>%
    dplyr::mutate("chr" = as.factor(chr), "start" = as.integer(as.integer(pos) + 1), "end" = as.integer(as.integer(pos) + 1)) %>%
    dplyr::mutate("snp_id" = paste(chr, start, sep = "_")) %>%
    select(snp_id, chr, start, end)

  colnames(ref) <- colnames(alt) <- cells$barcodes
  rownames(ref) <- rownames(alt) <- snps$snp_id

  # Load phasing if provided
  if (!is.null(phased_vcf)) {
    logger::log_info("Reading input vcf: {input_vcf}")
    vcf_df <- vcf_to_df(input_vcf)
    logger::log_info("Reading phased vcf: {phased_vcf}")
    phased_vcf_df <- vcf_to_df(phased_vcf)
    phased_vcf_df <- dplyr::rename(phased_vcf_df, "gt_phased" = "gt")

    vcf_df <- vcf_df %>%
      dplyr::left_join(phased_vcf_df, by = c("snp_id", "chr", "start", "end"), suffix = c("", "_phased")) %>%
      dplyr::filter(ref != "*", alt != "*")

    # Merge
    snps <- snps %>%
      filter(snp_id %in% vcf_df$snp_id) %>%
      dplyr::left_join(vcf_df)

    logger::log_success("VCF data loaded and merged")
  }

  # Reorder
  snps <- snps %>%
    dplyr::mutate(chr = factor(chr, levels = chr_reorder(unique(levels(chr))))) %>%
    dplyr::arrange(chr, start)

  # Filter ref and alt matrices and reorder at the same time
  ref <- ref[snps$snp_id, ]
  alt <- alt[snps$snp_id, ]

  sce <- SingleCellExperiment(list(ref = ref, alt = alt), rowRanges = GenomicRanges::makeGRangesFromDataFrame(snps, keep.extra.columns = TRUE), colData = cells)

  return(sce)
}

#' @export
vcf_to_df <- function(vcf, verbose = FALSE) {
  vcf <- vcfR::read.vcfR(file = vcf, verbose = verbose)
  gts <- vcfR::extract.gt(vcf, element = "GT", return.alleles = F)[, 1] %>% as.factor()
  alleles <- vcfR::extract.gt(vcf, element = "GT", return.alleles = T)[, 1] %>% stringr::str_split_fixed(string = ., pattern = "/|\\|", n = 2)
  stopifnot(all(names(gts) == names(alleles)))

  # Remove indels and non het SNPs
  keeps <- names(which((!vcfR::is.indel(vcf) & vcfR::is_het(as.matrix(gts)))[, 1]))
  vcf_df <- cbind.data.frame(gts, alleles)[keeps, ] %>%
    dplyr::mutate(across(where(is.character), as.factor))
  colnames(vcf_df) <- c("gt", "ref", "alt")
  vcf_df$snp_id <- rownames(vcf_df)

  pos <- stringr::str_split_fixed(vcf_df$snp_id, pattern = "_", n = 2)
  vcf_df$chr <- as.factor(pos[, 1])
  vcf_df$start <- vcf_df$end <- as.integer(pos[, 2])

  # Reordering columns and sorting
  vcf_df <- vcf_df %>%
    dplyr::select(snp_id, chr, start, end, ref, alt, gt) %>%
    dplyr::mutate(chr = factor(chr, levels = chr_reorder(unique(levels(chr))))) %>%
    dplyr::arrange(chr, start)

  return(vcf_df)
}
