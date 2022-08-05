#' Load atac binned depth data
#'
#' Loads binned atac reads, merges cell-wise and bin-wise metadata, and performs QC.
#'
#' @param directory Folder containing samples
#' @param ArchR_Proj Optional: ArchR project with matching cells and metadata
#' @param bins Optional: GRanges bins object
#' @param BPPARAM Options to pass to `bplapply` for data loading. Provides minor speedup if loading many samples
#' @param save_to File path in which to save the final output. Note: Will still return the sce object for downstream analysis.
#' @param save_as NOT IMPLEMENTED YET. Select file formats to save the object. Can provide multiple values
#' @param verbose Message verbosity (TRUE/FALSE)
#'
#' @inheritParams DropletUtils::read10xCounts
#'
#' @return A `SingleCellExperiment` object.
#' @export
#'
load_atac_bins <- function(samples,
                           sample.names,
                           ArchR_Proj = NULL,
                           bins = NULL,
                           BPPARAM = BiocParallel::bpparam(),
                           save_to = NULL,
                           verbose = FALSE,
                           save_as = c("sce", "adata", "seurat")) {
  if (verbose) {
    logger::log_info("Loading bin counts in {length(samples)} samples using {BPPARAM$workers} threads")
  }

  sce <- DropletUtils::read10xCounts(samples = samples, sample.names = sample.names, col.names = TRUE, BPPARAM = BPPARAM)

  # Save raw counts in a seperate slot
  assay(sce, "raw_counts") <- assay(sce, "counts")

  # Reset the barcodes to remove prepended index (in case of multi-sample loading)
  colnames(sce) <- sce$Barcode

  # Additional processing if matching ArchR project is provided
  if (!is.null(ArchR_Proj)) {
    if (verbose) {
      logger::log_info("ArchR project provided. Merging cell metadata")
    }
    # Filter down these cells to be those in the ArchR project
    # TODO: Allow for custom filtering
    # TODO: Enable more flexible metadata provision
    sce <- sce[, which(colnames(sce) %in% rownames(ArchR_Proj@cellColData))]

    if (!all(rownames(ArchR_Proj@cellColData[colnames(sce), ]) == colnames(sce))) {
      logger::log_error("Cells in ArchR project not matching with sce object")
    }
    # Line up and merge
    colData(sce) <- cbind(colData(sce), subset(ArchR_Proj@cellColData[colnames(sce), ], select = -c(Sample)))
  }

  # Merge bin level information if provided
  if (!is.null(bins)) {
    # Merge bin level information as GRanges
    rowRanges(sce) <- sort(GenomicRanges::makeGRangesFromDataFrame(merge(rowData(sce), as.data.frame(bins), by.x = "ID", by.y = "bin_id"), keep.extra.columns = TRUE))
    rownames(sce) <- rowData(sce)$ID

    # Grab these values from the bins df
    if (!is.null(bins$gc)) {
      gc <- bins$gc
    }
    if (!is.null(bins$n_freq)) {
      n_freq <- bins$n_freq
    }

    # TODO: add bin arm information
    # rowRanges(sce)$chr_arm <- with(rowRanges(sce), interaction(seqnames, arm, sep = "", lex.order = T))
  }

  if (verbose) {
    logger::log_info("Adding cellwise and binwise QC metrics")
  }
  sce <- scuttle::addPerCellQCMetrics(sce)
  sce <- scuttle::addPerFeatureQCMetrics(sce, subsets = get_f_idx(sce$Sample))

  # Throw error if row names are not matching
  stopifnot(all(rownames(sce) == rowData(sce)$ID))

  if (!is.null(save_to)) {
    .save_to(object = sce, save_to = save_to, verbose = verbose)
  }

  if (verbose) {
    logger::log_success("Fragments loaded successfully!")
  }

  return(sce)
}



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
read_vartrix <- function(dir_path = NULL, mtx_ref = NULL, mtx_alt = NULL, barcodes = NULL, variants = NULL, phased_vcf = NULL, verbose = TRUE) {
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
    dplyr::mutate("chr" = as.factor(chr), "start" = as.integer(as.integer(pos) + 1), "end" = as.integer(as.integer(pos) + 1))

  colnames(ref) <- cells$barcodes
  colnames(alt) <- cells$barcodes

  # Load phasing if provided
  if (!is.null(phased_vcf)) {
    if (verbose) {
      logger::log_info("Loading phasing information from: {phased_vcf}")
    }

    phased_vcfR <- vcfR::read.vcfR(phased_vcf, verbose = verbose)

    phased_df <- cbind(phased_vcfR@fix, gt = phased_vcfR@gt[, 2]) %>%
      dplyr::as_tibble() %>%
      dplyr::select(CHROM, POS, REF, ALT, gt) %>%
      dplyr::rename(chr = CHROM, ref = REF, alt = ALT) %>%
      dplyr::mutate(start = as.integer(POS), end = as.integer(POS)) %>%
      dplyr::mutate(across(where(is.character), as.factor))

    # TODO: Check that variants are in common

    snps <- snps %>%
      dplyr::left_join(phased_df, by = c("chr", "start", "end")) %>%
      select(chr, start, end, ref, alt, gt)

    if (verbose) {
      logger::log_success("Phased VCF data loaded")
    }
  }

  sce <- SingleCellExperiment(list(ref = ref, alt = alt), rowRanges = GenomicRanges::makeGRangesFromDataFrame(snps, keep.extra.columns = TRUE), colData = cells)

  return(sce)
}
