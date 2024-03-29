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
    # TODO: check this line for error:
    # Error in .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProject",  :
    # Input value for 'ArchRProj' is not a null,archrproject, (ArchRProj = function) please supply valid input!
    colData(sce) <- cbind(colData(sce), subset(ArchR_Proj@cellColData[colnames(sce), ], select = -c(Sample)))

    # Merge blacklist information if present
    try(metadata(sce)$blacklist <- ArchR::getBlacklist(ArchR_Proj))
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
    print(sce)
  }

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

#' Title
#'
#' @param snp_granges SNP granges object with GT column
#' @param binsize Size of bins
#' @param select_chrs Chromosomes to include
#' @param bins Optional override of bins
#'
#' @return SCE object with phased binned snps
#' @export
#'
bin_snp_data <- function(snp_sce, binsize = 500000, select_chrs = NULL, bins = NULL) {
  # THIS FUNCTION NEEDS WORK
  # TODO
  if (is.null(select_chrs)) {
    select_chrs <- paste("chr", c(1:22, "X"), sep = "")
  }

  # Create bins
  if (is.null(bins)) {
    bins <- get_tiled_bins(bs_genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, tilewidth = binsize, select_chrs = select_chrs)
  } else {
    bins <- bins[seqnames(bins) %in% select_chrs]
  }

  # Pull out snp granges
  snp_granges <- rowRanges(snp_sce)

  # Then overlap with findoverlaps
  hits <- GenomicRanges::findOverlaps(bins, snp_granges)

  # Add the bin index for aggregation
  snp_granges$bin_idx <- S4Vectors::queryHits(hits)

  # For each cell we want to aggregate the SNP depths per cell correcting for phasing
  # First get indices of the 0|1 vs 1|0 gts
  test <- snp_sce[, 1]

  colSums(assay(test, "ref"))

  which(snp_granges$gt == "0|1")
  which(snp_granges$gt == "1|0")

  # get the bin_ids
  bin_ids <- as.data.frame(bins) %>%
    select(seqnames, start, end) %>%
    unite("bin_id") %>%
    pull()

  # Have to aggregate per cell and end up with a matrix same shape as the bins
  snp_bins <- snp_granges %>%
    as_tibble() %>%
    filter(GT != "1|1") %>%
    # Adjust the phasing
    mutate(AD_phased = ifelse(GT == "1|0", DP - AD, AD)) %>%
    group_by(cell, bin_idx, GT) %>%
    summarise(
      DP = sum(DP),
      AD = sum(AD),
      AD_phased = sum(AD_phased),
      n_snps = n()
    ) %>%
    mutate(BAF = AD / DP)


  snp_bins$bin_id <- bin_ids[snp_bins$bin_idx]

  snp_bins <- snp_bins %>%
    separate(col = bin_id, into = c("chr", "start", "end"), sep = "_", remove = F)

  return(snp_bins)
}
