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
#' @return A `SingleCellExperiment` object.
#' @export
#'
load_atac_bins <- function(directory,
                           ArchR_Proj = NULL,
                           bins = NULL,
                           BPPARAM = BiocParallel::bpparam(),
                           save_to = NULL,
                           verbose = FALSE,
                           save_as = c("sce", "adata", "seurat")) {
  # Get sample names and directories
  sample_names <- dir(directory)
  samples <- dir(directory, full.names = TRUE)

  if (verbose) {
    logger::log_info("Loading bin counts")
  }
  sce <- DropletUtils::read10xCounts(samples = samples, sample.names = sample_names, col.names = TRUE, BPPARAM = BPPARAM)

  # Save raw counts in a seperate slot
  assay(sce, "raw_counts") <- assay(sce, "counts")

  # Reset the barcodes (since the above adds index to each sample)
  colnames(sce) <- sce$Barcode

  # Additional processing if matching ArchR project is provided
  if (!is.null(ArchR_Proj)) {
    if (verbose) {
      logger::log_info("ArchR project provided. Merging cell metadata")
    }
    # Filter down these cells to match the ArchR project (and in correct order)
    sce <- sce[, match(rownames(ArchR_Proj@cellColData), colnames(sce))]

    # Line up and merge
    colData(sce) <- cbind(colData(sce), subset(ArchR_Proj@cellColData[colnames(sce), ], select = -c(Sample)))

    # Throw error if the merging is imperfect
    stopifnot(all(sce$Barcode == rownames(ArchR_Proj@cellColData)))
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
    dplyr::mutate("chr" = as.factor(chr), "pos" = as.integer(as.integer(pos) + 1))

  colnames(ref) <- cells$barcodes
  colnames(alt) <- cells$barcodes

  sce <- SingleCellExperiment(list(ref = ref, alt = alt), rowData = snps, colData = cells)

  return(sce)
}
