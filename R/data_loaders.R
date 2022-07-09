#' Load atac binned depth data
#'
#' Loads binned atac reads, merges cell-wise and bin-wise metadata, and performs QC.
#'
#' @param directory Folder containing samples
#' @param ArchR_Proj Optional: ArchR project with matching cells and metadata
#' @param bins Optional: GRanges bins object
#' @param BPPARAM Options to pass to `bplapply` for data loading. Provides minor speedup if many samples
#' @param flag_ideal Logical. Whether or not to flag ideal bins for each cell. Can be slow if not parallelized. See [get_ideal_mat()] for more information.
#' @param save_to File path in which to save the final output. Note: Will still return the sce object for downstream analysis.
#' @param save_as NOT IMPLEMENTED YET. Select file formats to save the object. Can provide multiple values
#'
#' @inheritParams get_ideal_mat
#'
#' @return A `SingleCellExperiment` object.
#' @export
#'
load_atac_bins <- function(directory, ArchR_Proj = NULL, bins = NULL, BPPARAM = bpparam(), gc = NULL, n_freq = NULL, map = 1, min_reads = 1, max_N_freq = 0.05, reads_outlier = 0.01, gc_outlier = 0.001, min_map = 0.9, ncores = 1, flag_ideal = TRUE, save_to = NULL, save_as = c("sce", "adata", "seurat")) {
  sample_names <- dir(directory)
  samples <- dir(directory, full.names = TRUE)

  message("Loading bin counts")
  sce <- DropletUtils::read10xCounts(samples = samples, sample.names = sample_names, col.names = TRUE, BPPARAM = BPPARAM)

  # Save raw counts in a seperate slot
  assay(sce, "raw_counts") <- assay(sce, "counts")

  # Reset the barcodes (since the above adds index to each sample)
  colnames(sce) <- sce$Barcode

  # Additional processing if matching ArchR project is provided
  if (!is.null(ArchR_Proj)) {
    message("ArchR project provided. Merging cell metadata")
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

  # Throw error if row names are not matching
  stopifnot(all(rownames(sce) == rowData(sce)$ID))


  message("Performing bin-length normalization. Storing as assay(sce, 'counts_permb')")
  # Bin length normalize upfront
  # Mainly necessary for our chr arm analysis where bins are variable in size
  # This is effectively an reads per Megabase (RPBMb) calculation. For reminder: https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

  # Add binwidth if not provided
  if (is.null(rowData(sce)$binwidth)) {
    spl <- stringr::str_split(rowData(sce)$ID, "_", simplify = T)
    rowData(sce)$binwidth <- as.numeric(spl[, 3]) - as.numeric(spl[, 2]) + 1
  }

  assay(sce, "counts_permb") <- round(assay(sce, "raw_counts") / (rowData(sce)$binwidth / 1e6), 2)

  # Perform initial QC
  message("Adding cellwise and binwise QC metrics")
  sce <- scuttle::addPerCellQCMetrics(sce)
  sce <- scuttle::addPerFeatureQCMetrics(sce, subsets = get_f_idx(sce$Sample))

  # Add valid and ideal information to mirror HMMcopy's QC
  if (flag_ideal) {
    message("Flagging ideal bins using ", ncores, " threads")

    id_val_mats <- get_ideal_mat(
      mat = assay(sce, "counts_permb"),
      gc = gc,
      n_freq = n_freq,
      map = map,
      min_reads = min_reads,
      max_N_freq = max_N_freq,
      reads_outlier = reads_outlier,
      gc_outlier = gc_outlier,
      min_map = min_map,
      ncores = ncores
    )

    assay(sce, "ideal_bins") <- id_val_mats$ideal
    assay(sce, "valid_bins") <- id_val_mats$valid
  }

  message("Computing library size factors")
  sce <- scuttle::computeLibraryFactors(sce)

  cat("\n")
  print(sce)

  if (!is.null(save_to)) {
    dir.create(dirname(save_to), recursive = TRUE, showWarnings = FALSE)
    cat("\nSaving sce object to ", '"', save_to, '"', "\n", sep = "")
    save(sce, file = save_to)
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
