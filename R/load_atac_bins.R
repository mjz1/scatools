#' Write a sparse count matrix in 10x (v3) MatrixMarket format
#'
#' Base-R replacement for `DropletUtils::write10xCounts()` writing
#' `matrix.mtx.gz`, `barcodes.tsv.gz`, and `features.tsv.gz`.
#'
#' @noRd
write_10x_counts <- function(path, x, barcodes = colnames(x), gene.id = rownames(x),
                             gene.symbol = gene.id, gene.type = "Gene Expression",
                             overwrite = FALSE) {
  if (file.exists(file.path(path, "matrix.mtx.gz")) && !overwrite) {
    return(invisible(NULL))
  }
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  x <- methods::as(x, "CsparseMatrix")

  # matrix.mtx.gz (writeMM has no gzip option; write then gzip-copy)
  tmp <- tempfile(fileext = ".mtx")
  Matrix::writeMM(x, tmp)
  con <- gzfile(file.path(path, "matrix.mtx.gz"), "wb")
  writeBin(readBin(tmp, what = "raw", n = file.size(tmp)), con)
  close(con)
  unlink(tmp)

  bc <- gzfile(file.path(path, "barcodes.tsv.gz"), "wt")
  writeLines(as.character(barcodes), bc)
  close(bc)

  feats <- data.frame(gene.id, gene.symbol, gene.type, stringsAsFactors = FALSE)
  ft <- gzfile(file.path(path, "features.tsv.gz"), "wt")
  utils::write.table(feats, ft, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  close(ft)
  invisible(NULL)
}

#' Read a 10x (v3) MatrixMarket count directory into a SingleCellExperiment
#'
#' Base-R replacement for `DropletUtils::read10xCounts()`.
#'
#' @noRd
read_10x_counts <- function(path, sample.id = basename(path)) {
  pick <- function(stem) {
    gz <- file.path(path, paste0(stem, ".gz"))
    if (file.exists(gz)) gz else file.path(path, stem)
  }
  mfile <- pick("matrix.mtx")
  con <- if (grepl("\\.gz$", mfile)) gzfile(mfile) else mfile
  mat <- methods::as(Matrix::readMM(con), "CsparseMatrix")

  barcodes <- readLines(pick("barcodes.tsv")) # file() auto-decompresses .gz
  feats <- utils::read.delim(pick("features.tsv"), header = FALSE, stringsAsFactors = FALSE)
  while (ncol(feats) < 3) feats[[ncol(feats) + 1]] <- feats[[1]]

  rownames(mat) <- feats[[1]]
  colnames(mat) <- barcodes

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = mat),
    rowData = S4Vectors::DataFrame(ID = feats[[1]], Symbol = feats[[2]], Type = feats[[3]]),
    colData = S4Vectors::DataFrame(Sample = sample.id, Barcode = barcodes)
  )
  rownames(sce) <- feats[[1]]
  colnames(sce) <- barcodes
  sce
}

#' Load atac binned depth data
#'
#' Loads binned atac reads, merges cell-wise and bin-wise metadata, and performs QC.
#'
#' @param bin_dir Directory with the bin counts
#' @param sample_id Sample ID
#' @param bins Optional: GRanges bins object
#' @param save_to File path in which to save the final output. Note: Will still return the sce object for downstream analysis.
#' @param verbose Message verbosity
#'
#' @return A `SingleCellExperiment` object.
#' @export
#'
load_atac_bins <- function(bin_dir,
                           sample_id,
                           bins = NULL,
                           save_to = NULL,
                           verbose = TRUE) {
  sce <- read_10x_counts(path = bin_dir, sample.id = sample_id)

  # Save raw counts in a seperate slot
  assay(sce, "raw_counts") <- assay(sce, "counts")

  # Reset the barcodes to remove prepended index (in case of multi-sample loading)
  colnames(sce) <- sce$Barcode

  # Merge bin level information if provided
  if (!is.null(bins)) {
    # Merge bin level information as GRanges
    rowRanges(sce) <- sort(GenomicRanges::makeGRangesFromDataFrame(merge(rowData(sce), as.data.frame(bins), by.x = "ID", by.y = "bin_id"), keep.extra.columns = TRUE))
    rownames(sce) <- rowData(sce)$ID
  }

  if (verbose) {
    cli::cli_alert_info("Adding cellwise and binwise QC metrics")
  }

  sce <- scuttle::addPerCellQCMetrics(sce)
  sce <- scuttle::addPerFeatureQCMetrics(sce, subsets = get_f_idx(sce$Sample))

  if (!is.null(save_to)) {
    save_to(object = sce, save_to = save_to, verbose = verbose)
  }

  if (verbose) {
    cli::cli_alert_success("Fragments loaded successfully!")
    print(sce)
  }

  return(sce)
}


#' Bin SNP Data
#'
#' @param snp_sce SCE with snp data
#' @param binsize Size of bins
#' @param select_chrs Chromosomes to include
#' @param bins Optional override of bins
#'
#' @return SCE object with phased binned snps
#' @export
#'
bin_snp_data2 <- function(snp_sce, binsize = 500000, select_chrs = NULL, bins = NULL) {
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
