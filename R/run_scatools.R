#' Run scATAC Copy Number Analysis
#'
#' Convenience wrapper to run a sample
#'
#' @param sample_id Sample id
#' @param fragment_file Fragments file
#' @param cells A vector of cell barcodes. If not provided will use all barcodes
#' @param bins Genomic ranges object with bins to use
#' @param bin_name Name for the bins (e.g. 10Mb, 500Kb)
#' @param blacklist Genome blacklist regions to filter against
#' @param outdir Output directory
#' @param rmsb_size Remove small bins below a given size. Useful in cases where small bins are leftover at the ends of chromosomes. Defaults to 10% of the binwidth.
#' @param gc_range GC range for bins to keep. Removes large GC outliers
#' @param overwrite Whether to overwrite previous results
#' @param verbose Verbosity
#' @param save_h5ad Logical. Whether to save raw and processed h5ad files. Requires packages `zellkonverter` and `anndata`
#' @param ncores Number of cores
#'
#' @return `SingleCellExperiment` object
#' @export
#'
run_scatools <- function(sample_id,
                         fragment_file,
                         cells = NULL,
                         bins = bins_10mb,
                         bin_name = prettyMb(getmode(width(bins))),
                         blacklist = NULL,
                         outdir = sample_id,
                         rmsb_size = 0.1*getmode(width(bins)),
                         gc_range = c(0.3, 0.8),
                         overwrite = FALSE,
                         verbose = TRUE,
                         save_h5ad = TRUE,
                         ncores = 1) {
  # TODO INPUT VALIDATION
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  outdir <- normalizePath(outdir)

  logger::log_info("Output directory: {outdir}")

  # sample directories
  # cnv_out <- file.path(outdir)
  # bins_out <- file.path(cnv_out, "binned_counts")
  # processed <- file.path(cnv_out, "processed")

  # Check for final output and skip if present
  final_out <- file.path(outdir, glue::glue("{bin_name}_processed.sce"))

  if (file.exists(final_out) & !overwrite) {
    logger::log_info("Processed output already exists: {final_out}")
    return(get(load(final_out)))
  }

  if (is.null(fragment_file)) {
    stop(logger::log_error("Must provide fragments.tsv"), call. = FALSE)
  }

  # Scatools processing
  logger::log_info("BINNING FRAGMENTS")

  # Bin the fragments
  bin_atac_frags(
    sample_id = sample_id,
    fragment_file = fragment_file,
    cells = cells,
    blacklist = blacklist,
    bins = bins,
    bin_name = bin_name,
    outdir = outdir,
    ncores = ncores,
    return_mat = FALSE,
    overwrite = overwrite
  )

  # Load the binned fragments and process
  sce <- load_atac_bins(
    bin_dir = file.path(outdir, bin_name),
    sample_id = sample_id,
    bins = bins,
    verbose = verbose
  )

  save_to(object = sce, save_to = file.path(outdir, glue::glue("{bin_name}_raw.sce")))

  # Remove any small leftover bins less than 10% of the full bin length
  sce_processed <- sce[!width(rowRanges(sce)) <= rmsb_size, ]

  # Length normalize so the tail end bins are in the same scale
  sce_processed <- length_normalize(sce_processed, assay_name = "counts", assay_to = "counts")

  sce_processed <- sce_processed %>%
    filter_sce(gc_range = gc_range) %>%
    add_ideal_mat(verbose = TRUE, ncores = ncores) %>%
    add_gc_cor(method = "modal", verbose = TRUE, ncores = ncores) %>%
    smooth_counts(assay_name = "counts_gc_modal", ncores = ncores) %>%
    calc_ratios(assay = "counts_gc_modal_smoothed") %>%
    logNorm(assay = "counts_gc_modal_smoothed_ratios", name = "logr_modal") %>%
    cluster_seurat(assay_name = "counts_gc_modal_smoothed_ratios", resolution = 0.5, verbose = FALSE) %>%
    segment_cnv(assay_name = "counts_gc_modal_smoothed", ncores = ncores) %>%
    merge_segments(smooth_assay = "counts_gc_modal_smoothed", segment_assay = "counts_gc_modal_smoothed_segment", ncores = ncores) %>%
    identify_normal(assay_name = "segment_merged_logratios", group_by = "clusters", method = "gmm")
  save_to(object = sce_processed, save_to = final_out)

  # Save anndata
  if (save_h5ad == TRUE) {
    if (requireNamespace("zellkonverter", quietly = T) & requireNamespace("anndata")) {
      logger::log_info("Writing output to anndata")
      # Save raw anndata
      adata <- zellkonverter::SCE2AnnData(sce)
      anndata::write_h5ad(adata, filename = file.path(outdir, glue::glue("{bin_name}_raw.h5ad")))
      adata2 <- zellkonverter::SCE2AnnData(sce_processed, verbose = TRUE, reducedDims = FALSE) # Bug in reduced dims conversion nbd
      anndata::write_h5ad(adata2, filename = file.path(outdir, glue::glue("{bin_name}_processed.h5ad")))
    }
  }
  return(sce_processed)
}
