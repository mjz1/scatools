#' Run scATAC Copy Number Analysis
#'
#' Convenience wrapper to run a sample
#'
#' @param sample_id Sample id
#' @param fragment_file Fragments file
#' @param cells A vector of cell barcodes. If not provided will use all barcodes
#' @param bins Genomic ranges object with bins to use
#' @param bin_name bin_name Name of the bins (e.g. `'10Mb'`, `'500Kb'`, `'chr_arm'`). If not provided is automatically detected based on binwidth.
#' @param blacklist Genome blacklist regions to filter against
#' @param outdir Output directory
#' @param rmsb_size Remove small bins below a given size. Useful in cases where small bins are leftover at the ends of chromosomes. Defaults to 10% of the binwidth.
#' @param gc_range GC range for bins to keep. Removes large GC outliers
#' @param segment Logical. Determines whether to run segmentation or not
#' @param min_bin_counts A bin requires at least `min_bin_counts` across `min_bin_prop` proportion of cells to be kept
#' @param min_bin_prop Minimum proportion of cells with at least `min_bin_counts` per bin in order to keep a bin
#' @param min_cell_counts  A cell requires at least `min_cell_counts` across `min_cell_prop` proportion of bins to be kept
#' @param min_cell_prop Minimum proportion of bins with at least `min_cell_counts` per cell in order to keep a cell
#' @param overwrite Whether to overwrite previous binned counts
#' @param verbose Verbosity
#' @param save_h5ad Logical. Whether to save raw and processed h5ad files. Requires packages `zellkonverter` and `anndata`
#' @param ncores Number of cores
#' @param bpparam BiocParallel Parameters
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
                         segment = FALSE,
                         rmsb_size = 0.1 * getmode(width(bins)),
                         gc_range = c(0.3, 0.8),
                         overwrite = FALSE,
                         verbose = TRUE,
                         save_h5ad = TRUE,
                         min_bin_counts = 1,
                         min_bin_prop = 0.95,
                         min_cell_counts = min_bin_counts,
                         min_cell_prop = min_bin_prop,
                         ncores = 1,
                         bpparam = BiocParallel::SerialParam()) {
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
    bpparam = bpparam,
    return_mat = FALSE,
    overwrite = overwrite
  )

  # Load the binned fragments and process
  sce <- load_atac_bins(
    bin_dir = file.path(outdir, paste0(bin_name, "_counts")),
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
    filter_sce(
      gc_range = gc_range,
      min_bin_counts = min_bin_counts,
      min_bin_prop = min_bin_prop,
      min_cell_counts = min_cell_counts,
      min_cell_prop = min_cell_prop
    ) %>%
    add_ideal_mat(verbose = TRUE, ncores = ncores) %>%
    add_gc_cor(method = "modal", verbose = TRUE, bpparam = bpparam) %>%
    smooth_counts(assay_name = "counts_gc_modal", ncores = ncores) %>%
    calc_ratios(assay_name = "counts_gc_modal_smoothed") %>%
    logNorm(assay_name = "counts_gc_modal_smoothed_ratios", name = "logr_modal") %>%
    cluster_seurat(assay_name = "counts_gc_modal_smoothed_ratios", resolution = 0.5, verbose = FALSE)

  if (segment) {
    sce_processed <- segment_cnv(assay_name = "counts_gc_modal_smoothed", bpparam = bpparam) %>%
      merge_segments(smooth_assay = "counts_gc_modal_smoothed", segment_assay = "counts_gc_modal_smoothed_segment", bpparam = bpparam) %>%
      identify_normal(assay_name = "segment_merged_logratios", group_by = "clusters", method = "gmm")
  }

  save_to(object = sce_processed, save_to = final_out)
  logger::log_success("SCATools run completed!")

  # Save anndata
  if (save_h5ad == TRUE) {
    if (!requireNamespace("zellkonverter", quietly = T)) {
      logger::log_error("Package 'zellkonverter' must be installed to save as h5ad")
    } else {
      logger::log_info("Writing output to anndata")
      # Save raw anndata
      zellkonverter::writeH5AD(sce, file = file.path(outdir, glue::glue("{bin_name}_raw.h5ad")), compression = "gzip")
      zellkonverter::writeH5AD(sce_processed, file = file.path(outdir, glue::glue("{bin_name}_processed.h5ad")), compression = "gzip")
    }
  }

  return(sce_processed)
}
