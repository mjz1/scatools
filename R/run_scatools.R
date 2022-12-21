#' Run scATAC Copy Number Analysis
#'
#' Convenience wrapper to run a sample
#'
#' @param sample_id Sample id
#' @param fragments Fragments file
#' @param ArrowFile ArrowFile
#' @param outdir Outdir
#' @param overwrite Whether to overwrite previous results
#' @param force_arrow recreate arrow file
#' @param verbose Verbosity
#' @param bins Bins granges
#' @param ncores Number of cores
#'
#' @return `SingleCellExperiment` object
#' @export
#'
run_scatools <- function(sample_id,
                         fragments = NULL,
                         ArrowFile = NULL,
                         outdir = sample_id,
                         overwrite = FALSE,
                         force_arrow = FALSE,
                         verbose = TRUE,
                         bins,
                         ncores = 1) {

  # TODO INPUT VALIDATION

  bin_name <- prettyMb(getmode(width(bins)))

  # sample directories
  archr_dir <- file.path(outdir, "ArchR")
  cnv_out <- file.path(outdir, "cnv")
  bins_out <- file.path(cnv_out, "binned_counts")
  processed <- file.path(cnv_out, "processed")

  # Check for final output and skip if present
  final_out <- file.path(processed, glue("{bin_name}_processed.sce"))

  if (file.exists(final_out) & !overwrite) {
    logger::log_info("Processed output already exists: {final_out}")
    return(get(load(final_out)))
  }

  if (is.null(ArrowFile) & is.null(fragments)) {
    stop(logger::log_error("Must provide either fragments.tsv or ArrowFile"), call. = FALSE)
  }

  # Default arrow file
  if (is.null(ArrowFile)) {
    ArrowFile <- glue::glue("{archr_dir}/{sample_id}.arrow")
  }

  proj <- create_arrow_file(archr_dir = archr_dir,
                            ArrowFile = ArrowFile,
                            fragments = fragments,
                            sample_id = sample_id,
                            ncores = ncores,
                            force = force_arrow)

  # Binning and CNV processing
  dir.create(bins_out, recursive = TRUE, showWarnings = FALSE)
  dir.create(processed, recursive = TRUE, showWarnings = FALSE)

  # Scatools processing
  logger::log_info("BINNING FRAGMENTS")

  # Bin the fragments
  bin_atac_frags(
    ArrowFile = ArchR::getArrowFiles(proj),
    blacklist = ArchR::getBlacklist(proj),
    bins = bins,
    outdir = bins_out,
    ncores = ncores,
    return_mat = FALSE,
    overwrite = FALSE
  )

  # Load the binned fragments and process
  sce <- load_atac_bins(
    samples = file.path(bins_out, bin_name),
    sample.names = sample_id,
    ArchR_Proj = proj,
    bins = bins,
    BPPARAM = BiocParallel::bpparam(),
    verbose = verbose
  )

  sce <- flagDoublets(sce, filterRatio = 2, remove = FALSE)

  save_to(object = sce, save_to = file.path(processed, glue("{bin_name}_raw.sce")))

  # Remove any small leftover bins less than 10% of the full bin length
  sce_processed <- sce[rowRanges(sce)$binwidth > 0.1*getmode(width(bins)), ]

  # Length normalize so the tail end bins are in the same scale
  sce_processed <- length_normalize(sce_processed, assay_name = "counts", assay_to = "counts")

  # Flag and remove doublets
  sce_processed <- flagDoublets(sce_processed, filterRatio = 2, remove = TRUE)

  sce_processed <- sce_processed %>%
    filter_sce(gc_range = c(0.3, 1)) %>%
    add_ideal_mat(verbose = TRUE, ncores = ncores) %>%
    add_gc_cor(method = "modal", verbose = TRUE, ncores = ncores) %>%
    smooth_counts(assay_name = "counts_gc_modal", ncores = ncores) %>%
    calc_ratios(assay = "counts_gc_modal_smoothed") %>%
    copykit::logNorm(assay = "counts_gc_modal_smoothed_ratios", name = "logr_modal") %>%
    cluster_seurat(assay_name = "counts_gc_modal_smoothed_ratios", resolution = 0.5, verbose = FALSE) %>%
    segment_cnv(assay_name = "counts_gc_modal_smoothed", ncores = ncores) %>%
    merge_segments(smooth_assay = "counts_gc_modal_smoothed", segment_assay = "counts_gc_modal_smoothed_segment", ncores = ncores) %>%
    identify_normal(assay_name = "segment_merged_logratios", group_by = "clusters", method = "gmm")
  save_to(object = sce_processed, save_to = final_out)

  # Save anndata
  if (requireNamespace("zellkonverter", quietly = T) & requireNamespace("anndata")) {
    logger::log_info("Writing output to anndata")
    # Save raw anndata
    adata <- zellkonverter::SCE2AnnData(sce)
    anndata::write_h5ad(adata, filename = file.path(processed, glue("{bin_name}_raw.h5ad")))
    adata2 <- zellkonverter::SCE2AnnData(sce_processed, verbose = TRUE, reducedDims = FALSE) # Bug in reduced dims conversion nbd
    anndata::write_h5ad(adata2, filename = file.path(processed, glue("{bin_name}_processed.h5ad")))
  }

  return(sce_processed)
}


create_arrow_file <- function(fragments = NULL, ArrowFile, archr_dir, sample_id = "sample", genome = "hg38", ncores = 1, force = FALSE) {

  ArchR::addArchRThreads(ncores)
  ArchR::addArchRGenome(genome)

  dir.create(archr_dir, recursive = TRUE, showWarnings = FALSE)
  # Set working directory for ArchR
  # knitr::opts_current$set(root.dir = archr_dir)
  setwd(archr_dir)
  if (file.exists(ArrowFile) & !force) {
    logger::log_info("LOADING ARROWFILE")
    proj <- ArchR::ArchRProject(
      ArrowFile = ArrowFile,
      showLogo = FALSE,
      outputDirectory = archr_dir,
      copyArrows = FALSE # This is recommended so that you maintain an unaltered copy for later usage.
    )
    # Just in case no doublet calls
    if (is.null(proj$DoubletEnrichment)) {
      logger::log_info("CALCULATING DOUBLETS")
      doubScores <- ArchR::addDoubletScores(
        input = ArrowFile,
        k = 10, # Refers to how many cells near a "pseudo-doublet" to count.
        knnMethod = "LSI", # Refers to the embedding to use for nearest neighbor search.
        LSIMethod = 1,
        force = FALSE
      )
      # Reload
      proj <- ArchR::ArchRProject(
        ArrowFile = ArrowFile,
        outputDirectory = archr_dir,
        copyArrows = FALSE # This is recommended so that you maintain an unaltered copy for later usage.
      )
    }
  } else {
    logger::log_info("CREATING ARROWFILE")
    logger::log_info("Working in: ", getwd())
    ArrowFile <- ArchR::createArrowFiles(
      inputFiles = fragments,
      sampleNames = sample_id,
      minTSS = 4, # Dont set this too high because you can always increase later
      minFrags = 1000,
      addTileMat = TRUE,
      addGeneScoreMat = TRUE,
      force = force
    )
    logger::log_info("CALCULATING DOUBLETS")
    doubScores <- ArchR::addDoubletScores(
      input = ArrowFile,
      k = 10, # Refers to how many cells near a "pseudo-doublet" to count.
      knnMethod = "LSI", # Refers to the embedding to use for nearest neighbor search.
      LSIMethod = 1,
      force = FALSE
    )
    proj <- ArchR::ArchRProject(
      ArrowFile = ArrowFile,
      outputDirectory = archr_dir,
      copyArrows = FALSE # This is recommended so that you maintain an unaltered copy for later usage.
    )
  }
  return(proj)
}

