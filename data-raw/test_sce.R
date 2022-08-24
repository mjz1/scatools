## code to prepare `test_sce` dataset goes here

library(scatools)
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)

ncores <- parallelly::availableCores()[[1]]
overwrite = TRUE
verbose = TRUE

addArchRThreads(ncores)
addArchRGenome("hg38")

repo_dir <- getwd()

fragment_file <- file.path(repo_dir, "/data-raw/fragments.bed.gz")
names(fragment_file) <- "test_sample"

# Set up example output directories
arrow_dir = file.path(repo_dir, "/data-raw/example/ArrowFiles")
bindepth_dir = file.path(repo_dir, "/data-raw/example/binned_depth")
scatools_dir = file.path(repo_dir, "/data-raw/example/scatools_analysis")

invisible(lapply(list(arrow_dir, bindepth_dir, scatools_dir), dir.create, showWarnings = FALSE, recursive = TRUE))

setwd(arrow_dir)

ArrowFiles <- createArrowFiles(
  inputFiles = fragment_file,
  sampleNames = names(fragment_file),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = FALSE
)

# Calculate doublet scores
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "LSI", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1,
  force = FALSE
)

# Create an ArchR project file
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "./",
  copyArrows = FALSE #This is recommended so that you maintain an unaltered copy for later usage.
)

proj <- filterDoublets(proj)

# Generate bins
bins <- get_tiled_bins(bs_genome = BSgenome.Hsapiens.UCSC.hg38, tilewidth = 1e7)

bin_name <- prettyMb(getmode(width(bins)))

# Bin the fragments
bin_atac_frags(ArrowFiles = getArrowFiles(proj), bins = bins, outdir = bindepth_dir, ncores = ncores, overwrite = FALSE, return_mat = FALSE)

samples <- file.path(bindepth_dir, bin_name, names(fragment_file))

for (i in seq_along(samples)) {
  samp_dir <- samples[i]
  samp_name <- basename(samp_dir)

  samp_outdir <- file.path(scatools_dir, bin_name, samp_name)

  raw_out <- file.path(samp_outdir, "sce", "01_raw.sce")
  final_out <- file.path(samp_outdir, "sce", "02_hmm.sce")
  hmm_out <- file.path(samp_outdir, "hmm", "hmm_results.rda")

  logger::log_info("Processing sample {i} of {length(samples)}: {samp_name}")

  if (file.exists(final_out) & !overwrite) {
    logger::log_info("Final output exists! Skipping to next sample...")
    next
  }

  if (file.exists(raw_out) & !overwrite) {
    logger::log_info("Raw sce object found -- Loading...")
    sce <- get(load(raw_out))
  } else {
    sce <- load_atac_bins(
      samples = samp_dir,
      sample.names = samp_name,
      ArchR_Proj = proj,
      bins = bins,
      BPPARAM = BiocParallel::bpparam(),
      save_to = raw_out,
      verbose = verbose
    )
  }

  test_sce <- sce %>%
    add_ideal_mat(ncores = ncores, verbose = verbose) %>%
    add_gc_cor(method = "modal", verbose = verbose, ncores = ncores) %>%
    add_hmmcopy(verbose = verbose, ncores = ncores, save_raw_hmm = hmm_out)

  save_to(object = sce, save_to = final_out)
}

usethis::use_data(test_sce, overwrite = overwrite)
