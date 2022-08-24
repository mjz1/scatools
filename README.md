
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scatools

<!-- badges: start -->
<!-- badges: end -->

scatools is a package meant for performing copy number analysis in
scATAC data.

## Installation

You can install the development version of scatools from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mjz1/scatools")
```

## Example

We start by loading `scatools` and `ArchR`, as well as the filepath to
an example `fragments.bed.gz` file containing fragments from 100 normal
mammary cells. This bed file was generated using the
`reformatFragmentFiles()` function from the `ArchR` package.

Note that all steps in this vignette work with a list of fragment files
from multiple samples as well.

``` r
library(scatools)
library(ArchR, quietly = TRUE)
library(BSgenome.Hsapiens.UCSC.hg38, quietly = TRUE)
```

``` r
ncores <- 8 # Adjust accordingly
addArchRThreads(ncores)
addArchRGenome("hg38")

fragment_file <- system.file("extdata", "fragments.bed.gz", package = "scatools")
names(fragment_file) <- "test_sample"

# Set up example output directories
arrow_dir = "./example/ArrowFiles"
bindepth_dir = "./example/binned_depth"
scatools_dir = "./example/scatools_analysis"

invisible(lapply(list(arrow_dir, bindepth_dir, scatools_dir), dir.create, showWarnings = FALSE, recursive = TRUE))
```

We first create an ArrowFile from the fragments file using the `ArchR`
package. If you already have processed ArrowFiles you can skip to the
`scatools` processing steps.

``` r
knitr::opts_chunk$set(root.dir = arrow_dir)

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
```

Now we process this data using `scatools`. Helper functions help us to
create `GenomicRanges` bins, and compute GC content for downstream
usage. Here we demonstrate using 10Mb bins.

``` r
# Generate bins
bins <- get_tiled_bins(bs_genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, tilewidth = 1e7)

bin_name <- prettyMb(getmode(width(bins)))
message(bin_name)
#> 10Mb

head(bins)
#> GRanges object with 6 ranges and 4 metadata columns:
#>       seqnames            ranges strand |  binwidth                 bin_id
#>          <Rle>         <IRanges>  <Rle> | <integer>            <character>
#>   [1]     chr1        1-10000000      * |  10000000        chr1_1_10000000
#>   [2]     chr1 10000001-20000000      * |  10000000 chr1_10000001_20000000
#>   [3]     chr1 20000001-30000000      * |  10000000 chr1_20000001_30000000
#>   [4]     chr1 30000001-40000000      * |  10000000 chr1_30000001_40000000
#>   [5]     chr1 40000001-50000000      * |  10000000 chr1_40000001_50000000
#>   [6]     chr1 50000001-60000000      * |  10000000 chr1_50000001_60000000
#>              gc    n_freq
#>       <numeric> <numeric>
#>   [1]  0.497180 0.0203509
#>   [2]  0.475757 0.0100000
#>   [3]  0.474349 0.0001602
#>   [4]  0.462978 0.0000000
#>   [5]  0.440598 0.0000000
#>   [6]  0.423609 0.0000000
#>   -------
#>   seqinfo: 23 sequences from hg38 genome
```

Next we bin the atac fragments from the input ArrowFiles.

``` r
# Bin the fragments
bin_atac_frags(ArrowFiles = getArrowFiles(proj), bins = bins, outdir = bindepth_dir, ncores = ncores, overwrite = FALSE, return_mat = FALSE)
```

Then we set up to perform the analysis for all samples

``` r
samples <- file.path(bins_out, bin_name, sample_name)

for (i in seq_along(samples)) {
  samp_dir <- samples[i]
  samp_name <- sample_name[i]

  samp_outdir <- file.path(scatoolsdir, bin_name, samp_name)

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
  
  sce <- sce %>%
    add_ideal_mat(ncores = ncores, verbose = verbose) %>%
    add_gc_cor(method = "modal", verbose = verbose, ncores = ncores) %>%
    add_hmmcopy(verbose = verbose, ncores = ncores, save_raw_hmm = hmm_out)
  
  save_to(object = sce, save_to = final_out)
}
```

We can visualize the results as a heatmap.

``` r
cnaPlot <- cnaHeatmap(sce = sce, assay_name = "counts", log2 = TRUE, scale = "both", legend_name = "Log Normalized Counts", verbose = FALSE)
cnaPlot
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

Or plot individual cells

``` r
plot_cell_cna(sce = sce, cell_id = colnames(sce)[1:5], assay_name = "counts_gc_modal") +
  coord_cartesian(ylim = c(0, 10))
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />
