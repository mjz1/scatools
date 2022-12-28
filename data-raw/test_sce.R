## code to prepare `test_sce` dataset goes here
library(scatools)

ncores = parallelly::availableCores()

fragment_file <- system.file("extdata", "fragments.bed.gz", package = "scatools")
names(fragment_file) <- "test_sample"

# Generate bins
bins <- get_tiled_bins(bs_genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, tilewidth = 1e7)
bin_name <- prettyMb(getmode(width(bins)))

outdir <- tempdir()


test_sce <- run_scatools(sample_id = names(fragment_file),
                    fragments = fragment_file,
                    outdir = outdir,
                    verbose = TRUE,
                    overwrite = TRUE,
                    force_arrow = TRUE,
                    ncores = ncores,
                    bins = bins)

usethis::use_data(test_sce, overwrite = TRUE)
