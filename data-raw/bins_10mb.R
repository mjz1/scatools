## code to prepare `bins_10mb` dataset goes here

bins_10mb <- get_tiled_bins(bs_genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, tilewidth = 1e7)


usethis::use_data(bins_10mb, overwrite = TRUE)
