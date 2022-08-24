# Create test fragments file with 100 cells
library(tidyverse)

# frags_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5589390&format=file&file=GSM5589390%5Fovary%5FSM%2DIOBHR%5Frep1%5Ffragments%2Ebed%2Egz"

# download.file(url = frags_url, destfile = "data-raw/GSM5589390_ovary_SM-IOBHR_rep1_fragments.bed")

target_file <- "./data-raw/GSM5589390_ovary_SM-IOBHR_rep1_fragments.bed"

bedfile <- read_table(target_file, col_names = F)

# Take the first 100 unique cells
bedfilt <- bedfile[which(bedfile$X4 %in% unique(bedfile$X4)[1:100]), ]

write_delim(bedfilt, file = "data-raw/fragments.bed", col_names = F)

ArchR::reformatFragmentFiles("data-raw/fragments.bed")

file.copy(from = "data-raw/fragments.bed.gz", to = "inst/extdata/fragments.bed.gz", overwrite = TRUE)
